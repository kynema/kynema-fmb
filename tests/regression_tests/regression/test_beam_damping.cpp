#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "math/gl_quadrature.hpp"
#include "math/gll_quadrature.hpp"
#include "model/model.hpp"
#include "step/extract_system_matrices.hpp"
#include "step/step.hpp"
#include "test_utilities.hpp"

namespace kynema_fmb::tests {

class DynamicBeamTest : public ::testing::TestWithParam<size_t> {};

TEST_P(DynamicBeamTest, Damping) {
    // Mass matrix for uniform composite beam section
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},
        std::array{0., 8.538e-2, 0., 0., 0., 0.},
        std::array{0., 0., 8.538e-2, 0., 0., 0.},
        std::array{0., 0., 0., 1.4433e-2, 0., 0.},
        std::array{0., 0., 0., 0., 0.40972e-2, 0.},
        std::array{0., 0., 0., 0., 0., 1.0336e-2},
    };

    // Stiffness matrix for uniform composite beam section
    // Made diagonal so can test directional mu.
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 0., 0.},
        std::array{0., 0., 0., 0., 59.1200e3, 0.},
        std::array{0., 0., 0., 0., 0., 141.470e3},
    };
    

    // Node locations (GLL quadrature)
    const auto num_nodes = 5UL;
    const auto gll_locations = math::GetGllLocations(num_nodes - 1);
    std::vector<double> node_s(gll_locations.size());
    std::ranges::transform(gll_locations, node_s.begin(), [](auto xi) { return 0.5 * (xi + 1.0); });

    // Create model for managing nodes and constraints
    auto model = Model();

    // Set gravity in model
    model.SetGravity(0., 0., 0.);

    // Build vector of nodes (straight along x axis, no rotation)
    std::vector<size_t> beam_node_ids;
    std::ranges::transform(node_s, std::back_inserter(beam_node_ids), [&](auto s) {
        return model.AddNode()
            .SetElemLocation(s)
            .SetPosition(10 * s, 0., 0., 1., 0., 0., 0.)
            .Build();
    });

    const auto array_mu = std::array{0.0001, 0.0004, 0.0002,
                                     0.0003, 0.0002, 0.0004}; // 1/s

    const auto quad_order = 7UL;
    const auto gl_locations = math::GetGlLocations(quad_order);
    const auto gl_weights = math::GetGlWeights(quad_order);

    std::vector<std::array<double, 2>> quad_points(gl_locations.size());
    for (size_t i = 0; i < gl_locations.size(); ++i) {
        quad_points[i] = {gl_locations[i], gl_weights[i]};
    }

    // Add beam element
    model.AddBeamElement(
        beam_node_ids,
        std::array{
            BeamSection(0., mass_matrix, stiffness_matrix),
            BeamSection(1., mass_matrix, stiffness_matrix),
        },
        quad_points,
        array_mu
    );

    // Fix first node position
    model.AddFixedBC(beam_node_ids[0]);

    // Solution parameters
    const bool is_dynamic_solve(true);
    const size_t max_iter(5);
    constexpr size_t steps_per_cycle = 256;
    const double num_cycles = 1.5;
    const auto num_steps = static_cast<int>(steps_per_cycle*num_cycles);
    double step_size(0.0001);  // seconds, updated below based on natural frequency
    const double rho_inf(1.0); // no numerical damping, want accurate model damping

    // Create solver parameters (step_size is updated after the eigenanalysis).
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Mode to use for the initial velocity shape.
    const size_t mode_ind = GetParam();

    // Extract the system matrices at the initial (undeformed) state so we can
    // solve the generalized eigenvalue problem K v = lambda M v.
    auto matrices = step::ExtractSystemMatrices(parameters, solver, elements, state, constraints);

    const auto row_map =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrices.row_map);
    const auto col_ids =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrices.col_indices);
    const auto mass_vals =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrices.mass_matrix_values);
    const auto stiff_vals =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrices.stiffness_matrix_values);

    // Apply the fixed root constraint by eliminating the first node's DOFs.
    const auto num_reduced_dofs = 6 * (beam_node_ids.size() - 1);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(num_reduced_dofs, num_reduced_dofs);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(num_reduced_dofs, num_reduced_dofs);
    for (size_t row = 6; row < num_reduced_dofs + 6; ++row) {
        for (auto j = row_map(row); j < row_map(row + 1); ++j) {
            const auto col = static_cast<size_t>(col_ids(j));
            if (col >= 6 && col < num_reduced_dofs + 6) {
                M(row - 6, col - 6) = mass_vals(j);
                K(row - 6, col - 6) = stiff_vals(j);
            }
        }
    }

    const Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> eigensolver(K, M);
    const auto raw_eigenvalues = eigensolver.eigenvalues();
    const auto raw_eigenvectors = eigensolver.eigenvectors();

    // Sort eigenpairs by ascending real part of the eigenvalue.
    std::vector<Eigen::Index> sort_indices(raw_eigenvalues.size());
    std::iota(sort_indices.begin(), sort_indices.end(), 0);
    std::stable_sort(
        sort_indices.begin(), sort_indices.end(),
        [&](Eigen::Index i, Eigen::Index j) {
            return raw_eigenvalues(i).real() < raw_eigenvalues(j).real();
        }
    );

    ASSERT_LT(mode_ind, static_cast<size_t>(raw_eigenvalues.size()));
    const auto sorted_col = sort_indices[mode_ind];
    const double eigenvalue = raw_eigenvalues(sorted_col).real();
    Eigen::VectorXd eigvec = raw_eigenvectors.col(sorted_col).real();

    // Find direction of mode shape
    // Only check the first 4 of the last 6 components 
    // interested in translation + torsion (not bending rotations)
    Eigen::Index dominant_tip_dof = 0;
    eigvec.tail(6).head(4).cwiseAbs().maxCoeff(&dominant_tip_dof);
    const size_t dir_ind = static_cast<size_t>(dominant_tip_dof);

    // Determine which tip DOF corresponds to the current mode
    // 3 translations + 4 quaternions per node
    auto tip_dof = dir_ind;
    if (dir_ind >= 3) {
        // need to go past the first entry of the quaternion if it is a rotation.
        tip_dof++;
    }

    // Analytical modal properties from eigensolution / stiffness prop damping
    const double omega_eig = std::sqrt(eigenvalue);
    const double zeta_eig = array_mu[dir_ind] * omega_eig / 2.0;

    // Scale so (last node value in the dominant direction) / sqrt(eigenvalue) = 1e-3.
    const double tip_amplitude = eigvec.tail(6)(static_cast<Eigen::Index>(dir_ind));
    eigvec *= (1.0e-3 * omega_eig / tip_amplitude);

    // Set initial velocity from the (scaled) mode shape. The root node stays fixed.
    for (size_t i = 1; i < beam_node_ids.size(); ++i) {
        auto node_v = Kokkos::subview(state.v, beam_node_ids[i], Kokkos::ALL);
        for (int j = 0; j < 6; ++j) {
            node_v(j) = eigvec(static_cast<Eigen::Index>((i - 1) * 6 + j));
        }
    }

    // Set num_steps to cover exactly 1.5 cycles of the natural frequency
    // note, this is with the damped natural frequency.
    step_size = (2. * std::numbers::pi) / (omega_eig * std::sqrt(1 - zeta_eig * zeta_eig) 
                                            * static_cast<double>(steps_per_cycle));
    parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Time Stepping Loop for num_steps time steps saving the response at each step
    std::vector<double> tip_dof_history;
    tip_dof_history.reserve(num_steps);

    for ([[maybe_unused]] auto i : std::views::iota(0, num_steps)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);

        const auto q_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
        tip_dof_history.push_back(q_host(beam_node_ids.back(), tip_dof));

    }

    // Calculate the damping values based on log decrement
    // and verify against the analytical stiffness proportional damping
    // Also check the natural frequency

    // Only simulating two peaks
    std::vector<double> peak_times; // Peak times are just (peak_inds + 1) * step_size
    peak_times.reserve(2);
    std::vector<double> peak_vals;
    peak_vals.reserve(2);

    for (size_t i = 1; i + 1 < tip_dof_history.size(); ++i) {
        const auto prev = tip_dof_history[i - 1];
        const auto curr = tip_dof_history[i];
        const auto next = tip_dof_history[i + 1];
        if (curr > 0. && curr > prev && curr > next) {
            // Require peaks > 0 to avoid undefined log below.
            peak_times.push_back(step_size * (i+1));
            peak_vals.push_back(tip_dof_history[i]);
        }
    }

    // Should have exactly two peaks by doing 1.5 cycles
    // starting towards the positive direction.
    ASSERT_EQ(peak_times.size(), 2U);

    constexpr auto two_pi = 2. * std::numbers::pi;

    const auto delta = std::log(peak_vals[0] / peak_vals[1]);
    const auto zeta_num = delta / std::sqrt(two_pi * two_pi + delta * delta);
    const auto omega_num = (two_pi / (peak_times[1] - peak_times[0]))
                         / std::sqrt(1. - zeta_num * zeta_num);

    // _num = time integration, _eig = eigensolution
    // Tolerances that pass with 256 time steps/cycle
    ASSERT_NEAR(zeta_num, zeta_eig, 4e-4 * zeta_eig);
    ASSERT_NEAR(omega_num, omega_eig, 1e-6 * omega_eig);

    // // Tolerances that pass with 1024 time steps/cycle
    // ASSERT_NEAR(zeta_num, zeta_eig, 2e-5 * zeta_eig);
    // ASSERT_NEAR(omega_num, omega_eig, 1e-6 * omega_eig);

    // // These tolerances pass with 8192 steps/cycle for first 6 modes,
    // // helps verify convergence, but keeping test at faster/fewer steps
    // // with above tolerances.
    // ASSERT_NEAR(zeta_num, zeta_eig, 2e-6 * zeta_eig);
    // ASSERT_NEAR(omega_num, omega_eig, 1e-8 * omega_eig);
}

INSTANTIATE_TEST_SUITE_P(
    ModeIndices, DynamicBeamTest,
    ::testing::Values(0U, 1U, 2U, 3U, 4U, 8U), // 3U is torsion, 8U is axial
    [](const ::testing::TestParamInfo<DynamicBeamTest::ParamType>& info) {
        return "mode" + std::to_string(info.param);
    }
);

}  // namespace kynema_fmb::tests
