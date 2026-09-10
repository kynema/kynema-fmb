#include <array>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "math/gl_quadrature.hpp"
#include "math/gll_quadrature.hpp"
#include "model/model.hpp"
#include "step/extract_system_matrices.hpp"
#include "step/step.hpp"
#include "test_utilities.hpp"

namespace {
template <typename T>
void WriteMatrixToFile(const std::vector<std::vector<T>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << "\n";
        return;
    }
    for (const auto& innerVector : data) {
        for (const auto& element : innerVector) {
            file << element << ",";
        }
        file << "\n";
    }
    file.close();
}

}  // namespace

namespace kynema_fmb::tests {

TEST(DynamicBeamTest, Damping) {
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
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
        std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
        std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
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

    const double scalar_mu(0.0001);  // 1/s
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
    // step_size and num_steps are updated below once the mode's natural frequency is known.
    constexpr size_t steps_per_cycle = 128;
    const double num_cycles = 1.5;
    double step_size(0.0001);  // seconds
    const double rho_inf(1.0);

    // Create solver parameters (step_size is refined after the eigenanalysis).
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    // Eventually will want this to only track the tip displacement
    // // Get ID of last node
    // const auto last_node_id = beam_node_ids.back();

    // Ideally would pull from an eigenanalysis to get the mode shapes
    // and initial velocity, just prototyping here.
    // (x/L)^3 is a rough approx of first bending mode,
    // technically need rotational velocity as well

    /*
    OLD Initialization
    // Start the beam with an initial velocity in the first component
    for (size_t i = 0; i < beam_node_ids.size(); ++i) {
        auto node_v = Kokkos::subview(state.v, beam_node_ids[i], Kokkos::ALL);
        node_v(0) = node_s[i] * node_s[i] * node_s[i];
        node_v(1) = 0.;
        node_v(2) = 0.;
        node_v(3) = 0.;
        node_v(4) = 0.;
        node_v(5) = 0.;
    }
    */


    // Mode to use for the initial velocity shape.
    const size_t mode_ind = 0;

    // Extract the system matrices at the initial (undeformed) state so we can
    // solve the generalized eigenvalue problem K v = lambda M v.
    auto matrices = step::ExtractSystemMatrices(parameters, solver, elements, state, constraints);

    const auto row_map =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solver.A.graph.row_map);
    const auto col_ids =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solver.A.graph.entries);
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
    const double omega = std::sqrt(eigenvalue);
    Eigen::VectorXd eigvec = raw_eigenvectors.col(sorted_col).real();

    // Find direction of mode shape
    // Only check the first 4 of the last 6 components 
    // interested in translation + torsion (not bending rotations)
    Eigen::Index dominant_tip_dof = 0;
    eigvec.tail(6).head(4).cwiseAbs().maxCoeff(&dominant_tip_dof);
    const size_t dir_ind = static_cast<size_t>(dominant_tip_dof);

    // Scale so |last node value in the dominant direction| / sqrt(eigenvalue) = 1e-3.
    const double max_last_abs = eigvec.tail(6)(static_cast<Eigen::Index>(dir_ind));
    eigvec *= (1.0e-3 * omega / max_last_abs);

    // Set initial velocity from the (scaled) mode shape. The root node stays fixed.
    for (size_t i = 1; i < beam_node_ids.size(); ++i) {
        auto node_v = Kokkos::subview(state.v, beam_node_ids[i], Kokkos::ALL);
        for (int j = 0; j < 6; ++j) {
            node_v(j) = eigvec(static_cast<Eigen::Index>((i - 1) * 6 + j));
        }
    }

    // Set num_steps to cover exactly 1.5 cycles of the natural frequency
    // note, this is the undamped natural frequency so will be slightly off from
    // a perfect (damped) 1.5 cycles.
    step_size = (2. * std::numbers::pi) / (omega * static_cast<double>(steps_per_cycle));
    const auto num_steps = static_cast<int>(steps_per_cycle*num_cycles);
    parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Time Stepping Loop for num_steps time steps saving the response at each step

    std::vector<std::vector<double>> displacement_history;
    displacement_history.reserve(num_steps);
    std::vector<std::vector<double>> velocity_history;
    velocity_history.reserve(num_steps);

    for ([[maybe_unused]] auto i : std::views::iota(0, num_steps)) {
        auto converged = Step(parameters, solver, elements, state, constraints);
        EXPECT_TRUE(converged);

        const auto q_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.q);
        auto& displacement_row = displacement_history.emplace_back();
        displacement_row.reserve(1 + q_host.extent(0) * q_host.extent(1));
        displacement_row.push_back(static_cast<double>(state.time_step) * step_size);
        for (size_t node = 0; node < q_host.extent(0); ++node) {
            for (size_t dof = 0; dof < q_host.extent(1); ++dof) {
                displacement_row.push_back(q_host(node, dof));
            }
        }

        const auto v_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.v);
        auto& velocity_row = velocity_history.emplace_back();
        velocity_row.reserve(1 + v_host.extent(0) * v_host.extent(1));
        velocity_row.push_back(static_cast<double>(state.time_step) * step_size);
        for (size_t node = 0; node < v_host.extent(0); ++node) {
            for (size_t dof = 0; dof < v_host.extent(1); ++dof) {
                velocity_row.push_back(v_host(node, dof));
            }
        }
    }

    WriteMatrixToFile(displacement_history, "beam_damping_displacement_history.csv");
    WriteMatrixToFile(velocity_history, "beam_damping_velocity_history.csv");

    // Calculate the damping values based on log decrement
    // and verify against the analytical stiffness proportional damping


    // 3 translations + 4 quaternions per node
    const auto tip_displacement_col = 1 + dir_ind + 7*(beam_node_ids.size() - 1);

    std::vector<size_t> peak_inds;
    // Assume that should have at least 50 steps per cycle
    // for sizing the memory (100 is rule of thumb to be accurate)
    peak_inds.reserve(num_steps / 50);

    for (size_t i = 1; i + 1 < displacement_history.size(); ++i) {
        const auto prev = displacement_history[i - 1][tip_displacement_col];
        const auto curr = displacement_history[i][tip_displacement_col];
        const auto next = displacement_history[i + 1][tip_displacement_col];
        if (curr > 0. && curr > prev && curr > next) {
            // Require peaks > 0 to avoid undefined log below.
            peak_inds.push_back(i);
        }
    }

    // Need at least two peaks for log decrement
    ASSERT_GE(peak_inds.size(), 2U);

    std::vector<double> peak_vals;
    peak_vals.reserve(peak_inds.size());

    // Peak times are really just (peak_inds + 1) * step_size
    std::vector<double> peak_times;
    peak_times.reserve(peak_inds.size());
    
    for (const auto peak_ind : peak_inds) {
        peak_vals.push_back(displacement_history[peak_ind][tip_displacement_col]);
        peak_times.push_back(displacement_history[peak_ind][0]);
    }

    std::vector<double> zeta;
    zeta.reserve(peak_vals.size() - 1);

    std::vector<double> omega_n;
    omega_n.reserve(peak_vals.size() - 1);

    constexpr auto two_pi = 2. * std::numbers::pi;
    
    for (size_t i = 0; i + 1 < peak_vals.size(); ++i) {

        const auto delta_i = std::log(peak_vals[i] / peak_vals[i + 1]);
        
        const auto zeta_i = delta_i / std::sqrt(two_pi * two_pi +
                                    delta_i * delta_i);

        const auto omega_n_i = (two_pi / (peak_times[i + 1] - peak_times[i]))
                        / std::sqrt(1. - zeta_i * zeta_i);

        zeta.push_back(zeta_i);
        omega_n.push_back(omega_n_i);
    }

    // Check only the last few damping values starting from the end
    // The initialization is imperfect, so the early ones have other transients
    // If above gets updated to an eigensolution for the initialization
    // then should be able to run fewer cycles and immediately evaluate
    // the damping and at tighter tolerance.
    for (size_t eval_count = 0, idx = zeta.size() - 1;
            eval_count < 3 && idx > 4;
            ++eval_count, --idx) 
    {
        // Needs to actually refer to the array instead of the scalar here.
        ASSERT_NEAR(zeta[idx], scalar_mu * omega_n[idx] / 2.0, 1e-3*zeta[idx]);
    }
}

}  // namespace kynema_fmb::tests
