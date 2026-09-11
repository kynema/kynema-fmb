#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "assemble_constraints_matrix.hpp"
#include "assemble_system_matrix.hpp"
#include "constraints/constraints.hpp"
#include "elements/elements.hpp"
#include "reset_solver.hpp"
#include "solver/solver.hpp"
#include "state/state.hpp"
#include "step_parameters.hpp"
#include "update_constraint_variables.hpp"
#include "update_system_variables.hpp"
#include "update_tangent_operator.hpp"

namespace kynema_fmb {

template <typename DeviceType>
struct SystemMatrices {
    using CrsMatrixType = typename Solver<DeviceType>::CrsMatrixType;
    using ValuesType = typename CrsMatrixType::values_type::non_const_type;
    using RowMapType = typename CrsMatrixType::staticcrsgraph_type::row_map_type;
    using EntriesType = typename CrsMatrixType::staticcrsgraph_type::entries_type;

    // CRS graph shared by all assembled matrices
    RowMapType row_map;
    EntriesType col_indices;

    ValuesType mass_matrix_values;
    ValuesType stiffness_matrix_values;
    ValuesType damping_matrix_values;
    ValuesType constraint_matrix_values;
};

namespace step {

template <typename DeviceType>
inline SystemMatrices<DeviceType> ExtractSystemMatrices(
    const StepParameters& base_parameters, Solver<DeviceType>& solver,
    Elements<DeviceType>& elements, State<DeviceType>& state,
    Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("Extract System Matrices");

    using ValuesType = typename SystemMatrices<DeviceType>::ValuesType;
    SystemMatrices<DeviceType> result;
    const auto num_values = solver.A.values.extent(0);

    // Capture the shared sparsity pattern that is created with the solver.
    result.row_map = solver.A.graph.row_map;
    result.col_indices = solver.A.graph.entries;

    // Tangent depends only on state and base parameters — compute once
    auto params_for_tangent = base_parameters;
    params_for_tangent.h = 0.0; // effectively turn off the tangent operator
    step::UpdateTangentOperator(params_for_tangent, state);

    // --- Mass pass ---
    {
        auto params = base_parameters;
        params.h = 0.0;
        params.beta_prime = 1.0;
        params.gamma_prime = 0.0;
        params.conditioner = 1.0;
        params.include_stiffness = false;

        step::ResetSolver(solver);
        step::UpdateSystemVariables(params, elements, state);
        step::AssembleSystemMatrix(params, solver, elements);

        result.mass_matrix_values = ValuesType("mass_values", num_values);
        Kokkos::deep_copy(result.mass_matrix_values, solver.A.values);
    }

    // --- Stiffness pass ---
    {
        auto params = base_parameters;
        params.h = 0.0;
        params.beta_prime = 0.0;
        params.gamma_prime = 0.0;
        params.conditioner = 1.0;
        params.include_stiffness = true;

        step::ResetSolver(solver);
        step::UpdateSystemVariables(params, elements, state);
        step::AssembleSystemMatrix(params, solver, elements);

        result.stiffness_matrix_values = ValuesType("stiffness_values", num_values);
        Kokkos::deep_copy(result.stiffness_matrix_values, solver.A.values);
    }

    // --- Damping pass ---
    {
        auto params = base_parameters;
        params.h = 0.0;
        params.beta_prime = 0.0;
        params.gamma_prime = 1.0;
        params.conditioner = 1.0;
        params.include_stiffness = false;

        step::ResetSolver(solver);
        step::UpdateSystemVariables(params, elements, state);
        step::AssembleSystemMatrix(params, solver, elements);

        result.damping_matrix_values = ValuesType("damping_values", num_values);
        Kokkos::deep_copy(result.damping_matrix_values, solver.A.values);
    }

    // --- Constraint pass ---
    {
        step::ResetSolver(solver);
        step::UpdateConstraintVariables(state, constraints);
        step::AssembleConstraintsMatrix(solver, constraints);

        result.constraint_matrix_values = ValuesType("constraint_values", num_values);
        Kokkos::deep_copy(result.constraint_matrix_values, solver.A.values);
    }

    // --- Restore: re-run with original parameters so element
    //     storage is consistent for the next solve step ---
    {
        auto params = base_parameters;
        step::ResetSolver(solver);
        step::UpdateSystemVariables(params, elements, state);
        step::AssembleSystemMatrix(params, solver, elements);
    }

    return result;
}

}  // namespace step
}  // namespace kynema
