#pragma once

#include <array>
#include <fstream>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "interfaces/generalized/coupling_data.hpp"
#include "interfaces/host_constraints.hpp"
#include "interfaces/host_state.hpp"
#include "interfaces/outputs.hpp"
#include "interfaces/turbine/turbine_interface.hpp"
#include "model/model.hpp"
#include "step/step_parameters.hpp"

namespace kynema_fmb::interfaces {

/**
 * @brief Interface for structural components of FMB for a generalized multi-physics framework.
 *
 * This class represents the generalized interface for structural simulations
 * used coupled to the Kynema suite of codes.
 */
class GeneralStructureInterface {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    template <typename ValueType>
    using HostView = typename Kokkos::View<ValueType, DeviceType>::host_mirror_type;

    /**
     * @brief Constructs a GeneralStructureInterface from a yaml description
     * @param structural_description Configuration parameters for the structure
     * @param id Module instance identifier
     * @param enable_output Whether this instance writes diagnostic messages and configured CSV
     * output
     */
    explicit GeneralStructureInterface(
        const YAML::Node& structural_description, int id = 1, bool enable_output = true
    );

    void prepare();
    nlohmann::json initialization_data() const;
    void receive_inputs(const nlohmann::json& inputs);
    nlohmann::json outputs() const;
    std::array<double, 7> node_position(size_t node) const;
    size_t num_nodes() const { return this->state.x.extent(0); }
    void set_inputs(size_t iteration);
    int time_index() const { return this->module_step; }
    std::string identifier() const { return "KynemaFMB-" + std::to_string(this->module_id); }
    double get_time() const { return this->module_time; }
    double get_timestep_size() const { return this->module_dt; }
    void set_timestep_size(double dt);
    void additional_picard_iterations(int count);
    void dump_simulation_time();

    GeneralStructureInterface(GeneralStructureInterface& other) = delete;

    GeneralStructureInterface(GeneralStructureInterface&& other) = delete;

    GeneralStructureInterface& operator=(const GeneralStructureInterface&) = delete;

    GeneralStructureInterface& operator=(GeneralStructureInterface&&) = delete;

    ~GeneralStructureInterface() = default;

    /**
     * @brief Do any additional initialization and return a yaml of node ids etc.
     */
    YAML::Node initialize_module();

    /**
     * @brief Predict the inputs to the model at time n+1 from stored data.
     */
    void predict_inputs(size_t iteration = 0);

    /**
     * @brief Reset the state to the previous saved time.
     * If there is no substepping, this is just time n
     * If there is substepping, this is to the start of the largest time step.
     */
    void reset_states();

    /**
     * @brief Solve for the response at time n+1
     */
    void advance_states(size_t iteration = 0);

    /**
     * @brief Apply a host force view directly for standalone callers.
     */
    void receive_inputs(HostView<double* [6]> forces);

    /**
     * @brief Copy structural state to the host for output.
     */
    void calculate_outputs(size_t iteration = 0);

    /**
     * @brief Return the host state for standalone callers.
     */
    HostState<DeviceType>& get_outputs();

    /**
     * @brief Finalize the stored data for the time step and shift time step index.
     */
    void finalize_time_step();

    /*
    1. `initialize_module`
    2. `initialize_stage1` - have multiple extra stages based on information sharing at init.
    3. `predict_inputs`
    4. `advance_states`
    5. `calculate_outputs`
    6. `finalize_time_step`

    */

private:
    int module_id;
    bool output_enabled;
    bool prepared{false};
    int module_step{0};
    double module_time{0.};
    double module_dt{0.};
    StructuralInputs module_inputs;
    std::string output_filename;
    std::ofstream output_file;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;

    YAML::Node model_description;
    Model model;                    ///< Kynema-FMB class for model construction
    State<DeviceType> state;        ///< Kynema-FMB class for storing system state (solved for n+1)
    Elements<DeviceType> elements;  ///< Kynema-FMB class for model elements (beams, masses, springs)
    Constraints<DeviceType>
        constraints;                   ///< Kynema-FMB class for constraints tying elements together
    StepParameters parameters;         ///< Kynema-FMB class containing solution parameters
    Solver<DeviceType> solver;         ///< Kynema-FMB class for solving the dynamic system
    State<DeviceType> state_n;         ///< Kynema-FMB class state for storing data at time n
    View<double* [6]> f_nm1;           //< External forces at time n-1
    HostState<DeviceType> host_state;  ///< Host local copy of node state data
    HostConstraints<DeviceType> host_constraints;  ///< Host local copy of constraint data

    /**
     @brief Initializes everything for FMB
     @param structural_description Configuration parameters for the structure
     @param model_description Summary of information and indices of as built model
     */
    static Model BuildModel(const YAML::Node& structural_description, YAML::Node& model_description);

    static void Add6DOF(
        Model& model, const YAML::Node& dof6_description, YAML::Node& output_description
    );
};

}  // namespace kynema_fmb::interfaces
