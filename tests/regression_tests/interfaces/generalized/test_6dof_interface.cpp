#include <array>
#include <fstream>
#include <numbers>
#include <ranges>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "elements/beams/hollow_circle_properties.hpp"
#include "interfaces/generalized/generalized_structure_interface.hpp"

namespace kynema_fmb::tests {

/*
This is intended to test a single 6 DOF model run through
the general structural interface
*/
TEST(GeneralizedInterface, DOF6) {

    const YAML::Node input_yaml =
        YAML::LoadFile("generalized_test_files/general_6dof_inputs.yaml");

    const auto duration =
        input_yaml["simulation"]["duration"].as<double>();  // Simulation duration in seconds
    const auto time_step =
        input_yaml["simulation"]["step_size"].as<double>();  // Time step for the simulation
    const auto write_output{true};                           // Write output file

    auto structure_interface = interfaces::GeneralStructureInterface(input_yaml);

    const int num_steps(40);

    // Open a CSV to record the first row of host_state.x (node 0) at each time step
    std::ofstream x_output;
    if (write_output) {
        x_output.open("general_6dof_host_state_x.csv");
        x_output << std::scientific;
        x_output << "time_step,x,y,z,qw,qx,qy,qz\n";
    }

    YAML::Node fmb_model_description = structure_interface.initialize_module();


    auto num_nodes = fmb_model_description["model"]["num_nodes"].as<size_t>();
    auto mass_node_id = fmb_model_description["components"][0]["mass_node_id"].as<size_t>();

    // Right now, I just know that there are 2 nodes because the 6 dof gets a node
    // and then each spring gets an additional node for the fixed end of the spring.
    // This should probably be passed out from the initialization at some point.
    auto inputs = interfaces::StructuralInputs{std::vector<std::array<double, 6>>(num_nodes)};
    nlohmann::json invalid = inputs;
    invalid["version"] = interfaces::structural_schema_version + 1;
    EXPECT_THROW(structure_interface.receive_inputs(invalid), std::invalid_argument);
    invalid = inputs;
    invalid["nodal_forces"][0].push_back(0.);
    EXPECT_THROW(structure_interface.receive_inputs(invalid), std::invalid_argument);
    invalid = inputs;
    invalid["nodal_forces"].erase(0);
    EXPECT_THROW(structure_interface.receive_inputs(invalid), std::invalid_argument);

    // TODO : Response looks a bit funny on this, maybe some steps aren't converging?

    for (auto i : std::views::iota(0, num_steps)) {

        structure_interface.predict_inputs();

        structure_interface.advance_states();

        if (i > 20){
            // Put a vertical force on the mass.
            inputs.nodal_forces[mass_node_id][2] = 50.;
        }

        structure_interface.receive_inputs(nlohmann::json(inputs));
        structure_interface.set_inputs(1);

        // calling a second advance states recalculates n+1, but start again from n
        // These two steps mimic what would happen with nonlinear outer iterations.
        structure_interface.reset_states();
        structure_interface.advance_states();

        structure_interface.calculate_outputs();
        const auto output = structure_interface.outputs();
        const auto poses = output.at("nodal_poses").get<std::vector<std::array<double, 7>>>();
        ASSERT_EQ(poses.size(), num_nodes);
        EXPECT_EQ(poses[0], structure_interface.node_position(0));

        structure_interface.finalize_time_step();

        // Write out the first node's pose at each time step to a CSV.
        if (write_output) {
            x_output << i;
            for (const auto component : poses[0]) {
                x_output << "," << component;
            }
            x_output << "\n";
        }

        EXPECT_EQ(structure_interface.time_index(), i + 1);
        EXPECT_DOUBLE_EQ(structure_interface.get_time(), (i + 1) * time_step);

    }

}

}  // namespace kynema_fmb::tests
