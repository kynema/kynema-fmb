#include <stdexcept>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "interfaces/generalized/generalized_structure_interface.hpp"

namespace kynema_fmb::tests {

TEST(GeneralizedInterface, RejectsUnsupportedTurbine) {
    const YAML::Node input_yaml =
        YAML::LoadFile("generalized_test_files/general_turbine_inputs.yaml");

    EXPECT_THROW((void)interfaces::GeneralStructureInterface(input_yaml), std::invalid_argument);
}

}  // namespace kynema_fmb::tests
