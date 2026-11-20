#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <nlohmann/json.hpp>

namespace kynema_fmb::interfaces {

inline constexpr int structural_schema_version = 1;

struct StructuralInputs {
    std::vector<std::array<double, 6>> nodal_forces;
};

inline void to_json(nlohmann::json& json, const StructuralInputs& inputs) {
    json = {{"version", structural_schema_version}, {"nodal_forces", inputs.nodal_forces}};
}

inline void from_json(const nlohmann::json& json, StructuralInputs& inputs) {
    const auto& version = json.at("version");
    if (!version.is_number_integer() || version != structural_schema_version) {
        throw std::invalid_argument("Unsupported FMB JSON schema version");
    }
    const auto& forces = json.at("nodal_forces");
    if (!forces.is_array()) {
        throw std::invalid_argument("FMB nodal_forces must be an array");
    }
    for (const auto& force : forces) {
        if (!force.is_array() || force.size() != 6) {
            throw std::invalid_argument("Each FMB nodal force must have six components");
        }
        for (const auto& component : force) {
            if (!component.is_number() || !std::isfinite(component.get<double>())) {
                throw std::invalid_argument("FMB nodal forces must be finite numbers");
            }
        }
    }
    forces.get_to(inputs.nodal_forces);
}

struct StructuralOutputs {
    std::vector<std::array<double, 7>> nodal_poses;
};

inline void to_json(nlohmann::json& json, const StructuralOutputs& outputs) {
    json = {{"version", structural_schema_version}, {"nodal_poses", outputs.nodal_poses}};
}

}  // namespace kynema_fmb::interfaces
