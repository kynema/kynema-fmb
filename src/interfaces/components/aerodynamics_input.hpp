#pragma once

#include <vector>

namespace kynema_fmb::interfaces::components {

struct AerodynamicSection {
    size_t id{};
    double s{};
    double chord{};
    double section_offset_x{};
    double section_offset_y{};
    double aerodynamic_center{};
    double twist{};
    std::vector<double> aoa;
    std::vector<double> cl;
    std::vector<double> cd;
    std::vector<double> cm;

    AerodynamicSection() = default;

    AerodynamicSection(
        size_t id_, double s_, double chord_, double section_offset_x_, double section_offset_y_,
        double aerodynamic_center_, double twist_, const std::vector<double>& aoa_,
        const std::vector<double>& cl_, const std::vector<double>& cd_,
        const std::vector<double>& cm_
    )
        : id(id_),
          s(s_),
          chord(chord_),
          section_offset_x(section_offset_x_),
          section_offset_y(section_offset_y_),
          aerodynamic_center(aerodynamic_center_),
          twist(twist_),
          aoa(aoa_),
          cl(cl_),
          cd(cd_),
          cm(cm_) {}
};

struct AerodynamicBodyInput {
    size_t id{};
    std::vector<size_t> beam_node_ids;
    std::vector<AerodynamicSection> aero_sections;
    AerodynamicBodyInput(
        size_t id_val, const std::vector<size_t>& node_ids,
        const std::vector<AerodynamicSection>& sections
    )
        : id(id_val), beam_node_ids(node_ids), aero_sections(sections) {}
};

class AerodynamicsInput {
public:
    bool is_enabled = false;
    std::vector<std::vector<AerodynamicSection>> aero_inputs;
    std::vector<size_t> airfoil_map;
};
}  // namespace kynema_fmb::interfaces::components
