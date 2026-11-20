target_sources(kynema_fmb_library
  PRIVATE
  generalized_structure_interface.cpp
)

install(FILES
        coupling_data.hpp
        generalized_structure_interface.hpp
        DESTINATION include/KynemaFMB/interfaces/generalized/
)
