target_sources(kynema_fmb_library
  PRIVATE
  generalized_structure_interface.cpp
)

configure_file("${KYNEMA_FMB_DRIVER_MODULE_DIR}/Module.h"
  "${CMAKE_BINARY_DIR}/driver-interface/driver/Module.h" COPYONLY)
install(FILES "${CMAKE_BINARY_DIR}/driver-interface/driver/Module.h"
  DESTINATION include/KynemaFMB/driver)
target_include_directories(kynema_fmb_library SYSTEM PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/driver-interface>)

install(FILES
        generalized_structure_interface.hpp
        DESTINATION include/KynemaFMB/interfaces/generalized/
)
