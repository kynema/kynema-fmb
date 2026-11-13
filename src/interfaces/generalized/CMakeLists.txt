target_sources(kynema_fmb_library
  PRIVATE
  generalized_structure_interface.cpp
)

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/driver-interface/driver")
foreach(header Module.h Timers.h ParallelPrinter.h)
  configure_file("${KYNEMA_FMB_DRIVER_MODULE_DIR}/${header}"
    "${CMAKE_BINARY_DIR}/driver-interface/driver/${header}" COPYONLY)
  install(FILES "${CMAKE_BINARY_DIR}/driver-interface/driver/${header}"
    DESTINATION include/KynemaFMB/driver)
endforeach()
target_include_directories(kynema_fmb_library SYSTEM PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/driver-interface>)
target_link_libraries(kynema_fmb_library PUBLIC MPI::MPI_CXX)

install(FILES
        generalized_structure_interface.hpp
        DESTINATION include/KynemaFMB/interfaces/generalized/
)
