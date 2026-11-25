# Specify the source files for the regression test executable
target_sources(
    kynema_fmb_regression_tests
    PRIVATE
    test_6dof_interface.cpp
    test_turbine_interface.cpp
)

add_custom_target(
    copy_generalized_test_files ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/generalized_test_files
            $<TARGET_FILE_DIR:kynema_fmb_regression_tests>/generalized_test_files)