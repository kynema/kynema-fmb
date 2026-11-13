set(CTEST_CUSTOM_WARNING_EXCEPTION ".*: warning: could not create compact unwind*"
                                   ".*has no symbols*"
                                   ".*submake: resetting jobserver mode*"
                                   ".*[Ww]arning #20014-D.*"  # calling a __host__ function from a __host__ __device__ function
                                   ".*[Ww]arning #20015-D.*"  # related host/device call diagnostic
)

set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "131072")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "131072")
