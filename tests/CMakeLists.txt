add_executable(PionCloudTests test_all.cc)
target_link_libraries(PionCloudTests PRIVATE PionCloudLib)
add_test(NAME PionCloudTests COMMAND PionCloudTests)

# v1.4 -- why these properties:
#  * WORKING_DIRECTORY: the CTEQ6 Fortran reads cteq6m.tbl from the cwd; the
#    tables are copied to the top build directory, not to build/tests.
#  * PASS_REGULAR_EXPRESSION: in v1.3 a missing table made the Fortran STOP
#    with exit status 0 after ~80 of the tests, and ctest reported "Passed".
#    Now the test passes only if the suite reaches its summary line with
#    zero failures; any early exit or crash is a failure.
set_tests_properties(PionCloudTests PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    PASS_REGULAR_EXPRESSION "Results: [0-9]+ passed, 0 failed"
    FAIL_REGULAR_EXPRESSION "\\[FAIL\\]")

# v1.4: regression against the PUBLISHED TDIS curves (docs/validation/ref_data).
# Fails if PionCloud stops reproducing PR12-15-006 Figs. 5, 6, 37, 46, 47
# within the tolerances documented in tests/validation.cc.
add_executable(PionCloudValidation validation.cc)
target_link_libraries(PionCloudValidation PRIVATE PionCloudLib)
target_compile_definitions(PionCloudValidation PRIVATE
    PIONCLOUD_REF_DIR="${PROJECT_SOURCE_DIR}/docs/validation/ref_data")
add_test(NAME PionCloudValidation COMMAND PionCloudValidation)
set_tests_properties(PionCloudValidation PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    PASS_REGULAR_EXPRESSION "Results: [0-9]+ passed, 0 failed"
    FAIL_REGULAR_EXPRESSION "\\[FAIL\\]")
