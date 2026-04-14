add_executable(PionCloudTests test_all.cc)
target_link_libraries(PionCloudTests PRIVATE PionCloudLib)
add_test(NAME PionCloudTests COMMAND PionCloudTests)
