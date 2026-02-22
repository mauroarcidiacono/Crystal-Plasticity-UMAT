# Create executable for testing
add_executable(test_ludcmp test_ludcmp.f90)
add_executable(test_lubksb test_lubksb.f90)
add_executable(test_cross test_cross.f90)

# Link against umatlib
target_link_libraries(test_ludcmp umatlib)
target_link_libraries(test_lubksb umatlib)
target_link_libraries(test_cross umatlib)

# If using Intel Fortran, link the runtime statically so tests are portable
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(test_ludcmp PRIVATE /libs:static /threads)
    target_compile_options(test_lubksb PRIVATE /libs:static /threads)
    target_compile_options(test_cross PRIVATE /libs:static /threads)
endif()

# Test with CTest
add_test(NAME LUDCMP_Test COMMAND test_ludcmp)
add_test(NAME LUBKSB_Test COMMAND test_lubksb)
add_test(NAME CROSS_Test COMMAND test_cross)
