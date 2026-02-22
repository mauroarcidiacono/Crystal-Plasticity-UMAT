cmake_minimum_required(VERSION 3.10)
project(CrystalPlasticityUMAT Fortran)

enable_testing()

# Only compile specific routines for testing
set(SRC_FILES
    ludcmp.for
    lubksb.for
    utils.for
    ${CMAKE_SOURCE_DIR}/tests/xit_stub.f90
    # Add more!
)

add_library(umatlib ${SRC_FILES})

# Add the tests folder
add_subdirectory(tests)