cmake_minimum_required(VERSION 3.16)

project(euler1d-fvm LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(euler1d_core STATIC
    src/eos.cpp
    src/flux.cpp
    src/reconstruction.cpp
    src/mesh.cpp
    src/problems.cpp
    src/boundary.cpp
    src/time_integration.cpp
    src/solver.cpp
    src/io.cpp
)

target_include_directories(euler1d_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_executable(euler1d src/main.cpp)
target_link_libraries(euler1d PRIVATE euler1d_core)

foreach(target euler1d_core euler1d)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
    )
endforeach()

option(EULER1D_BUILD_TESTS "Build unit tests" ON)

if(EULER1D_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    enable_testing()
    include(GoogleTest)

    add_executable(test_eos tests/test_eos.cpp)
    target_link_libraries(test_eos PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_eos)

    add_executable(test_flux tests/test_flux.cpp)
    target_link_libraries(test_flux PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_flux)

    add_executable(test_mesh_boundary tests/test_mesh_boundary.cpp)
    target_link_libraries(test_mesh_boundary PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_mesh_boundary)

    add_executable(test_reconstruction tests/test_reconstruction.cpp)
    target_link_libraries(test_reconstruction PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_reconstruction)

    add_executable(test_solver tests/test_solver.cpp)
    target_link_libraries(test_solver PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_solver)

    add_executable(test_problems tests/test_problems.cpp)
    target_link_libraries(test_problems PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_problems)

    add_executable(test_io tests/test_io.cpp)
    target_link_libraries(test_io PRIVATE euler1d_core GTest::gtest_main)
    gtest_discover_tests(test_io)

    # End-to-end CLI runs on the shipped case files.
    foreach(case density_wave sod lax)
        add_test(NAME cli_${case}
                 COMMAND euler1d ${CMAKE_CURRENT_SOURCE_DIR}/cases/${case}.cfg
                         --n-cells 50 --quiet
                         --output-dir ${CMAKE_CURRENT_BINARY_DIR}/cli_out/${case})
    endforeach()
    add_test(NAME cli_scheme_overrides
             COMMAND euler1d ${CMAKE_CURRENT_SOURCE_DIR}/cases/density_wave.cfg
                     --n-cells 50 --quiet --flux rusanov --reconstruction MUSCL
                     --limiter van_leer
                     --output-dir ${CMAKE_CURRENT_BINARY_DIR}/cli_out/overrides)
    add_test(NAME cli_missing_config
             COMMAND euler1d ${CMAKE_CURRENT_BINARY_DIR}/does_not_exist.cfg --quiet)
    set_tests_properties(cli_missing_config PROPERTIES WILL_FAIL TRUE)
endif()
