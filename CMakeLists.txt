cmake_minimum_required(VERSION 3.16)

# Note: Using local FLINT 3.3.1 
# vcpkg is optional - only used if explicitly provided or if system libraries are missing
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "Vcpkg toolchain file")
    message(STATUS "Found vcpkg, using it for dependency management")
endif()

project(rur-cpp)

# Use C++23 standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ensure we're using the system's default GCC (version 14.2.0)
# This is typically already the case, but we can be explicit
if(NOT CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER gcc)
endif()
if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER g++)
endif()

# Add warning flags to ALL build configurations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic")

# Add debug flags for memory issues and undefined behavior
# ASAN disabled for now - will re-enable when ready to fix memory management
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined -fno-omit-frame-pointer -g")
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address,undefined -fno-omit-frame-pointer -g")
# set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address,undefined")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")

# PRIORITY: Force use of local FLINT 3.3.1 and ignore system FLINT headers
# This prevents conflicts between FLINT 3.0.1 (system) and FLINT 3.3.1 (local)

# Explicitly set local FLINT paths to bypass any interference from vcpkg/system
set(LOCAL_FLINT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/external/flint-install")
set(FLINT_INCLUDE_DIR "${LOCAL_FLINT_ROOT}/include")
set(FLINT_LIBRARY_DIR "${LOCAL_FLINT_ROOT}/lib")

# Verify local FLINT installation exists
if(NOT EXISTS "${FLINT_INCLUDE_DIR}/flint/flint.h")
    message(FATAL_ERROR "FLINT headers not found at ${FLINT_INCLUDE_DIR}/flint/flint.h. Please run ./setup-deps.sh to build and install FLINT.")
endif()

# Find the FLINT library - explicitly use local version
set(FLINT_LIBRARY "${FLINT_LIBRARY_DIR}/libflint.a")
if(NOT EXISTS "${FLINT_LIBRARY}")
    # Try shared library if static not found
    set(FLINT_LIBRARY "${FLINT_LIBRARY_DIR}/libflint.so")
    if(NOT EXISTS "${FLINT_LIBRARY}")
        message(FATAL_ERROR "FLINT 3.3.1 library not found in ${FLINT_LIBRARY_DIR}. Make sure FLINT has been built and installed.")
    endif()
endif()

message(STATUS "Using local FLINT installation: ${FLINT_INCLUDE_DIR}")
message(STATUS "FLINT library: ${FLINT_LIBRARY}")

# Add global compiler flags to prioritize local FLINT over system/vcpkg FLINT
# Note: Prepending to ensure FLINT headers take priority
set(CMAKE_CXX_FLAGS "-I${FLINT_INCLUDE_DIR} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS "-I${FLINT_INCLUDE_DIR} ${CMAKE_C_FLAGS}")

# Check FLINT version
file(READ "${FLINT_INCLUDE_DIR}/flint/flint.h" FLINT_VERSION_FILE)
string(REGEX MATCH "#define __FLINT_VERSION ([0-9]+)" FLINT_VERSION_MAJOR_MATCH "${FLINT_VERSION_FILE}")
if(FLINT_VERSION_MAJOR_MATCH)
    set(FLINT_VERSION_MAJOR ${CMAKE_MATCH_1})
endif()
string(REGEX MATCH "#define __FLINT_VERSION_MINOR ([0-9]+)" FLINT_VERSION_MINOR_MATCH "${FLINT_VERSION_FILE}")
if(FLINT_VERSION_MINOR_MATCH)
    set(FLINT_VERSION_MINOR ${CMAKE_MATCH_1})
endif()
message(STATUS "Found FLINT ${FLINT_VERSION_MAJOR}.${FLINT_VERSION_MINOR} with integrated ARB support")

# FLINT 3.x with ARB support is required
set(USE_FLINT_ARB ON CACHE BOOL "FLINT ARB support is required" FORCE)

# Find GMP, MPFR (FLINT dependencies) - use system versions for compatibility with FLINT 3.3.1
# Force system libraries to avoid version conflicts between vcpkg and local FLINT
find_library(GMP_LIBRARY NAMES gmp PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib NO_DEFAULT_PATH)
find_library(GMPXX_LIBRARY NAMES gmpxx PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib NO_DEFAULT_PATH)  
find_library(MPFR_LIBRARY NAMES mpfr PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib NO_DEFAULT_PATH)

if(NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY OR NOT MPFR_LIBRARY)
    message(STATUS "System libraries not found, falling back to vcpkg GMP/MPFR libraries")
    find_library(GMP_LIBRARY NAMES gmp PATHS ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib NO_DEFAULT_PATH)
    find_library(GMPXX_LIBRARY NAMES gmpxx PATHS ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib NO_DEFAULT_PATH)
    find_library(MPFR_LIBRARY NAMES mpfr PATHS ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib NO_DEFAULT_PATH)
endif()

message(STATUS "GMP library: ${GMP_LIBRARY}")
message(STATUS "GMPXX library: ${GMPXX_LIBRARY}")
message(STATUS "MPFR library: ${MPFR_LIBRARY}")

# ARB support (integrated in FLINT 3.x)
if(USE_FLINT_ARB)
    # In FLINT 3.x, ARB functionality is integrated into libflint
    message(STATUS "ARB support enabled - using FLINT 3.x integrated ARB for complex polynomial root finding")
    set(ARB_LIBRARIES "")  # ARB is included in libflint for FLINT 3.x
    add_compile_definitions(HAVE_FLINT_ARB)
else()
    message(STATUS "ARB support disabled - using Eigen fallback for complex polynomial root finding")
    set(ARB_LIBRARIES "")
endif()

# Create FLINT interface library with priority include paths
add_library(flint::flint INTERFACE IMPORTED)
# Fix library order: MPFR needs GMP, so GMP must come after MPFR
target_link_libraries(flint::flint INTERFACE ${FLINT_LIBRARY} ${ARB_LIBRARIES} ${MPFR_LIBRARY} ${GMPXX_LIBRARY} ${GMP_LIBRARY})
# Use BEFORE to ensure local FLINT headers have highest priority
target_include_directories(flint::flint BEFORE INTERFACE ${FLINT_INCLUDE_DIR})
# Add compile options to prioritize local FLINT over system FLINT
target_compile_options(flint::flint INTERFACE "-I${FLINT_INCLUDE_DIR}")

# Find Eigen3 using pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(EIGEN3 REQUIRED eigen3)

# Find nlohmann/json (header-only, just check include path)
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp PATHS /usr/include /usr/local/include)
if(NOT NLOHMANN_JSON_INCLUDE_DIR)
    message(FATAL_ERROR "nlohmann/json not found. Please install nlohmann-json3-dev package.")
endif()

# Find Google Test for modern testing framework
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    message(STATUS "Google Test not found via find_package, trying vcpkg...")
    # Try to use vcpkg's gtest
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib/libgtest.a")
        set(GTEST_LIBRARIES 
            "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib/libgtest.a"
            "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/lib/libgtest_main.a"
        )
        set(GTEST_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-linux/include")
        set(GTEST_FOUND TRUE)
        message(STATUS "Found Google Test via vcpkg: ${GTEST_LIBRARIES}")
    else()
        message(STATUS "Google Test not found. Install via: vcpkg install gtest")
        set(GTEST_FOUND FALSE)
    endif()
else()
    message(STATUS "Found Google Test via system installation")
endif()

# Add source files
set(SOURCES
    src/monomial.cpp
    src/polynomial.cpp
    src/axf4_wrapper.c
    src/f4_solver.cpp
    src/axf4_lib.c
    src/quotient_ring.cpp
    src/bivariate_lex.cpp
    src/rur_solver.cpp
    src/direct_rur_solver.cpp
    src/flint_reducer.cpp
    src/flint_linear_solver.cpp
    src/polynomial_basis.cpp
    src/flint_reducer_v2.cpp
    src/flint_mpoly_wrappers.cpp
    src/flint_linear_algebra.cpp
    src/flint_rur_solver.cpp
    src/julia_rur/data_structures.cpp
    src/julia_rur/multiplication_tables.cpp
    src/julia_rur/bivariate_algorithm.cpp
    src/julia_rur/multi_modular.cpp
    src/julia_rur/quotient_basis.cpp
    src/julia_rur/f4_monomial_decoder.cpp
    src/julia_rur/f4_integration.cpp
    src/julia_rur/univariate_parameterization.cpp
    src/julia_rur/separating_element_search.cpp
    src/julia_rur/separating_element_systematic.cpp
    src/julia_rur/rur_main_algorithm.cpp
    src/julia_rur/numerical_roots_eigen.cpp
    src/julia_rur/polynomial_solver.cpp
    src/julia_rur/numerical_roots_flint.cpp
    src/julia_rur/polynomial_solver_enhanced.cpp
    src/julia_rur/f4_polynomial_formatter.cpp
    src/julia_rur/prime_utils.cpp
    src/julia_rur/hyperplane_sections.cpp
    src/julia_rur/polynomial_evaluator.cpp
    src/julia_rur/polynomial_expander.cpp
)

# Create library
add_library(rur STATIC ${SOURCES})

# CRITICAL: Use BEFORE to ensure local FLINT headers take absolute priority
# This prevents any conflicts with system FLINT 3.0.1 headers
target_include_directories(rur BEFORE PUBLIC 
    ${FLINT_INCLUDE_DIR}  # Local FLINT 3.3.1 - MUST be first
)

target_include_directories(rur PUBLIC 
    src/
    ${EIGEN3_INCLUDE_DIRS}
    ${NLOHMANN_JSON_INCLUDE_DIR}
)

# Add compiler-specific flags to completely ignore system FLINT
target_compile_options(rur PRIVATE "-I${FLINT_INCLUDE_DIR}")
target_compile_definitions(rur PUBLIC "FLINT_LOCAL_INSTALL=1")

target_link_libraries(rur PUBLIC 
    flint::flint
    ${EIGEN3_LIBRARIES}
    m  # Math library for axf4.c
)

# Test executables
add_executable(test_basic test/test_basic.cpp)
target_link_libraries(test_basic rur)

add_executable(test_f4_wrapper test/test_f4_wrapper.cpp)
target_link_libraries(test_f4_wrapper rur)

add_executable(test_f4_simple test/test_f4_simple.cpp)
target_link_libraries(test_f4_simple rur)

add_executable(test_polynomial_formatter_fractions test/test_polynomial_formatter_fractions.cpp)
target_link_libraries(test_polynomial_formatter_fractions rur)

add_executable(test_circle_line_fix test/test_circle_line_fix.cpp)
target_link_libraries(test_circle_line_fix rur)

add_executable(test_parabola_line_debug test/test_parabola_line_debug.cpp)
target_link_libraries(test_parabola_line_debug rur)

add_executable(test_biv_lex_bug test/test_biv_lex_bug.cpp)
target_link_libraries(test_biv_lex_bug rur)

add_executable(test_minpoly_comparison test/test_minpoly_comparison.cpp)
target_link_libraries(test_minpoly_comparison rur)

add_executable(test_rur_solver test/test_rur_solver.cpp)
target_link_libraries(test_rur_solver rur)

add_executable(test_debug_gb test/test_debug_gb.cpp)
target_link_libraries(test_debug_gb rur)

add_executable(test_raw_axf4 test/test_raw_axf4.cpp)
target_link_libraries(test_raw_axf4 rur)

add_executable(test_json_input test/test_json_input.cpp)
target_link_libraries(test_json_input rur)

add_executable(test_direct_rur test/test_direct_rur.cpp)
target_link_libraries(test_direct_rur rur)

add_executable(test_flint_linear_algebra test/test_flint_linear_algebra.cpp)
target_link_libraries(test_flint_linear_algebra rur)

add_executable(test_flint_rur_solver test/test_flint_rur_solver.cpp)
target_link_libraries(test_flint_rur_solver rur)

add_executable(test_julia_data_structures test/test_julia_data_structures.cpp)
target_link_libraries(test_julia_data_structures rur)

add_executable(test_multiplication_tables test/test_multiplication_tables.cpp)
target_link_libraries(test_multiplication_tables rur)

add_executable(test_bivariate_algorithm test/test_bivariate_algorithm.cpp)
target_link_libraries(test_bivariate_algorithm rur)

add_executable(test_multi_modular test/test_multi_modular.cpp)
target_link_libraries(test_multi_modular rur)

add_executable(test_f4_structured_api test/test_f4_structured_api.cpp)
target_link_libraries(test_f4_structured_api rur)

add_executable(test_f4_simple_structured test/test_f4_simple_structured.cpp)
target_link_libraries(test_f4_simple_structured rur)

add_executable(test_f4_debug_structured test/test_f4_debug_structured.cpp)
target_link_libraries(test_f4_debug_structured rur)

add_executable(test_f4_format_debug test/test_f4_format_debug.cpp)
target_link_libraries(test_f4_format_debug rur)

add_executable(test_quotient_basis test/test_quotient_basis.cpp)
target_link_libraries(test_quotient_basis rur)

add_executable(test_f4_quotient_basis_integration test/test_f4_quotient_basis_integration.cpp)
target_link_libraries(test_f4_quotient_basis_integration rur)

add_executable(test_f4_multiplication_tables test/test_f4_multiplication_tables.cpp)
target_link_libraries(test_f4_multiplication_tables rur)

add_executable(test_debug_f4_output test/test_debug_f4_output.cpp)
target_link_libraries(test_debug_f4_output rur)

add_executable(demo_rur demo/demo_rur.cpp)
target_link_libraries(demo_rur rur)

add_executable(test_rational_reconstruction test/test_rational_reconstruction.cpp)
target_link_libraries(test_rational_reconstruction rur)

add_executable(test_simple_rational_reconstruction test/test_simple_rational_reconstruction.cpp)
target_link_libraries(test_simple_rational_reconstruction rur)

add_executable(test_simple_rur test/test_simple_rur.cpp)
target_link_libraries(test_simple_rur rur)

add_executable(test_rational_no_verify test/test_rational_no_verify.cpp)
target_link_libraries(test_rational_no_verify rur)

add_executable(test_numerical_roots_eigen test/test_numerical_roots_eigen.cpp)
target_link_libraries(test_numerical_roots_eigen rur)

add_executable(test_polynomial_solver test/test_polynomial_solver.cpp)
target_link_libraries(test_polynomial_solver rur)

add_executable(test_multivariate_debug test/test_multivariate_debug.cpp)
target_link_libraries(test_multivariate_debug rur)

add_executable(test_multivariate_rational test/test_multivariate_rational.cpp)
target_link_libraries(test_multivariate_rational rur)

add_executable(test_f4_multivariate test/test_f4_multivariate.cpp)
target_link_libraries(test_f4_multivariate rur)

add_executable(test_multivariate_primes test/test_multivariate_primes.cpp)
target_link_libraries(test_multivariate_primes rur)

add_executable(test_simple_multivariate test/test_simple_multivariate.cpp)
target_link_libraries(test_simple_multivariate rur)

add_executable(test_linear_debug test/test_linear_debug.cpp)
target_link_libraries(test_linear_debug rur)

add_executable(test_linear_rur test/test_linear_rur.cpp)
target_link_libraries(test_linear_rur rur)

add_executable(test_linear_debug2 test/test_linear_debug2.cpp)
target_link_libraries(test_linear_debug2 rur)

add_executable(test_linear_simple test/test_linear_simple.cpp)
target_link_libraries(test_linear_simple rur)

add_executable(test_linear_verbose test/test_linear_verbose.cpp)
target_link_libraries(test_linear_verbose rur)

add_executable(test_f4_crash_debug test/test_f4_crash_debug.cpp)
target_link_libraries(test_f4_crash_debug rur)

add_executable(test_f4_single_thread test/test_f4_single_thread.cpp)
target_link_libraries(test_f4_single_thread rur)

add_executable(test_x2_minus_2_debug test/test_x2_minus_2_debug.cpp)
target_link_libraries(test_x2_minus_2_debug rur)

add_executable(test_f4_parsing test/test_f4_parsing.cpp)
target_link_libraries(test_f4_parsing rur)

add_executable(test_x2_minus_2_formatter test/test_x2_minus_2_formatter.cpp)
target_link_libraries(test_x2_minus_2_formatter rur)

add_executable(test_f4_negative_parsing test/test_f4_negative_parsing.cpp)
target_link_libraries(test_f4_negative_parsing rur)

add_executable(test_minimal_poly_debug test/test_minimal_poly_debug.cpp)
target_link_libraries(test_minimal_poly_debug rur)

add_executable(test_crt_reconstruction_debug test/test_crt_reconstruction_debug.cpp)
target_link_libraries(test_crt_reconstruction_debug rur)

add_executable(test_reconstruction_issue test/test_reconstruction_issue.cpp)
target_link_libraries(test_reconstruction_issue rur)

add_executable(test_x2_small_prime test/test_x2_small_prime.cpp)
target_link_libraries(test_x2_small_prime rur)

add_executable(test_debug_rational_coeff test/test_debug_rational_coeff.cpp)
target_link_libraries(test_debug_rational_coeff rur)

add_executable(test_polynomial_printing test/test_polynomial_printing.cpp)
target_link_libraries(test_polynomial_printing rur)

add_executable(test_poly_solver_debug test/test_poly_solver_debug.cpp)
target_link_libraries(test_poly_solver_debug rur)

add_executable(test_complete_x2_trace test/test_complete_x2_trace.cpp)
target_link_libraries(test_complete_x2_trace rur)

add_executable(test_eigen_roots_debug test/test_eigen_roots_debug.cpp)
target_link_libraries(test_eigen_roots_debug rur)

add_executable(test_solver_roots_debug test/test_solver_roots_debug.cpp)
target_link_libraries(test_solver_roots_debug rur)

add_executable(test_cubic_debug test/test_cubic_debug.cpp)
target_link_libraries(test_cubic_debug rur)

add_executable(test_cubic_rur_debug test/test_cubic_rur_debug.cpp)
target_link_libraries(test_cubic_rur_debug rur)

add_executable(test_multivariate_simple_new test/test_multivariate_simple_new.cpp)
target_link_libraries(test_multivariate_simple_new rur)

add_executable(test_simple_3var test/test_simple_3var.cpp)
target_link_libraries(test_simple_3var rur)

add_executable(test_linear_3x3 test/test_linear_3x3.cpp)
target_link_libraries(test_linear_3x3 rur)



# add_executable(test_f4_gb_consistency test/test_f4_gb_consistency.cpp)
# target_link_libraries(test_f4_gb_consistency rur)

add_executable(test_f4_gb_consistency_simple test/test_f4_gb_consistency_simple.cpp)
target_link_libraries(test_f4_gb_consistency_simple rur)

add_executable(test_multivariate_complex test/test_multivariate_complex.cpp)
target_link_libraries(test_multivariate_complex rur)

add_executable(test_simple_quadratic test/test_simple_quadratic.cpp)
target_link_libraries(test_simple_quadratic rur)

add_executable(test_separating_debug test/test_separating_debug.cpp)
target_link_libraries(test_separating_debug rur)

add_executable(test_circle_parabola_debug test/test_circle_parabola_debug.cpp)
target_link_libraries(test_circle_parabola_debug rur)

add_executable(test_variable_indexing test/test_variable_indexing.cpp)
target_link_libraries(test_variable_indexing rur)

add_executable(test_indexing_edge_cases test/test_indexing_edge_cases.cpp)
target_link_libraries(test_indexing_edge_cases rur)

add_executable(test_tensor_product_debug test/test_tensor_product_debug.cpp)
target_link_libraries(test_tensor_product_debug rur)

add_executable(test_polynomial_solver_enhanced test/test_polynomial_solver_enhanced.cpp)
target_link_libraries(test_polynomial_solver_enhanced rur)

add_executable(test_hyperplane_sections test/test_hyperplane_sections.cpp)
target_link_libraries(test_hyperplane_sections rur)

add_executable(test_flint_arb_support test/test_flint_arb_support.cpp)
target_link_libraries(test_flint_arb_support rur)

add_executable(test_comprehensive test/test_comprehensive.cpp)
target_link_libraries(test_comprehensive rur)

add_executable(test_suite test/test_suite.cpp)
target_link_libraries(test_suite rur)

add_executable(run_tests test/run_tests.cpp)
target_link_libraries(run_tests rur)

add_executable(test_polynomial_evaluator test/test_polynomial_evaluator.cpp)
target_link_libraries(test_polynomial_evaluator rur)

add_executable(test_multiplicities test/test_multiplicities.cpp)
target_link_libraries(test_multiplicities rur)

add_executable(test_simple_multiplicities test/test_simple_multiplicities.cpp)
target_link_libraries(test_simple_multiplicities rur)

add_executable(test_benchmark_systems test/test_benchmark_systems.cpp)
target_link_libraries(test_benchmark_systems rur)

add_executable(test_cas_examples test/test_cas_examples.cpp)
target_link_libraries(test_cas_examples rur)

add_executable(test_cyclic4_zerodim test/test_cyclic4_zerodim.cpp)
target_link_libraries(test_cyclic4_zerodim rur)

# Modern unified test suite using Google Test
if(GTEST_FOUND)
    # Create unified test executable
    add_executable(rur_tests
        test/gtest/test_main.cpp
        test/gtest/test_basic_operations.cpp
        test/gtest/test_univariate_systems.cpp
        test/gtest/test_linear_systems.cpp
        test/gtest/test_nonlinear_systems.cpp
        test/gtest/test_numerical_roots.cpp
        test/gtest/test_integration.cpp
        test/gtest/test_benchmark_systems.cpp
        test/gtest/test_positive_dimensional.cpp
        test/gtest/test_underdetermined_systems.cpp
    )
    
    target_link_libraries(rur_tests 
        rur 
        ${GTEST_LIBRARIES}
        pthread  # Required by gtest on Linux
    )
    
    target_include_directories(rur_tests PRIVATE 
        ${GTEST_INCLUDE_DIRS}
        test/gtest/include
    )
    
    # Enable CTest integration
    enable_testing()
    include(GoogleTest)
    gtest_discover_tests(rur_tests)
    
    message(STATUS "Modern test suite 'rur_tests' configured successfully")
else()
    message(WARNING "Google Test not available - unified test suite disabled")
    message(STATUS "To enable: run 'vcpkg install gtest' then reconfigure")
endif()

# F4 wrapper (will be added later)
# add_subdirectory(axf4)
