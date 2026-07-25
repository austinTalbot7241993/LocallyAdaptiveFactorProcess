cmake_minimum_required(VERSION 3.16)

project(lafp
    VERSION 1.0.0
    DESCRIPTION "Locally Adaptive Factor Process MCMC sampler"
    LANGUAGES C
)

# ── Build type default ─────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# ── Dependencies ───────────────────────────────────────────────────────────────
find_package(GSL REQUIRED)       # provides GSL::gsl and GSL::gslcblas
find_package(Threads REQUIRED)   # provides Threads::Threads (pthreads on POSIX)

# ── Library sources ────────────────────────────────────────────────────────────
# ── SSsimulate2 variant (serial vs pthreads) ──────────────────────────────────
option(LAFP_USE_PTHREADS
    "Use the pthreads-parallel SSsimulate2 implementation instead of the serial one"
    OFF)

if(LAFP_USE_PTHREADS)
    set(SSSIM_SRC src/SSsimulate2_pthreads.c)
else()
    set(SSSIM_SRC src/SSsimulate2.c)
endif()

# NOTE: mdarray_complex.c, KalmanFilter2c.c, and fastStateSmoother2c.c are
# excluded from the build. They originate from the experimental Code/Complex/
# subdirectory, contain pre-existing declaration/definition mismatches, and
# are not referenced by any MCMC code. They require dedicated fixes before
# being added back here.
set(LAFP_SOURCES
    src/gmlib.c
    src/mdarray.c
    src/matrix_utils.c
    # src/mvnorm.c is an exact duplicate of distributions.c (same functions,
    # same signatures). Include only one to avoid duplicate symbol errors.
    src/distributions.c
    src/KalmanFilter2.c
    src/fastStateSmoother2.c
    ${SSSIM_SRC}
    src/NGPtools2.c
    src/NGPsample2.c
    src/NGPmcmc.c
    src/NGPlaf2.c
)

# ── Shared library ─────────────────────────────────────────────────────────────
add_library(lafp_shared SHARED ${LAFP_SOURCES})
set_target_properties(lafp_shared PROPERTIES OUTPUT_NAME "lafp")

target_include_directories(lafp_shared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(lafp_shared
    PUBLIC GSL::gsl GSL::gslcblas Threads::Threads m
)

target_compile_options(lafp_shared PRIVATE -O3 -Wall -Wextra)

# ── Static library ─────────────────────────────────────────────────────────────
add_library(lafp_static STATIC ${LAFP_SOURCES})
set_target_properties(lafp_static PROPERTIES OUTPUT_NAME "lafp")

target_include_directories(lafp_static
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(lafp_static
    PUBLIC GSL::gsl GSL::gslcblas Threads::Threads m
)

target_compile_options(lafp_static PRIVATE -O3 -Wall -Wextra)

# ── Subdirectories ─────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
add_subdirectory(examples)

# ── Install rules ──────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS lafp_shared lafp_static
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/lafp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# ── CMake package config (for find_package(lafp)) ─────────────────────────────
include(CMakePackageConfigHelpers)

install(TARGETS lafp_shared lafp_static
    EXPORT  lafpTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(EXPORT lafpTargets
    FILE      lafpTargets.cmake
    NAMESPACE lafp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lafp
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lafpConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/lafpConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lafp
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/lafpConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/lafpConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/lafpConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lafp
)
