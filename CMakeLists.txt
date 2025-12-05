cmake_minimum_required(VERSION 3.28)
project(OTPQ-MATH LANGUAGES CXX)

# ============================================================
# C++23 Standard
# ============================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================
# project_options (warnings, features, etc.)
# ============================================================
add_library(project_options INTERFACE)

target_compile_features(project_options INTERFACE cxx_std_23)

target_compile_options(project_options INTERFACE
        -Wall -Wextra -Wpedantic -Wconversion
)

# ============================================================
# Choose underlying modular integer type (default: int)
# ============================================================
set(OTPQ_ALS_UNDERLYING_TYPE "int" CACHE STRING
        "Underlying variable type for modular arithmetic")

message(STATUS "OTPQ modular integer type = ${OTPQ_ALS_UNDERLYING_TYPE}")



option(OTPQ_SEC_MEMORY_BOUND "Enable secure memory bound checking" OFF)

message(STATUS "OTPQ secure memory bound = ${OTPQ_SEC_MEMORY_BOUND}")

if(OTPQ_SEC_MEMORY_BOUND)
    add_compile_definitions(OTPQ_SEC_MEMORY_BOUND)
endif()


# ============================================================
# otpqmath (header-only library)
# ============================================================
add_library(libotpqmath INTERFACE)

target_include_directories(libotpqmath INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
)

# Inject macro: OTPQ_ALS_UNDERLYING_TYPE=int or whatever user chooses
target_compile_definitions(libotpqmath INTERFACE
        OTPQ_ALS_UNDERLYING_TYPE=${OTPQ_ALS_UNDERLYING_TYPE}
)

target_link_libraries(libotpqmath INTERFACE project_options)

# ============================================================
# Test executable
# ============================================================

# ============================================================
# Examples
# ============================================================
add_subdirectory(examples)
