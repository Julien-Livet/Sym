cmake_minimum_required(VERSION 3.10)
project(Sym)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

find_package(Boost REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(include)

enable_testing()
add_executable(test_sym tests/test_sym.cpp)
target_link_libraries(test_sym
#  asan ubsan
  gtest gtest_main)

math(EXPR TIMEOUT "3 * 60 * 60")
add_test(NAME TestSym COMMAND test_sym)
set_tests_properties(TestSym PROPERTIES TIMEOUT ${TIMEOUT})

include(GoogleTest)
gtest_discover_tests(test_sym)

install(DIRECTORY include/Sym DESTINATION include)

add_library(Sym INTERFACE)
target_include_directories(Sym INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SymConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SymConfig.cmake
    INSTALL_DESTINATION lib/cmake/Sym
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SymConfigVersion.cmake
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SymConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SymConfigVersion.cmake
    DESTINATION lib/cmake/Sym
)

install(TARGETS Sym EXPORT SymTargets)

install(EXPORT SymTargets
    NAMESPACE Sym::
    DESTINATION lib/cmake/Sym
)
