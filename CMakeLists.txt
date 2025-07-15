cmake_minimum_required(VERSION 3.10)
project(Sym)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 23)
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

find_package(PkgConfig REQUIRED)
pkg_check_modules(GINAC REQUIRED ginac)

include_directories(include)

add_compile_definitions(EIGEN_NO_DEBUG)

enable_testing()
add_executable(test_sym tests/test_sym.cpp)
target_link_libraries(test_sym
#  asan ubsan
  gtest gtest_main ${GINAC_LIBRARIES})

math(EXPR TIMEOUT "3 * 60 * 60")
add_test(NAME TestSym COMMAND test_sym)
set_tests_properties(TestSym PROPERTIES TIMEOUT ${TIMEOUT})

include(GoogleTest)
gtest_discover_tests(test_sym)
