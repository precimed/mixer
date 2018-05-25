cmake_minimum_required(VERSION 2.8)

project(bgmg)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE "Release")
message("-- Build type: ${CMAKE_BUILD_TYPE}")

set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0" "1.56" "1.56.0" "1.55" "1.55.0" "1.54" "1.54.0" "1.53" "1.53.0" "1.52" "1.52.0" "1.51" "1.51.0" "1.50" "1.50.0" "1.49" "1.49.0" "1.48" "1.48.0" "1.47" "1.47.0" "1.46" "1.46.0" "1.45" "1.45.0" "1.44" "1.44.0" "1.42" "1.42.0" "1.41.0" "1.41" "1.40.0" "1.40")

# find boost
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find boost libraries.")
endif (NOT Boost_FOUND)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_LIBS   ON)
set(BGMG_BOOST_COMPONENTS program_options filesystem system)
find_package(Boost COMPONENTS REQUIRED ${L2RTOOL_BOOST_COMPONENTS})
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find required boost libraries.")
  return()
else (NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif (NOT Boost_FOUND)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else (COMPILER_SUPPORTS_CXX11)
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif (COMPILER_SUPPORTS_CXX11)

#SET(BUILD_GMOCK OFF CACHE BOOL "Builds the googlemock subproject")
#SET(BUILD_GTEST ON CACHE BOOL "Builds the googletest subproject")
#add_subdirectory( googletest )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(SRC_LIST
	source.cc
	bgmg_calculator.cc
	bgmg.h
	bgmg_calculator.h
)

add_library(bgmg SHARED ${SRC_LIST})
# set_target_properties(bgmg PROPERTIES LINK_FLAGS "-static" )
target_link_libraries(bgmg ${Boost_LIBRARIES})

set(SRC_TEST_LIST
	test/bgmg_main_unittest.cc
	bgmg_calculator.cc
	${PROJECT_SOURCE_DIR}/googletest/googletest/src/gtest_main.cc
	${PROJECT_SOURCE_DIR}/googletest/googletest/src/gtest-all.cc
)

add_executable(bgmg-test ${SRC_TEST_LIST})

target_include_directories(bgmg PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(bgmg PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest/include)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR})

include(GenerateExportHeader)
generate_export_header(bgmg
		EXPORT_MACRO_NAME DLL_PUBLIC)

install( TARGETS bgmg DESTINATION bin )
