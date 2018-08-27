cmake_minimum_required(VERSION 2.8)

project(bgmg)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE "Release")
message("-- Build type: ${CMAKE_BUILD_TYPE}")

set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0" "1.56" "1.56.0" "1.55" "1.55.0" "1.54" "1.54.0" "1.53" "1.53.0" "1.52" "1.52.0" "1.51" "1.51.0" "1.50" "1.50.0" "1.49" "1.49.0" "1.48" "1.48.0" "1.47" "1.47.0" "1.46" "1.46.0" "1.45" "1.45.0" "1.44" "1.44.0" "1.42" "1.42.0" "1.41.0" "1.41" "1.40.0" "1.40")
#set(BGMG_BOOST_COMPONENTS program_options filesystem system date_time iostreams zlib)
set(BGMG_BOOST_COMPONENTS program_options filesystem system date_time)

# find boost (static runtime version of BOOST for BGMG command-line interface)
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_STATIC_RUNTIME_LIBS   ON)
find_package(Boost COMPONENTS REQUIRED ${BGMG_BOOST_COMPONENTS})
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find required boost libraries.")
  return()
endif (NOT Boost_FOUND)
set(Boost_LIBRARIES_STATIC_RUNTIME_LIBS ${Boost_LIBRARIES})
message("-- Boost_LIBRARIES_STATIC_RUNTIME_LIBS: ${Boost_LIBRARIES_STATIC_RUNTIME_LIBS}")

# find boost (static version of BOOST for libbgmg.so / libbgmg.dynlib / bgmglib.dll )
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find boost libraries.")
endif (NOT Boost_FOUND)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_STATIC_RUNTIME_LIBS   OFF)
find_package(Boost COMPONENTS REQUIRED ${BGMG_BOOST_COMPONENTS})
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find required boost libraries.")
  return()
endif (NOT Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
set(Boost_LIBRARIES_STATIC_LIBS ${Boost_LIBRARIES})
message("-- Boost_LIBRARIES_STATIC_LIBS: ${Boost_LIBRARIES_STATIC_LIBS}")


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
	bgmg_calculator.cc
	ld_matrix_csr.cc
	bgmg.h
	bgmg_calculator.h
	ld_matrix_csr.h
	bgmg_log.cc
	bgmg_log.h
	fmath.hpp
	semt/semt/VectorExpr.cpp
	TurboPFor/vsimple.c
	TurboPFor/vint.c
	TurboPFor/bitutil.c
	nlopt/nldrmd.c
	nlopt/sbplx.c
	nlopt/redblack.c
	nlopt/stop.c
	nlopt/qsort_r.c
	nlopt/timer.c
	FastDifferentialCoding/fastdelta.c
)

add_library(bgmg-static STATIC ${SRC_LIST})
set_target_properties(bgmg-static PROPERTIES LINK_FLAGS "-static" )
target_include_directories(bgmg-static PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(bgmg-static PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(bgmg-static PUBLIC ${PROJECT_SOURCE_DIR}/semt)
target_compile_definitions(bgmg-static PRIVATE USE_SSE)

add_library(bgmg SHARED source.cc)
target_link_libraries(bgmg bgmg-static ${Boost_LIBRARIES_STATIC_LIBS})

add_executable(bgmg-cli source.cc bgmg_main.cc)
if (MSVC)
	target_link_libraries(bgmg-cli bgmg-static ${Boost_LIBRARIES_STATIC_LIBS})
	target_compile_definitions(bgmg-cli PRIVATE BGMG_STATIC_DEFINE)
else (MSVC)
	set_target_properties(bgmg-cli PROPERTIES LINK_FLAGS "-static" )
	target_link_libraries(bgmg-cli bgmg-static ${Boost_LIBRARIES_STATIC_RUNTIME_LIBS})
endif (MSVC)

set(SRC_TEST_LIST
	source.cc
	test/bgmg_main_unittest.cc
	test/semt_test.cc
	test/streamvbyte_test.cc
	test/nlopt_test.cc
	${PROJECT_SOURCE_DIR}/googletest/googletest/src/gtest_main.cc
	${PROJECT_SOURCE_DIR}/googletest/googletest/src/gtest-all.cc
)
add_executable(bgmg-test ${SRC_TEST_LIST})
target_link_libraries(bgmg-test bgmg-static ${Boost_LIBRARIES_STATIC_LIBS}) 
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest/include)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/semt)
target_compile_definitions(bgmg-test PRIVATE BGMG_STATIC_DEFINE)

include(GenerateExportHeader)
generate_export_header(bgmg
		EXPORT_MACRO_NAME DLL_PUBLIC)

install( TARGETS bgmg DESTINATION bin )
install( TARGETS bgmg-cli DESTINATION bin )
