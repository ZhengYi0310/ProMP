cmake_minimum_required(VERSION 2.6)

project(ProMP)
set(PROJECT_VERSION 1.0)
set(PROJECT_DESCRIPTION "Probabilistic Movement Primitives.")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_path(EIGEN3_INCLUDE_DIRS Eigen/Dense
  ${CMAKE_INSTALL_PREFIX}/include/eigen3
  /usr/include/eigen3
  /opt/local/include/eigen3
  DOC "Eigen 3 include directory")

set(EIGEN3_FOUND ${EIGEN3_INCLUDE_DIRS} CACHE BOOL "" FORCE)
find_package(Eigen3)

include_directories(${EIGEN3_INCLUDE_DIRS} 
                    include)

set(sources src/PhaseSystem.cpp)
set(headers include/ProMP/PhaseSystem.hpp)

add_library(ProMP SHARED ${sources} ${headers})
target_link_libraries(ProMP ${EIGEN3_LIBRARIES})
#BUILD_DMP_BEHAVIOR_STANDALONE enables the behavior macros for the DmpBehavior
set_target_properties(ProMP PROPERTIES COMPILE_FLAGS "-DBUILD_PROMP_BEHAVIOR_STANDALONE")

install (TARGETS ProMP DESTINATION lib)
install (FILES ${headers} DESTINATION include/${PROJECT_NAME})

enable_testing()
add_subdirectory(test)
