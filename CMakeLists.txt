#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(KNN_simple_polarimeter)

#----------------------------------------------------------------------------
# Setup compile options
#
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 -Wall")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# include Root
#
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake)
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIR})
include_directories($ENV{ROOTSYS}/include)
#

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
  ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(execute-KNN_simple_polarimeter KNN_simple_polarimeter.cc ${sources} ${headers})
target_link_libraries(execute-KNN_simple_polarimeter ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lexpat)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory.
#
set(KNN_simple_polarimeter_SCRIPTS
  icons.mac
  gui.mac
  init.mac 
  init_vis.mac 
  vis.mac
  run.mac 
  run.png
  )

foreach(_script ${KNN_simple_polarimeter_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Copy all macros to the build directory.
#
set(KNN_simple_polarimeter_MACRO
  Analysis.C 
  DrawTools.h 
  )

foreach(_script ${KNN_simple_polarimeter_MACRO})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macro/${_script}
    ${PROJECT_BINARY_DIR}/macro/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(KNN_simple_polarimeter DEPENDS execute-KNN_simple_polarimeter)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS execute-KNN_simple_polarimeter DESTINATION bin)

