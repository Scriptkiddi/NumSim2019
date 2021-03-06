cmake_minimum_required(VERSION 3.10)
project(numsim_fluid)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME}
        main.cpp
        Computation/Computation.cpp Computation/Computation.h
        Settings.cpp Settings.h
        Util/Utils.cpp Util/Utils.h
        StaggeredGrid/StaggeredGrid.h StaggeredGrid/StaggeredGrid.cpp
        StaggeredGrid/Discretization.cpp StaggeredGrid/Discretization.h
        StaggeredGrid/DonorCell.cpp StaggeredGrid/DonorCell.h
        StaggeredGrid/CentralDifferences.h StaggeredGrid/CentralDifferences.cpp
        PressureSolver/SOR.h PressureSolver/SOR.cpp
        PressureSolver/PressureSolver.cpp PressureSolver/PressureSolver.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        output_writer/output_writer.cpp output_writer/output_writer.h
        Array2D/Array2D.cpp Array2D/Array2D.h
        Array2D/FieldVariable.cpp Array2D/FieldVariable.h
        Geometry.cpp Geometry.h Util/GeometryParser.cpp Util/GeometryParser.h)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
message("Including ${PROJECT_SOURCE_DIR}")

# Search for the external package "VTK"
find_package(VTK REQUIRED)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})  # add the libraries for the linker # old: target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
endif(VTK_FOUND)


install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Add preCICE to project
find_package(precice REQUIRED CONFIG)
message("${PRECICE_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} PRIVATE ${PRECICE_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PRIVATE precice::precice)

project(numsim_solid)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(WITHPRECICE "Compile with preCICE activated." ON)

if(WITHPRECICE)
    # Add preCICE to project
    #add_definitions(-DWITHPRECICE)
endif()

add_executable(${PROJECT_NAME}
        main_solid.cpp
        Computation/Computation_solid.cpp Computation/Computation_solid.h
        TemperatureSolver/TemperatureSolver.h TemperatureSolver/TemperatureSolver.cpp
        TemperatureSolver/GaussSeidel.h TemperatureSolver/GaussSeidel.cpp
        Settings.cpp Settings.h
        Util/Utils.cpp Util/Utils.h
        StaggeredGrid/StaggeredGrid.h StaggeredGrid/StaggeredGrid.cpp
        StaggeredGrid/Discretization.cpp StaggeredGrid/Discretization.h
        StaggeredGrid/DonorCell.cpp StaggeredGrid/DonorCell.h
        StaggeredGrid/CentralDifferences.h StaggeredGrid/CentralDifferences.cpp
        PressureSolver/SOR.h PressureSolver/SOR.cpp
        PressureSolver/PressureSolver.cpp PressureSolver/PressureSolver.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        output_writer/output_writer.cpp output_writer/output_writer.h
        Array2D/Array2D.cpp Array2D/Array2D.h
        Array2D/FieldVariable.cpp Array2D/FieldVariable.h
        Geometry.cpp Geometry.h Util/GeometryParser.cpp Util/GeometryParser.h)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
message("Including ${PROJECT_SOURCE_DIR}")

# Search for the external package "VTK"
find_package(VTK REQUIRED)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})  # add the libraries for the linker # old: target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
endif(VTK_FOUND)


install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Add preCICE to project
find_package(precice REQUIRED CONFIG)
target_link_libraries(${PROJECT_NAME} PRIVATE precice::precice)
