cmake_minimum_required(VERSION 3.10)
project(numsim)
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
        PressureSolver/GaussSeidel.cpp PressureSolver/GaussSeidel.h
        output_writer/output_writer_paraview.cpp output_writer/output_writer_paraview.h
        output_writer/output_writer_text.cpp output_writer/output_writer_text.h
        output_writer/output_writer.cpp output_writer/output_writer.h
        output_writer/output_writer_paraview_parallel.cpp output_writer/output_writer_paraview_parallel.h
        output_writer/output_writer_text_parallel.cpp output_writer/output_writer_text_parallel.h
        Array2D/Array2D.cpp Array2D/Array2D.h
        Array2D/FieldVariable.cpp Array2D/FieldVariable.h
        Computation/ComputationParallel.cpp Computation/ComputationParallel.h Partitioning/Partitioning.cpp Partitioning/Partitioning.h Util/InputParser.cpp Util/InputParser.h Communication.cpp Communication.h)


target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
message("Including ${PROJECT_SOURCE_DIR}")

# Search for the external package "VTK"
find_package(VTK REQUIRED HINTS C:/Program Files/VTK 7.1.1/bin)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

find_package(MPI REQUIRED)
message("If MPI was found on the system: MPI_FOUND: ${MPI_FOUND}")

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
