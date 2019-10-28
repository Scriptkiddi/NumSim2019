cmake_minimum_required(VERSION 3.15)
project(NumSim2019)

set(CMAKE_CXX_STANDARD 14)

#add_executable(NumSim2019 Computation.cpp Settings.cpp Array2D/Array2D.cpp Array2D/FieldVariable.cpp output_writer/output_writer.cpp output_writer/output_writer_paraview.cpp output_writer/output_writer_text.cpp PressureSolver/PressureSolver.cpp PressureSolver/GaussSeidel.cpp PressureSolver/SOR.cpp StaggeredGrid/StaggeredGrid.cpp StaggeredGrid/CentralDifferences.cpp StaggeredGrid/Discretization.cpp StaggeredGrid/DonorCell.cpp  main.cpp)

add_executable(NumSim2019 Array2D/Array2D.cpp Array2D/FieldVariable.cpp main.cpp)