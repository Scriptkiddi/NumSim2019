#pragma once

#include "Partitioning/Partitioning.h"

#include <memory>
#include <StaggeredGrid/Discretization.h>

/** Inteface class for writing simulation data output.
 */
class OutputWriter
{
public:
  //! constructor
  OutputWriter(std::shared_ptr<Discretization> discretization, Partitioning &partitioning);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime) = 0;

protected:

  std::shared_ptr<Discretization> discretization_;  //< a shared pointer to the discretization which contains all data that will be written to the file
  Partitioning partitioning_;                 //< the partitioning object that knowns about the domain decomposition, only significant when executing in parallel
  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};
