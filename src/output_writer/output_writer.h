#pragma once

#include "../Grid/Grid.h"

#include <memory>

/** Inteface class for writing simulation data output.
 */
class OutputWriter
{
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that will contain all the data to be written to the file
  OutputWriter(std::shared_ptr<Grid> staggeredGrid);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime) = 0;

protected:

  std::shared_ptr<Grid> staggeredGrid_;  //< a shared pointer to the discretization which contains all data that will be written to the file
  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};
