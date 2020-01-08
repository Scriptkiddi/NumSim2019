#include "output_writer.h"

#include <iostream>
#include <Settings.h>

OutputWriter::OutputWriter(shared_ptr<Discretization> discretization, Settings settings)
 : discretization_(discretization), settings_(settings), fileNo_(0)
{
  // create "out" subdirectory if it does not yet exist
  int returnValue = system("mkdir -p out");
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;
}

