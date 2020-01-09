#include "output_writer.h"

#include <iostream>
#include <Settings.h>

OutputWriter::OutputWriter(shared_ptr<Discretization> discretization, Settings settings)
 : discretization_(discretization), settings_(settings), fileNo_(0)
{
  // create "out" subdirectory if it does not yet exist
  std::string out = "mkdir -p out_" + settings_.participantName;
  int returnValue = system(out.c_str());
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;
}

