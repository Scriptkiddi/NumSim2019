#include "output_writer/output_writer_text.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ios>

void OutputWriterText::writeFile(double currentTime)
{
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << std::setfill('0') << fileNo_ << ".txt";

  // increment file no.
  fileNo_++;

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write time
  file << "t: " << currentTime << std::endl;

  // write mesh width
  file << "nCells: " << staggeredGrid_->nCells()[0] << "x" << staggeredGrid_->nCells()[1]
    << ", dx: " << staggeredGrid_->dx() << ", dy: " << staggeredGrid_->dy() << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write u
  // ---------
  // write header lines
  file << "u (" << staggeredGrid_->u().size()[0] << "x" << staggeredGrid_->u().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ')  << "  "<< "|";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+2; i++)
  {
    file << std::setw(fieldWidth) << i  << "  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->u().size()[0]+2)+1, '-') << std::endl;

  // write u values
  for (int j = staggeredGrid_->jEnd()+1; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file << "  "<< std::setw(fieldWidth) << j << "|" <<  "  " ;
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+2; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->u(i,j) << "  ";

      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->u(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write v
  // ---------
  // write header lines
  file << "v (" << staggeredGrid_->v().size()[0] << "x" << staggeredGrid_->v().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "  "<< "|";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i <<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->v().size()[0]+2)+1, '-') << std::endl;

  // write v values
  for (int j = staggeredGrid_->jEnd()+2; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file <<"  "<< std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->v(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->v(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

  // write p
  // ---------
  // write header lines
  file << "p (" << staggeredGrid_->p().size()[0] << "x" << staggeredGrid_->p().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') <<"  "<< "|";
  for (int i = staggeredGrid_->iBegin()-1; i < staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i<<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = staggeredGrid_->jEnd()+1; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file <<"  " << std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->p(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->p(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

  // write T
  // ---------
  // write header lines
  file << "T (" << staggeredGrid_->t().size()[0] << "x" << staggeredGrid_->t().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|"<<"  ";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i<<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->t().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = staggeredGrid_->jEnd()+1; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file <<"  "<< std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->T(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->t(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

  // write f
  // ---------
  // write header lines
  file << "F (" << staggeredGrid_->u().size()[0] << "x" << staggeredGrid_->u().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|"<<"  ";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+2; i++)
  {
    file << std::setw(fieldWidth) << i<<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->u().size()[0]+2)+1, '-') << std::endl;

  // write f values
  for (int j = staggeredGrid_->jEnd()+1; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file <<"  "<< std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+2; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->f(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->f(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

  // write g
  // ---------
  // write header lines
  file << "G (" << staggeredGrid_->v().size()[0] << "x" << staggeredGrid_->v().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|"<<"  ";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i <<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->v().size()[0]+2)+1, '-') << std::endl;

  // write g values
  for (int j = staggeredGrid_->jEnd()+2; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file<<"  " << std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->g(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->g(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

  // write rhs
  // ---------
  // write header lines
  file << "rhs (" << staggeredGrid_->p().size()[0] << "x" << staggeredGrid_->p().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|"<<"  ";
  for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i<<"  ";
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->p().size()[0]+2)+1, '-') << std::endl;

  // write rhs values
  for (int j = staggeredGrid_->jEnd()+1; j >= staggeredGrid_->jBegin()-1; j--)
  {
    file <<"  "<< std::setw(fieldWidth) << j << "|"<<"  ";
    for (int i = staggeredGrid_->iBegin()-1; i <= staggeredGrid_->iEnd()+1; i++)
    {
      // file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->rhs(i,j);
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << std::scientific << staggeredGrid_->rhs(i,j)<< "  ";
    }
    file << std::endl;
  }
  file << std::endl;

}

void OutputWriterText::writePressureFile()
{
  // counter for files, counter value is part of the file name
  static int pressurefileNo = 0;

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/pressure_" << std::setw(4) << std::setfill('0') << pressurefileNo++ << ".txt";

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write mesh width
  file << "nCells: " << staggeredGrid_->nCells()[0] << "x" << staggeredGrid_->nCells()[1]
    << ", dx: " << staggeredGrid_->dx() << ", dy: " << staggeredGrid_->dy() << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write p
  // ---------
  // write header lines
  file << "p (" << staggeredGrid_->p().size()[0] << "x" << staggeredGrid_->p().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|";
  for (int i = staggeredGrid_->iBegin()-1; i < staggeredGrid_->iEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(staggeredGrid_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = staggeredGrid_->jEnd(); j >= staggeredGrid_->jBegin()-1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = staggeredGrid_->iBegin()-1; i < staggeredGrid_->iEnd()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << staggeredGrid_->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

}
