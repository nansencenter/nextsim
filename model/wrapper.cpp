// Simple C-functions to instantiate a FinitElement class, call the run
// function and delete the class. These are used to interface with fortran
// through finiteelement_mod.f90

#include <finiteelement.hpp>

extern "C" {
// Interface to instantiate a new class
    Nextsim::FiniteElement *CFE__new () {
    return new Nextsim::FiniteElement();
  }

// Interface to call the function 'run'
  void CFE__run(Nextsim::FiniteElement *This) {
    This->run();
  }

// Interface to delete an instance of the FiniteElement class
  void CFE__delete (Nextsim::FiniteElement *This) {
    delete This;
  }
}

