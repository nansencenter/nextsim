// Simple C-functions to instantiate instances of the Environment and
// FinitElement classes, call the run, init, step, and finalise functions, and
// delete the class instances. These are used to interface with fortran through
// finiteelement_mod.f90

#include <finiteelement.hpp>

namespace Nextsim
{
    po::options_description descrOptions();
}

extern "C" {
// Interface to instantiate a new object of the FiniteElement class
    Nextsim::FiniteElement *FiniteElementNew () {
        return new Nextsim::FiniteElement();
  }

// Interface to instantiate a new object of the Environment class
    Nextsim::Environment *EnvironmentNew () {
        // Still need to work out how to pass argv from Fortran
        int argc = 0;
        char **argv;
        return new Nextsim::Environment(argc, argv, Nextsim::descrOptions() );
  }

// Interface to call the function 'run'
  void FiniteElementRun(Nextsim::FiniteElement *This) {
        This->run();
  }

// Interface to call the function 'init'
  void FiniteElementInit(Nextsim::FiniteElement *This, int *pcpt) {
        int pcpt_tmp = This->init();
        pcpt = &pcpt_tmp;
  }

// Interface to call the function 'step'
  void FiniteElementStep(Nextsim::FiniteElement *This, int pcpt) {
        This->step(pcpt);
  }

// Interface to call the function 'finalise'
  void FiniteElementFinalise(Nextsim::FiniteElement *This) {
        This->finalise();
  }

// Interface to delete an instance of the FiniteElement class
  void FiniteElementDelete (Nextsim::FiniteElement *This) {
        delete This;
  }

// Interface to delete an instance of the Environment class
  void EnvironmentDelete (Nextsim::Environment *This) {
        delete This;
  }
}

