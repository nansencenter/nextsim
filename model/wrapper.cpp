// Simple C-functions to instantiate a FinitElement class, call the run
// function and delete the class. These are used to interface with fortran
// through finiteelement_mod.f90

#include <finiteelement.hpp>

namespace Nextsim
{
    po::options_description descrOptions();
}

extern "C" {
// Interface to instantiate a new object of the FiniteElement class
    Nextsim::FiniteElement *CFE__new () {
        return new Nextsim::FiniteElement();
  }

// Interface to instantiate a new object of the Environment class
    Nextsim::Environment *Cenv__new () {
        int argc = 0;
        char **argv;
        return new Nextsim::Environment(argc, argv, Nextsim::descrOptions() );
  }

// Interface to call the function 'run'
  void CFE__run(Nextsim::FiniteElement *This) {
        This->run();
  }

// Interface to call the function 'init'
  void CFE__init(Nextsim::FiniteElement *This, int *pcpt) {
        int pcpt_tmp = This->init();
        pcpt = &pcpt_tmp;
  }

// Interface to call the function 'step'
  void CFE__step(Nextsim::FiniteElement *This, int pcpt) {
        This->step(pcpt);
  }

// Interface to call the function 'finalise'
  void CFE__finalise(Nextsim::FiniteElement *This) {
        This->finalise();
  }

// Interface to delete an instance of the FiniteElement class
  void CFE__delete (Nextsim::FiniteElement *This) {
        delete This;
  }

// Interface to delete an instance of the Environment class
  void Cenv__delete (Nextsim::Environment *This) {
        delete This;
  }
}

