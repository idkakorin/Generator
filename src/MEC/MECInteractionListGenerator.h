//____________________________________________________________________________
/*!

\class    genie::MECInteractionListGenerator

\brief    Concrete implementations of the InteractionListGeneratorI interface.
          Generates a list of all the interactions that can be generated by the 
          MEC EventGenerator.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 22, 2008

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_INTERACTION_GENERATOR_H_
#define _MEC_INTERACTION_GENERATOR_H_

#include "EVGCore/InteractionListGeneratorI.h"

namespace genie {

class MECInteractionListGenerator : public InteractionListGeneratorI {

public :
  MECInteractionListGenerator();
  MECInteractionListGenerator(string config);
 ~MECInteractionListGenerator();

  //-- implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _MEC_INTERACTION_GENERATOR_H_