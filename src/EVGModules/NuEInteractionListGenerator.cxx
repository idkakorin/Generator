//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/NuEInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
NuEInteractionListGenerator::NuEInteractionListGenerator() :
InteractionListGeneratorI("genie::NuEInteractionListGenerator")
{

}
//___________________________________________________________________________
NuEInteractionListGenerator::NuEInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::NuEInteractionListGenerator", config)
{

}
//___________________________________________________________________________
NuEInteractionListGenerator::~NuEInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  if(fIsIMD)  return this -> IMDInteractionList   (init_state);
  else        return this -> NuEELInteractionList (init_state);
}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::IMDInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// numu + e- -> mu- + nu_e [CC] -- 'inverse muon decay'

  if(init_state.ProbePdg() != kPdgNuMu) {
     LOG("IntLst", pDEBUG) 
          << "Return *null* interaction list (non nu_mu probe in IMD thread)";
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  ProcessInfo   proc_info(kScInverseMuDecay, kIntWeakCC);
  Interaction * interaction = new Interaction(init, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
InteractionList * NuEInteractionListGenerator::NuEELInteractionList(
                                       const InitialState & init_state) const
{
// channels:
// nue      + e- -> nue      + e-   [CC + NC + interference]
// nuebar   + e- -> nuebar   + e-   [CC + NC + interference]
// numu     + e- -> numu     + e-   [NC]
// numu     + e- -> mu-      + nu_e [CC] -- handled by the IMD thread
// nutau    + e- -> nutau    + e-   [NC]
// nutau    + e- -> tau    - + nu_e [CC] -- neglected
// numubar  + e- -> numubar  + e-   [NC]
// nutaubar + e- -> nutaubar + e-   [NC]

  int nupdg = init_state.ProbePdg();
  InteractionList * intlist = new InteractionList;

  // clone init state and de-activate the struck nucleon info
  InitialState init(init_state);
  init_state.TgtPtr()->SetHitNucPdg(0);

  // NC
  if(nupdg == kPdgNuMu  || nupdg == kPdgAntiNuMu || 
     nupdg == kPdgNuTau || nupdg == kPdgAntiNuTau) {
     ProcessInfo   proc_info(kScNuElectronElastic,  kIntWeakNC);
     Interaction * interaction = new Interaction(init, proc_info);
     intlist->push_back(interaction);
  }

  // CC+NC+interference
  if(nupdg == kPdgNuE  || nupdg == kPdgAntiNuE) { 
     ProcessInfo   proc_info(kScNuElectronElastic,  kIntWeakMix); 
     Interaction * interaction = new Interaction(init, proc_info);
     intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
void NuEInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuEInteractionListGenerator::LoadConfig(void)
{
  fIsIMD = fConfig->GetBoolDef("is-IMD", false);
}
//____________________________________________________________________________
