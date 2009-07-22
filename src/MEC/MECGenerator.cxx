//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Sep 22, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 22, 2008 - CA
   This event generation modules was first added in version 2.5.1 as part of
   the new event generation thread handling MEC interactions. 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "MEC/MECGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
MECGenerator::MECGenerator() :
EventRecordVisitorI("genie::MECGenerator")
{

}
//___________________________________________________________________________
MECGenerator::MECGenerator(string config) :
EventRecordVisitorI("genie::MECGenerator", config)
{

}
//___________________________________________________________________________
MECGenerator::~MECGenerator()
{

}
//___________________________________________________________________________
void MECGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this -> SelectKinematics    (evrec);
  this -> AddNucleonCluster   (evrec);
  this -> AddTargetRemnant    (evrec);
  this -> AddFinalStateLepton (evrec);
  this -> DecayNucleonCluster (evrec);

  LOG("MEC", pNOTICE) << *evrec;
}
//___________________________________________________________________________
void MECGenerator::SelectKinematics(GHepRecord * evrec) const
{
  evrec->Summary()->KinePtr()->SetQ2(1.0, true);
  evrec->Summary()->KinePtr()->Sety (0.5, true);
}
//___________________________________________________________________________
void MECGenerator::AddFinalStateLepton(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  // Look-up selected kinematics
  double Q2 = interaction->Kine().Q2(true);
  double y  = interaction->Kine().y(true);

  // Auxiliary params
  double Ev  = init_state.ProbeE(kRfLab);
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  // Compute the final state primary lepton energy and momentum components
  // along and perpendicular the neutrino direction 
  double El  = (1-y)*Ev;
  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("LeptonicVertex", pNOTICE)
          << "fsl: E = " << El << ", |p//| = " << plp << "[pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  // Take a unit vector along the neutrino direction
  TVector3 unit_nudir = evrec->Probe()->P4()->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 p3l(pltx,plty,plp);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  TLorentzVector p4l(p3l,El);

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Lepton 4-position (= interacton vtx)
  TLorentzVector v4(*evrec->Probe()->X4());

  evrec->AddParticle(pdgc, kIStStableFinalState, -1, 0, -1, -1, p4l, v4);
}
//___________________________________________________________________________
void MECGenerator::AddNucleonCluster(GHepRecord * evrec) const
{
  const int nc = 3;

  int    clusters [nc] = { kPdgClusterNN, kPdgClusterNP, kPdgClusterPP };
  double prob     [nc] = { 0.25,          0.50,          0.25          };

  GHepParticle * target = evrec->TargetNucleus();
  int Z = target->Z();
  int A = target->A();
  int np = Z;
  int nn = A-Z;

  RandomGen * rnd = RandomGen::Instance();

  bool selected = false;
  int  cluster  = 0;

  while(!selected) {
     double prob_gen = rnd->RndGen().Rndm();
     double prob_sum = 0;
     for(int ic = 0; ic < nc; ic++) {
        prob_sum += prob[ic];
        if(prob_gen < prob_sum) {
           cluster = clusters[ic];
           break;
        }
     }
     if ( cluster==kPdgClusterNN && nn>=2          ) selected=true;
     if ( cluster==kPdgClusterNP && nn>=1 && np>=1 ) selected=true;
     if ( cluster==kPdgClusterPP && np>=2          ) selected=true;
  }

  double mc = PDGLibrary::Instance()->Find(cluster)->Mass();

  const TLorentzVector p4(0.,0.,0., mc);
  const TLorentzVector v4(0.,0.,0., 0.);

  LOG("MEC", pINFO) << "Adding nucleon cluster [pdgc = " << cluster << "]";

  evrec->AddParticle(cluster, kIStNucleonTarget, -1, -1, -1, -1, p4, v4);
}
//___________________________________________________________________________
void MECGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
  GHepParticle * target  = evrec->TargetNucleus();
  GHepParticle * cluster = evrec->Particle(2);

  int Z = target->Z();
  int A = target->A();

  if(cluster->Pdg() == kPdgClusterNN) { A-=2; ;     }
  if(cluster->Pdg() == kPdgClusterNP) { A-=2; Z-=1; }
  if(cluster->Pdg() == kPdgClusterPP) { A-=2; Z-=2; }

  int ipdgc = pdg::IonPdgCode(A, Z);

  const TLorentzVector & p4cluster = *(cluster->P4());
  const TLorentzVector & p4tgt     = *(target ->P4());

  const TLorentzVector p4 = p4tgt - p4cluster;
  const TLorentzVector v4(0.,0.,0., 0.);

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(ipdgc,kIStStableFinalState, mom,-1,-1,-1, p4,v4);  
}
//___________________________________________________________________________
void MECGenerator::DecayNucleonCluster(GHepRecord * evrec) const
{
  int cluster_pos = 2;

  GHepParticle * cluster = evrec->Particle(cluster_pos);

  int nuc1=0, nuc2=0;

  if(cluster->Pdg() == kPdgClusterNN) { nuc1=kPdgNeutron; nuc2=kPdgNeutron; }
  if(cluster->Pdg() == kPdgClusterNP) { nuc1=kPdgNeutron; nuc2=kPdgProton;  }
  if(cluster->Pdg() == kPdgClusterPP) { nuc1=kPdgProton;  nuc2=kPdgProton;  }

  int mom = cluster_pos;

  double m1 = PDGLibrary::Instance()->Find(nuc1)->Mass();
  double m2 = PDGLibrary::Instance()->Find(nuc2)->Mass();

  const TLorentzVector p4_1(0.,0., 1., TMath::Sqrt(1.+m1*m1));
  const TLorentzVector p4_2(0.,0.,-1., TMath::Sqrt(1.+m2*m2));

  const TLorentzVector v4(0.,0.,0.,0.);

  evrec->AddParticle(nuc1, kIStHadronInTheNucleus, mom,-1,-1,-1, p4_1, v4);  
  evrec->AddParticle(nuc2, kIStHadronInTheNucleus, mom,-1,-1,-1, p4_2, v4);  
}
//___________________________________________________________________________

