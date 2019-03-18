//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <vector>

#include <TMath.h>

// For the minimizer
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/BrentMinimizer1D.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/QuasiElastic/EventGen/QELEventGeneratorSuSAMstar.h"


#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

namespace { // anonymous namespace (file only visibility)
  const double eps = std::numeric_limits<double>::epsilon();
}
//___________________________________________________________________________
QELEventGeneratorSuSAMstar::QELEventGeneratorSuSAMstar() :
KineGeneratorWithCache("genie::QELEventGeneratorSuSAMstar")
{

}
//___________________________________________________________________________
QELEventGeneratorSuSAMstar::QELEventGeneratorSuSAMstar(string config) :
KineGeneratorWithCache("genie::QELEventGeneratorSuSAMstar", config)
{

}
//___________________________________________________________________________
QELEventGeneratorSuSAMstar::~QELEventGeneratorSuSAMstar()
{

}
//___________________________________________________________________________
void QELEventGeneratorSuSAMstar::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("QELEvent", pINFO) << "Generating QE event kinematics...";

  if(fGenerateUniformly) {
    LOG("QELEvent", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  // Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  Kinematics * kinematics = interaction->KinePtr();
  const InitialState & init_state = interaction -> InitState();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);
  
  // Skip if no hit nucleon is set
  if(! evrec->HitNucleon())
  {
    LOG("QELEvent", pFATAL) << "No hit nucleon was set";
    gAbortingInErr = true;
    exit(1);
  }

  // Access the target from the interaction summary
  Target * tgt = init_state.TgtPtr();
  GHepParticle * nucleon = evrec->HitNucleon();
  // Store position of nucleon
  double hitNucPos = nucleon->X4()->Vect().Mag();
  tgt->SetHitNucPosition( hitNucPos );

  // Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  
  // heavy nucleus is nucleus that heavier than hydrogen and deuterium
  fIsHeavyNucleus = tgt->A()>=3;

  // phase space for heavy nucleus is different from light one
  fkps = fIsHeavyNucleus?kPSTlctl:kPSQ2fE;
  // Try to calculate the maximum cross-section in kinematical limits
  // if not pre-computed already
  
  double Enu = interaction->InitState().ProbeE(kRfLab);
  double M = init_state.Tgt().HitNucP4().M();
  
  if (fIsHeavyNucleus)
  {
     sm_utils->SetInteraction(interaction);
     vector<vector<double>> Coords(2);
     Coords[0] = vector<double>(3);
     Coords[1] = vector<double>(3);
     vector<double> maxs(4);
     MaxsSubAreas(evrec, Coords, maxs);
         
     double probs[4];
     probs[0] =            maxs[0]*(Coords[0][1]-Coords[0][0])*(Coords[1][1]-Coords[1][0]);
     probs[1] = probs[0] + maxs[1]*(Coords[0][1]-Coords[0][0])*(Coords[1][2]-Coords[1][1]);
     probs[2] = probs[1] + maxs[2]*(Coords[0][2]-Coords[0][1])*(Coords[1][1]-Coords[1][0]);
     probs[3] = probs[2] + maxs[3]*(Coords[0][2]-Coords[0][1])*(Coords[1][2]-Coords[1][1]);
     
     probs[0] /= probs[3];
     probs[1] /= probs[3];
     probs[2] /= probs[3];
     probs[3] /= probs[3];
     
     
     unsigned int iter = 0;
     bool accept = false;
     double El, cost;
     double xsec   = -1., xsec_max;
     while(1)
     {
        LOG("QELEvent", pINFO) << "Attempt #: " << iter;
        if(iter > 100*kRjMaxIterations)
        {
           LOG("QELEvent", pWARN)
             << "Couldn't select a valid kinematics after " << iter << " iterations";
           evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
           genie::exceptions::EVGThreadException exception;
           exception.SetReason("Couldn't select kinematics");
           exception.SwitchOnFastForward();
           throw exception;
        }
        
        double ml  = interaction->FSPrimLepton()->Mass();
        if (Enu < ml)
        {
           evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
           genie::exceptions::EVGThreadException exception;
           exception.SetReason("Couldn't select kinematics");
           exception.SwitchOnFastForward();
           throw exception;
        }
  
        double r = rnd->RndKine().Rndm();

        if (r < probs[0])
        {
           El = (Coords[0][1]-Coords[0][0])*rnd->RndKine().Rndm() + Coords[0][0];
           cost = (Coords[1][1]-Coords[1][0])*rnd->RndKine().Rndm() + Coords[1][0];
           xsec_max = maxs[0];
        }
        else if (r < probs[1])
        {
           El = (Coords[0][1]-Coords[0][0])*rnd->RndKine().Rndm() + Coords[0][0];
           cost = (Coords[1][2]-Coords[1][1])*rnd->RndKine().Rndm() + Coords[1][1];
           xsec_max = maxs[1];
        }
        else if (r < probs[2])
        {
           El = (Coords[0][2]-Coords[0][1])*rnd->RndKine().Rndm() + Coords[0][1];
           cost = (Coords[1][1]-Coords[1][0])*rnd->RndKine().Rndm() + Coords[1][0];
           xsec_max = maxs[2];
        }
        else
        {
           El = (Coords[0][2]-Coords[0][1])*rnd->RndKine().Rndm() + Coords[0][1];
           cost = (Coords[1][2]-Coords[1][1])*rnd->RndKine().Rndm() + Coords[1][1];
           xsec_max = maxs[3];
        }
        

        double Pl = TMath::Sqrt(El*El - ml*ml);
        double sint = TMath::Sqrt(1 - cost*cost);
        double phi = 2*TMath::Pi()*rnd->RndKine().Rndm();
        kinematics->SetFSLeptonP4(Pl*sint*TMath::Cos(phi), Pl*sint*TMath::Sin(phi), Pl*cost, El);
        xsec = fXSecModel->XSec(interaction, kPSTlctl);
        
     
         //-- Decide whether to accept the current kinematics
        if(!fGenerateUniformly)
        {
          this->AssertXSecLimits(interaction, xsec, xsec_max);
          double t = xsec_max * rnd->RndKine().Rndm();
          accept = (t < xsec);
        }
        else 
        {
           accept = (xsec>0);
        }
     
        // If the generated kinematics are accepted, finish-up module's job
        if(accept)
        {
          interaction->ResetBit(kISkipProcessChk);
          interaction->ResetBit(kISkipKinematicChk);
          break;
        }
        iter++;
     }
     
     
     // Momentum of initial neutrino in LAB frame
     TLorentzVector * tempTLV = evrec->Probe()->GetP4();
     TLorentzVector neutrinoMom = *tempTLV;
     delete tempTLV;
     
     
     // 4-momentum of final lepton in LAB frame
     const TLorentzVector outLeptonMom = kinematics->FSLeptonP4();
     
     TLorentzVector transferMom = neutrinoMom - outLeptonMom;
     double gQ2 = -transferMom.Mag2();
     
     int rpdgc = interaction->RecoilNucleonPdg();
     assert(rpdgc);
     double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();
     LOG("QELEvent", pNOTICE) << "Selected: W = "<< gW;
     

     double costheta = -1. + 2. * rnd->RndKine().Rndm();
     double sintheta = TMath::Sqrt(1.-costheta*costheta);
     double fi       = 2 * kPi * rnd->RndKine().Rndm();
     double cosfi    = TMath::Cos(fi);
     double sinfi    = TMath::Sin(fi);
     
     double pF = sm_utils->GetFermiMomentum();
     
     double px = pF*sintheta*cosfi;
     double py = pF*sintheta*sinfi;
     double pz = pF*costheta;
     
     TVector3 pV(px, py, pz);
     TVector3 temp = pV + transferMom.Vect();
     
     
     TLorentzVector inNucleonMom(pV, TMath::Sqrt(gW*gW + temp.Mag2()) - transferMom.E());
     
     TLorentzVector outNucleonMom = inNucleonMom + transferMom;
     
     // set 4-momentum of struck nucleon
     TLorentzVector * p4 = tgt->HitNucP4Ptr();
     p4->SetPx( inNucleonMom.Px() );
     p4->SetPy( inNucleonMom.Py() );
     p4->SetPz( inNucleonMom.Pz() );
     p4->SetE ( inNucleonMom.E() );
          
     // (W,Q2) -> (x,y)
     double gx=0, gy=0;
     kinematics::WQ2toXY(Enu,M,gW,gQ2,gx,gy);
     
     // lock selected kinematics & clear running values
     interaction->KinePtr()->SetQ2(gQ2, true);
     interaction->KinePtr()->SetW (gW,  true);
     interaction->KinePtr()->Setx (gx,  true);
     interaction->KinePtr()->Sety (gy,  true);
     interaction->KinePtr()->SetKV(kKVSelTl, El);    // lepton kinetic energy
     interaction->KinePtr()->SetKV(kKVSelctl, cost); // cosine lepton theta
     interaction->KinePtr()->ClearRunningValues();
     
     // set the cross section for the selected kinematics
     evrec->SetDiffXSec(xsec,fkps);
     if(fGenerateUniformly) 
     {
        double vol     = sm_utils->PhaseSpaceVolume(fkps);
        double totxsec = evrec->XSec();
        double wght    = (vol/totxsec)*xsec;
        LOG("QELEvent", pNOTICE)  << "Kinematics wght = "<< wght;
        
        // apply computed weight to the current event weight
        wght *= evrec->Weight();
        LOG("QELEvent", pNOTICE) << "Current event wght = " << wght;
        evrec->SetWeight(wght);
     }
     TLorentzVector x4l(*(evrec->Probe())->X4());
     
     evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, outLeptonMom, x4l);
     
     GHepParticle outNucleon(interaction->RecoilNucleonPdg(), kIStHadronInTheNucleus, evrec->HitNucleonPosition(),-1,-1,-1, outNucleonMom , x4l);
     evrec->AddParticle(outNucleon);
     
     // Store struck nucleon momentum
     LOG("QELEvent",pNOTICE) << "pn: " << inNucleonMom.X() << ", " <<inNucleonMom.Y() << ", " <<inNucleonMom.Z() << ", " <<inNucleonMom.E();
     nucleon->SetMomentum(inNucleonMom);
     
     //nucleon->SetRemovalEnergy(0.);
     
     // skip if not a nuclear target
     if(evrec->Summary()->InitState().Tgt().IsNucleus())
       // add a recoiled nucleus remnant
       this->AddTargetNucleusRemnant(evrec); 
  }
  else
  {
    qel_kin_gen->ProcessEventRecord(evrec);
    qel_primlep_gen->ProcessEventRecord(evrec);  
    qel_hadrsys_gen->ProcessEventRecord(evrec); 
  }
  LOG("QELEvent", pINFO) << "Done generating QE event kinematics!";
  return; 
}
//___________________________________________________________________________
void QELEventGeneratorSuSAMstar::AddTargetNucleusRemnant(GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("QELEvent", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  for(int id = fd; id <= ld; id++) 
  {

     // compute A,Z for final state nucleus & get its PDG code and its mass
     GHepParticle * particle = evrec->Particle(id);
     assert(particle);
     int  pdgc = particle->Pdg();
     bool is_p  = pdg::IsProton (pdgc);
     bool is_n  = pdg::IsNeutron(pdgc);
     
     if (is_p) Z--;
     if (is_p || is_n) A--;
     
     Px += particle->Px();
     Py += particle->Py();
     Pz += particle->Pz();
     E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn)
  {
    LOG("QELEvent", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("QELEvent", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

  LOG("QELEvent", pINFO) << "Done";
  LOG("QELEvent", pINFO) << *evrec;
  
}
//___________________________________________________________________________
void QELEventGeneratorSuSAMstar::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorSuSAMstar::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r("QELEventGeneratorSuSAMstar_specific", false ) ;
  r.Set("sm_utils_algo", RgAlg("genie::SuSAMstarUtils","Default") ) ;
  r.Set("qel_kin_gen_algo", RgAlg("genie::QELKinematicsGenerator","CC-Default") ) ;
  r.Set("qel_primlep_gen_algo", RgAlg("genie::QELPrimaryLeptonGenerator","Default") ) ;
  r.Set("qel_hadrsys_gen_algo", RgAlg("genie::QELHadronicSystemGenerator","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorSuSAMstar::LoadConfig(void)
{

  // Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,   1.2) ;

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  GetParamDef( "Cache-MinEnergy", fEMin,   1.00) ;
 //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false);
     
  GetParamDef( "NumOfKnots_x", Nx, 20);
  GetParamDef( "NumOfKnots_y", Ny, 20);
  
  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999.);
  assert(fMaxXSecDiffTolerance>=0);

  
  sm_utils = const_cast<genie::SuSAMstarUtils *>(
             dynamic_cast<const genie::SuSAMstarUtils *>(
             this -> SubAlg( "sm_utils_algo" ) ) ) ;
             
  qel_kin_gen = const_cast<genie::QELKinematicsGenerator *>(
                dynamic_cast<const genie::QELKinematicsGenerator *>(
                this -> SubAlg( "qel_kin_gen_algo" ) ) ) ;
             
  qel_primlep_gen = const_cast<genie::QELPrimaryLeptonGenerator *>(
                    dynamic_cast<const genie::QELPrimaryLeptonGenerator *>(
                    this -> SubAlg( "qel_primlep_gen_algo" ) ) ) ;
             
  qel_hadrsys_gen = const_cast<genie::QELHadronicSystemGenerator *>(
                    dynamic_cast<const genie::QELHadronicSystemGenerator *>(
                    this -> SubAlg( "qel_hadrsys_gen_algo" ) ) ) ;
  
}
//____________________________________________________________________________
double QELEventGeneratorSuSAMstar::ComputeMaxXSec(const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small dQ2 step.

  double max_xsec = 0.0;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  if(rQ2.min <=0 || rQ2.max <= rQ2.min) return 0.;

  const double logQ2min = TMath::Log(rQ2.min + kASmallNum);
  const double logQ2max = TMath::Log(rQ2.max - kASmallNum);

  const int N  = 15;
  const int Nb = 10;

  double dlogQ2   = (logQ2max - logQ2min) /(N-1);
  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->KinePtr()->SetQ2(Q2);
     double xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("QELEvent", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
     max_xsec = TMath::Max(xsec, max_xsec);
     increasing = xsec-xseclast>=0;
     xseclast   = xsec;

     // once the cross section stops increasing, I reduce the step size and
     // step backwards a little bit to handle cases that the max cross section
     // is grossly underestimated (very peaky distribution & large step)
     if(!increasing) {
       dlogQ2/=(Nb+1);
       for(int ib=0; ib<Nb; ib++) {
	 Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
         if(Q2 < rQ2.min) continue;
         interaction->KinePtr()->SetQ2(Q2);
         xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("QELEvent", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
     }
  }//Q^2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("QELEvent", pDEBUG) << interaction->AsString();
  SLOG("QELEvent", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("QELEvent", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//____________________________________________________________________________
void QELEventGeneratorSuSAMstar::MaxsSubAreas(GHepRecord * event_rec, vector<vector<double>> & Coords, vector<double> & maxs) const
{
  
  Interaction * interaction = event_rec->Summary();

  // Find a cached max xsec for the specified xsec algorithm & interaction and
  // close to the specified energy
  
  double m_lep = interaction->FSPrimLepton()->Mass();
  // get neutrino energy
  double E_nu = interaction->InitState().ProbeE(kRfLab);
  Coords[0][0] = m_lep;
  Coords[0][2] = E_nu;
  Coords[1][0] = -1.;
  Coords[1][2] = 1.;

  if(E_nu < fEMin)
  {
     LOG("QELEvent", pINFO) << "Below minimum energy - Forcing explicit calculation";
     this->ComputeMaxsSubAreas(interaction, Coords, maxs);
     return;
  }

  vector<CacheBranchFx*> branches(6);
  
  // access the the cache branches
  this->AccessCacheBranches(interaction, branches);
  double spl_val;
  bool success;
  bool computed = false;
  for (int i = 0; i < 6; i++)
  {
    success = false;
    // access the the cache branch
    CacheBranchFx * cb = branches[i];
    // if there are enough points stored in the cache buffer to build a
    // spline, then intepolate
    if( cb->Spl() )
    {
       if( E_nu >= cb->Spl()->XMin() && E_nu <= cb->Spl()->XMax())
       {
         spl_val = cb->Spl()->Evaluate(E_nu);
         if (spl_val > 0) success = true;
       }
    }
    else
    {
      // if there are not enough points at the cache buffer to have a spline,
      // look whether there is another point that is sufficiently close
      double dE = TMath::Min(0.25, 0.05*E_nu);
      const map<double,double> & fmap = cb->Map();
      map<double,double>::const_iterator iter = fmap.lower_bound(E_nu);
      if(iter != fmap.end()) 
      {
         if(TMath::Abs(E_nu - iter->first) < dE)
         {
           spl_val = iter->second;
           success = true;
         }
      }
    
    }
    if (success)
    {
      if (i == 0)
        Coords[0][1] = spl_val;
      else if (i == 1)
        Coords[1][1] = spl_val;
      else if (i == 2)
        maxs[0] = spl_val;
      else if (i == 3)
        maxs[1] = spl_val;
      else if (i == 4)
        maxs[2] = spl_val;
      else if (i == 5)
        maxs[3] = spl_val;
    }
    else
    {
      if (!computed)
      {
        this->ComputeMaxsSubAreas(interaction, Coords, maxs);
        if (maxs[0]<0 || maxs[1]<0 || maxs[2]<0 || maxs[3]<0)
        {
          // xsec for selected kinematics = 0
          event_rec->SetDiffXSec(0,kPSNull);
          // switch on error flag
          event_rec->EventFlags()->SetBitNumber(kKineGenErr, true);
          // reset 'trust' bits
          interaction->ResetBit(kISkipProcessChk);
          interaction->ResetBit(kISkipKinematicChk);
          // throw exception
          genie::exceptions::EVGThreadException exception;
          exception.SetReason("kinematics generation: max_xsec({K};E)<=0");
          exception.SwitchOnFastForward();
          throw exception;
        }
        computed = true;
      }
      if (i == 0)
        cb->AddValues(E_nu,Coords[0][1]);
      else if (i == 1)
        cb->AddValues(E_nu,Coords[1][1]);
      else if (i == 2 && maxs[0]>0)
        cb->AddValues(E_nu,maxs[0]);
      else if (i == 3 && maxs[1]>0)
        cb->AddValues(E_nu,maxs[1]);
      else if (i == 4 && maxs[2]>0)
        cb->AddValues(E_nu,maxs[2]);
      else if (i == 5 && maxs[3]>0)
        cb->AddValues(E_nu,maxs[3]);
                
      if(! cb->Spl() ) {
        if( cb->Map().size() > 40 ) 
          cb->CreateSpline();
      }
      
      if( cb->Spl() ) {
         if( E_nu < cb->Spl()->XMin() || E_nu > cb->Spl()->XMax() ) {
            cb->CreateSpline();
         }
      }
      
    }
    
  }

}
//____________________________________________________________________________
void QELEventGeneratorSuSAMstar::AccessCacheBranches(const Interaction * interaction, vector<CacheBranchFx*> & branches) const
{
  // Returns the cache branches for this algorithm and this interaction. If no
  // branches is found then one is created.

  Cache * cache = Cache::Instance();

  // build the cache branch key as: namespace::algorithm/config/interaction
  string algkey = this->Id().Key();
  string intkey = interaction->AsString();
  string key;
  CacheBranchFx * cache_branch;
  
  for (int i = 0; i < 6; i++)
  {
    if (i == 0)
      key = cache->CacheBranchKey(algkey, intkey, "X1");
    else if (i == 1)
      key = cache->CacheBranchKey(algkey, intkey, "Y1");
    else if (i == 2)
      key = cache->CacheBranchKey(algkey, intkey, "max0");
    else if (i == 3)
      key = cache->CacheBranchKey(algkey, intkey, "max1");
    else if (i == 4)
      key = cache->CacheBranchKey(algkey, intkey, "max2");
    else if (i == 5)
      key = cache->CacheBranchKey(algkey, intkey, "max3");
    
    cache_branch = dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
    
    if(!cache_branch) 
    {
      //-- create the cache branch at the first pass
      LOG("QELEvent", pINFO) << "Creating cache branch - key = " << key;
      char branch_name[30];
      sprintf(branch_name, "branch#%d for SuSAM* model", i);
      cache_branch = new CacheBranchFx(branch_name);
      cache->AddCacheBranch(key, cache_branch);
      assert(cache_branch);
    }
    
    branches[i] = cache_branch;
    
  }

  
}
//____________________________________________________________________________
void QELEventGeneratorSuSAMstar::ComputeMaxsSubAreas(const Interaction * interaction, vector<vector<double>> & Coords, vector<double> & maxs) const
{
  /*
     Areas:
     Y:
     2 ------------
       | 1  |  3  |
     1 ------------
       | 0  |  2  |
     0 ------------
   X:0      1     2
     
     The point 3 lies in area 3 and corresponds to global maxiumum


  */
  
  
  
  Interaction * in = new Interaction (*interaction);
  Kinematics * kinematics = in->KinePtr();
  vector<double> Xs(3);
  vector<double> Ys(3);
  double m_lep = interaction->FSPrimLepton()->Mass();
  double E_nu = interaction->InitState().ProbeE(kRfLab);
  Xs[0] = m_lep;
  Xs[2] = E_nu;
  Ys[0] = -1.;
  Ys[2] = 1.;
  double Xs3, Ys3 = 1.;
  // fast maximum search
  double Fmax = -1.;
  int ixmax = -1;
  for (int ix = 0; ix <= Nx; ix++)
  {
    double El  = E_nu + m_lep - m_lep*TMath::Power(E_nu/m_lep, 1.*ix/Nx);
    double Pl2  = El*El-m_lep*m_lep;
    double Pl = Pl2>0?TMath::Sqrt(Pl2):0;
    kinematics->SetFSLeptonP4(0, 0, Pl, El);
    double F = fXSecModel->XSec(in, kPSTlctl);
    //std::cout << "f1(" << El << "," << cost << "): "<< F << std::endl;
    if (F > Fmax)
    {
      Fmax = F;
      ixmax = ix;
    }
  }
  int ixleft = TMath::Max(ixmax-1,0);
  int ixright = TMath::Min(ixmax+1,Nx);
  
  //std::cout << "Minimum rude: f(" << E_nu + m_lep - m_lep*TMath::Power(E_nu/m_lep, 1.*ixmax/Nx) << "," << 2. - TMath::Power(3., 1.*iymax/Ny) << "): "<< Fmax << std::endl;
  
  // The feature of SuSAM* model, that it has a maximum at cos_theta = 1
  ROOT::Math::BrentMinimizer1D bm;
  ROOT::Math::IBaseFunctionOneDim * func = new genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed(fXSecModel, interaction);
  bm.SetFunction(*func, E_nu + m_lep - m_lep*TMath::Power(E_nu/m_lep, 1.*ixleft/Nx), E_nu + m_lep - m_lep*TMath::Power(E_nu/m_lep, 1.*ixright/Nx));
  bm.SetLogScan(true);
  bm.Minimize(100000);
  //std::cout << "Minimum: f(" << bm.XMinimum() << ") = " <<bm.FValMinimum() << std::endl;
  double max_xsec = -bm.FValMinimum();
  
  if (max_xsec>Fmax)
  {
    Xs3 = bm.XMinimum();
  }
  else
  {
    max_xsec = Fmax;
    Xs3 = E_nu + m_lep - m_lep*TMath::Power(E_nu/m_lep, 1.*ixmax/Nx);
  }
  maxs[3] = fSafetyFactor*max_xsec;  //global maximum in area 3
  
  int status = bm.Status();
  //std::cout << "Status: " << status << std::endl;

  delete func;
  
  ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
  
  // Empirical coefficient, which determines the area of the peak
  double MaxCoeff = 1e13;
  if (E_nu < 0.4)
    MaxCoeff = 1e1;
  else if (E_nu < 0.5)
    MaxCoeff = 1e2;
  else if (E_nu < 0.9)
    MaxCoeff = 1e3;
  else if (E_nu < 2)
    MaxCoeff = 1e4;
  else if (E_nu < 3)
    MaxCoeff = 1e6;
  else if (E_nu < 5)
    MaxCoeff = 1e7;
  else if (E_nu < 7)
    MaxCoeff = 1e8;
  else if (E_nu < 10)
    MaxCoeff = 1e9;
  else if (E_nu < 20)
    MaxCoeff = 1e10;
  else if (E_nu < 30)
    MaxCoeff = 1e11;
  else if (E_nu < 40)
    MaxCoeff = 1e12;  
  double level = max_xsec/MaxCoeff;
  double xrf, yrf;
  
  ROOT::Math::IBaseFunctionOneDim * fEl    = new genie::utils::gsl::d2Xsec_dEldCosThetal_fixed(fXSecModel, interaction, 1, Ys3, level);
  ROOT::Math::IBaseFunctionOneDim * fcostl = new genie::utils::gsl::d2Xsec_dEldCosThetal_fixed(fXSecModel, interaction, 0, Xs3, level);
  
  rfgb.Solve(*fEl, Xs[0], Xs3, 1000, 0., 1.0e-08);
  xrf = rfgb.Root();
  status = rfgb.Status();
  if (status == 0)
    Xs[1] = xrf;
  else
    Xs[1] = Xs[0];
      
  rfgb.Solve(*fcostl, Ys[0], Ys3, 1000, 0., 1.0e-08);
  yrf = rfgb.Root();
  status = rfgb.Status();
  if (status == 0)
    Ys[1] = yrf;
  else
    Ys[1] = Ys[0];
      
  Coords[0] = Xs;
  Coords[1] = Ys;
  
  delete fcostl;
  delete fEl;
  
  
  ROOT::Math::Minimizer * max = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
  ROOT::Math::IBaseFunctionMultiDim * f = new genie::utils::gsl::d2Xsec_dEldCosThetal4Min(fXSecModel, interaction);
  max->SetFunction( *f );
  max->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
  max->SetMaxIterations(10000);     // for GSL
  max->SetTolerance(0.001);
  max->SetPrintLevel(0);
  double step[2] = {1e-7,1e-7};
  double variable[2];
  
  // fast maximum search in area 0
  Fmax = -1.;
  ixmax = -1;
  int iymax = -1;
  for (int ix = 0; ix <= Nx; ix++)
  {
    double El  = Xs[0] + ix*(Xs[1]-Xs[0])/Nx;
    double Pl2  = El*El-m_lep*m_lep;
    double Pl = Pl2>0?TMath::Sqrt(Pl2):0;
    for (int iy = 0; iy <= Ny; iy++)
    {
      double cost = Ys[0] + iy*(Ys[1]-Ys[0])/Ny;
      double sint   = TMath::Sqrt(1. - cost*cost);
      kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
      double F = fXSecModel->XSec(in, kPSTlctl);
      if (F > Fmax)
      {
        Fmax = F;
        ixmax = ix;
        iymax = iy;
      }
    }
  }
  ixleft = TMath::Max(ixmax-1,0);
  ixright = TMath::Min(ixmax+1,Nx);
  int iydown = TMath::Max(iymax-1,0);
  int iyup = TMath::Min(iymax+1,Ny);
    
  variable[0] = Xs[0] + ixmax*(Xs[1]-Xs[0])/Nx;
  variable[1] = Ys[0] + iymax*(Ys[1]-Ys[0])/Ny;
  max->SetVariable(0,"El",variable[0], step[0]);
  max->SetVariable(1,"cost",variable[1], step[1]);
  max->SetVariableLimits(0, Xs[0] + ixleft*(Xs[1]-Xs[0])/Nx, Xs[0] + ixright*(Xs[1]-Xs[0])/Nx);
  max->SetVariableLimits(1, Ys[0] + iydown*(Ys[1]-Ys[0])/Ny, Ys[0] + iyup*(Ys[1]-Ys[0])/Ny);
  max->Minimize();
  max_xsec = -max->MinValue();
  
  if (max_xsec < Fmax)
  {
    max_xsec = Fmax;
  }
  maxs[0] = fSafetyFactor*max_xsec;  //global maximum in area 0
  
  // fast maximum search in area 1
  Fmax = -1.;
  ixmax = -1;
  iymax = -1;
  for (int ix = 0; ix <= Nx; ix++)
  {
    double El  = Xs[0] + ix*(Xs[1]-Xs[0])/Nx;
    double Pl2  = El*El-m_lep*m_lep;
    double Pl = Pl2>0?TMath::Sqrt(Pl2):0;
    for (int iy = 0; iy <= Ny; iy++)
    {
      double cost = Ys[1] + iy*(Ys[2]-Ys[1])/Ny;
      double sint   = TMath::Sqrt(1. - cost*cost);
      kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
      double F = fXSecModel->XSec(in, kPSTlctl);
      //std::cout << "f1(" << El << "," << cost << "): "<< F << std::endl;
      if (F > Fmax)
      {
        Fmax = F;
        ixmax = ix;
        iymax = iy;
      }
    }
  }
  ixleft = TMath::Max(ixmax-1,0);
  ixright = TMath::Min(ixmax+1,Nx);
  iydown = TMath::Max(iymax-1,0);
  iyup = TMath::Min(iymax+1,Ny);
  
  variable[0] = Xs[0] + ixmax*(Xs[1]-Xs[0])/Nx;
  variable[1] = Ys[1] + iymax*(Ys[2]-Ys[1])/Ny;
  max->SetVariable(0,"El",variable[0], step[0]);
  max->SetVariable(1,"cost",variable[1], step[1]);
  max->SetVariableLimits(0, Xs[0] + ixleft*(Xs[1]-Xs[0])/Nx, Xs[0] + ixright*(Xs[1]-Xs[0])/Nx);
  max->SetVariableLimits(1, Ys[1] + iydown*(Ys[2]-Ys[1])/Ny, Ys[1] + iyup*(Ys[2]-Ys[1])/Ny);
  max->Minimize();
  max_xsec = -max->MinValue();
  
  if (max_xsec < Fmax)
  {
    max_xsec = Fmax;
  }
  maxs[1] = fSafetyFactor*max_xsec;  //global maximum in area 1
  
  // fast maximum search in area 2
  Fmax = -1.;
  ixmax = -1;
  iymax = -1;
  for (int ix = 0; ix <= Nx; ix++)
  {
    double El  = Xs[1] + ix*(Xs[2]-Xs[1])/Nx;
    double Pl2  = El*El-m_lep*m_lep;
    double Pl = Pl2>0?TMath::Sqrt(Pl2):0;
    for (int iy = 0; iy <= Ny; iy++)
    {
      double cost = Ys[0] + iy*(Ys[1]-Ys[0])/Ny;
      double sint   = TMath::Sqrt(1. - cost*cost);
      kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
      double F = fXSecModel->XSec(in, kPSTlctl);
      //std::cout << "f1(" << El << "," << cost << "): "<< F << std::endl;
      if (F > Fmax)
      {
        Fmax = F;
        ixmax = ix;
        iymax = iy;
      }
    }
  }
  ixleft = TMath::Max(ixmax-1,0);
  ixright = TMath::Min(ixmax+1,Nx);
  iydown = TMath::Max(iymax-1,0);
  iyup = TMath::Min(iymax+1,Ny);
  
  variable[0] = Xs[1] + ixmax*(Xs[2]-Xs[1])/Nx;
  variable[1] = Ys[0] + iymax*(Ys[1]-Ys[0])/Ny;
  max->SetVariable(0,"El",variable[0], step[0]);
  max->SetVariable(1,"cost",variable[1], step[1]);
  max->SetVariableLimits(0, Xs[1] + ixleft*(Xs[2]-Xs[1])/Nx, Xs[1] + ixright*(Xs[2]-Xs[1])/Nx);
  max->SetVariableLimits(1, Ys[0] + iydown*(Ys[1]-Ys[0])/Ny, Ys[0] + iyup*(Ys[1]-Ys[0])/Ny);
  max->Minimize();
  max_xsec = -max->MinValue();
  
  if (max_xsec < Fmax)
  {
    max_xsec = Fmax;
  }
  maxs[2] = fSafetyFactor*max_xsec;  //global maximum in area 2
  
  delete max;
  delete f;
  delete in;
  
  
  
}
//_____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal4Min::d2Xsec_dEldCosThetal4Min(const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(interaction)
{

}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal4Min::~d2Xsec_dEldCosThetal4Min()
{
  
}   
//____________________________________________________________________________
unsigned int genie::utils::gsl::d2Xsec_dEldCosThetal4Min::NDim(void) const
{
  return 2;
}
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dEldCosThetal4Min::DoEval(const double * xin) const
{
// inputs:
//    lepton energy
//    cosine of lepton angle
// outputs:
//   differential cross section
//
 
  Kinematics * kinematics = fInteraction->KinePtr();
  double ml   = fInteraction->FSPrimLepton()->Mass();
  double El   = xin[0];
  double Pl = TMath::Sqrt(El*El - ml*ml);
  double cost = xin[1];
  double sint = TMath::Sqrt(1 - cost*cost);
  kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
  double xsec=fModel->XSec(fInteraction, kPSTlctl);
  //std::cout << "f2(" << El << "," << cost << "): "<< xsec << std::endl;
  return -xsec;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2Xsec_dEldCosThetal4Min::Clone() const
{
  return new genie::utils::gsl::d2Xsec_dEldCosThetal4Min(fModel, fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal_fixed::d2Xsec_dEldCosThetal_fixed(const XSecAlgorithmI * m, const Interaction * interaction, int numvar, double valvar, double root) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(interaction),
fNumvar(numvar),
fValvar(valvar),
fRoot(root)
{

}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal_fixed::~d2Xsec_dEldCosThetal_fixed()
{
  
}   
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dEldCosThetal_fixed::DoEval(const double xin) const
{
  double El, cost;
  if (fNumvar == 0)
  {
    El = fValvar;
    cost = xin;
  }
  else
  {
    El = xin;
    cost = fValvar;
  }
    
  Kinematics * kinematics = fInteraction->KinePtr();
  double ml   = fInteraction->FSPrimLepton()->Mass();
  double Pl   = TMath::Sqrt(El*El - ml*ml);
  double sint = TMath::Sqrt(1 - cost*cost);
  kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
  double xsec=fModel->XSec(fInteraction, kPSTlctl);
  
  return xsec-fRoot;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2Xsec_dEldCosThetal_fixed::Clone() const
{
  return new genie::utils::gsl::d2Xsec_dEldCosThetal_fixed(fModel, fInteraction, fNumvar, fValvar, fRoot);
}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed::d2Xsec_dEldCosThetal4Min_fixed(const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(interaction)
{

}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed::~d2Xsec_dEldCosThetal4Min_fixed()
{
  
}   
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed::DoEval(const double xin) const
{
// inputs:
//    normalized El from 0 to 1
// outputs:
//   differential cross section
//
 
    
  Kinematics * kinematics = fInteraction->KinePtr();
  double ml   = fInteraction->FSPrimLepton()->Mass();
  double El     = xin;
  double Pl = TMath::Sqrt(El*El - ml*ml);
  kinematics->SetFSLeptonP4(0, 0, Pl, El);
  double xsec=fModel->XSec(fInteraction, kPSTlctl);
  //std::cout << "f2(" << El << "," << 1.0 << "): "<< xsec << std::endl;
  return -xsec;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed::Clone() const
{
  return new genie::utils::gsl::d2Xsec_dEldCosThetal4Min_fixed(fModel, fInteraction);
}
