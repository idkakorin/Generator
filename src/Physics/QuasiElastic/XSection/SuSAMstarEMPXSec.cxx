//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________
#include <TLorentzVector.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarEMPXSec.h"


using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
SuSAMstarEMPXSec::SuSAMstarEMPXSec() :
XSecAlgorithmI("genie::SuSAMstarEMPXSec")
{

}
//____________________________________________________________________________
SuSAMstarEMPXSec::SuSAMstarEMPXSec(string config) :
XSecAlgorithmI("genie::SuSAMstarEMPXSec", config)
{

}
//____________________________________________________________________________
SuSAMstarEMPXSec::~SuSAMstarEMPXSec()
{

}
//____________________________________________________________________________
double SuSAMstarEMPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{
  
  if(! this -> ValidProcess (interaction) ) 
  {
     LOG("SuSAM*",pWARN) << "not a valid process"; 
     return 0.;
  }
  
  // Get kinematics and init-state parameters
  Kinematics *  kinematics = interaction -> KinePtr();
  sm_utils->SetInteraction(interaction);
  
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  int N = target.N();
  int Z = target.Z();
  int A = target.A();
  
  
  const TLorentzVector leptonMom = kinematics->FSLeptonP4();
  
  // The calculation of double differential cross section dsigma/dE_e*dcos(theta_e)
  // according to Ref.1
    
  // mass of isoscalar nucleon equal to masses of ingoing and outgoing nucleons
  double mN = kNucleonMass;
  
  // bellow all calculated values correspond to star values in the paper
  double meff = mN*sm_utils->GetEffectiveMassRatio();
  // Fermi momentum
  double fm      = sm_utils->GetFermiMomentum();
  
  // neutrino energy in Lab system of reference
  double Ei = init_state.ProbeE(kRfLab);
  // outgoing lepton energy
  double Ef = leptonMom.E();
  // energy transfer
  double omega = Ei - Ef;
  
  
  //  mass of outgoing charged lepton
  double m_lep = interaction->FSPrimLepton()->Mass();
  double mm_lep = m_lep*m_lep;
  // square of speed of initial electron
  double beta2 = 1 - mm_lep/Ei/Ei;
  // magnitude of the 3-momentum of outgoing lepton
  double Pl  = leptonMom.P();
  // cosine of angle of outgoing charged lepton
  double cost = leptonMom.CosTheta();
  // magnitude of the 3-momentum transfer
  double q2 = Ei*Ei + Pl*Pl - 2.*Ei*Pl*cost;
  double q = TMath::Sqrt(q2);
  // magnitude of the 4-momentum transfer
  double Q2  = q2 - omega*omega;
  
  if (Q2 <= 0)
    return 0.;
    
  kinematics->SetQ2(Q2);
    
  double sin_halftheta2  = (1 - cost)/2;
  double cos_halftheta2  = (1 + cost)/2;
  
  // Values defined underneath Eq.(5)
  double vl = Q2*Q2/q2/q2;
  double vt = sin_halftheta2/cos_halftheta2 + Q2/q2/2;
  
  // adimensional variables
  double kappa  = q/meff/2;
  double lambda = omega/meff/2;
  double kappa2 = kappa*kappa;
  double tau    = kappa2 - lambda*lambda;
  double tau1   = tau + 1;
  double k2tau  = kappa2/tau;
  
  
  // Calculate the EM form factors and response functions for proton and neutron
  fFormFactors.Calculate(interaction);
  double Ge, Gm, fmN;
  int xN;
  double R = 0;
  for (int t = 0; t < 2; t++)
  {
    if (t == 0)
    {
      fmN     = fm*TMath::Power(2.*Z/A, 1./3);
      xN = N;
      double M0 = PDGLibrary::Instance()->Find(kPdgProton)->Mass();
      double tau0 = -Q2/4/M0/M0;
      double T0 = 1/(1 - tau0);
      double F1p = T0*(fFormFactors.Gep() - fFormFactors.Gmp()*tau0);
      double F2p = T0*(fFormFactors.Gmp() - fFormFactors.Gep());
      Ge = F1p - meff*F2p*tau/mN;                                     //  Eq. (14)
      Gm = F1p + meff*F2p/mN;                                         //  Eq. (15)
    }
    else
    {
      fmN     = fm*TMath::Power(2.*N/A, 1./3);
      xN = Z;
      double M0 = PDGLibrary::Instance()->Find(kPdgNeutron)->Mass();
      double tau0 = -Q2/4/M0/M0;
      double T0 = 1/(1 - tau0);
      double F1n = T0*(fFormFactors.Gen() - fFormFactors.Gmn()*tau0);
      double F2n = T0*(fFormFactors.Gmn() - fFormFactors.Gen());
      Ge = F1n - meff*F2n*tau/mN;                                     //  Eq. (14)
      Gm = F1n + meff*F2n/mN;                                         //  Eq. (15)
    }
  
  
    double etaf    = fmN/meff;
    double epsif   = TMath::Sqrt(1 + etaf*etaf);                   // 1. + etaf*etaf is always greater than 0
    double sqt1t   = TMath::Sqrt(tau1/tau);                        // tau1/tau is always greater than 0 if Q2>0
    double xif     = epsif - 1;
    double e1      = kappa*sqt1t - lambda;
    double e2      = epsif - 2*lambda;
    double e0      = TMath::Max(e1, e2);                           // for non-Pauli blocking e0 = e1
    double xi0     = e0 - 1;
    
    
    // scaling variable psi^2
    double psi2 = xi0/xif;   // Eqs. (16) and (17)
    
    // small functions Delta
    double d1    = (1 + lambda)*(1 + lambda) - k2tau*tau1;
    double d2    = xif*(1 + lambda)*(1 + psi2);
    double d3    = xif*xif*(1. + psi2 + psi2*psi2)/3.;
    double delta = (d1 + d2 + d3)/k2tau;
    
    // single nucleon invariant functions
    double w1vv = tau*Gm*Gm;
    double w2vv = (Ge*Ge + tau*Gm*Gm)/tau1;
    double Ul = k2tau*(tau1*w2vv - w1vv + w2vv*delta);                  // Eq. (11)
    double Ut = 2*w1vv + w2vv*delta;                                    // Eq. (12)
    

    // scaling functions
    /*
    // well-known function of Fermi gas
    double scalfun = 0.;
    if(psi2 <= 1.)
      scalfun = 3./4.*(1. - psi2);
    */ 
  
    // Phenomenological scaling function  central
    // calculation of scaling function by Eq. (23) of Ref.2
    if (psi2 < 0.)  psi2 = 0.;
    double psi = TMath::Sqrt(psi2);   // there is strong evidence that psi2>=0
    if(lambda <= tau) 
       psi = -psi;
    double x = psi;
    if ( TMath::Abs((x - fa1)/fa2) > 5. && TMath::Abs((x - fb1)/fb2) > 5. ) return 0.;
    double f1 = fa3*TMath::Exp(-(x - fa1)*(x - fa1)/(2.*fa2*fa2));
    double f2 = fb3*TMath::Exp(-(x - fb1)*(x - fb1)/(2.*fb2*fb2));
    double f3 = 1 + fc3*TMath::Exp(-(x - fc1)/fc2);
    double scalfun = (f1 + f2)/f3; 
    
    // Eqs. (8)-(10) - the resulting nuclear response function 
    double r0 = xN*xif/(meff*TMath::Power(etaf, 3)*kappa)*scalfun;
    R += (Ul*vl + Ut*vt)*r0;                                            

  }
  
 
//  Eq. (10a) of Ref. 3
  double sMott           = kAem2/4/Ei/Ei/sin_halftheta2/sin_halftheta2*(1 - beta2*sin_halftheta2);
// Eq. (3.65) of Ref. 4, the same as Eq. (10a) of Ref. 3 in the ultrarelativistic approach
//  double sMott = kAem2*cos_halftheta2/4./Ei/Ei/sin_halftheta2/sin_halftheta2; 

  double xsec            = 2*kPi*sMott*R;                                                                // Eq. (5)
  
  // The algorithm computes d^2xsec/dEldctl
  // Check whether variable tranformation is needed
  if (kps != kPSElctl ) 
  {
     double J = utils::kinematics::Jacobian(interaction, kPSElctl, kps);
     xsec *= J;
  }

  return fXSecScale*xsec;
  
}
//____________________________________________________________________________
double SuSAMstarEMPXSec::Integral(const Interaction * in) const
{
  return fXSecIntegrator->Integrate(this,in);

}
//____________________________________________________________________________
bool SuSAMstarEMPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsEM()) return false;
  
  const Target & target = init_state.Tgt();
  if (target.A() < 4) return false;

  int  nuc = target.HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP         = pdg::IsProton(nuc);
  bool isN         = pdg::IsNeutron(nuc);
  bool isel        = pdg::IsElectron(nu);
  bool ispos       = pdg::IsPositron(nu);

  bool prcok = (isel || ispos) && (isP || isN);
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void SuSAMstarEMPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarEMPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r("SuSAMstarEMPXSec_specific", false ) ;
  r.Set("sm_utils_algo", RgAlg("genie::SuSAMstarUtils","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarEMPXSec::LoadConfig(void)
{
  
  // Cross section scaling factor
  GetParam( "EM-XSecScale", fXSecScale ) ;

  // parameters of SuSAM* scaling function
  // Cross section scaling factor
  GetParam( "SuSAM-a1", fa1 );
  GetParam( "SuSAM-a2", fa2 );
  GetParam( "SuSAM-a3", fa3 );
  GetParam( "SuSAM-b1", fb1 );
  GetParam( "SuSAM-b2", fb2 );
  GetParam( "SuSAM-b3", fb3 );
  GetParam( "SuSAM-b1", fc1 );
  GetParam( "SuSAM-b2", fc2 );
  
  // load EM form factors model
  fFormFactorsModel = dynamic_cast<const ELFormFactorsModelI *> (this->SubAlg("ElasticFormFactorsModel"));
  assert(fFormFactorsModel);
  
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  sm_utils = const_cast<genie::SuSAMstarUtils *>(
               dynamic_cast<const genie::SuSAMstarUtils *>(
                 this -> SubAlg( "sm_utils_algo" ) ) ) ;

}


