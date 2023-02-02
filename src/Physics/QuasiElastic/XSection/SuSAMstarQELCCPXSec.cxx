//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: I. Kakorin <kakorin@jinr.ru>,                Joint Institute for Nuclear Research
          V. Naumov  <vnaumov@theor.jinr.ru>,          Joint Institute for Nuclear Research
          adapted from  fortran code provided by 
          I. Ruiz Simo <ruizsig@ugr.es>,               University of Granada
          V.L. Martinez-Consentino <victormc@ugr.es>,  University of Granada
          J.E. Amaro <amaro@ugr.es>,                   University of Granada
          E. Ruiz Arriola1 <earriola@ugr.es>,          University of Granada
 
 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________
#include <TLorentzVector.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarQELCCPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
SuSAMstarQELCCPXSec::SuSAMstarQELCCPXSec() :
XSecAlgorithmI("genie::SuSAMstarQELCCPXSec")
{

}
//____________________________________________________________________________
SuSAMstarQELCCPXSec::SuSAMstarQELCCPXSec(string config) :
XSecAlgorithmI("genie::SuSAMstarQELCCPXSec", config)
{

}
//____________________________________________________________________________
SuSAMstarQELCCPXSec::~SuSAMstarQELCCPXSec()
{

}
//____________________________________________________________________________
double SuSAMstarQELCCPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{
  
  double xsec = 0.;
  // dimension of kine phase space
  std::string s = KinePhaseSpace::AsString(kps);
  int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
  
  if(! this -> ValidProcess (interaction) ) 
  {
     LOG("SuSAM*",pWARN) << "not a valid process"; 
     return 0.;
  }

  if(kpsdim == 1)
  {
     if(! this -> ValidKinematics (interaction) )
     {
        LOG("SuSAM*",pWARN) << "not valid kinematics"; 
        return 0.;
     }
     xsec = this->dsQES_dQ2(interaction);
  }

  if(kpsdim == 2) 
  {
    xsec = this->d2sQES_dEldCosThetal(interaction);
  }
  
  // The algorithm computes d^1xsec/dQ2, d^2xsec/dEldctl
  // Check whether variable tranformation is needed
  if ( kps != kPSQ2fE && kps != kPSElctl ) 
  {
     double J = 1.;
     if (kpsdim == 1)
       J = utils::kinematics::Jacobian(interaction, kPSQ2fE, kps);
     else if (kpsdim == 2)
       J = utils::kinematics::Jacobian(interaction, kPSElctl, kps);
     xsec *= J;
  }

  return fXSecScale*xsec;
  
}
//____________________________________________________________________________
double SuSAMstarQELCCPXSec::d2sQES_dEldCosThetal(const Interaction * interaction) const
{
  
  // Get kinematics and init-state parameters
  Kinematics *  kinematics = interaction -> KinePtr();
  sm_utils->SetInteraction(interaction);
  
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  
  
  const TLorentzVector leptonMom = kinematics->FSLeptonP4();
  
  // The calculation of double differential cross section dsigma/dE_mu*dcos(theta_mu)
  // according to Ref.1
  
  // one of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : +1;
  int N = target.N();
  int Z = target.Z();
  int A = target.A();
  int XN = (is_neutrino) ? N : Z;
  
  // mass of isoscalar nucleon equal to masses of ingoing and outgoing nucleons
  double mN = kNucleonMass;
  
  // neutrino energy in Lab system of reference
  double enu = init_state.ProbeE(kRfLab);
  // outgoing lepton energy
  double emu = leptonMom.E();
  // energy transfer
  double omega = enu - emu;
  
  
  //  mass of outgoing charged lepton
  double m_lep = interaction->FSPrimLepton()->Mass();
  double mm_lep = m_lep*m_lep;
  // magnitude of the 3-momentum of outgoing lepton
  double xkmu  = leptonMom.P();
  // cosine of angle of outgoing charged lepton
  double cost = leptonMom.CosTheta();
  // magnitude of the 3-momentum transfer
  double q2 = enu*enu + xkmu*xkmu - 2.*enu*xkmu*cost;
  double q = TMath::Sqrt(q2);
  // magnitude of the 4-momentum transfer
  double fourq2  = omega*omega - q2;
  
  if (fourq2 >= 0)
    return 0.;
  
  //kinematics->Setq2(fourq2);
  
  // kinematic factor
  double v0      = (enu+emu)*(enu+emu) - q2;
  
  // dimensionless factors
  double delta2  = -mm_lep/fourq2;
  double rho    = -fourq2/q2;
  double rhop   = q/(enu + emu);
  double xnu    = omega/q;
  
  // auxiliary factors
  double tant2   = -fourq2/v0;
  
  // the v-coefficients to which response functions are multiplied
  double vcc = 1. - tant2*delta2;                                              // Eq. (3)
  double vcl = xnu + tant2*delta2/rhop;                                        // Eq. (4)
  double vll = xnu*xnu + tant2*(1. + 2.*xnu/rhop + rho*delta2)*delta2;         // Eq. (5)
  double vt  = rho/2. + tant2 - tant2*(xnu + rho*rhop*delta2/2.)*delta2/rhop;  // Eq. (6)
  double vtp = tant2*(1. - xnu*rhop*delta2)/rhop;                              // Eq. (7)
  
  // bellow all calculated values correspond to star values in the paper
  double am = mN*sm_utils->GetEffectiveMassRatio(); 

  // adimensional variables
  double xkappa  = q/2./am;
  double xlambda = omega/2./am;
  double xkappa2 = xkappa*xkappa;
  double tau     = xkappa2 - xlambda*xlambda;
  double tau1    = tau + 1.;
  double xk2tau  = xkappa2/tau;
  double sqttau1 = TMath::Sqrt(tau*tau1); // tau*tau1 is always greater than 0 if Q2>0
  
  double fm      = sm_utils->GetFermiMomentum();
  double fmp     = fm*TMath::Power(2.*Z/A, 1./3);
  double fmn     = fm*TMath::Power(2.*N/A, 1./3);
  
  double etafini   = fmp/am;   // antineutrino case
  double etaffin   = fmn/am;   // antineutrino case
  if (is_neutrino)
  {
      etafini   = fmn/am;      // neutrino case
      etaffin   = fmp/am;      // neutrino case
  }
  
  double epsiffin  = TMath::Sqrt(1. + etaffin*etaffin);                 // 1. + etaffin*etaffin is always greater than 0
  double epsifini  = TMath::Sqrt(1. + etafini*etafini);                 // 1. + etafini*etafini is always greater than 0
  double sqt1t     = TMath::Sqrt(tau1/tau);                             // tau1/tau is always greater than 0 if Q2>0
  double e1        = xkappa*sqt1t - xlambda;                            
  double e2        = epsiffin - 2.*xlambda;                             
  double e0        = TMath::Max(e1, e2);                                // for non-Pauli blocking e0 = e1
  double xif       = epsifini - 1.;
  double xi0       = e0 - 1.;
  

  // scaling variable psi^2
  double psi2 = xi0/xif;   // Eq. (19)
  
  
  // scaling functions
  /*
  // well-known function of Fermi gas
  double scalfun = 0.;    // Eq. (18)
  if(psi2 <= 1.)
    scalfun = 3./4.*(1. - psi2);
  */ 
  
  // Phenomenological scaling function  central
  // calculation of scaling function by Eq. (37)
  if (psi2 < 0.)  psi2 = 0.;
  double psi = TMath::Sqrt(psi2);   // there is strong evidence that psi2>=0
  if(xlambda <= tau) 
     psi = -psi;
  double x = psi;
  if ( TMath::Abs((x - fa1)/fa2) > 5. && TMath::Abs((x - fb1)/fb2) > 5. ) return 0.;
  double f1 = fa3*TMath::Exp(-(x - fa1)*(x - fa1)/(2.*fa2*fa2));
  double f2 = fb3*TMath::Exp(-(x - fb1)*(x - fb1)/(2.*fb2*fb2));
  double f3 = 1 + fc3*TMath::Exp(-(x - fc1)/fc2)
  double scalfun = (f1 + f2)/f3; 
  
  // Linhard function
  double r0 = XN*xif/(am*TMath::Power(etafini, 3)*xkappa)*scalfun;      // Eq. (17) - the resulting nuclear response function 
  
  // non relativistic linhard function:
  //double r0 = XN/(2.d0*xkappa*fmn)*scalfun;

  // small functions Delta and Delta~
  double d1 = (1. + xlambda)*(1. + xlambda) - xk2tau*tau1;
  double d2 = xif*(1. + xlambda)*(1. + psi2);
  double d3 = xif*xif*(1. + psi2 + psi2*psi2)/3.;
  double delta = (d1 + d2 + d3)/xk2tau;                            

  /*
  // alternative definition of delta (non-Pauli blocking)
  double d1 = xkappa*TMath::Sqrt(tau1/tau) + xif*(1.-psi2)/3.;
  double delta = d1*xif*(1. - psi2)/xk2tau;                             // Eq. (24)
  */

  double d4 = tau/xkappa*(xlambda + 1.) - sqttau1;
  double d5 = tau/xkappa*xif*(1. + psi2)/2.;
  double deltat = (d4 + d5)/sqttau1;
  
  /*
  // alternative definition of tilde delta (non-Pauli blocking)
  double deltat = xif*(1. - psi2)/xkappa/sqt1t/2.;                      // Eq. (36)
  */
  
  
  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);
  
  // The form factors to use in the vector current are twice the vector ones
  double ge1 = fFormFactors.F1V() - tau*am*fFormFactors.xiF2V()/mN;     //  Eq. (22)
  double gm1 = fFormFactors.F1V() + am*fFormFactors.xiF2V()/mN;         //  Eq. (23)
  
  //axial form factor is not modified in the medium because in its definition in the axial current, Eq. (16), 
  // the nucleon mass does not appear explicit 
  double ga = fFormFactors.FA();
  double gp = 2*am*fFormFactors.Fp()/mN;                                //  Eq. (27)
  double gap = ga - tau*gp;

  // single nucleon invariant functions
  double w1vv = tau*gm1*gm1;
  double w2vv = (ge1*ge1 + tau*gm1*gm1)/tau1;
  double w1aa = tau1*ga*ga;
  double w2aa = ga*ga;
  double w4aa = gap*gap/(4.*tau);
  double w3va = gm1*ga;
  double w3 = w3va;
  double w4 = w4aa;

  // VV responses  
  double ucv = xk2tau*(tau1*w2vv - w1vv + w2vv*delta);                  //  Eq. (21), the first term of Eq. (20)
  double uclv = -xlambda/xkappa*ucv;                                    //  the first term of Eq. (28)
  double ulv = (xlambda/xkappa)*(xlambda/xkappa)*ucv;                   //  the first term of Eq. (29)
  
  // AA responses conserved part
  double ucac = xk2tau*(tau1*w2aa - w1aa + w2aa*delta);                 // Eq. (25), the second term of Eq. (20) (see (276) from Notes)
  // double ucac = xk2tau*w2aa*delta;                                   // Eq. (25) for non-Pauli blocking case
  double uclac = -xlambda/xkappa*ucac;                                  // the second term of Eq. (28)
  double ulac = (xlambda/xkappa)*(xlambda/xkappa)*ucac;                 // the second term of Eq. (29)
  
  // AA responses non conserved paer
  double ucanc = 4.*xlambda*xlambda*w4;                                 // Eq. (26), the third term of Eq. (20)
  double uclanc = -4.*xlambda*xkappa*w4;                                // Eq. (30), the third term of Eq. (28)
  double ulanc = 4.*xkappa*xkappa*w4;                                   // Eq. (31), the third term of Eq. (29)
  double uca = ucac + ucanc;                                            // Eq. (20), the sum of the last two terms: Eq. (25) and Eq. (26)
  double ucla = uclac + uclanc;                                         // Eq. (28), the sum of the last two terms: Eq. (30)
  double ula = ulac + ulanc;                                            // Eq. (29), the sum of the last two terms: Eq. (31)
  //double w2vv = f11*f11 + tau*f21*f21;
  double utv = 2.*w1vv + w2vv*delta;                                    // Eq. (33), the first  term of Eq. (32)
  double uta = 2.*w1aa + w2aa*delta;                                    // Eq. (34), the second term of Eq. (32)
  double utp = 2.*w3*sqttau1*(1. + deltat);                             // Eq. (35)


  // total responses multiplied by R0
  double rccv = ucv;                                                    // see Eq. (8), (20)
  double rcca = uca;                                                    // see Eq. (8), (25), (26)
  double rcc = rccv + rcca;                                             // Eq. (8)
                                                                        
  double rclv = uclv;                                                   // see Eq. (9), (28)
  double rcla = ucla;                                                   // see Eq. (9), (28), (30)
  double rcl = rclv + rcla;                                             // Eq. (9)
                                                                        
  double rllv = ulv;                                                    // see Eq. (10), (29)
  double rlla = ula;                                                    // see Eq. (10), (29), (31)
  double rll = rllv + rlla;                                             // Eq. (10)
                                                                        
  double rtv  = utv;                                                    // see Eq. (11), (32), (33)
  double rta  = uta;                                                    // see Eq. (11), (32), (34)
  double rt  = rtv + rta;                                               // Eq. (11)
                                                                          
  double rtp  = utp;                                                    // Eq. (12), see Eq. (35)
  
  double f = (rcc*vcc + 2.*rcl*vcl + rll*vll + rt*vt + 2.*sign*rtp*vtp)*r0;  // See Eq. (1)
  double se0 = TMath::Power(kGF*fVud, 2)*xkmu*v0/4./kPi/enu;            // Eq. (2)
  double dsigde = se0*f;                                                // Eq. (1)
  

  return dsigde;
  
}
//____________________________________________________________________________
double SuSAMstarQELCCPXSec::dsQES_dQ2(const Interaction * interaction) const
{
  // Get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double E  = init_state.ProbeE(kRfHitNucRest);
  double E2 = TMath::Power(E,2);
  double ml = interaction->FSPrimLepton()->Mass();
  double M  = target.HitNucMass();
  double q2 = kinematics.q2();

  // One of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : 1;

  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);    

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  double Fp    = fFormFactors.Fp();


  // Calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);
  double M4      = TMath::Power(M2,    2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double Gfactor = M2*TMath::Power(kGF*fVud, 2)*(kMw2/(kMw2-q2))*(kMw2/(kMw2-q2)) / (8*kPi*E2);
  double s_u     = 4*E*M + q2 - ml2;
  double q2_M2   = q2/M2;

  // Compute free nucleon differential cross section
  double A = (0.25*(ml2-q2)/M2) * (
      (4-q2_M2)*FA2 - (4+q2_M2)*F1V2 - q2_M2*xiF2V2*(1+0.25*q2_M2)
              -4*q2_M2*F1V*xiF2V - (ml2/M2)*( 
               (F1V2+xiF2V2+2*F1V*xiF2V)+(FA2+4*Fp2+4*FA*Fp)+(q2_M2-4)*Fp2));
  double B = -1 * q2_M2 * FA*(F1V+xiF2V);
  double C = 0.25*(FA2 + F1V2 - 0.25*q2_M2*xiF2V2);

  double xsec = Gfactor * (A + sign*B*s_u/M2 + C*s_u*s_u/M4);
  
  // Apply given scaling factor
  xsec *= fXSecScale;
  
  // Deuterium and tritium is a special case
  if (target.A()>1 && target.A()<4)
  {
    double Q2 = -q2;
    double fQES_Pauli = 1.-0.529*TMath::Exp((Q2*(228.-531.*Q2)-48.)*Q2);
    xsec *= fQES_Pauli;
  }

  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 

  xsec *= NNucl; // nuclear xsec
  
  // Apply radiative correction to the cross section for IBD processes
  // Refs: 
  // 1) I.S. Towner, Phys. Rev. C 58 (1998) 1288;
  // 2) J.F. Beacom, S.J. Parke, Phys. Rev. D 64 (2001) 091302;
  // 3) A. Kurylov, M.J. Ramsey-Musolf, P. Vogel, Phys. Rev. C 65 (2002) 055501;
  // 4) A. Kurylov, M.J. Ramsey-Musolf, P. Vogel, Phys. Rev. C 67 (2003) 035502.
  double rc = 1.;
  if ( (target.IsProton() && pdg::IsAntiNuE(init_state.ProbePdg())) || (target.IsNeutron() && pdg::IsNuE(init_state.ProbePdg()) ))
  {
    const double mp  = kProtonMass;
    const double mp2 = kProtonMass2;
    const double mn2 = kNeutronMass2;
    const double Ee  = E + ( (q2 - mn2 + mp2) / 2. / mp ); 
    assert(Ee > 0.); // must be non-zero and positive
    rc  = 6. + (1.5 * TMath::Log(kProtonMass / 2. / Ee));
    rc += 1.2 * TMath::Power((kElectronMass / Ee), 1.5);
    rc *= kAem / kPi;
    rc += 1.;
  }

  xsec *= rc;
  return xsec;
}
//____________________________________________________________________________
double SuSAMstarQELCCPXSec::Integral(const Interaction * in) const
{
  return fXSecIntegrator->Integrate(this,in);

}
//____________________________________________________________________________
bool SuSAMstarQELCCPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP     = pdg::IsProton(nuc);
  bool isN     = pdg::IsNeutron(nuc);
  bool isnumu  = pdg::IsNuMu(nu);
  bool isnumub = pdg::IsAntiNuMu(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnumub) || (isN&&isnumu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void SuSAMstarQELCCPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarQELCCPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r("SuSAMstarQELCCPXSec_specific", false ) ;
  r.Set("sm_utils_algo", RgAlg("genie::SuSAMstarUtils","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarQELCCPXSec::LoadConfig(void)
{
  
  GetParam( "CKM-Vud", fVud ) ;

  // Cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

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

   // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
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


