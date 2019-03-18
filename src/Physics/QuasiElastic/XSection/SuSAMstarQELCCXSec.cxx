//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          
 For the class documentation see the corresponding header file.

 
*/
//____________________________________________________________________________
#include <limits>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarQELCCXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using std::ostringstream;

//____________________________________________________________________________
SuSAMstarQELCCXSec::SuSAMstarQELCCXSec() :
XSecIntegratorI("genie::SuSAMstarQELCCXSec")
{

}
//____________________________________________________________________________
SuSAMstarQELCCXSec::SuSAMstarQELCCXSec(string config) :
XSecIntegratorI("genie::SuSAMstarQELCCXSec", config)
{

}
//____________________________________________________________________________
SuSAMstarQELCCXSec::~SuSAMstarQELCCXSec()
{

}
//____________________________________________________________________________
double SuSAMstarQELCCXSec::Integrate(
                  const XSecAlgorithmI * model, const Interaction * in) const
{
  LOG("SuSAM*QELXSec",pDEBUG) << "Beginning integrate";
  if(! model->ValidProcess(in)) return 0.;
  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold())
  {
    LOG("SuSAM*QELXSec", pDEBUG)  << "*** Below energy threshold";
    return 0;
  }

  const InitialState & init_state = in -> InitState();
  if (init_state.ProbeE(kRfLab)>fEmax)
  {
     if (fXSecAtEmax<0)
      in->InitStatePtr()->SetProbeE(fEmax);
    else
      return fXSecAtEmax;
  }

    
  
  const Target & target = init_state.Tgt();
  if (target.A()<3)
  {
     Range1D_t rQ2 = kps.Limits(kKVQ2);
     if(rQ2.min<0 || rQ2.max<0) return 0;
     Interaction * interaction = new Interaction(*in);
     interaction->SetBit(kISkipProcessChk);
     interaction->SetBit(kISkipKinematicChk);
     ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::dXSec_dQ2_E(model, interaction);
     ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
     double abstol = 0; //We mostly care about relative tolerance
     ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxSizeOfSubintervals, fGSLRule);
     double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);
     delete func;
     delete interaction;
     return xsec;
  }
  else
  {   
     Interaction * interaction = new Interaction(*in);
     interaction->SetBit(kISkipProcessChk);
     interaction->SetBit(kISkipKinematicChk);
     double xsec = 0;
     
     double kine_min[2] = {0, 0}; 
     double kine_max[2] = {1, 1}; 
     double abstol = 0; //We mostly care about relative tolerance.
     
     ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType2D);
     ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d2Xsec_dEldCosThetal(model, interaction);
     ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol2D, fGSLMaxEval);
     xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
     delete func;
     delete interaction;
     if (init_state.ProbeE(kRfLab)>=fEmax)
       fXSecAtEmax = xsec;
     return xsec;
  }
  
}
//____________________________________________________________________________
void SuSAMstarQELCCXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarQELCCXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarQELCCXSec::LoadConfig(void)
{
  
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("gauss") );
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-3 );
  int max_size_of_subintervals;
  GetParamDef( "gsl-max-size-of-subintervals", max_size_of_subintervals, 40000);
  fGSLMaxSizeOfSubintervals = (unsigned int) max_size_of_subintervals;
  int rule;
  GetParamDef( "gsl-rule", rule, 3);
  fGSLRule = (unsigned int) rule;
  if (fGSLRule>6) fGSLRule=3;
  GetParamDef( "gsl-integration-type-2D", fGSLIntgType2D, string("adaptive") );
  GetParamDef( "gsl-relative-tolerance-2D", fGSLRelTol2D, 1e-7);
  GetParamDef( "gsl-max-eval", fGSLMaxEval, 1000000000);
  GetParamDef( "Emax", fEmax, 45.0);
  fXSecAtEmax = -1.0;
  
}


//_____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal::d2Xsec_dEldCosThetal(const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(interaction)
{

}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dEldCosThetal::~d2Xsec_dEldCosThetal()
{
  
}   
//____________________________________________________________________________
unsigned int genie::utils::gsl::d2Xsec_dEldCosThetal::NDim(void) const
{
  return 2;
}
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dEldCosThetal::DoEval(const double * xin) const
{
// inputs:
//    normalized El from 0 to 1
//    normalized cosine of lepton angle from 0 to 1 which transformed into  -1,1
// outputs:
//   differential cross section [10^-38 cm^2]
//
 
    
  Kinematics * kinematics = fInteraction->KinePtr();
  double Enu = fInteraction->InitState().ProbeE(kRfLab);
  double ml   = fInteraction->FSPrimLepton()->Mass();
  if (Enu < ml) return 0.0;
  
  double El     = (Enu - ml)*xin[0] + ml;
  double Pl = TMath::Sqrt(El*El - ml*ml);
  double cost   = 2*xin[1] - 1.0;
  double sint = TMath::Sqrt(1 - cost*cost);
  kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
  double xsec=fModel->XSec(fInteraction, kPSTlctl);
  double J  = 2*(Enu - ml); // Jacobian for transformation 
  xsec *= J;
  
  return xsec/(1E-38 * units::cm2);
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2Xsec_dEldCosThetal::Clone() const
{
  return new genie::utils::gsl::d2Xsec_dEldCosThetal(fModel, fInteraction);
}

