//____________________________________________________________________________
/*!

\class    genie::SuSAMstarEMPXSec

\brief    Computes differential cross section of EM scatering electron neutrino on a heavy nucleus (mass greater than 3) within SuSAM* model.
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref
             1. Amaro J. E., Ruiz Arriola E., Ruiz Simo I., Phys.Rev.C 92 (2015) 054607
             2. Amaro J. E., Ruiz Arriola E., Ruiz Simo I., Phys.Rev.C 98 (2018) 024627
             3. de Forest T. Jr., Walecka J.D., Adv.Phys. 15(57) (1966) 1-109
      

\created  January 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SUSAMSTAR_EM_CROSS_SECTION_H_
#define _SUSAMSTAR_EM_CROSS_SECTION_H_

#include <vector>
#include <map>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarUtils.h"

namespace genie {

class ELFormFactorsModelI;
class XSecIntegratorI;

class SuSAMstarEMPXSec : public XSecAlgorithmI {

public:
  SuSAMstarEMPXSec();
  SuSAMstarEMPXSec(string config);
  virtual ~SuSAMstarEMPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  mutable SuSAMstarUtils * sm_utils;
  
  void   LoadConfig (void);
  mutable ELFormFactors        fFormFactors;      ///<
  const ELFormFactorsModelI  * fFormFactorsModel; ///<
  const XSecIntegratorI *      fXSecIntegrator;   ///<

  double                       fXSecScale;        ///< external xsec scaling factor
  
  double fa1, fa2, fa3, fb1, fb2, fb3, fc1, fc2, fc3;            ///< parameters of the first gaussian f1 for SuSAM* scaling function 
                                                                 ///< f = (a3*exp[-(psi*-a1)^2/(2*a2^2)]+b3*exp[-(psi*-b1)^2/(2*b2^2)])/(1+c3*exp[-(psi-c1)/c2])

};
}       // genie namespace
#endif  // _SUSAMSTAR_EM_CROSS_SECTION_H_
