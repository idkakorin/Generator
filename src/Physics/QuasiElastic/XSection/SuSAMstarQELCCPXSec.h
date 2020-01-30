//____________________________________________________________________________
/*!

\class    genie::SuSAMstarQELCCPXSec

\brief    Computes neutrino-nucleus QELCC differential cross section
          by SuSAM* model
          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      Physical Review D 97, 116006 (2018)
          Physical Review C 98, 024627 (2018)

\author   I. Kakorin <kakorin@jinr.ru>,                Joint Institute for Nuclear Research \n
          V. Naumov  <vnaumov@theor.jinr.ru>,          Joint Institute for Nuclear Research \n
          adapted from fortran code provided by  \n
          I. Ruiz Simo <ruizsig@ugr.es>,               University of Granada \n
          V.L. Martinez-Consentino <victormc@ugr.es>,  University of Granada \n
          J.E. Amaro <amaro@ugr.es>,                   University of Granada \n
          E. Ruiz Arriola1 <earriola@ugr.es>,          University of Granada
      

\created  June 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SUSAMSTAR_QELCC_CROSS_SECTION_H_
#define _SUSAMSTAR_QELCC_CROSS_SECTION_H_

#include <vector>
#include <map>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarUtils.h"

namespace genie {

class QELFormFactorsModelI;
class XSecIntegratorI;

class SuSAMstarQELCCPXSec : public XSecAlgorithmI {

public:
  SuSAMstarQELCCPXSec();
  SuSAMstarQELCCPXSec(string config);
  virtual ~SuSAMstarQELCCPXSec();

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
  mutable QELFormFactors       fFormFactors;      ///<
  const QELFormFactorsModelI * fFormFactorsModel; ///<
  const XSecIntegratorI *      fXSecIntegrator;   ///<
  double                       fVud;             ///< |Vud|- magnitude of ud-element of CKM-matrix

  double                       fXSecScale;        ///< external xsec scaling factor
  
  double fa1, fa2, fa3, fb1, fb2, fb3;            ///< parameters of the first gaussian f1 for SuSAM* scaling function 
                                                  ///< f=a3*exp[-(psi*-a1)^2/(2*a2^2)]+b3*exp[-(psi*-b1)^2/(2*b2^2)]
                                                                  
  //Functions needed to calculate XSec:                           

  double d2sQES_dEldCosThetal(const Interaction * i) const;
  double dsQES_dQ2(const Interaction * i) const;

};
}       // genie namespace
#endif  // _SUSAMSTAR_QELCC_CROSS_SECTION_H_
