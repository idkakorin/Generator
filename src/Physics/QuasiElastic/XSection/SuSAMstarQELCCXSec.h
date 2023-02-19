//____________________________________________________________________________
/*!

\class    genie::SuSAMstarQELCCXSec

\brief    Computes the Quasi Elastic (QEL) cross section by Smith Moniz model. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          
\created  June 2019

\cpright  Copyright (c) 2003-2019, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SUSAM_STAR_QEL_XSEC_H_
#define _SUSAM_STAR_QEL_XSEC_H_

#include <Math/IFunction.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarUtils.h"

namespace genie {



class SuSAMstarQELCCXSec : public XSecIntegratorI {

public:
  SuSAMstarQELCCXSec();
  SuSAMstarQELCCXSec(string config);
  virtual ~SuSAMstarQELCCXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);
  
protected:
  string fGSLIntgType2D;  ///< name of GSL 2D numerical integrator 
  double fGSLRelTol2D;    ///< required relative tolerance (error) for 2D integrator

private:
  void LoadConfig (void);
  double fEmax;
  double ftheta_min;
mutable double fXSecAtEmax;
  
};

//_____________________________________________________________________________________
// 
// GSL wrappers 
//
//_____________________________________________________________________________________

 namespace utils {
  namespace gsl   {

   class d2Xsec_dEldCosThetal: public ROOT::Math::IBaseFunctionMultiDim
   {
    public:
      d2Xsec_dEldCosThetal(const XSecAlgorithmI * m, const Interaction * i, double cos_theta_max);
     ~d2Xsec_dEldCosThetal();
      // ROOT::Math::IBaseFunctionMultiDim interface
      unsigned int                        NDim   (void)               const;
      double                              DoEval (const double * xin) const;
      ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
    
    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
      double fcos_theta_max;
   };
   

  } // gsl   namespace
 } // utils namespace
 
}       // genie namespace
#endif  // _SUSAM_STAR_QEL_XSEC_H_
