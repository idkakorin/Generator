//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorSuSAMstar

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events for SuSAM* model.
          Is a concrete implementation of the EventRecordVisitorI interface.


\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          

\created  June, 2019

\cpright  Copyright (c) 2003-2019, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _QEL_EVENT_GENERATORSUSAMSTAR_H_
#define _QEL_EVENT_GENERATORSUSAMSTAR_H_

#include <Math/IFunction.h>
#include <Math/RootFinder.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarUtils.h"
#include "Physics/QuasiElastic/EventGen/QELKinematicsGenerator.h"
#include "Physics/QuasiElastic/EventGen/QELKinematicsGenerator.h"
#include "Physics/QuasiElastic/EventGen/QELPrimaryLeptonGenerator.h"
#include "Physics/QuasiElastic/EventGen/QELHadronicSystemGenerator.h"



namespace genie {

class QELEventGeneratorSuSAMstar: public KineGeneratorWithCache {

public :
  QELEventGeneratorSuSAMstar();
  QELEventGeneratorSuSAMstar(string config);
 ~QELEventGeneratorSuSAMstar();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction * in) const;
  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant
  void MaxsSubAreas(GHepRecord * event_rec, vector<vector<double>> & Coords, vector<double> & maxs) const;
  void ComputeMaxsSubAreas(const Interaction * in, vector<vector<double>> & Coords, vector<double> & maxs) const;
  void AccessCacheBranches(const Interaction * interaction, vector<CacheBranchFx*> & branches) const;

  
  mutable KinePhaseSpace_t fkps;
  
  mutable bool fIsHeavyNucleus;             ///< heavy nucleus is nucleus that heavier than hydrogen and deuterium
  mutable SuSAMstarUtils * sm_utils;
  mutable QELKinematicsGenerator * qel_kin_gen;
  mutable QELPrimaryLeptonGenerator * qel_primlep_gen;
  mutable QELHadronicSystemGenerator * qel_hadrsys_gen;
  int Nx;             ///<  Number of x-knots to search maximum on the grid
  int Ny;             ///<  Number of y-knots to search maximum on the grid
 



}; // class definition

//_____________________________________________________________________________________
// 
// GSL wrappers 
//
//_____________________________________________________________________________________

namespace utils {
  namespace gsl   {

   class d2Xsec_dEldCosThetal4Min: public ROOT::Math::IBaseFunctionMultiDim
   {
    public:
      d2Xsec_dEldCosThetal4Min(const XSecAlgorithmI * m, const Interaction * i);
     ~d2Xsec_dEldCosThetal4Min();
      // ROOT::Math::IBaseFunctionMultiDim interface
      unsigned int                        NDim   (void)               const;
      double                              DoEval (const double * xin) const;
      ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
    
    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
   };
   
   class d2Xsec_dEldCosThetal_fixed: public ROOT::Math::IBaseFunctionOneDim
   {
    public:
      d2Xsec_dEldCosThetal_fixed(const XSecAlgorithmI * m, const Interaction * i, int numvar, double valvar, double root);
     ~d2Xsec_dEldCosThetal_fixed();
      // ROOT::Math::IBaseFunctionMultiDim interface
      double                            DoEval (const double xin) const;
      ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;
    
    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
      const int fNumvar;
      const double fValvar;
      const double fRoot;
   };
   
   class d2Xsec_dEldCosThetal4Min_fixed: public ROOT::Math::IBaseFunctionOneDim
   {
    public:
      d2Xsec_dEldCosThetal4Min_fixed(const XSecAlgorithmI * m, const Interaction * i);
     ~d2Xsec_dEldCosThetal4Min_fixed();
      // ROOT::Math::IBaseFunctionMultiDim interface
      double                            DoEval (const double xin) const;
      ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;
    
    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
   };
   

  } // gsl   namespace
 } // utils namespace  
 
}       // genie namespace

#endif // _QEL_EVENT_GENERATORSUSAMSTAR_H_
