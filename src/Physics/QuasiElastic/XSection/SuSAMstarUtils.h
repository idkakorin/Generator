//____________________________________________________________________________
/*!

\class    genie::SuSAMstarUtils

\brief    Contains auxiliary functions for SuSAM* model. \n

\ref      Physical Review D 97, 116006 (2018)
          Physical Review C 98, 024627 (2018)

\author   I. Kakorin <kakorin@jinr.ru>,                Joint Institute for Nuclear Research \n
          V. Naumov  <vnaumov@theor.jinr.ru>,          Joint Institute for Nuclear Research \n
          adapted from  fortran code provided by  \n
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

#ifndef _SUSAM_STAR_UTILS_H_
#define _SUSAM_STAR_UTILS_H_

#include <vector>

#include <Math/IFunction.h>
#include <TLorentzVector.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"

namespace genie {



class SuSAMstarUtils : public Algorithm {

public:
        SuSAMstarUtils();
        SuSAMstarUtils(string config);
        virtual ~SuSAMstarUtils();
        void SetInteraction(const Interaction * i);
        double GetEffectiveMassRatio(void) const;
        double GetFermiMomentum(void) const;
        double PhaseSpaceVolume(KinePhaseSpace_t ps) const;
        
        //! methods overloading the Algorithm() interface implementation
        //! to build the fragmentation function from configuration data
        void Configure(const Registry & config);
        void Configure(string config);

private:
        mutable SuSAMstarUtils * sm_utils; 
        void   LoadConfig (void);
        double EffectiveMassRatio(const Target & t) const;
        double FermiMomentumSuSAMstar(const Target & t) const;


       
        std::vector<double>   fSuSAMpFvsA;       ///<  values of Fermi momentum obtained in the global fit to electron scattering data for some
                                                 ///<  nuclei and presented in PRCC98(2018)024627 ,for other nuclei cubic interpolation is used
   
        map <int,double>      fSuSAMEffectsvsA;  ///< ratios of effective mass to nucleon one, obtained in the global fit 
                                                 ///< to electron scattering data for some nuclei and presented in PRCC98(2018)024627
        
        const Interaction *  fInteraction;
        
        // Some often used variables of class.
        // To not calculate them again and again and for speed increase
        // they are initialized at once for multiple use
        double  E_nu;       ///<  Neutrino energy (GeV)
        double  m_lep;      ///<  Mass of final charged lepton (GeV)
        double  mm_lep;     ///<  Squared mass of final charged lepton (GeV)
        double  P_Fermi;    ///<  Maximum value of Fermi momentum of target nucleon (GeV)
        double  eff_rat;    ///<  Effectve mass ration to ordinary mass
};



}       // genie namespace
#endif  // _SUSAM_STAR_UTILS_H_
