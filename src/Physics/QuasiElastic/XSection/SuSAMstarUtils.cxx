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
#include <TMath.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/SuSAMstarUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
//____________________________________________________________________________
SuSAMstarUtils::SuSAMstarUtils() :
Algorithm("genie::SuSAMstarUtils")
{

}
//____________________________________________________________________________
SuSAMstarUtils::SuSAMstarUtils(string config) :
Algorithm("genie::SuSAMstarUtils", config)
{

}
//____________________________________________________________________________
SuSAMstarUtils::~SuSAMstarUtils()
{

}
//____________________________________________________________________________
void SuSAMstarUtils::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarUtils::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SuSAMstarUtils::LoadConfig(void)
{
  
  // hardcoded values of Fermi momenta in MeV (for A from 1 to 250), 
  // which are obtained in the global fit to electron scattering 
  // data for some nuclei, for others the interpolated values are used. 
  double pFvalues[] = {
         65.08351, 81.00000, 143.0835, 159.0000, 162.3472, 164.0000, 167.1633, 170.8081, 176.9159,
         187.9457, 201.4385, 212.0000, 218.2456, 222.6905, 225.7902, 228.0000, 229.7475, 231.1788, 
         232.3275, 233.2269, 233.9105, 234.4115, 234.7636, 235.0000, 233.8940, 231.5697, 230.0872, 
         229.7955, 229.6502, 229.6207, 229.6766, 229.7875, 229.9230, 230.0525, 230.1457, 230.1720, 
         230.1010, 229.9023, 229.5455, 229.0000, 227.9275, 226.1855, 224.0377, 221.7477, 219.5791, 
         217.7956, 216.6608, 216.4383, 217.3107, 219.1054, 221.5526, 224.3830, 227.3268, 230.1147, 
         232.4769, 234.1440, 233.6822, 231.4719, 229.9647, 229.6296, 229.3011, 228.9793, 228.6642, 
         228.3559, 228.0544, 227.7598, 227.4721, 227.1913, 226.9175, 226.6507, 226.3910, 226.1384, 
         225.8930, 225.6548, 225.4238, 225.2002, 224.9838, 224.7749, 224.5733, 224.3793, 224.1927, 
         224.0137, 223.8423, 223.6785, 223.5224, 223.3740, 223.2334, 223.1006, 222.9757, 222.8629, 
         222.7659, 222.6837, 222.6155, 222.5603, 222.5172, 222.4852, 222.4635, 222.4511, 222.4471, 
         222.4506, 222.4606, 222.4762, 222.4965, 222.5206, 222.5475, 222.5763, 222.6061, 222.6360, 
         222.6651, 222.6923, 222.7169, 222.7378, 222.7541, 222.7650, 222.7695, 222.7666, 222.7555, 
         222.7352, 222.7049, 222.6700, 222.6367, 222.6050, 222.5748, 222.5459, 222.5182, 222.4917, 
         222.4662, 222.4415, 222.4178, 222.3947, 222.3722, 222.3502, 222.3286, 222.3073, 222.2861, 
         222.2651, 222.2440, 222.2227, 222.2013, 222.1794, 222.1572, 222.1343, 222.1108, 222.0866, 
         222.0615, 222.0353, 222.0081, 221.9798, 221.9501, 221.9190, 221.8864, 221.8522, 221.8163, 
         221.7786, 221.7389, 221.6972, 221.6534, 221.6073, 221.5589, 221.5080, 221.4546, 221.3985, 
         221.3396, 221.2779, 221.2132, 221.1454, 221.0744, 221.0001, 220.9224, 220.8412, 220.7564, 
         220.6679, 220.5755, 220.4793, 220.3790, 220.2745, 220.1658, 220.0528, 219.9353, 219.8132, 
         219.6865, 219.0319, 217.6522, 216.0369, 214.6752, 214.0564, 214.0855, 214.3132, 214.6973, 
         215.1956, 215.7659, 216.3661, 216.9540, 217.4874, 217.9242, 218.2221, 218.3391, 218.2243, 
         217.8923, 217.3988, 216.7997, 216.1506, 215.5073, 214.9256, 214.4611, 214.1698, 214.1072, 
         214.3293, 214.7301, 215.1657, 215.6342, 216.1338, 216.6625, 217.2187, 217.8005, 218.4060, 
         219.0334, 219.6809, 220.3467, 221.0290, 221.7258, 222.4355, 223.1561, 223.8859, 224.6229, 
         225.3655, 226.1117, 226.8597, 227.6077, 228.3539, 229.0964, 229.8335, 230.5632, 231.2838, 
         231.9935, 232.6904, 233.3726, 234.0384, 234.6859, 235.3134, 235.9189, 236.5006, 237.0568, 
         237.5856, 238.0851, 238.5536, 238.9892, 239.3901, 239.7544, 240.0804
      };
  
  fSuSAMpFvsA = std::vector<double>(pFvalues, pFvalues + sizeof(pFvalues) / sizeof(pFvalues[0]) );
  
  
  // hardcoded ratios of effective mass to nucleon one, which are 
  // obtained in the global fit to electron scattering for some nuclei
  fSuSAMEffectsvsA[2]   = 0.99;
  fSuSAMEffectsvsA[3]   = 0.96;
  fSuSAMEffectsvsA[4]   = 0.87;
  fSuSAMEffectsvsA[6]   = 0.80;
  fSuSAMEffectsvsA[9]   = 0.80;
  fSuSAMEffectsvsA[12]  = 0.83;
  fSuSAMEffectsvsA[16]  = 0.80;
  fSuSAMEffectsvsA[24]  = 0.75;
  fSuSAMEffectsvsA[27]  = 0.80;
  fSuSAMEffectsvsA[40]  = 0.74;
  fSuSAMEffectsvsA[48]  = 0.75;
  fSuSAMEffectsvsA[56]  = 0.72;
  fSuSAMEffectsvsA[59]  = 0.68;
  fSuSAMEffectsvsA[89]  = 0.65;
  fSuSAMEffectsvsA[119] = 0.66;
  fSuSAMEffectsvsA[181] = 0.66;
  fSuSAMEffectsvsA[186] = 0.80;
  fSuSAMEffectsvsA[197] = 0.74;
  fSuSAMEffectsvsA[208] = 0.63;
  fSuSAMEffectsvsA[238] = 0.65;

}
//____________________________________________________________________________
// Set the variables necessary for further calculations
void SuSAMstarUtils::SetInteraction(const Interaction * interaction)
{
  
  fInteraction = interaction;
  // Get kinematics & init-state parameters
  // unused // const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  
  eff_rat = EffectiveMassRatio(target);
  P_Fermi = FermiMomentumSuSAMstar(target);
  E_nu = interaction->InitState().ProbeE(kRfLab);
  m_lep = interaction->FSPrimLepton()->Mass();
  mm_lep     = TMath::Power(m_lep,    2);

  return;   


}
//____________________________________________________________________________
double SuSAMstarUtils::GetFermiMomentum(void) const
{
  return P_Fermi;
}
//____________________________________________________________________________
double SuSAMstarUtils::GetEffectiveMassRatio(void) const
{
   return eff_rat;
}
//____________________________________________________________________________
double SuSAMstarUtils::EffectiveMassRatio(const Target & target) const
{
  int Z = target.Z();
  int A = target.A();
  
  map<int, double >::const_iterator miter;
  
  if( (miter=fSuSAMEffectsvsA.find(A)) != fSuSAMEffectsvsA.end()) 
  {
     // tritium is a special case
     if (A == 3 && Z == 1)
       return 0.97;
     return miter->second;
  }
  // if we can't find predefined values we use
  // interpolation of theoretical dependence of effective mass from
  // Fermi momentum, obtained in Walecka model
  double pF = FermiMomentumSuSAMstar(target);
  double effect = 0.91891 + 0.00181*pF - (0.003570*pF)*(0.003570*pF);
  return effect;
  
}
//____________________________________________________________________________
double SuSAMstarUtils::FermiMomentumSuSAMstar(const Target & target) const
{
    int Z = target.Z();
    int A = target.A();
    if (A > 250) A = 250;
    if (A < 1) A = 1;
    double pF = fSuSAMpFvsA[A-1];
    
    // tritium is a special case
    if ( A == 3 && Z == 1)
        pF = 0.120;
        
    // converse from MeV to GeV
    return pF/1000.0;
    
}
//____________________________________________________________________________
double SuSAMstarUtils::PhaseSpaceVolume(KinePhaseSpace_t ps) const
{
   double vol = 0.0;
   if (ps == kPSQ2fE)
   {
     vol = kinematics::PhaseSpaceVolume(fInteraction,kPSQ2fE);
   }
   else if (ps == kPSTlctl)
   {
     vol = 2*(E_nu - m_lep);
   }
   return vol;
}
