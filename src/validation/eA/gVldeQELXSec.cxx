//____________________________________________________________________________________________
/*!

\program gvld_e_qel_xsec

\brief   Compares GENIE with quasi-elastic electron scattering data.

         The data come from the Quasielastic Electron Nucleus Scattering Data
         Archive maintained by Donal Day:
         http://http://faculty.virginia.edu/qes-archive/
         See also O. Benhar, D. Day and I. Sick, Rev.Mod.Phys.80, 189, 2008
         A local copy of the data may be found in:
         $GENIE/data/validation/eA/xsec/differential/qe
         The archive contains ~19k data points.

         Syntax:
           gvld_e_qel_xsec 
                 [-g input_file_list] 
                 [-d data_archive]
                 [-s datasets_to_plot]

         Options:

           [] Denotes an optional argument.

           -d Full path to the electron QE archive.
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/qe/eQE.root

           -s Info on which datasets to plot. 
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/qe/datasets.txt

           -g An XML file with GENIE inputs (cross sections and event samples).
              If not set, only data -no GENIE predictions- will be displayed.
              Multiple models can be included in the input file, each identified
              by a "name" (all model predictions will be overlayed).
              For info on the XML file format see the GSimFiles class documentation.
              Notes:
              - A script for preparing inputs for this benchmark test can be found in: 
                $GENIE/src/scripts/production/batch/batch/submit-vld_e_semi_incl_xsec.pl
              - The input event files are `gst' summary ntuples generated by
                GENIE gntpc utility.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Oct 16, 2009 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>
#include <TBox.h>

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/GSimFiles.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"

using std::ostringstream;
using std::ifstream;
using std::string;
using std::vector;

using namespace genie;
using namespace genie::constants;

//.................................................................................
// Utility class to hold info on plotted datasets
//
class eQEDataSetDescription
{
public:
  eQEDataSetDescription(
     int tgtpdg, string citation, double E, double theta) :
    fTgtPdg   (tgtpdg), 
    fCitation (citation), 
    fE        (E), 
    fTheta    (theta)
  {   
  }
  eQEDataSetDescription() 
  {   
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  string TgtName  (void) const { 
    if(fTgtPdg == 1000010010) return "1H";
    if(fTgtPdg == 1000010020) return "2D";
    if(fTgtPdg == 1000010030) return "3H";
    if(fTgtPdg == 1000020030) return "3He";
    if(fTgtPdg == 1000020040) return "4He";
    if(fTgtPdg == 1000060120) return "12C";
    if(fTgtPdg == 1000080160) return "16O";
    if(fTgtPdg == 1000130270) return "27Al";
    if(fTgtPdg == 1000200400) return "40Ca";
    if(fTgtPdg == 1000200480) return "48Ca";
    if(fTgtPdg == 1000260560) return "56Fe";
    if(fTgtPdg == 1000791970) return "197Au";
    if(fTgtPdg == 1000822080) return "208Pb";
    if(fTgtPdg == 1000922380) return "238U";
    return "Other";
  }
  string Citation (void) const { return fCitation; }
  double E        (void) const { return fE; }
  double Theta    (void) const { return fTheta; }
  string LabelTeX (void) const { 
    ostringstream label;
    label << this->TgtName() << "(e,e^{'}) ";
    label << "E = " << fE << " GeV, ";
    label << "#theta = " << fTheta << "^{o}";
    label << " [" << fCitation << "]";
    return label.str();
  }
private:
  int    fTgtPdg;   //
  string fCitation; //
  double fE;        //
  double fTheta;    //
};
//.................................................................................

//
// constants
//

// defaults
const char * kDefDataArchiveFilename = "data/validation/eA/xsec/differential/qe/eQE.root";  
const char * kDefDataSetsFilename    = "data/validation/eA/xsec/differential/qe/datasets.txt";  

// plot config
const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

// model line styles
const int kNMaxNumModels = 5;
const int kLStyle[kNMaxNumModels] = 
{
   1, 2,  3,  5, 6
};
string kLStyleTxt[kNMaxNumModels] = 
{
  "solid", "dashed", "dotted", "dot-dashed", "dot-dot-dashed"
};

//
// globals
//

string gOptDataArchiveFilename = ""; // -d command-line argument
string gOptDataSetsFilename    = ""; // -s command-line argument

GSimFiles gOptGenieInputs;

TFile *        gQEDataFile  = 0;
TTree *        gQEDataTree  = 0;
TPostScript *  gPS          = 0;
TCanvas *      gC           = 0;
bool           gShowModel   = false;

vector<eQEDataSetDescription *> gDataSets; // list of plotted datasets

//
// function prototypes
//

void           Init               (void);
void           End                (void);
TGraphErrors * Data               (int iset);
TGraph *       Model              (int iset, int imodel);
void           Draw               (int iset);
void           GetCommandLineArgs (int argc, char ** argv);
void           PrintSyntax        (void);


//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();

  // loop over data sets and plot data and corresponding GENIE predictions
  for(unsigned int iset = 0; iset < gDataSets.size(); iset++) 
  {
    LOG("gvldtest", pNOTICE) 
      << "Producing plots for: " << gDataSets[iset]->LabelTeX();
    Draw(iset);
  }

  End();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";

  // Set GENIE style
  utils::style::SetDefaultStyle();

  //
  // Get TTree with electron-scattering data
  //
  if( ! utils::system::FileExists(gOptDataArchiveFilename) ) {
      LOG("gvldtest", pFATAL) 
         << "Can not find file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gQEDataFile = new TFile(gOptDataArchiveFilename.c_str(),"read");  
  gQEDataTree = (TTree *) gQEDataFile->Get("qent");
  if(!gQEDataTree) {
      LOG("gvldtest", pFATAL) 
         << "Can not find TTree `qent' in file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }

  //
  // Read information on which data-sets to plot
  //
  LOG("gvldtest", pDEBUG) 
    << "Reading dataset summary info from: " << gOptDataSetsFilename;
  ifstream summary_file(gOptDataSetsFilename.c_str());
  if (!summary_file.good() ) {
      LOG("gvldtest", pFATAL) 
         << "Can't open data summary file: " << gOptDataSetsFilename;
      gAbortingInErr = true;
      exit(1);
  }
  while(1) {
      // skip header lines staring with #
      if(summary_file.peek() == '#') {
         summary_file.ignore(1000, '\n');
       } else {
         int    target   = 0;
         string citation = "";
         double E        = 0;
         double theta    = 0;
         summary_file >> target >> citation >> E >> theta;
         summary_file.ignore(1000, '\n');
         if(summary_file.eof()) break;            

         LOG("gvldtest", pNOTICE) 
            << "target: " << target << ", citation: " << citation
            << ", E = " << E << " GeV, theta = " << theta << " deg";

         eQEDataSetDescription * dataset = 
                 new eQEDataSetDescription(target,citation,E,theta);
         gDataSets.push_back(dataset);
       }
  }
  summary_file.close();
  LOG("gvldtest", pNOTICE) 
     << "Read "  << gDataSets.size() << " datasets";



  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // Create output postscript file
  string localtime = utils::system::LocalTimeAsString("%d.%d.%d_%d.%d.%d");
  string filename  = Form("genie-e_qe_data_comp-%s.ps",localtime.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with e QE data");
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;

  gQEDataFile->Close();
}
//_________________________________________________________________________________
TGraph * Model(int iset, int imodel)
{
// Corresponding GENIE prediction for the `iset' data set 

  if(!gShowModel) return 0;

  LOG("gvldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  return 0;
/*

  TFile * xsec_file = gOptGenieInputs.XSecFile(imodel);
  if(!xsec_file) {
     LOG("vldtest", pWARN)
        << "No corresponding cross section file";
     return 0;
  }

  TChain * event_chain = gOptGenieInputs.EvtChain(imodel);
  if(!event_chain) {
     LOG("vldtest", pWARN)
        << "No corresponding event chain.";
     return 0;
  }
  
  // electron energy and scattering angle

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();

  double dE       = 0.01; // GeV
  double costheta = TMath::Cos(kPi*theta/180.);

  // total cross section
  TGraph * xsec_em_gr = 0; // get it from the input xsec_file
  if(!xsec_em_gr) {
     LOG("vldtest", pWARN)
        << "Null E/M cross section graph";
     return 0;
  }
  double xsec_em = xsec_em_gr->Eval(E);
  if(xsec_em <= 0) {
     LOG("vldtest", pWARN)
        << "Null E/M cross section graph";
     return 0;
  }
  
  //
  // book histograms
  //

  // E', final state primary lepton energy
  double Ep_min   = 0;
  double Ep_max   = E;
  double Ep_bin   = 0.01;
  double Ep_nbins = (Ep_max-Ep_min)/Ep_bin;

  // theta, scattering angle
  double costh_min   = -1;
  double costh_max   =  1;
  double costh_bin   =  0.01;
  double costh_nbins = (costh_max-costh_min)/costh_bin;
   
  TH2D * h2EpOmega =
     new TH2D("h2EpOmega",  "N=N(E',Omega)|{fixed:E}",
         Ep_nbins, Ep_min, Ep_max, costh_nbins, costh_min, costh_max);
     
  //
  // estimate d^2 sigma / dE' dOmega at the current incoming lepton energy E0
  //

  char cut[500];
  Form(cut,"em&&(fabs(Ev-%d)<%d)", E, dE);

  event_chain->Draw("((pxv*pxl+pyv*pyl+pzv*pzl)/(Ev*El)):El>>h2EpOmega", cut, "GOFF");

  double integral = h2EpOmega->Integral("width");
  if(integral <= 0) {
     LOG("vldtest", pWARN)
        << "Non-positive d^2N / dEp dOmega integral";
     return 0;
  }
  double normalization = 2*kPi*xsec_em/integral;
     
  h2EpOmega->Scale(normalization); // units: 1E-38 cm^2 / GeV /sterad
  
  //
  // now pick a slice at selected theta and return
  // d^2 sigma / dE' dOmega (fixed: E, theta) = f(v = E-E')
  // in the same units as the expt data (nbar/GeV/sterad)
  //

  int N = Ep_nbins;
  
  double * x = new double[N]; // v
  double * y = new double[N]; // d^2 sigma / dE' dOmega
  
  int costheta_bin = h2EpOmega->GetYaxis()->FindBin(costheta);
  
  for(int i = 0; i < h2EpOmega->GetNbinsX(); i++) {
    int Ep_bin = i+1;
    double Ep  = h2EpOmega->GetXaxis()->GetBinCenter(Ep_bin);

    double v      = E - Ep;
    double d2xsec = h2EpOmega->GetBinContent(Ep_bin, costheta_bin);
  
    x[i] = v;
    y[i] = d2xsec;
  }
     
  TGraph * gr = new TGraph(N,x,y);

  delete [] x;
  delete [] y;
    
  return gr;
*/
}
//_________________________________________________________________________________
TGraphErrors * Data(int iset)
{
  const double dE      = 2.0E-3;
  const double dtheta  = 2.5E-2;

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();

  const char * selection = 
    Form("E > %f && E < %f && theta > %f && theta < %f && Z == %d && A == %d",
         E     - dE,
         E     + dE,
         theta - dtheta,
         theta + dtheta,
         Z,A);

  gQEDataTree->Draw("v:xsec:xsec_err", selection, "goff");
  
  int n = gQEDataTree->GetSelectedRows();

  LOG("gvldtest", pNOTICE) 
    << "Found " << n << " data points in the xsec archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in W
  // Do the ordering here before building the graph
  int    *  idx = new int   [n];
  double *  xv  = new double[n];
  double *  yv  = new double[n];
  double *  dyv = new double[n];

  TMath::Sort(n,gQEDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     xv [i] = (gQEDataTree->GetV1())[ii];
     yv [i] = (gQEDataTree->GetV2())[ii];
     dyv[i] = (gQEDataTree->GetV3())[ii];
  }

  TGraphErrors * gr = new TGraphErrors(n,xv,yv,0,dyv);

  delete [] idx;
  delete [] xv;
  delete [] yv;
  delete [] dyv;

  return gr;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  TGraphErrors * data = Data(iset);

  // get the corresponding GENIE model prediction
  TGraph * model = Model(iset,0);

  if(!model && !data) return;

  int plots_per_page = kNCx * kNCy;
  int iplot = 1 + iset % plots_per_page;

  if(iplot == 1) {
     gPS->NewPage();
     gC -> Clear();
     gC -> Divide(kNCx,kNCy);
  }

  gC -> GetPad(iplot) -> Range(0,0,100,100);
  gC -> GetPad(iplot) -> SetFillColor(0);
  gC -> GetPad(iplot) -> SetBorderMode(0);
  gC -> GetPad(iplot) -> cd();

  double xmin = 0.0, scale_xmin = 0.5;
  double xmax = 0.0, scale_xmax = 1.2;
  double ymin = 0.0, scale_ymin = 0.4;
  double ymax = 0.0, scale_ymax = 1.2;

  TH1F * hframe = 0;
  bool have_frame = false;

  if(data) {
    xmin  = ( data->GetX() )[TMath::LocMin(data->GetN(),data->GetX())];
    xmax  = ( data->GetX() )[TMath::LocMax(data->GetN(),data->GetX())];
    ymin  = ( data->GetY() )[TMath::LocMin(data->GetN(),data->GetY())];
    ymax  = ( data->GetY() )[TMath::LocMax(data->GetN(),data->GetY())];
    if(model) {
       ymin  = TMath::Min(
        ymin, ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())]);
       ymax  = TMath::Max(
        ymax, ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())]);
    }
    hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    have_frame = true;
    utils::style::Format(data, 1,1,1,1,8,0.8);
    data->Draw("P");
  }//data?

  if(model) {
    if(!have_frame) {
       xmin  = ( model->GetX() )[TMath::LocMin(model->GetN(),model->GetX())];
       xmax  = ( model->GetX() )[TMath::LocMax(model->GetN(),model->GetX())];
       ymin  = ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())];
       ymax  = ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())];
       hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
         scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    }
    utils::style::Format(model, 1,1,1,1,1,1);
    model->Draw("L");
  }

  //hframe->Draw();
  hframe->GetXaxis()->SetTitle("#nu = E-E^{'} (GeV)");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");
  hframe->GetXaxis()->SetLabelFont(62);
  hframe->GetYaxis()->SetLabelFont(62);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.65);

/*
  // scaling region
  TBox * scaling_region = 0;
  if(kDrawHatchcedScalingRegion) {
    double W2c = kWcut*kWcut;
    if(W2c > scale_xmin*xmin && W2c < scale_xmax*xmax) {
       scaling_region = new TBox(
           W2c, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
       scaling_region->SetFillColor(kRed);
       scaling_region->SetFillStyle(3005);
       scaling_region->Draw();
    }
  }
*/

/*
  // some data show the elastic peak - mark the are to avoid confusion
  if(xmin < 1) {
    double Wm2 = 1.21; // between the QE and Delta peaks
    TBox * qe_peak = new TBox(
       scale_xmin*xmin, scale_ymin*ymin, Wm2, scale_ymax*ymax);
     qe_peak->SetFillColor(kBlue);
     qe_peak->SetFillStyle(3005);
     qe_peak->Draw();
  }
*/

  // title
  TLatex * title = new TLatex(
     scale_xmin*xmin + 0.2*(scale_xmax*xmax-scale_xmin*xmin),
    1.01*scale_ymax*ymax,gDataSets[iset]->LabelTeX().c_str());
  title->SetTextSize(0.027);
  title->Draw();

  gC->GetPad(iplot)->Update();
  gC->Update();
}
//_________________________________________________________________________________
// Parsing command-line arguments, check/form filenames, etc
//.................................................................................
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataArchiveFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataArchiveFilename;
        gOptDataArchiveFilename = filename;
     } else {
        LOG("gvldtest", pFATAL)
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  if(parser.OptionExists('s')){
     string filename = parser.ArgAsString('s');
     gOptDataSetsFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataSetsFilename;
        gOptDataSetsFilename = filename;
     } else {
        LOG("gvldtest", pFATAL)
          << "\n Please make sure that $GENIE is defined, or use the -s option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // get GENIE inputs
  gShowModel = true;
  if(parser.OptionExists('g')) {
     string inputs = parser.ArgAsString('g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) {
        LOG("gvldtest", pFATAL)
           << "Could not read validation program inputs: " << inputs;
        gAbortingInErr = true;
        exit(1);
     }
  } else {
    gShowModel = false;
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "  gvld_e_qel_xsec -g inputs [-d data_archive_location]\n";
}
//_________________________________________________________________________________