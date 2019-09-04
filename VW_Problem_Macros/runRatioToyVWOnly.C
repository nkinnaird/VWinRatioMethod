R__LOAD_LIBRARY(/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_16_00/slf6.x86_64.e15.prof/lib/libgm2util_blinders.so)

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TVectorD.h>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <Math/Minimizer.h>

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
// #include "ratioAnalysisConfig.hh"
#include "ratioToyAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "TmethodFit.hh"
#include "ratioFit.hh"
#include "ratioCBOFit.hh"
#include "residualPlots.hh"
#include "copyFile.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

class ThreeParamRatioClass {

public:
  ThreeParamRatioClass(double xMin, double xMax) {
    threeParamRatioFunc = new TF1("threeParamRatioFunc", this, &ThreeParamRatioClass::fullFunc, xMin, xMax, 3);
    threeParamRatioFunc->SetNpx(10000);
  }

  double fullFunc(double* x, double* par){
    double freq = par[0]*1e-3;

    double fitVal = par[1] * cos(freq*x[0] + par[2]);
    return fitVal;
  };

  double Evaluate(double* x, double* p) {
    return threeParamRatioFunc->EvalPar(x, p);
    //        threeParamRatioFunc->SetParameters(p);
    // return threeParamRatioFunc->Integral(x[0] - 0.5*binWidth, x[0] + 0.5*binWidth) / binWidth;
  }

private:
  TF1* threeParamRatioFunc;
  double halfPeriodGuess = (1/blindingFa)/2.;

}; // end ThreeParamRatioClass class


/////////////////////////////////////////////////////////////////////////////////////

class RatioStyleFuncClass {

public:
  RatioStyleFuncClass(double xMin, double xMax) {
    ratioStyleFunc = new TF1("ratioStyleFunc", this, &RatioStyleFuncClass::fullFunc, xMin, xMax, 3);
    ratioStyleFunc->SetNpx(10000);
  }

  double VWPart(double t, double* par){
    double freq = par[0]*1e-3;
    double VW_part = (1 + par[1] * cos(freq*t + par[2]));
    return VW_part;
  };

  double fullFunc(double* x, double* par){
    double time = x[0];

    double f0 = VWPart(time, par);
    double fplus  = VWPart(time+halfPeriodGuess, par);
    double fminus = VWPart(time-halfPeriodGuess, par);

    double fitVal = (2*f0 - fplus - fminus) / (2*f0 + fplus + fminus);
    return fitVal;
  };

  double Evaluate(double* x, double* p) {
    return ratioStyleFunc->EvalPar(x, p);
    //        ratioStyleFunc->SetParameters(p);
    // return ratioStyleFunc->Integral(x[0] - 0.5*binWidth, x[0] + 0.5*binWidth) / binWidth;
  }

private:
  TF1* ratioStyleFunc;
  double halfPeriodGuess = (1/blindingFa)/2.;

}; // end RatioStyleFuncClass class

/*
class RatioStyleFuncClass {

public:
  RatioStyleFuncClass() {
  }

  double VWPart(double t, double* par){
    double freq = par[0]*1e-3;
    double VW_part = (1 + par[1] * cos(freq*t + par[2]));
    return VW_part;
  };

  double Evaluate(double* x, double* par){
    double time = x[0];

    double f0 = VWPart(time, par);
    double fplus  = VWPart(time+halfPeriodGuess, par);
    double fminus = VWPart(time-halfPeriodGuess, par);

    double fitVal = (2*f0 - fplus - fminus) / (2*f0 + fplus + fminus);
    return fitVal;
  };

private:
  double halfPeriodGuess = (1/blindingFa)/2.;

}; // end RatioStyleFuncClass class
*/

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int runRatioToyVWOnly(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

  // cout << "Remember to comment out VW lifetime stuff in fit classes again." << endl;
  // return 0;

  cout << "Make sure starting VW fit frequency is right for the root file I'm fitting." << endl;

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputVW.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  double toyFitStart = 30000;
  double toyFitEnd = 650000;

  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)

/////////////////////////////////////////////////////////////////////////////////////
 
  auto toyMCDir = topDir->mkdir("ToyMC");

/////////////////////////////////////////////////////////////////////////////////////

  ResidualPlots residualPlotsClass;

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC
  
  RatioFit ratioFitToyClass(myBlinder);

/////////////////////////////////////////////////////////////////////////////////////

// 5 parameter histogram fit function

    TF1* hist_fitFunc = new TF1("hist_fitFunc", "[0]*exp(-x/([1]*1000))*(1+[3]*cos(([2]/1000)*x+[4]))", toyFitStart, toyFitEnd);

    hist_fitFunc->SetLineColor(2);
    hist_fitFunc->SetNpx(10000);
    hist_fitFunc->SetParNames(N_string.c_str(), tau_string.c_str(), vw_w_string.c_str(), vw_amp_string.c_str(), vw_phi_string.c_str());

// 3 parameter ratio graph fit function

    // TF1* graph_fitFunc = new TF1("graph_fitFunc", "([1]*cos(([0]/1000)*x+[2]))", toyFitStart, toyFitEnd);
    auto myThreeParamRatioClass = new ThreeParamRatioClass(toyFitStart, toyFitEnd); // use this one for an integral fit
    auto graph_fitFunc = new TF1("graph_fitFunc", myThreeParamRatioClass, &ThreeParamRatioClass::Evaluate, toyFitStart, toyFitEnd, 3);

    graph_fitFunc->SetLineColor(2);
    graph_fitFunc->SetNpx(10000);
    graph_fitFunc->SetParNames(vw_w_string.c_str(), vw_amp_string.c_str(), vw_phi_string.c_str());


// fit function in the ratio style

    auto myRatioStyleFuncClass = new RatioStyleFuncClass(toyFitStart, toyFitEnd); // use this one for an integral fit
    // auto myRatioStyleFuncClass = new RatioStyleFuncClass();
    auto ratioStyleFunction = new TF1("ratioStyleFunction", myRatioStyleFuncClass, &RatioStyleFuncClass::Evaluate, toyFitStart, toyFitEnd, 3);
    ratioStyleFunction->SetParNames(vw_w_string.c_str(), vw_amp_string.c_str(), vw_phi_string.c_str());
    ratioStyleFunction->SetNpx(10000);
    ratioStyleFunction->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFunc = (TF1*) inputFile->Get("truthFunc"); // could maybe put this outside the file loop, but could cause problems if truth func changes even though it shouldn't
  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  // for (int iter = 0; iter < 1; ++iter)
  {
    cout << "Iter: " << iter << endl;
    auto iterDir = toyMCDir->mkdir(Form("Iter%i",iter));
    iterDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  TF1* iter_truthFunc = (TF1*) inputFile->Get(Form("topDir/Iter%d/SavedParameters/truthFunc", iter));
  for (int i = 0; i < iter_truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << iter_truthFunc->GetParameter(i) << endl;
  iter_truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////

      TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist", iter)))->Clone("5FitHist");

          int firstBin = 0;
          double firstBinTime = 0;

          for (int bin = 1; bin <= toyFiveParamHist->GetNbinsX(); ++bin){
            if(toyFiveParamHist->GetBinContent(bin) > 0){
              firstBin = bin;
              firstBinTime = toyFiveParamHist->GetBinCenter(bin);
              break;
            } 
          }
          double startingN0Guess = toyFiveParamHist->Integral("WIDTH") / (defaultLifetime * exp(-firstBinTime/defaultLifetime));

          double startingVWFreqGuess = iter_truthFunc->GetParameter(2);
          double startingVWPhaseGuess = iter_truthFunc->GetParameter(4);

          hist_fitFunc->SetParameters(startingN0Guess, defaultLifetime/1000., startingVWFreqGuess, startingVWAmp, startingVWPhaseGuess); 
          for (int parNum = 0; parNum < hist_fitFunc->GetNpar(); ++parNum) hist_fitFunc->FixParameter(parNum, hist_fitFunc->GetParameter(parNum));

          hist_fitFunc->ReleaseParameter(0);
          hist_fitFunc->SetParLimits(1, 64, 65);
    
              toyFiveParamHist->Fit(hist_fitFunc, "QMR");

          hist_fitFunc->SetParLimits(2, 0.9*startingVWFreqGuess, 1.1*startingVWFreqGuess);
          hist_fitFunc->SetParLimits(3, 0, 10);
          hist_fitFunc->SetParLimits(4, -3*pi, 3*pi);

              toyFiveParamHist->Fit(hist_fitFunc, "MR");
              cout << "Hist fit p value: " << hist_fitFunc->GetProb() << endl;

              residualPlotsClass.makeResidualPlots("histogramFit", toyFiveParamHist, "histogramFit_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist", iter)))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist", iter)))->Clone();

      TH1F* toyNumHist = (TH1F*) toyVHist->Clone("numeratorHist");
      toyNumHist->Add(toyUHist, -1);

      TH1F* toyDenomHist = (TH1F*) toyVHist->Clone("denominatorHist");
      toyDenomHist->Add(toyUHist);

      TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);

// fit parameters from starting guess

          graph_fitFunc->SetParameters(startingVWFreqGuess, startingVWAmp, startingVWPhase); 
          for (int parNum = 0; parNum < graph_fitFunc->GetNpar(); ++parNum) graph_fitFunc->FixParameter(parNum, graph_fitFunc->GetParameter(parNum));

          graph_fitFunc->SetParLimits(0, 0.9*startingVWFreqGuess, 1.1*startingVWFreqGuess);
          graph_fitFunc->SetParLimits(1, 0, 10);
          graph_fitFunc->SetParLimits(2, -3*pi, 3*pi);

// fit parameters from histogram fit

            graph_fitFunc->SetParameter(0, hist_fitFunc->GetParameter(2));
            graph_fitFunc->SetParameter(1, hist_fitFunc->GetParameter(3));
            graph_fitFunc->SetParameter(2, hist_fitFunc->GetParameter(4));

            graph_fitFunc->FixParameter(0, hist_fitFunc->GetParameter(2));
            // graph_fitFunc->FixParameter(1, hist_fitFunc->GetParameter(3));
            graph_fitFunc->FixParameter(2, hist_fitFunc->GetParameter(4));


              toyRatioGraph->Fit(graph_fitFunc, "MR");
              cout << "Three param fit func p value: " << graph_fitFunc->GetProb() << endl;

              toyRatioGraph->SetTitle("Toy Ratio Graph; Time (ns); R (unitless)");
              toyRatioGraph->SetName("Toy_Ratio_Graph");
              toyRatioGraph->Write();

              residualPlotsClass.makeResidualPlots("ratioGraphFit", toyRatioGraph, "ratioGraphFit_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


      TGraphErrors* ratioStyleFitFuncGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);

// fit parameters from starting guess

          ratioStyleFunction->SetParameters(startingVWFreqGuess, startingVWAmp, startingVWPhase); 
          for (int parNum = 0; parNum < ratioStyleFunction->GetNpar(); ++parNum) ratioStyleFunction->FixParameter(parNum, ratioStyleFunction->GetParameter(parNum));

          ratioStyleFunction->SetParLimits(0, 0.9*startingVWFreqGuess, 1.1*startingVWFreqGuess);
          ratioStyleFunction->SetParLimits(1, 0, 10);
          ratioStyleFunction->SetParLimits(2, -3*pi, 3*pi);

// fit parameters from histogram fit

            ratioStyleFunction->SetParameter(0, hist_fitFunc->GetParameter(2));
            ratioStyleFunction->SetParameter(1, hist_fitFunc->GetParameter(3));
            ratioStyleFunction->SetParameter(2, hist_fitFunc->GetParameter(4));

            ratioStyleFunction->FixParameter(0, hist_fitFunc->GetParameter(2));
            // ratioStyleFunction->FixParameter(1, hist_fitFunc->GetParameter(3));
            ratioStyleFunction->FixParameter(2, hist_fitFunc->GetParameter(4));

            // ratioStyleFunction->FixParameter(1, graph_fitFunc->GetParameter(1));

              ratioStyleFitFuncGraph->Fit(ratioStyleFunction, "R");
              cout << "Full ratio fit func p value: " << ratioStyleFunction->GetProb() << endl;

              ratioStyleFitFuncGraph->SetTitle("Toy Ratio Graph Full Func; Time (ns); R (unitless)");
              ratioStyleFitFuncGraph->SetName("Toy_Ratio_Graph_Full_Func");
              ratioStyleFitFuncGraph->Write();

              residualPlotsClass.makeResidualPlots("ratioGraphRatioStyleFit", ratioStyleFitFuncGraph, "ratioGraphRatioStyleFit_residual", toyFitStart, toyFitEnd, true);


} // end loop over iterations

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
