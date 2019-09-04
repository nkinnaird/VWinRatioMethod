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

#include "../../ratioMacroHeaders/ratioAnalysisDefs.hh" // include paths for potential future grid jobs
#include "../../ratioMacroHeaders/ratioToyAnalysisConfig.hh"

using namespace std;

class TruthFuncClass{

public:
  TruthFuncClass()
  {}

  double VWPart(double t, double* par){
    double VW_part;
    // double lifetime = par[] * 1000.;  // convert from us to ns for the fit function
    double frequency;

    frequency = par[2] * 1e-3; // convert to rad/ns for fit function from rad/us, which is the units the parameter is in

    // VW_part = (1 + exp(-t/lifetime) * par[1] * cos(frequency*t + par[2]));
    VW_part = (1 + par[3] * cos(frequency*t + par[4]));
    
    return VW_part;
  };

  double Evaluate(double* x, double* par){
      double time = x[0];
      double lifetime = par[1] * 1000.; // convert from us to ns for fit function

      double value = VWPart(time, par) * par[0] * exp(-time/lifetime);
      return value;
  };

private:

}; // end TruthFuncClass


int makeRatioToyHistsVWOnly()
{
  clock_t t; // for time profiling - put after input
  t = clock();

/////////////////////////////////////////////////////////////////////////////////////

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(0);//444); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0);//43210);// for number of entries
  TRandom3* randGen_UV = new TRandom3(0);//112233);

  int func_iterRandSeed = gRandom->Integer(1e8);
  int UV_iterRandSeed = randGen_UV->Integer(1e8);


  double nPts = 1.5e9; // 1.5e9 for about the same precision as the 60 hr dataset

  // create output file that will hold plots
  TFile* outputFile = new TFile("ratioToyHists.root","RECREATE");

  // make top directory for output file
  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double functionStart = 10000; // can't be 0 because the frequency model blows up there

  int numParameters = 5;

  auto myTruthFuncClass = new TruthFuncClass();
  auto truthFunc = new TF1("truthFunc", myTruthFuncClass, &TruthFuncClass::Evaluate, functionStart, histMaxTime, numParameters);
  truthFunc->SetNpx(100000);

  truthFunc->SetParameter(0, 1);
  truthFunc->SetParameter(1, defaultLifetime/1000.);

  truthFunc->SetParameter(2, startingVWFreq);
  truthFunc->SetParameter(3, startingVWAmp);
  truthFunc->SetParameter(4, startingVWPhase);

    TF1* temp_expFunc = new TF1("temp_expFunc", "[0]*exp(-x/[1])", 0, histMaxTime);
    temp_expFunc->SetParameter(0, 1); 
    temp_expFunc->SetParameter(1, defaultLifetime);
    double Ntoy = nPts / ( temp_expFunc->Integral(0, histMaxTime) / binWidth );

  truthFunc->SetParameter(0, Ntoy);


  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

  truthFunc->Write();



/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 1;

  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");

/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFuncCopies[totalIters];
  TRandom3* truthFuncGRandoms[totalIters];
  TRandom3* randGen_UV_iters[totalIters];

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double gm2PeriodGuess = 1/blindingFa;

/////////////////////////////////////////////////////////////////////////////////////

for (int iter = 0; iter < totalIters; ++iter) // loop over totalIters
{

  histIterDirs[iter] = topDir->mkdir(Form("Iter%d",iter));
  histIterDirs[iter]->cd();

  toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);

  auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
  savDirs->cd();

/////////////////////////////////////////////////////////////////////////////////////

  truthFuncCopies[iter] = new TF1();
  truthFunc->Copy(*truthFuncCopies[iter]);

  // truthFuncCopies[iter]->SetParameter(2, (1 + 0.25*iter)*blindingWa*1e3);
  truthFuncCopies[iter]->SetParameter(4, 1 + pi/2);
  truthFuncCopies[iter]->SetNpx(100000);
  truthFuncCopies[iter]->Write();

  // cout << "Iter: " << iter << " freq value: " << truthFuncCopies[iter]->GetParameter(2) << endl;
  cout << "Iter: " << iter << " phase value: " << truthFuncCopies[iter]->GetParameter(4) << endl;

  truthFuncGRandoms[iter] = new TRandom3(func_iterRandSeed);
  randGen_UV_iters[iter] = new TRandom3(UV_iterRandSeed);

}


    double randEntries = randGen_positrons->PoissonD(nPts);
    double tenIncrement = 0.1 * randEntries;
    double countdown = tenIncrement;
    int tenPercent = 0;

    double savedTime = 0;

        for (double entry = 0; entry < randEntries; ++entry)
        {

          if(--countdown <= 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }

          // double time_original = truthFunc->GetRandom();

/////////////////////////////////////////////////////////////////////////////////////
         
              for (int iter = 0; iter < totalIters; ++iter)
              {

                gRandom = truthFuncGRandoms[iter];
                double time_original = truthFuncCopies[iter]->GetRandom();

/////////////////////////////////////////////////////////////////////////////////////

                double time = time_original;

                toyFiveParamHists[iter]->Fill(time);

                double randNum = randGen_UV_iters[iter]->Uniform();

                double halfPeriodGuess = gm2PeriodGuess/2;

                double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
                else if(randNum < 1) toyVHists[iter]->Fill(time);

              }
        }


/////////////////////////////////////////////////////////////////////////////////////

    t = clock() - t;
    printf ("It took me %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);


/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
