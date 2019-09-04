#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

void verticalWaistCaloByCalo(){

  // Params
  long long nEntries = 1e9;
  double maxTime = 650064.4;
  double binWidth = 149.2;
  int nBins = maxTime/binWidth;

  // Detector / fast rotation parameters
  const int ndet = 24;
  double pmumean = 3.094349631;
  double pmusd = 0.001 * pmumean; //0.1%
  double rMagic = 711.214761;
  double trevav = TMath::TwoPi()*rMagic/(0.99942*29.979246);
  double freqav = 1.0/trevav;

  // Input information for realistic pulses
  Double_t inputT0Content[26] = { 0, 1033.210, 2878.228, 10258.30, 12250.92, 12841.32, 13284.13, 15350.55, 20590.40, 24575.64, 32693.72, 37859.77, 41697.41, 38892.98, 35129.15, 26937.26, 21107.01, 15276.75, 10996.30, 9667.896, 8634.686, 7011.070, 4428.044, 1623.616, 811.80, 0};
  Double_t inputRadBins[38] = { 708.0275, 708.2379, 708.4483, 708.6589, 708.8696, 709.0804, 709.2913, 709.5024, 709.7136, 709.9249, 710.1363, 710.3479, 710.5596, 710.7714, 710.9834, 711.1954, 711.4076, 711.6200, 711.8324, 712.0450, 712.2577, 712.4706, 712.6835, 712.8966, 713.1098, 713.3232, 713.5367, 713.7503, 713.9640, 714.1778, 714.3918, 714.6059, 714.8202, 715.0346, 715.2491, 715.4637, 715.6784, 715.8933};
  Double_t inputRadContent[38] = { 0.003010, 0.001786, -0.001229, -0.001294, 0.002887, 0.003798, 0.006922, 0.007412, 0.025320, 0.080977, 0.208575, 0.401118, 0.627794, 0.819469, 0.933363, 0.996938, 1.000000, 0.989440, 0.964585, 0.931940, 0.881744, 0.817763, 0.762939, 0.687833, 0.583787, 0.434679, 0.240783, 0.096240, 0.031282, 0.009576, 0.004438, 0.005673, 0.001774, -0.000256, -0.000667, -0.000568, 0.000937, 0.005438};

  // Input Spectra
  TH1D* hInputT0 = new TH1D("hInputT0","Input T_{0} Shape; Time [ns]",26,-65,65);
  for(int i = 1; i <= 26; i++) hInputT0->SetBinContent(i,inputT0Content[i-1]);
  TH1D* hInputRad = new TH1D("hInputRad","Input Radius; x_{e} [cm]", 37, inputRadBins);  // Values taken from Antoine's 60h FT analysis - could easily be off by half a bin
  for(int i = 1; i <= 37; i++){
    if(inputRadContent[i-1] > 0.008) hInputRad->SetBinContent(i,inputRadContent[i-1]);
  }

  // Truth parameters
  double tau = 64440;
  double A = 0; //0.3714;
  double w_a = 0.00144;
  double phi = 2.081;
  
  double A_vw = 0.2; // Exagerrated to make effect clearer
  double tau_vw = 250800;
  // double w_vw = 8.913*w_a; //10.06*w_a;
  double w_vw = 10.0*w_a;
  double phi_vw = 0;
  
  // Get N0 from integral of function (need to do numeric as we truncate exponential slightly and wiggle introduces error too)
  auto normalisedFunc = new TF1("normalisedFunc","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",0,maxTime);
  normalisedFunc->SetParameters(1, tau, 0, w_a, phi);
  double N0 = (binWidth*nEntries)/normalisedFunc->Integral(0,maxTime);

  // Declare truth function
  auto truthFunc = new TF1("truthFunc","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",0,maxTime);
  truthFunc->SetParameters(N0, tau, 0, w_a, phi);
  truthFunc->SetParNames("N0", "tau", "A", "w_a", "phi");
  truthFunc->SetLineColor(2);
  truthFunc->SetLineWidth(1);
  truthFunc->SetNpx(150000);

  // Do some pseudo data trials
  TRandom3* rng = new TRandom3(54321);
  double half_t_a = TMath::Pi()/w_a;

  // Normalisation for weighted parts
  double denom = 2+2*cosh(half_t_a/tau);
  double VThresh = 2./denom;
  double UThresh_p = VThresh + exp(-half_t_a/tau)/denom;
  double UThresh_m = UThresh_p + exp(+half_t_a/tau)/denom;
  std::cout << "Cut offs : " << VThresh << "\t" << UThresh_p << "\t" << UThresh_m << endl;

  // Hit times - all calos together
  auto hData = new TH1D("hData",";Time [ns]; Entries",nBins,0,maxTime);
  auto hData_U = new TH1D("hData_U",";Time [ns]; Entries",nBins,0,maxTime);
  auto hData_V = new TH1D("hData_V",";Time [ns]; Entries",nBins,0,maxTime);
  auto hData_num = new TH1D("hData_num",";Time [ns]; Entries",nBins,0,maxTime);
  auto hData_den = new TH1D("hData_den",";Time [ns]; Entries",nBins,0,maxTime);
  auto hDataFR = new TH1D("hDataFR","Raw Signal;Time [#mus]", 450000, -0.0005, 449999.5);

  // Hit times - individual calos
  TH1D* hDataCalo[ndet];
  TH1D* hDataCalo_U[ndet];
  TH1D* hDataCalo_V[ndet];
  TH1D* hDataCalo_num[ndet];
  TH1D* hDataCalo_den[ndet];
  for(int i = 0; i < ndet; i++){
    hDataCalo[i] = new TH1D(Form("hDataCalo_%02d",i+1),";Time [ns]; Entries",nBins,0,maxTime);
    hDataCalo_U[i] = new TH1D(Form("hDataCalo_U_%02d",i+1),";Time [ns]; Entries",nBins,0,maxTime);
    hDataCalo_V[i] = new TH1D(Form("hDataCalo_V_%02d",i+1),";Time [ns]; Entries",nBins,0,maxTime);
    hDataCalo_num[i] = new TH1D(Form("hDataCalo_num_%02d",i+1),";Time [ns]; Entries",nBins,0,maxTime);
    hDataCalo_den[i] = new TH1D(Form("hDataCalo_den_%02d",i+1),";Time [ns]; Entries",nBins,0,maxTime);
  }

  // Fast rotation histograms
  TH1D* hMomSpec = new TH1D("hMomSpec","Momentum Spectrum (GeV/c)",101,0,-1);
  TH1D* hT0Pulse = new TH1D("hT0Pulse","T0Pulse (ns) ",201,-100.5,100.5);
  TH1D* hDataCaloFR[ndet];
  for(int idet = 0; idet < ndet; idet++){
    hDataCaloFR[idet] = new TH1D(Form("hDataCaloFR_%02d",idet),"Raw Signal;Time [#mus]", 450000, -0.0005, 449999.5);
  }
  TH1D* hTimeDecay = new TH1D("hTimeDecay","Muon Decay Time ",150000,0.0,150000.0);
  TH1D* hDetector = new TH1D("hDetector","Detector Number ",24,0.5,24.5);
  TH1D* hnTurns = new TH1D("hnTurns","Number of turns",10000,0.0,6000.0);
  TH1D* hTTT = new TH1D("hTTT","Decay point (in number of turns)",10000,0.0,1.0);

  // Set seed for all the GetRandom() calls
  gRandom->SetSeed(12345);

  // Run pseudo experiment!
  for(long long i = 0; i < nEntries; i++){

    if(i % 1000000 == 0) cout << "Muon number " << i << endl;

    /////////////////////////////////////////////////////////////////////////////
    // Add in wiggle
    /////////////////////////////////////////////////////////////////////////////

    //Now choose decay time from wiggle plot - this is relative to injection (when the precession starts)
    double tdec = truthFunc->GetRandom();

    // Modulate decay to get wiggle
    double wiggle = (1-A)*(1+A*cos(w_a*tdec + phi));
    if(rng->Uniform() > wiggle) continue;

    /////////////////////////////////////////////////////////////////////////////
    // Do fast rotation magic to figure out which calorimeter we're filling and what time we should use
    /////////////////////////////////////////////////////////////////////////////

    //Select a T0 for the muon from histogram of t0 pulse
    double t0 = -hInputT0->GetRandom();  // minus as smaller t0 shows up later
    // double t0 = -hInputT0->GetMean();  // minus as smaller t0 shows up later
    
    //Now select "momentum" (actually a revolution time) from input radial distribution
    double rmu = hInputRad->GetRandom();
    double pmu = pmumean*rmu/rMagic;
    // double pmu = pmumean; // just magic momentum

    // Convert to time
    double trev = trevav*pmu/pmumean;

    //Calculate number of turns, detector
    double nturns = tdec/trev;    // Number of turns at which decay occurs
    double ttt = nturns-(int)nturns;   // This is distance past detector 0 in turns
    int caloNum = (int)(ndet*ttt) + 1; // This is the calorimeter we're going to fill (numbered 1 - 24) - just after decay point

    //Calculate time that particle would hit detector (assuming that drift time is just muon travel time)
    double tfill = tdec + (float(caloNum-1)/ndet - ttt)*trev + t0; 

    // // OVERWRITE FAST ROTATION - ALSO COMMENT OUT CALO-ACCEPT AND FR TIME-SMEARING
    // ttt = rng->Uniform();
    // caloNum = (int)(ndet*ttt) + 1;
    // tfill = tdec;

    /////////////////////////////////////////////////////////////////////////////
    // Modulate acceptance based on calorimeter number (so that it doesn't completely cancel)
    /////////////////////////////////////////////////////////////////////////////
    vector<double> caloAccept = {0.993, 0.665, 0.846, 0.845, 0.763, 0.830, 0.867, 0.846, 0.776, 0.874, 0.944, 1.000, 0.936, 0.926, 0.873, 0.832, 0.892, 0.996, 0.906, 0.907, 0.806, 0.802, 0.897, 0.867};
    if(rng->Uniform() > caloAccept[caloNum-1]) continue;

    /////////////////////////////////////////////////////////////////////////////
    // Modulate acceptance based on turn and time (for vertical waist)
    /////////////////////////////////////////////////////////////////////////////
    double VW = (1.0 - A_vw) * (1 + A_vw*exp(-tdec/tau_vw)*cos(w_vw*tdec - ttt*TMath::TwoPi()));
    if(rng->Uniform() > VW) continue;

    /////////////////////////////////////////////////////////////////////////////
    // Randomize time to get rid of fast rotation
    /////////////////////////////////////////////////////////////////////////////
    tfill += binWidth*(rng->Uniform() - 0.5);

    /////////////////////////////////////////////////////////////////////////////
    // Fill histograms
    /////////////////////////////////////////////////////////////////////////////

    hT0Pulse->Fill(t0);
    hMomSpec->Fill(pmu);
    hTimeDecay->Fill(tdec);
    hnTurns->Fill(nturns);
    hTTT->Fill(ttt);
    hDetector->Fill(caloNum);

    hData->Fill(tfill);
    hDataFR->Fill(tfill);

    hDataCalo[caloNum-1]->Fill(tfill);
    hDataCaloFR[caloNum-1]->Fill(tfill); // Array runs 0 - 23

    double r = rng->Uniform();
    if(r < VThresh){
      hData_V->Fill(tfill);
      hDataCalo_V[caloNum-1]->Fill(tfill);
    } else if (r < UThresh_p){
      hData_U->Fill(tfill+half_t_a);
      hDataCalo_U[caloNum-1]->Fill(tfill+half_t_a);
    } else {
      hData_U->Fill(tfill-half_t_a);
      hDataCalo_U[caloNum-1]->Fill(tfill-half_t_a);
    }
  }

  // Make ratio fit graph - all calos
  auto gData_ratio = new TGraphErrors();
  gData_ratio->SetName("Ratio");
  hData_num->Add(hData_U,hData_V,-1);
  hData_den->Add(hData_U,hData_V);
  for(int i = 1; i <= hData_num->GetNbinsX(); i++){
    int nPt = gData_ratio->GetN();
    if(hData_den->GetBinContent(i) > 0 && hData_num->GetBinCenter(i) > half_t_a){
      double ratio = hData_num->GetBinContent(i)/hData_den->GetBinContent(i);
      double ratioErr = sqrt((1. - ratio*ratio) / hData_den->GetBinContent(i));
      gData_ratio->SetPoint(nPt,hData_num->GetBinCenter(i),ratio);
      gData_ratio->SetPointError(nPt,0,ratioErr);
    }
  }

  // Make ratio fit graphs - individual calos
  TGraphErrors* gDataCalo_ratio[ndet];
  for(int iCalo = 0; iCalo < ndet; iCalo++){
    gDataCalo_ratio[iCalo] = new TGraphErrors();  
    gDataCalo_ratio[iCalo]->SetName(Form("RatioCalo_%02d",iCalo+1));
    hDataCalo_num[iCalo]->Add(hDataCalo_U[iCalo],hDataCalo_V[iCalo],-1);
    hDataCalo_den[iCalo]->Add(hDataCalo_U[iCalo],hDataCalo_V[iCalo]);
    for(int i = 1; i <= hDataCalo_num[iCalo]->GetNbinsX(); i++){
      int nPt = gDataCalo_ratio[iCalo]->GetN();
      if(hDataCalo_den[iCalo]->GetBinContent(i) > 0 && hDataCalo_num[iCalo]->GetBinCenter(i) > half_t_a){
	double ratio = hDataCalo_num[iCalo]->GetBinContent(i)/hDataCalo_den[iCalo]->GetBinContent(i);
	double ratioErr = sqrt((1. - ratio*ratio) / hDataCalo_den[iCalo]->GetBinContent(i));
	gDataCalo_ratio[iCalo]->SetPoint(nPt,hDataCalo_num[iCalo]->GetBinCenter(i),ratio);
	gDataCalo_ratio[iCalo]->SetPointError(nPt,0,ratioErr);
      }
    }
  }


  // Write all the plots to an output file
  TFile* outFile = new TFile("verticalWaistHists.root","RECREATE");
  hMomSpec->Write();
  hT0Pulse->Write();
  hTimeDecay->Write();
  hDetector->Write();
  hnTurns->Write();
  hTTT->Write();
  hData->Write();
  gData_ratio->Write();
  hDataFR->Write();
  for(int i = 0; i < ndet; i++){
    hDataCalo[i]->Write();
    gDataCalo_ratio[i]->Write();
    hDataCaloFR[i]->Write();
  }
  outFile->Close();
 
}

