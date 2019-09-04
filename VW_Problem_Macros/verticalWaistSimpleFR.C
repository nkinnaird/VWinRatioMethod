#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

void verticalWaistSimpleFR(){

  // Params
  long long nEntries = 1e9;
  double cycBinWidth = 149.2;
  double binWidth = cycBinWidth;

  double approxMaxTime = 150000; 
  int nBins = int(approxMaxTime/binWidth);
  double histMaxTime = nBins*binWidth;
  double binEdgeShift = 0.;

  // Truth parameters
  double tau = 64440;
  double w_a = 0.00144;  // Only needed for forming ratio
  
  double A_vw = 0.2; // Exagerrated to make effect clearer
  double tau_vw = 250000000;
  // double tau_vw = 2500000;
  double w_vw = 10*w_a; // Even multiple - we shouldn't see anything in the ratio
  double phi_vw = 0.5;


  double half_t_a = TMath::Pi()/w_a;


  cout << "Ratio of Ta/2 to bin width: " << half_t_a << " " << binWidth << " " << half_t_a/binWidth << endl;

  
  // Declare truth function

  // no FR
  // auto truthFunc = new TF1("truthFunc","[0]*exp(-x/[1])*(1+[2]*exp(-x/[3])*cos([4]*x+[5]))",0,histMaxTime);
  // truthFunc->SetParameters(1, tau, A_vw, tau_vw, w_vw, phi_vw);
  // truthFunc->SetParNames("N0", "tau", "A_vw", "tau_vw", "w_vw", "phi_vw");

  // FR
  auto truthFunc = new TF1("truthFunc","[0]*exp(-x/[1])*(1+[2]*exp(-x/[3])*cos([4]*x+[5]))*pow(sin([6]*x+[7]),16)",0,histMaxTime);
  truthFunc->SetParameters(1, tau, A_vw, tau_vw, w_vw, phi_vw, TMath::Pi()/cycBinWidth, 0);
  truthFunc->SetParNames("N0", "tau", "A_vw", "tau_vw", "w_vw", "phi_vw", "w_FR", "phi_FR");

  truthFunc->SetLineColor(2);
  truthFunc->SetLineWidth(1);
  truthFunc->SetNpx(150000);

  // Do some pseudo data trials
  TRandom3* rng = new TRandom3(54321);

  // Normalisation for weighted parts
  double denom = 2+2*cosh(half_t_a/tau);
  double VThresh = 2./denom;
  double UThresh_p = VThresh + exp(-half_t_a/tau)/denom;
  double UThresh_m = UThresh_p + exp(+half_t_a/tau)/denom;
  // std::cout << "Cut offs : " << VThresh << "\t" << UThresh_p << "\t" << UThresh_m << endl;

  // Hit times - all calos together
  auto hData = new TH1D("hData",";Time [ns]; Entries",nBins,0+binEdgeShift,histMaxTime+binEdgeShift);
  auto hData_U = new TH1D("hData_U",";Time [ns]; Entries",nBins,0+binEdgeShift,histMaxTime+binEdgeShift);
  auto hData_V = new TH1D("hData_V",";Time [ns]; Entries",nBins,0+binEdgeShift,histMaxTime+binEdgeShift);
  auto hData_num = new TH1D("hData_num",";Time [ns]; Entries",nBins,0+binEdgeShift,histMaxTime+binEdgeShift);
  auto hData_den = new TH1D("hData_den",";Time [ns]; Entries",nBins,0+binEdgeShift,histMaxTime+binEdgeShift);

  auto hData_FR = new TH1D("hData_FR","Raw Signal;Time [#mus]", 450000, -0.0005, 449999.5);
  auto hData_U_FR = new TH1D("hData_U_FR",";Time [ns]; Entries",450000, -0.0005, 449999.5);
  auto hData_V_FR = new TH1D("hData_V_FR",";Time [ns]; Entries",450000, -0.0005, 449999.5);
  auto hData_num_FR = new TH1D("hData_num_FR",";Time [ns]; Entries",450000, -0.0005, 449999.5);
  auto hData_den_FR = new TH1D("hData_den_FR",";Time [ns]; Entries",450000, -0.0005, 449999.5);


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

    /////////////////////////////////////////////////////////////////////////////
    // Randomize time to get rid of fast rotation
    /////////////////////////////////////////////////////////////////////////////
    // tdec += cycBinWidth*(rng->Uniform() - 0.5);

    /////////////////////////////////////////////////////////////////////////////
    // Fill histograms
    /////////////////////////////////////////////////////////////////////////////

    hData->Fill(tdec);
    hData_FR->Fill(tdec);

    double r = rng->Uniform();
    if(r < VThresh){
      hData_V->Fill(tdec);
      hData_V_FR->Fill(tdec);
    } else if (r < UThresh_p){
      hData_U->Fill(tdec+half_t_a);
      hData_U_FR->Fill(tdec+half_t_a);
    } else {
      hData_U->Fill(tdec-half_t_a);
      hData_U_FR->Fill(tdec-half_t_a);
    }
  }

  // Make ratio fit graph - all calos
  auto gData_ratio = new TGraphErrors();
  gData_ratio->SetName("Ratio");
  hData_num->Add(hData_U,hData_V,-1);
  hData_den->Add(hData_U,hData_V);
  for(int i = 1; i <= hData_num->GetNbinsX(); i++){
    int nPt = gData_ratio->GetN();
    if(hData_den->GetBinContent(i) > 0 && hData_num->GetBinCenter(i) > half_t_a && hData_num->GetBinCenter(i) < histMaxTime - half_t_a){
      double ratio = hData_num->GetBinContent(i)/hData_den->GetBinContent(i);
      double ratioErr = sqrt((1. - ratio*ratio) / hData_den->GetBinContent(i));
      if(abs(ratio) <= 1 && ratioErr != 0){ // both U and V should be non-zero and postive, which this takes care of (they can go negative sometimes when subtracting off pileup leading to a ratio greater than 1)
        gData_ratio->SetPoint(nPt,hData_num->GetBinCenter(i),ratio);
        gData_ratio->SetPointError(nPt,0,ratioErr);
      }
    }
  }

  auto gData_ratio_FR = new TGraphErrors();
  gData_ratio_FR->SetName("Ratio_FR");
  hData_num_FR->Add(hData_U_FR,hData_V_FR,-1);
  hData_den_FR->Add(hData_U_FR,hData_V_FR);
  for(int i = 1; i <= hData_num_FR->GetNbinsX(); i++){
    int nPt = gData_ratio_FR->GetN();
    if(hData_den_FR->GetBinContent(i) > 0 && hData_num_FR->GetBinCenter(i) > half_t_a && hData_num_FR->GetBinCenter(i) < histMaxTime - half_t_a){
      double ratio = hData_num_FR->GetBinContent(i)/hData_den_FR->GetBinContent(i);
      double ratioErr = sqrt((1. - ratio*ratio) / hData_den_FR->GetBinContent(i));
      if(abs(ratio) <= 1 && ratioErr != 0){ // both U and V should be non-zero and postive, which this takes care of (they can go negative sometimes when subtracting off pileup leading to a ratio greater than 1)
        gData_ratio_FR->SetPoint(nPt,hData_num_FR->GetBinCenter(i),ratio);
        gData_ratio_FR->SetPointError(nPt,0,ratioErr);
      }
    }
  }

  // Write all the plots to an output file
  TFile* outFile = new TFile("verticalWaistSimpleFR.root","RECREATE");
  hData->Write();
  hData_U->Write();
  hData_V->Write();
  hData_num->Write();
  hData_den->Write();
  gData_ratio->Write();

  hData_FR->Write();
  hData_U_FR->Write();
  hData_V_FR->Write();
  hData_num_FR->Write();
  hData_den_FR->Write();
  gData_ratio_FR->Write();

  outFile->Close();
 
}

