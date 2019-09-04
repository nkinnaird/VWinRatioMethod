#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TStyle.h>
#include <TLegend.h>

class BinIntegralFunc {

public:
  BinIntegralFunc(TF1* inputFunc, double inputBinWidth)
    : f(inputFunc)
    , binWidth(inputBinWidth)
  {}

  double Evaluate(double *x, double *p) {
    f->SetParameters(p);
    return f->Eval( x[0] );
    //return f->Integral( x[0] - 0.5*binWidth, x[0] + 0.5*binWidth) / binWidth;
  }

private:
  TF1* f;
  double binWidth;
};

double ratioForm(double* x, double* p){

  double tp = TMath::Pi()/0.00144;

  double f_0 = p[0]*exp(-x[0]/p[1])*(1+p[2]*cos(p[3]*x[0]+p[4]))*(1+p[5]*exp(-x[0]/p[6])*cos(p[7]*x[0]+p[8]));
  double f_p = p[0]*exp(-(x[0]+tp)/p[1])*(1+p[2]*cos(p[3]*(x[0]+tp)+p[4]))*(1+p[5]*exp(-(x[0]+tp)/p[6])*cos(p[7]*(x[0]+tp)+p[8]));
  double f_m = p[0]*exp(-(x[0]-tp)/p[1])*(1+p[2]*cos(p[3]*(x[0]-tp)+p[4]))*(1+p[5]*exp(-(x[0]-tp)/p[6])*cos(p[7]*(x[0]-tp)+p[8]));

  return (2*f_0 - f_p - f_m)/(2*f_0 + f_p + f_m);
}

void fitCaloByCalo(){


  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0000);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);


  // Params
  double maxTime = 650064.4;
  double binWidth = 149.2;
  int nBins = maxTime/binWidth;
  const int ndet = 24;

  // Truth parameters
  double tau = 64440;
  double A = 0.3714;
  double w_a = 0.00144;
  double phi = 2.081;
  
  double A_vw = 0.13; // Effect is damped compared to the generated value of 0.2 (not 100% sure why)
  double tau_vw = 250800;
  double w_vw = 8.913*w_a; //10.06*w_a; // 1.6*w_a;
  // double w_vw = 10.0*w_a; //10.06*w_a; // 1.6*w_a;
  double phi_vw = 0;

  double half_t_a = TMath::Pi()/w_a;

  // Open input  
  TFile* inFile = new TFile("verticalWaistHists_1e9.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_5e9.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_NoFRRandomisation_1e9.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_NoFR_1e9.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_60h_1e9.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_25us.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_magicP.root","READ");
  // TFile* inFile = new TFile("verticalWaistHists_magicPmeanT0.root","READ");

  // Get histograms
  TH1D* hData = (TH1D*)inFile->Get("hData");
  TH1D* hDataCalo[ndet];
  for(int idet = 0; idet < ndet; idet++){
    hDataCalo[idet] = (TH1D*)inFile->Get(Form("hDataCalo_%02d",idet+1));
  }

  // T-Method Fit
  auto fitFuncTMeth = new TF1("fitFuncTMeth","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))*(1+[5]*exp(-x/[6])*cos([7]*x+[8]))", half_t_a, 300000);
  fitFuncTMeth->SetParameters(1, tau, A, w_a, phi, A_vw, tau_vw, w_vw, phi_vw);
  double N0 = (binWidth*hData->GetEntries())/fitFuncTMeth->Integral(half_t_a, 300000);
  fitFuncTMeth->SetParameter(0,N0);
  fitFuncTMeth->SetParameter(5,A_vw*0.1);
  fitFuncTMeth->SetParNames("N0", "tau", "A", "w_a", "phi", "A_vw", "tau_vw", "w_vw", "phi_vw");
  fitFuncTMeth->SetNpx(10000);

  // Fit function with integral over bin width
  auto binIntegralFuncTMeth = new BinIntegralFunc(fitFuncTMeth,binWidth);
  auto fitFuncTMethInt = new TF1("f", binIntegralFuncTMeth, &BinIntegralFunc::Evaluate, fitFuncTMeth->GetXmin(), fitFuncTMeth->GetXmax(), fitFuncTMeth->GetNpar()); 
  fitFuncTMethInt->SetParameters(fitFuncTMeth->GetParameters());
  fitFuncTMethInt->FixParameter(2,0);
  fitFuncTMethInt->FixParameter(3,0);
  fitFuncTMethInt->FixParameter(4,0);
  fitFuncTMethInt->FixParameter(6,tau_vw);
  fitFuncTMethInt->FixParameter(7,w_vw);
  fitFuncTMethInt->SetNpx(10000);
  TFitResultPtr fitResultTMeth = hData->Fit(fitFuncTMethInt,"RSNQ");
  fitResultTMeth->Print();
  
  // Get results from integral fit and set fitFuncTMeth parameter errors based on these
  fitFuncTMeth->SetParErrors(fitFuncTMethInt->GetParErrors());

  cout << "TMETHOD:" << endl;
  for(int iPar = 0; iPar < 9; iPar++){
    cout << fitFuncTMeth->GetParName(iPar) << ":\t" << fitFuncTMeth->GetParameter(iPar) << " +/- " << fitFuncTMeth->GetParError(iPar) << endl;
  }
  cout << "PValue : " << TMath::Prob(fitResultTMeth->Chi2(),fitResultTMeth->Ndf()) << endl;
  
/////////////////////////////////////////////////////////////////////////////////////
  // T-Method Fit - per calo
  TGraphErrors* tgTMethFitResult[9];
  TGraphErrors* tgDiffFitResult[9];
  double N0Calo[ndet];
  for(int iPar = 0; iPar < 9; iPar++) {
    tgTMethFitResult[iPar] = new TGraphErrors();
    tgTMethFitResult[iPar]->SetTitle(Form(";Calo Num;%s",fitFuncTMeth->GetParName(iPar)));
    tgDiffFitResult[iPar] = new TGraphErrors();
    tgDiffFitResult[iPar]->SetTitle(Form(";Calo Num;#Delta %s (R-T)",fitFuncTMeth->GetParName(iPar)));
  }
  for(int iCalo = 0; iCalo < ndet; iCalo++){
    TF1* fitFuncTMethCalo = new TF1("fitFuncTMethCalo","[0]*exp(-x/[1])*(1+[2]*cos([3]*x+[4]))*(1+[5]*exp(-x/[6])*cos([7]*x+[8]))", half_t_a, 300000);
    fitFuncTMethCalo->SetParameters(1, tau, A, w_a, TMath::TwoPi()*(1.-float(iCalo)/24), A_vw, tau_vw, w_vw, phi_vw);
    N0Calo[iCalo] = (binWidth*hDataCalo[iCalo]->GetEntries())/fitFuncTMethCalo->Integral(half_t_a, 300000);
    fitFuncTMethCalo->SetParameter(0,N0Calo[iCalo]);
    fitFuncTMethCalo->FixParameter(2,0);
    fitFuncTMethCalo->FixParameter(3,0);
    fitFuncTMethCalo->FixParameter(4,0);
    fitFuncTMethCalo->SetParLimits(5,0,5*A_vw);
    // fitFuncTMethCalo->SetParLimits(6,0.5*tau_vw,2*tau_vw);
    fitFuncTMethCalo->FixParameter(6,tau_vw);
    fitFuncTMethCalo->FixParameter(7,w_vw);
    TFitResultPtr fitResultTMethCalo = hDataCalo[iCalo]->Fit(fitFuncTMethCalo,"RSNQ");
    for(int iPar = 0; iPar < 9; iPar++){
      if(iPar == 4 or iPar == 8){
	double phiTmp = fitFuncTMethCalo->GetParameter(iPar);
	while (phiTmp > TMath::TwoPi()) phiTmp -= TMath::TwoPi();
	while (phiTmp < 0)              phiTmp += TMath::TwoPi();
	tgTMethFitResult[iPar]->SetPoint(tgTMethFitResult[iPar]->GetN(), iCalo+1, phiTmp);
	tgTMethFitResult[iPar]->SetPointError(tgTMethFitResult[iPar]->GetN()-1, 0, fitFuncTMethCalo->GetParError(iPar));
      } else {
	tgTMethFitResult[iPar]->SetPoint(tgTMethFitResult[iPar]->GetN(), iCalo+1, fitFuncTMethCalo->GetParameter(iPar));
	tgTMethFitResult[iPar]->SetPointError(tgTMethFitResult[iPar]->GetN()-1, 0, fitFuncTMethCalo->GetParError(iPar));
      }
    }
    if(int(fitResultTMethCalo) > 0){
      cout << "FAILED FIT: TMETHOD iCalo = " << iCalo << endl;
      hDataCalo[iCalo]->Draw();
      fitFuncTMethCalo->SetLineColor(2);
      fitFuncTMethCalo->Draw("SAME");
      for(int iPar = 0; iPar < 9; iPar++){
       	cout << fitFuncTMethCalo->GetParName(iPar) << ":\t" << fitFuncTMethCalo->GetParameter(iPar) << " +/- " << fitFuncTMethCalo->GetParError(iPar) << endl;
      }
      cout << "PValue : " << TMath::Prob(fitResultTMethCalo->Chi2(),fitResultTMethCalo->Ndf()) << endl;
      return;
    }
    delete fitFuncTMethCalo;
  }
  for(int iPar = 0; iPar < 9; iPar++){
    if(iPar == 4 or iPar == 8){
      double phiTmp = fitFuncTMeth->GetParameter(iPar);
      while (phiTmp > TMath::TwoPi()) phiTmp -= TMath::TwoPi();
      while (phiTmp < 0)              phiTmp += TMath::TwoPi();
      tgTMethFitResult[iPar]->SetPoint(tgTMethFitResult[iPar]->GetN(), 30, phiTmp);
      tgTMethFitResult[iPar]->SetPointError(tgTMethFitResult[iPar]->GetN()-1, 0, fitFuncTMeth->GetParError(iPar));
    } else {
      tgTMethFitResult[iPar]->SetPoint(tgTMethFitResult[iPar]->GetN(), 30, fitFuncTMeth->GetParameter(iPar));
      tgTMethFitResult[iPar]->SetPointError(tgTMethFitResult[iPar]->GetN()-1, 0, fitFuncTMeth->GetParError(iPar));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Get histograms
  TGraphErrors* gData_ratio = (TGraphErrors*)inFile->Get("Ratio");
  TGraphErrors* gDataCalo_ratio[ndet];
  for(int idet = 0; idet < ndet; idet++){
    gDataCalo_ratio[idet] = (TGraphErrors*)inFile->Get(Form("RatioCalo_%02d",idet+1));
  }

  // R-Method Fit 
  auto fitFuncRMeth = new TF1("fitFuncRMeth",ratioForm, half_t_a, 300000, 9);
  fitFuncRMeth->SetParameters(N0, tau, A, w_a, phi, A_vw, tau_vw, w_vw, phi_vw);
  fitFuncRMeth->SetParNames("N0", "tau", "A", "w_a", "phi", "A_vw", "tau_vw", "w_vw", "phi_vw");
  fitFuncRMeth->SetNpx(10000);

  // Fit function with integral over bin width
  auto binIntegralFuncRMeth = new BinIntegralFunc(fitFuncRMeth,binWidth);
  auto fitFuncRMethInt = new TF1("f", binIntegralFuncRMeth, &BinIntegralFunc::Evaluate, fitFuncRMeth->GetXmin(), fitFuncRMeth->GetXmax(), fitFuncRMeth->GetNpar()); 
  fitFuncRMethInt->SetParameters(fitFuncRMeth->GetParameters());
  fitFuncRMethInt->FixParameter(0, fitFuncRMethInt->GetParameter(0));
  fitFuncRMethInt->FixParameter(1, fitFuncRMethInt->GetParameter(1));
  fitFuncRMethInt->FixParameter(2, 0);
  fitFuncRMethInt->FixParameter(3, 0);
  fitFuncRMethInt->FixParameter(4, 0);
  fitFuncRMethInt->SetParLimits(5,0,2*A_vw);
  // fitFuncRMethInt->SetParLimits(6,0.5*tau_vw,2*tau_vw);
  fitFuncRMethInt->FixParameter(6, tau_vw);
  fitFuncRMethInt->FixParameter(7, w_vw);
  fitFuncRMethInt->SetNpx(10000);
  TFitResultPtr fitResultRMeth = gData_ratio->Fit(fitFuncRMethInt,"RNSQ");
  fitResultRMeth->Print();
  
  // Get results from integral fit and set fitFuncRMeth parameter errors based on these
  fitFuncRMeth->SetParErrors(fitFuncRMethInt->GetParErrors());

  cout << "RMETHOD:" << endl;
  for(int iPar = 0; iPar < 9; iPar++){
    cout << fitFuncRMeth->GetParName(iPar) << ":\t" << fitFuncRMeth->GetParameter(iPar) << " +/- " << fitFuncRMeth->GetParError(iPar) << endl;
  }
  cout << "PValue : " << TMath::Prob(fitResultRMeth->Chi2(),fitResultRMeth->Ndf()) << endl;

  cout << endl; 
  cout << "T-Method - R-Method (% diff):" << endl;
  for(int iPar = 2; iPar < 9; iPar++){
    cout << fitFuncTMeth->GetParName(iPar) << ":\t" << 100*(fitFuncTMeth->GetParameter(iPar) - fitFuncRMeth->GetParameter(iPar))/fitFuncTMeth->GetParameter(iPar) << endl;
  }

  // R-Method Fit - per calo
  TGraphErrors* tgRMethFitResult[9];
  for(int iPar = 0; iPar < 9; iPar++){
    tgRMethFitResult[iPar] = new TGraphErrors();
    tgRMethFitResult[iPar]->SetTitle(Form(";Calo Num;%s",fitFuncRMeth->GetParName(iPar)));
  }
  for(int iCalo = 0; iCalo < ndet; iCalo++){

    auto fitFuncRMethCalo = new TF1("fitFuncRMethCalo",ratioForm, half_t_a, 300000, 9);
    fitFuncRMethCalo->SetParameters(N0Calo[iCalo], tau, A, w_a, TMath::TwoPi()*(1.-float(iCalo)/24), A_vw, tau_vw, w_vw, phi_vw);
    fitFuncRMethCalo->SetParNames("N0", "tau", "A", "w_a", "phi", "A_vw", "tau_vw", "w_vw", "phi_vw");
    fitFuncRMethCalo->SetNpx(10000);

    auto binIntegralFuncRMethCalo = new BinIntegralFunc(fitFuncRMethCalo,binWidth);
    auto fitFuncRMethCaloInt = new TF1("f", binIntegralFuncRMethCalo, &BinIntegralFunc::Evaluate, fitFuncRMethCalo->GetXmin(), fitFuncRMethCalo->GetXmax(), fitFuncRMethCalo->GetNpar()); 
    fitFuncRMethCaloInt->SetParameters(fitFuncRMethCalo->GetParameters());
    fitFuncRMethCaloInt->FixParameter(0, fitFuncRMethCaloInt->GetParameter(0));
    fitFuncRMethCaloInt->FixParameter(1, fitFuncRMethCaloInt->GetParameter(1));
    fitFuncRMethCaloInt->FixParameter(2, 0);
    fitFuncRMethCaloInt->FixParameter(3, 0);
    fitFuncRMethCaloInt->FixParameter(4, 0);
    fitFuncRMethCaloInt->SetParLimits(5,0,5*A_vw);
    // fitFuncRMethCaloInt->SetParLimits(6,0.5*tau_vw,2*tau_vw);
    fitFuncRMethCaloInt->FixParameter(6, tau_vw);
    fitFuncRMethCaloInt->FixParameter(7, w_vw);
    fitFuncRMethCaloInt->SetNpx(10000);
    gDataCalo_ratio[iCalo]->Fit(fitFuncRMethCaloInt,"RSNQ");
    fitFuncRMethCalo->SetParErrors(fitFuncRMethCaloInt->GetParErrors());
    for(int iPar = 0; iPar < 9; iPar++){
      if(iPar == 4 or iPar == 8){
	double phiTmp = fitFuncRMethCalo->GetParameter(iPar);
	while (phiTmp > TMath::TwoPi()) phiTmp -= TMath::TwoPi();
	while (phiTmp < 0)              phiTmp += TMath::TwoPi();
	tgRMethFitResult[iPar]->SetPoint(tgRMethFitResult[iPar]->GetN(), iCalo+1, phiTmp);
	tgRMethFitResult[iPar]->SetPointError(tgRMethFitResult[iPar]->GetN()-1, 0, fitFuncRMethCalo->GetParError(iPar));
	double tMethValX, tMethValY;
	tgTMethFitResult[iPar]->GetPoint(tgRMethFitResult[iPar]->GetN()-1, tMethValX, tMethValY);
	tgDiffFitResult[iPar]->SetPoint(tgDiffFitResult[iPar]->GetN(), iCalo+1, phiTmp - tMethValY);
      } else {
	tgRMethFitResult[iPar]->SetPoint(tgRMethFitResult[iPar]->GetN(), iCalo+1, fitFuncRMethCalo->GetParameter(iPar));
	tgRMethFitResult[iPar]->SetPointError(tgRMethFitResult[iPar]->GetN()-1, 0, fitFuncRMethCalo->GetParError(iPar));
	double tMethValX, tMethValY;
	tgTMethFitResult[iPar]->GetPoint(tgRMethFitResult[iPar]->GetN()-1, tMethValX, tMethValY);
	tgDiffFitResult[iPar]->SetPoint(tgDiffFitResult[iPar]->GetN(), iCalo+1, fitFuncRMethCalo->GetParameter(iPar) - tMethValY);
      }
    }
    delete fitFuncRMethCalo;
  }
  for(int iPar = 0; iPar < 9; iPar++){
    if(iPar == 4 or iPar == 8){
      double phiTmp = fitFuncRMeth->GetParameter(iPar);
      while (phiTmp > TMath::TwoPi()) phiTmp -= TMath::TwoPi();
      while (phiTmp < 0)              phiTmp += TMath::TwoPi();
      tgRMethFitResult[iPar]->SetPoint(tgRMethFitResult[iPar]->GetN(), 30, phiTmp);
      tgRMethFitResult[iPar]->SetPointError(tgRMethFitResult[iPar]->GetN()-1, 0, fitFuncRMeth->GetParError(iPar));
    } else {
      tgRMethFitResult[iPar]->SetPoint(tgRMethFitResult[iPar]->GetN(), 30, fitFuncRMeth->GetParameter(iPar));
      tgRMethFitResult[iPar]->SetPointError(tgRMethFitResult[iPar]->GetN()-1, 0, fitFuncRMeth->GetParError(iPar));
    }
  }

  // Draw results for complete fit
  auto cTMeth = new TCanvas("cTMeth","cTMeth",200,10,800,600);
  hData->Draw();
  fitFuncTMeth->SetLineColor(2);
  fitFuncTMeth->SetLineWidth(1);
  fitFuncTMeth->Draw("SAME");
  
  auto cRMeth = new TCanvas("cRMeth","cRMeth",200,10,800,600);
  gData_ratio->SetMarkerStyle(2);
  gData_ratio->Draw("AP");
  fitFuncRMeth->SetLineColor(2);
  fitFuncRMeth->SetLineWidth(1);
  fitFuncRMeth->Draw("SAME");

  TCanvas* cAVW = new TCanvas("cAVW","cAVW");
  tgTMethFitResult[5]->SetMarkerStyle(20);
  tgTMethFitResult[5]->GetYaxis()->SetRangeUser(0, 0.16);
  tgTMethFitResult[5]->Draw("AP");
  tgRMethFitResult[5]->SetMarkerStyle(20);
  tgRMethFitResult[5]->SetMarkerColor(2);
  tgRMethFitResult[5]->SetLineColor(2);
  tgRMethFitResult[5]->Draw("PSAME");
  
  auto legend_amp = new TLegend(0.2,0.2,.5,0.4);
  legend_amp->AddEntry(tgTMethFitResult[5], "T method", "p");
  legend_amp->AddEntry(tgRMethFitResult[5], "R method", "p");
  legend_amp->SetBorderSize(0);
  legend_amp->Draw("SAME");

  cAVW->SaveAs("Images/AVW_canv.png");

  TCanvas* cTauVW = new TCanvas("cTauVW","cTauVW");
  tgTMethFitResult[6]->SetMarkerStyle(20);
  tgTMethFitResult[6]->Draw("AP");
  tgRMethFitResult[6]->SetMarkerStyle(20);
  tgRMethFitResult[6]->SetMarkerColor(2);
  tgRMethFitResult[6]->SetLineColor(2);
  tgRMethFitResult[6]->Draw("PSAME");

  TCanvas* cwVW = new TCanvas("cwVW","cwVW");
  tgTMethFitResult[7]->SetMarkerStyle(20);
  tgTMethFitResult[7]->Draw("AP");
  tgRMethFitResult[7]->SetMarkerStyle(20);
  tgRMethFitResult[7]->SetMarkerColor(2);
  tgRMethFitResult[7]->SetLineColor(2);
  tgRMethFitResult[7]->Draw("PSAME");

  TCanvas* cPhiVW = new TCanvas("cPhiVW","cPhiVW");
  tgTMethFitResult[8]->SetMarkerStyle(20);
  tgTMethFitResult[8]->Draw("AP");
  tgRMethFitResult[8]->SetMarkerStyle(20);
  tgRMethFitResult[8]->SetMarkerColor(2);
  tgRMethFitResult[8]->SetLineColor(2);
  tgRMethFitResult[8]->Draw("PSAME");

  auto legend_phi = new TLegend(0.2,0.2,.5,0.4);
  legend_phi->AddEntry(tgTMethFitResult[8], "T method", "p");
  legend_phi->AddEntry(tgRMethFitResult[8], "R method", "p");
  legend_phi->SetBorderSize(0);
  legend_phi->Draw("SAME");

  cPhiVW->SaveAs("Images/PhiVW_canv.png");

  TCanvas* cAVWDiff = new TCanvas("cAVWDiff","cAVWDiff");
  tgDiffFitResult[5]->SetMarkerStyle(20);
  tgDiffFitResult[5]->Draw("AP");
  cAVWDiff->SaveAs("Images/AVW_diff_canv.png");

  TCanvas* cTauVWDiff = new TCanvas("cTauVWDiff","cTauVWDiff");
  tgDiffFitResult[6]->SetMarkerStyle(20);
  tgDiffFitResult[6]->Draw("AP");

  TCanvas* cwVWDiff = new TCanvas("cwVWDiff","cwVWDiff");
  tgDiffFitResult[7]->SetMarkerStyle(20);
  tgDiffFitResult[7]->Draw("AP");

  TCanvas* cPhiVWDiff = new TCanvas("cPhiVWDiff","cPhiVWDiff");
  tgDiffFitResult[8]->SetMarkerStyle(20);
  tgDiffFitResult[8]->Draw("AP");
  cPhiVWDiff->SaveAs("Images/PhiVW_diff_canv.png");
  
}

