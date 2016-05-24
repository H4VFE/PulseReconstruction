#define Capillary100GeV_cxx
#include "Capillary100GeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/SpecFunc.h>
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TProfile.h"
#include "TVirtualFFT.h"
#include "TString.h"
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

void Capillary100GeV::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  //defining the histograms to be used

  TH1F *TempHisto = new TH1F("TempHisto", "TempHisto", 1000, 0, 200);
  TH1F *FiberEnergySum = new TH1F("FiberEnergySum", "FiberEnergySum", 100, 5, 30);

  TH1F *Fiber1Energy = new TH1F("Fiber1Energy", "Fiber1Energy", 100, 0, 40);
  TH1F *Fiber2Energy = new TH1F("Fiber2Energy", "Fiber2Energy", 100, 0, 150);
  TH1F *Fiber3Energy = new TH1F("Fiber3Energy", "Fiber3Energy", 100, 0, 40);
  TH1F *Fiber4Energy = new TH1F("Fiber4Energy", "Fiber4Energy", 100, 0, 150);

  TH1F *Fiber1ProfileHisto = new TH1F("Fiber1ProfileHisto", "Fiber1ProfileHisto", 100, 5, 30);
  TH1F *Fiber2ProfileHisto = new TH1F("Fiber2ProfileHisto", "Fiber2ProfileHisto", 100, 5, 30);
  TH1F *Fiber3ProfileHisto = new TH1F("Fiber3ProfileHisto", "Fiber3ProfileHisto", 100, 5, 30);
  TH1F *Fiber4ProfileHisto = new TH1F("Fiber4ProfileHisto", "Fiber4ProfileHisto", 100, 5, 30);
  
  char name[20];
  char title[100];

  //defining the array of histograms that will store the data for each bin

  TH1F *Fiber1Prof[101];
  for (Int_t z=0;z<100;z++) {
    sprintf(name,"Fiber1Prof%d",z+1);
    sprintf(title,"bin%d histo", z+1);
    Fiber1Prof[z+1] = new TH1F(name,title,200, 0, 200);
  }

  TH1F *Fiber2Prof[101];
  for (Int_t z=0;z<100;z++) {
    sprintf(name,"Fiber2Prof%d",z+1);
    sprintf(title,"bin%d histo", z+1);
    Fiber2Prof[z+1] = new TH1F(name,title,200, 0, 200);
  }

  TH1F *Fiber3Prof[101];
  for (Int_t z=0;z<100;z++) {
    sprintf(name,"Fiber3Prof%d",z+1);
    sprintf(title,"bin%d histo", z+1);
    Fiber3Prof[z+1] = new TH1F(name,title,200, 0, 200);
  }

  TH1F *Fiber4Prof[101];
  for (Int_t z=0;z<100;z++) {
    sprintf(name,"Fiber4Prof%d",z+1);
    sprintf(title,"bin%d histo", z+1);
    Fiber4Prof[z+1] = new TH1F(name,title,200, 0, 200);
  }

  //defining the canvases

  TCanvas *Energy = new TCanvas("Energy", "Energy",25,145,1210,910);
  Energy->Divide(2,2);

  TCanvas *XProfile = new TCanvas("XProfile", "XProfile",25,145,900,1500);

  TCanvas *Temp = new TCanvas("Temp", "Temp",25,145,900,1500);

  //defining various variables used throughout the code in one place
 
  Long64_t nbytes = 0, nb = 0;
  Float_t Fiber1Mean, Fiber2Mean, Fiber3Mean, Fiber4Mean;
  Double_t rms1, rms2, rms3, rms4, baselineloc1, baseline1, baselineloc2, baseline2, baselineloc3, baseline3, baselineloc4, baseline4;
  Double_t mean, stddev, rms, Xmin, Xmax, mean1, mean2, mean3, mean4, max;
  TString profile, cut;

  Int_t histobin;

  //obtaining the values for the baseline to subtract and center the data at 0

  fChain->Draw("max_first>>TempHisto");
  TempHisto->GetXaxis()->SetRange(0,35);
  baselineloc1 = TempHisto->GetBinCenter(TempHisto->GetMaximumBin());
  stddev = TempHisto->GetStdDev();
  TempHisto->Fit("gaus","q","", baselineloc1-stddev, baselineloc1+stddev);
  baseline1 = TempHisto->GetFunction("gaus")->GetParameter(1);

  fChain->Draw("max_second>>TempHisto");
  TempHisto->GetXaxis()->SetRange(0,35);
  baselineloc2 = TempHisto->GetBinCenter(TempHisto->GetMaximumBin());
  stddev = TempHisto->GetStdDev();
  TempHisto->Fit("gaus","q","", baselineloc2-stddev, baselineloc2+stddev);
  baseline2 = TempHisto->GetFunction("gaus")->GetParameter(1);

  fChain->Draw("max_third>>TempHisto");
  TempHisto->GetXaxis()->SetRange(0,35);
  baselineloc3 = TempHisto->GetBinCenter(TempHisto->GetMaximumBin());
  stddev = TempHisto->GetStdDev();
  TempHisto->Fit("gaus","q","", baselineloc3-stddev, baselineloc3+stddev);
  baseline3 = TempHisto->GetFunction("gaus")->GetParameter(1);

  fChain->Draw("max_fourth>>TempHisto");
  TempHisto->GetXaxis()->SetRange(0,35);
  baselineloc4 = TempHisto->GetBinCenter(TempHisto->GetMaximumBin());
  stddev = TempHisto->GetStdDev();
  TempHisto->Fit("gaus","q","", baselineloc4-stddev, baselineloc4+stddev);
  baseline4 = TempHisto->GetFunction("gaus")->GetParameter(1);

  //looping over the events

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if ((Hodo_y1+Hodo_y2)*0.25<18 && (Hodo_y1+Hodo_y2)*0.25>14 && (Hodo_x1+Hodo_x2)*0.25>15 && (Hodo_x1+Hodo_x2)*0.25<18 && max_fourth>5 && max_first>2.5) {
       
    //drawing the Energy plots

      Energy->cd(1);
      Fiber1Energy->Fill(max_first-baseline1);
      mean = Fiber1Energy->GetMean();
      rms = Fiber1Energy->GetRMS();
      Xmin = mean-rms;
      Xmax = mean+2*rms;
      Fiber1Energy->GetYaxis()->SetTitle("Energy");
      Fiber1Energy->Draw();
  
      Energy->cd(2);
      Fiber2Energy->Fill(max_second-baseline2);
      mean = Fiber2Energy->GetMean();
      rms = Fiber2Energy->GetRMS();
      Xmin = 0.9*mean;               
      Xmax = mean+2*rms;
      Fiber2Energy->GetYaxis()->SetTitle("Energy");
      Fiber2Energy->Draw();
  
      Energy->cd(3);
      Fiber3Energy->Fill(max_third-baseline3);
      mean = Fiber3Energy->GetMean();
      rms = Fiber3Energy->GetRMS();
      Xmin = 0.9*mean;
      Xmax = mean+2*rms;
      Fiber3Energy->GetYaxis()->SetTitle("Energy");
      Fiber3Energy->Draw();
  
      Energy->cd(4);
      Fiber4Energy->Fill(max_fourth-baseline4);
      mean = Fiber4Energy-> GetMean();
      rms = Fiber4Energy-> GetRMS();
      Xmin = 0.9*mean;
      Xmax = mean+2*rms;
      Fiber4Energy->GetYaxis()->SetTitle("Energy");
      Fiber4Energy->Draw();
      }
  }

  //fitting the energy plots with a Gaussian

  Energy->cd(1);
  mean = Fiber1Energy->GetMean();
  stddev = Fiber1Energy->GetStdDev();
  Fiber1Energy->Fit("gaus","q","", mean-stddev , mean+stddev);
  mean1=Fiber1Energy->GetFunction("gaus")->GetParameter(1);

  Energy->cd(2);
  mean = Fiber2Energy->GetMean();
  stddev = Fiber2Energy->GetStdDev();
  Fiber2Energy->Fit("gaus","q","", mean-stddev , mean+stddev);
  mean2=Fiber2Energy->GetFunction("gaus")->GetParameter(1);

  Energy->cd(3);
  mean = Fiber3Energy->GetMean();
  stddev = Fiber3Energy->GetStdDev();
  Fiber3Energy->Fit("gaus","q","", mean-stddev , mean+stddev);
  mean3=Fiber3Energy->GetFunction("gaus")->GetParameter(1);

  Energy->cd(4);
  mean = Fiber4Energy->GetMean();
  stddev = Fiber4Energy->GetStdDev();
  Fiber4Energy->Fit("gaus","q","", mean-stddev , mean+stddev);
  mean4=Fiber4Energy->GetFunction("gaus")->GetParameter(1);
  
  Temp->cd();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);

    //looping over the events again to fill 100 histograms, with each histogram representing 0.25 mm in the x-direction (total 25 mm between 5 mm and 30 mm)

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ((Hodo_y1+Hodo_y2)*0.25>10.0 && (Hodo_y1+Hodo_y2)*0.25<15.0 && NFibers[0]>1 &&  NFibers[1]>1) {
      for (histobin=0;histobin<100;histobin++) {
        if ((Hodo_x1+Hodo_x2)*0.25<(5+((histobin+1)*25./100)) && (Hodo_x1+Hodo_x2)*0.25>(5+(histobin*25./100))) {
          Fiber1Prof[histobin+1]->Fill(max_first-baseline1);
          Fiber2Prof[histobin+1]->Fill(max_second-baseline2);
          Fiber3Prof[histobin+1]->Fill(max_third-baseline3);
          Fiber4Prof[histobin+1]->Fill(max_fourth-baseline4);
        }
      }
    }
  }

  //normalizing the profiles to have the max energy be 1

  for (Int_t s=0;s<100;s++) {
    Fiber1Prof[s+1]->Scale(1/mean1);
    Fiber2Prof[s+1]->Scale(1/mean2);
    Fiber3Prof[s+1]->Scale(1/mean3);
    Fiber4Prof[s+1]->Scale(1/mean4);

    //fitting the histograms

    Fiber1Prof[s+1]->GetXaxis()->SetRange(2,200);
    max = Fiber1Prof[s+1]->GetBinCenter(Fiber1Prof[s+1]->GetMaximumBin());
    Fiber1Prof[s+1]->GetXaxis()->SetRange(0,200);
    stddev = Fiber1Prof[s+1]->GetStdDev();
    Fiber1Prof[s+1]->Fit("gaus","q","", max-2*stddev, max+2*stddev);

    Fiber2Prof[s+1]->GetXaxis()->SetRange(2,200);
    max = Fiber2Prof[s+1]->GetBinCenter(Fiber2Prof[s+1]->GetMaximumBin());
    Fiber2Prof[s+1]->GetXaxis()->SetRange(0,200);
    stddev = Fiber2Prof[s+1]->GetStdDev();
    Fiber2Prof[s+1]->Fit("gaus","q","", max-2*stddev, max+2*stddev);

    Fiber3Prof[s+1]->GetXaxis()->SetRange(2,200);
    max = Fiber3Prof[s+1]->GetBinCenter(Fiber3Prof[s+1]->GetMaximumBin());
    Fiber3Prof[s+1]->GetXaxis()->SetRange(0,200);
    stddev = Fiber3Prof[s+1]->GetStdDev();
    Fiber3Prof[s+1]->Fit("gaus","q","", max-2*stddev, max+2*stddev);

    Fiber4Prof[s+1]->GetXaxis()->SetRange(2,200);
    max = Fiber4Prof[s+1]->GetBinCenter(Fiber4Prof[s+1]->GetMaximumBin());
    Fiber4Prof[s+1]->GetXaxis()->SetRange(0,200);
    stddev = Fiber4Prof[s+1]->GetStdDev();
    Fiber4Prof[s+1]->Fit("gaus","q","", max-2*stddev, max+2*stddev);

    Int_t Entry1 = Fiber1Prof[s+1]->GetEntries();
    Int_t Entry2 = Fiber2Prof[s+1]->GetEntries();
    Int_t Entry3 = Fiber3Prof[s+1]->GetEntries();
    Int_t Entry4 = Fiber4Prof[s+1]->GetEntries();

    //setting cuts on which histograms use the peak value (from the fit) and which use the mean value of the histogram based on number of entries (more entries = better fit = peak value used)

    if (Entry1 >= 50) {
      Fiber1Mean = Fiber1Prof[s+1]->GetFunction("gaus")->GetParameter(1);
      if (Fiber1Mean/mean1 > 1.15 || Fiber1Mean/mean1 < 0) Fiber1ProfileHisto->SetBinContent(s+1, Fiber1Prof[s+1]->GetMean());
      else Fiber1ProfileHisto->SetBinContent(s+1, Fiber1Mean);
    }
    if (Entry1 < 50 && Entry1 > 5) Fiber1ProfileHisto->SetBinContent(s+1, Fiber1Prof[s+1]->GetMean());
    if (Entry1 <= 5) Fiber1ProfileHisto->SetBinContent(s+1, 0);

    if (Entry2 >= 50) {
      Fiber2Mean = Fiber2Prof[s+1]->GetFunction("gaus")->GetParameter(1);
      if (Fiber2Mean/mean2 > 1.15 || Fiber2Mean/mean2 < 0) Fiber2ProfileHisto->SetBinContent(s+1, Fiber2Prof[s+1]->GetMean());
      else Fiber2ProfileHisto->SetBinContent(s+1, Fiber2Mean);
    }
    if (Entry2 < 50 && Entry2 > 5) Fiber2ProfileHisto->SetBinContent(s+1, Fiber2Prof[s+1]->GetMean());
    if (Entry2 <= 5) Fiber2ProfileHisto->SetBinContent(s+1, 0);

    if (Entry3 >= 50) {
      Fiber3Mean = Fiber3Prof[s+1]->GetFunction("gaus")->GetParameter(1);
      if (Fiber3Mean/mean3 > 1.15 || Fiber3Mean/mean3 < 0) Fiber3ProfileHisto->SetBinContent(s+1, Fiber3Prof[s+1]->GetMean());
      else Fiber3ProfileHisto->SetBinContent(s+1, Fiber3Mean);
    }
    if (Entry3 < 50 && Entry3 > 5) Fiber3ProfileHisto->SetBinContent(s+1, Fiber3Prof[s+1]->GetMean());
    if (Entry3 <= 5) Fiber3ProfileHisto->SetBinContent(s+1, 0);

    if (Entry4 >= 50) {
      Fiber4Mean = Fiber4Prof[s+1]->GetFunction("gaus")->GetParameter(1);
      if (Fiber4Mean/mean4 > 1.15 || Fiber4Mean/mean4 < 0) Fiber4ProfileHisto->SetBinContent(s+1, Fiber4Prof[s+1]->GetMean());
      else Fiber4ProfileHisto->SetBinContent(s+1, Fiber4Mean);
    }
    if (Entry4 < 50 && Entry4 > 5) Fiber4ProfileHisto->SetBinContent(s+1, Fiber4Prof[s+1]->GetMean());
    if (Entry4 <= 5) Fiber4ProfileHisto->SetBinContent(s+1, 0);
  }
  
  Fiber1ProfileHisto->SetLineColor(kRed);
  Fiber1ProfileHisto->Scale(1/mean1);
  Fiber2ProfileHisto->GetXaxis()->SetTitle("Position (mm)");
  Fiber2ProfileHisto->GetYaxis()->SetTitle("Energy Distribution");
  cout << "Fiber 1 = Red" << endl;
  Fiber2ProfileHisto->SetLineColor(kOrange);
  Fiber2ProfileHisto->Scale(1/mean2);
  cout << "Fiber 2 = Orange" << endl;
  Fiber3ProfileHisto->SetLineColor(kBlue);
  Fiber3ProfileHisto->Scale(1/mean3);
  cout << "Fiber 3 = Blue" << endl;
  Fiber4ProfileHisto->SetLineColor(kViolet);
  Fiber4ProfileHisto->Scale(1/mean4);
  cout << "Fiber 4 = Purple" << endl;

  XProfile->Divide(1,2);
  XProfile->cd(1);
  Fiber2ProfileHisto->Draw();
  Fiber1ProfileHisto->Draw("same");
  Fiber3ProfileHisto->Draw("same");
  Fiber4ProfileHisto->Draw("same");

  //making a plot that shows the normalized sum of all the fibers

  for (Int_t i=0;i<100;i++) {
    FiberEnergySum->SetBinContent(i+1, (Fiber1ProfileHisto->GetBinContent(i+1) + Fiber2ProfileHisto->GetBinContent(i+1) + Fiber3ProfileHisto->GetBinContent(i+1) + Fiber4ProfileHisto->GetBinContent(i+1))/4);
  }

  XProfile->cd(2);
  FiberEnergySum->GetXaxis()->SetTitle("Position (mm)");
  FiberEnergySum->GetYaxis()->SetTitle("Energy Distribution");
  FiberEnergySum->Draw();

  //writing the histograms to a root file

  TFile f("Capillary100GeVXProfile.root", "recreate");
  Fiber1ProfileHisto->Write();
  Fiber2ProfileHisto->Write();
  Fiber3ProfileHisto->Write();
  Fiber4ProfileHisto->Write();
  FiberEnergySum->Write();

  for (Int_t i=0;i<100;i++) {
    Fiber1Prof[i+1]->Write();
    Fiber2Prof[i+1]->Write();
    Fiber3Prof[i+1]->Write();
    Fiber4Prof[i+1]->Write();
  }
  for (Int_t i=0;i<100;i++) {
    delete Fiber1Prof[i+1];
    delete Fiber2Prof[i+1];
    delete Fiber3Prof[i+1];
    delete Fiber4Prof[i+1];
  }
  delete Temp;
}