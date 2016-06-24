#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TVirtualFFT.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
using namespace std;

TH1F* transform2Dto1D (TH2F *h2) {
    Int_t ybins = h2->GetNbinsY();
    Int_t xbins = h2->GetNbinsX();
    Float_t ymax = h2->GetYaxis()->GetBinUpEdge(h2->GetNbinsY());
    Float_t ymin = h2->GetYaxis()->GetBinLowEdge(1);
    Float_t xmax = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());
    Float_t xmin = h2->GetXaxis()->GetBinLowEdge(1);
    Float_t yrange = ymax - ymin;
    TString name = h2->GetName();
    name += "_Transformed";
    TH1F *h1 = new TH1F (name, "Pre- & Post-Filter Wave Pulses", xbins, xmin, xmax);
    for (Int_t xbin=0; xbin<xbins; xbin++) {
        for (Int_t ybin=0; ybin<ybins; ybin++) {
            if (h2->GetBinContent(xbin+1, ybin) != 0) {
                h1->SetBinContent(xbin+1,ybin*(yrange/ybins)+ymin);
            }
        }
    }
    return h1;
}

int main() {
	Int_t nbins = 800, count = 0;
	TFile *input = new TFile("FilterRMSComparison.root", "read");
    TFile *input2 = new TFile("Alignment.root", "read");
    TFile *templatefile = new TFile("/home/marko/Desktop/H4Analysis/ntuples/Templates_APDs.old.root", "read");
    TFile *newtemplatefile = new TFile("/home/marko/Desktop/H4Analysis/ntuples/Templates_APDs.root", "recreate");
	TTree *MyTree = (TTree*) input->Get("RMS");    
    MyTree->SetEntryList(0);
    TString listcut = "abs(unfilteredbslope)<6 && unfilteredampfit>400";
    MyTree->Draw(">>myList", listcut, "entrylist");
    TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    MyTree->SetEntryList(myList);
    Int_t nevents = myList->GetN();
	MyTree->Draw("unfilteredevent", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "goff");
    Double_t *vTemp = MyTree->GetV1();
    Int_t *vEvent = new Int_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vEvent[iEntry] = vTemp[iEntry];
    }
	//TString plot, plot2, cut;
	char name[50];
    TH1F *originaltemplate = (TH1F*) templatefile->Get("APD2_E50_G50_prof");
    TH1F *histoave = (TH1F*) input2->Get("histoave");
    TH1F *extendedave = new TH1F("extendedave", "Extended Ave Histo", 1000, -40, 160);
    TH1F *newtemplate = new TH1F("APD2_E50_G50_prof", "APD2_E50_G50_prof", 16000, -40, 160);
    for (int i=0;i<nbins;i++) {
        extendedave->SetBinContent(i+1, histoave->GetBinContent(i+1));
    }
    for (int i=nbins;i<1000;i++) {
        extendedave->SetBinContent(i+1, histoave->GetBinContent(nbins));
    }
    extendedave->Draw();
    Float_t x, y, x1, x2, y1, y2;
    for (int i=0; i<(1000)*16;i++) {
        if (i<=8) newtemplate->SetBinContent(i+1, extendedave->GetBinContent(1));
        else if (15992<=i) newtemplate->SetBinContent(i+1, extendedave->GetBinContent(1000));
        else { 
            x1 = extendedave->GetBinCenter(i/16+1);
            x2 = extendedave->GetBinCenter(i/16+2);
            y1 = extendedave->GetBinContent(i/16+1);
            y2 = extendedave->GetBinContent(i/16+2);
            x = newtemplate->GetBinCenter(i+1);
            y = (x-x1)*(y2-y1)/(x2-x1) + y1;
            newtemplate->SetBinContent(i+1, y);
        }
    }

    TH2F *APD2_E200_G50 = (TH2F*) templatefile->Get("APD2_E200_G50");     
    TH1F *APD2_E200_G50_prof = (TH1F*) templatefile->Get("APD2_E200_G50_prof");         
    TH2F *APD2_E150_G50 = (TH2F*) templatefile->Get("APD2_E150_G50");     
    TH1F *APD2_E150_G50_prof = (TH1F*) templatefile->Get("APD2_E150_G50_prof");
    TH2F *APD2_E100_G50 = (TH2F*) templatefile->Get("APD2_E100_G50");     
    TH1F *APD2_E100_G50_prof = (TH1F*) templatefile->Get("APD2_E100_G50_prof");
    TH2F *APD2_E50_G50 = (TH2F*) templatefile->Get("APD2_E50_G50");
    
    
    TH1F *APD2_E50_G50_prof = newtemplate; 
    
    
    TH2F *APD2_E50_G100 = (TH2F*) templatefile->Get("APD2_E50_G100");     
    TH1F *APD2_E50_G100_prof = (TH1F*) templatefile->Get("APD2_E50_G100_prof");
    TH2F *APD2_E50_G200 = (TH2F*) templatefile->Get("APD2_E50_G200");     
    TH1F *APD2_E50_G200_prof = (TH1F*) templatefile->Get("APD2_E50_G200_prof");
    TH2F *APD1_E50_G50 = (TH2F*) templatefile->Get("APD1_E50_G50");      
    TH1F *APD1_E50_G50_prof = (TH1F*) templatefile->Get("APD1_E50_G50_prof"); 
    TH2F *APD1_E50_G100 = (TH2F*) templatefile->Get("APD1_E50_G100");     
    TH1F *APD1_E50_G100_prof = (TH1F*) templatefile->Get("APD1_E50_G100_prof");
    TH2F *APD1_E50_G200 = (TH2F*) templatefile->Get("APD1_E50_G200");     
    TH1F *APD1_E50_G200_prof = (TH1F*) templatefile->Get("APD1_E50_G200_prof");
    TH2F *APD1_E100_G50 = (TH2F*) templatefile->Get("APD1_E100_G50");     
    TH1F *APD1_E100_G50_prof = (TH1F*) templatefile->Get("APD1_E100_G50_prof");
    TH2F *APD1_E150_G50 = (TH2F*) templatefile->Get("APD1_E150_G50");     
    TH1F *APD1_E150_G50_prof = (TH1F*) templatefile->Get("APD1_E150_G50_prof");
    TH2F *APD1_E200_G50 = (TH2F*) templatefile->Get("APD1_E200_G50");     
    TH1F *APD1_E200_G50_prof = (TH1F*) templatefile->Get("APD1_E200_G50_prof");

    newtemplatefile->cd();

    APD2_E150_G50->Write();
    APD2_E150_G50_prof->Write();
    APD2_E100_G50->Write();
    APD2_E100_G50_prof->Write();
    APD2_E50_G50->Write();
    APD2_E50_G50_prof->Write();
    APD2_E50_G100->Write();
    APD2_E50_G100_prof->Write();
    APD2_E50_G200->Write();
    APD2_E50_G200_prof->Write();
    APD1_E50_G50->Write();
    APD1_E50_G50_prof->Write();
    APD1_E50_G100->Write();
    APD1_E50_G100_prof->Write();
    APD1_E50_G200->Write();
    APD1_E50_G200_prof->Write();
    APD1_E100_G50->Write();
    APD1_E100_G50_prof->Write();
    APD1_E150_G50->Write();
    APD1_E150_G50_prof->Write();
    APD1_E200_G50->Write();
    APD1_E200_G50_prof->Write();

    newtemplate->SetLineWidth(4);
    newtemplate->SetLineColor(kGreen);
    extendedave->SetLineWidth(4);
    extendedave->SetLineColor(kRed);
    newtemplate->Draw();
    extendedave->Draw("same");
    gPad->SaveAs("NewTemplate.root");
    gPad->SaveAs("NewTemplate.png");

    input->Close();
    input2->Close();
    templatefile->Close();
    newtemplatefile->Close();
}