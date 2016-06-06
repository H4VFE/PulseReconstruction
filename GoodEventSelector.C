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
	Int_t nbins = 800;
    Int_t count = 0;
	TFile *input = new TFile("FilterRMSComparison.root");
    TFile *templatefile = new TFile("/home/marko/Desktop/H4Analysis/ntuples/Templates_APDs.root");
	TTree *MyTree = (TTree*) input->Get("RMS");
    TCanvas *can1 = new TCanvas("can1", "canvas", 1200,1200);
    TCanvas *can2 = new TCanvas("can2", "canvas", 1200,1200);
    TCanvas *can3 = new TCanvas("can3", "canvas", 1200,1200);
    can1->Divide(1,3);
    can2->Divide(1,2);
    can3->Divide(1,2);
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
	TString plot, plot2, cut;
	char name[50];
    TH1F *histoave = new TH1F("histoave","Wave Pulse Average", nbins, -42, 118);
    TH1F *originaltemplate = (TH1F*) templatefile->Get("APD2_E50_G50_prof");
    originaltemplate->Rebin(16);
    TH1F *templatehisto = new TH1F("templatehisto", "Template = Green, Average = Red", 16000*800/1024, -42, 118); //12500 bins
    for (Int_t i=0;i<16000*800/1024;i++) {
        templatehisto->SetBinContent(i+1, originaltemplate->GetBinContent(i+1));
    }
    templatehisto->SetLineColor(kGreen+3);
    templatehisto->SetLineWidth(4);
    can2->cd(2);
    histoave->GetYaxis()->SetRangeUser(-.120,1.2);
    //templatehisto->Draw();
    histoave->Draw();
	for (Int_t i=0;i<nevents;i++) {
    //for (Int_t i=0;i<10;i++) {
		TString histoname = "TempHisto_";
        histoname += i;
        TString histoname2 = "TempHisto2_";
        histoname2 += i;
        TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -42, 118, 1000, -120, 800); //nanoseconds
        TH2F* TempHisto2 = new TH2F (histoname2, "Temp Histo", nbins, -42, 118, 1000, -120, 800); //nanoseconds
        TString h1name = "h1001_";
        h1name += i;
        TString h1name2 = "h1002_";
        h1name2 += i;
        TString h1name3 = "h1003_";
        h1name3 += i;
        TString h1name4 = "h1004_";
        h1name4 += i;
        TH1F *h1001 = new TH1F(h1name,"Red = Unfiltered, Blue = Filtered", nbins, -42, 118);
        TH1F *h1002 = new TH1F(h1name2,"h1002", nbins, -42, 118);
        TH1F *h1003 = new TH1F(h1name3,"Filtered WF - Unfiltered WF", nbins, -42, 118);
        TH1F *h1004 = new TH1F(h1name4,"Percent Change in Filtered WF - Unfiltered WF", nbins, -42, 118);
        plot = "unfilteredwfval:(unfilteredwftime-unfilteredtimeref)>>";
        plot += histoname;
        plot2 = "filteredwfval:(filteredwftime-filteredtimeref)>>";
        plot2 += histoname2;
        cut = "abs(unfilteredbslope)<6 && unfilteredampfit>400 && unfilteredevent==";
        cut += vEvent[i];
        MyTree->Draw(plot, cut, "goff");
        TempHisto = (TH2F*) gDirectory->Get(histoname);
        h1001 = transform2Dto1D(TempHisto);
        MyTree->Draw(plot2, cut, "goff");
        TempHisto2 = (TH2F*) gDirectory->Get(histoname2);
        h1002 = transform2Dto1D(TempHisto2);
        sprintf(name, "Good Events/Event%d", vEvent[i]);
        strcat(name, ".png");
        h1001->SetLineColor(kRed);
        h1002->SetLineColor(kBlue);
        h1001->SetStats(0);
        h1002->SetStats(0);
        can1->cd(1);
        h1001->Draw();
        h1002->Draw("same");
        for (Int_t i=0;i<nbins;i++) {
            h1003->SetBinContent(i+1, (h1002->GetBinContent(i+1))-(h1001->GetBinContent(i+1)));
            if (h1001->GetBinContent(i+1) != 0) h1004->SetBinContent(i+1, ((h1002->GetBinContent(i+1))-(h1001->GetBinContent(i+1)))/(h1001->GetBinContent(i+1)));
            else h1004->SetBinContent(i+1, 0);
        }
        can1->cd(2);
        h1003->Draw();
        can1->cd(3);
        h1004->Draw();
        h1004->GetYaxis()->SetRangeUser(-0.1,0.1);
        gPad->SetGrid();
        can1->SaveAs(name);
        histoave->Add(histoave, h1002);
        can2->cd(2);
        h1002->Scale(1./(h1002->GetMaximum()));
        h1002->Draw("same");
        count++;
        delete TempHisto, TempHisto2, h1001, h1002, h1003, h1004, histoname, histoname2, h1name, h1name2, h1name3, h1name4;
        gDirectory->Clear();
	}
    can2->cd(1);
    histoave->Scale(1./count);
    histoave->SetLineColor(kRed);
    histoave->SetLineWidth(4);
    histoave->Scale(1./(histoave->GetMaximum()));
    templatehisto->Draw();
    histoave->Draw("same");
    can2->cd(2);
    //templatehisto->Draw("same");
    histoave->Draw("same");
    can2->SaveAs("AmpSpread.png");
    can2->SaveAs("AmpSpreadRoot.root");
    can3->cd(1);

}