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
    Double_t integral;
	TFile *input = new TFile("FilterRMSComparison.root");
    TFile *templatefile = new TFile("/home/marko/Desktop/H4Analysis/ntuples/Templates_APDs.old.root");
	TTree *MyTree = (TTree*) input->Get("RMS");
    TCanvas *can1 = new TCanvas("can1", "canvas", 1200,900);
    TCanvas *can2 = new TCanvas("can2", "canvas", 1200,1200);
    TCanvas *can3 = new TCanvas("can3", "canvas", 1200,1200);
    can1->Divide(1,3);
    can2->Divide(1,2);
    can3->Divide(1,2);
    MyTree->SetEntryList(0);
    TString listcut = "abs(unfilteredbslope)<6 && unfilteredampfit>700";
    MyTree->Draw(">>myList", listcut, "entrylist");
    TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    MyTree->SetEntryList(myList);
    Int_t nevents = myList->GetN();
	MyTree->Draw("unfilteredevent", "abs(unfilteredbslope)<6 && unfilteredampfit>700", "goff");
    Double_t *vTemp = MyTree->GetV1();
    Int_t *vEvent = new Int_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vEvent[iEntry] = vTemp[iEntry];
    }
	TString plot, plot2, cut;
	char name[50];
    TH1F *histoave = new TH1F("histoave","Wave Pulse Average", nbins, -40, 120);
    TH1F *histoavefft = new TH1F("histoavefft","Wave Pulse Average FFT", nbins, 0, 5);
    TH1F *histoaveph = new TH1F("histoaveph","Wave Pulse Average Phase", nbins, 0, 800);
    TH1F *originaltemplate = (TH1F*) templatefile->Get("APD2_E50_G50_prof");
    originaltemplate->Rebin(16);
    TH1F *templatehisto = new TH1F("templatehisto", "Template = Green, Average = Red", nbins, -40, 120); 
    TH1F *difference = new TH1F("difference","Template vs. Average Difference", nbins, -40, 120);
    TH1F *percentdifference = new TH1F("percentdifference","Template vs. Average Percent Difference", nbins, -40, 120);
    for (Int_t i=0;i<nbins;i++) {
        templatehisto->SetBinContent(i+1, originaltemplate->GetBinContent(i+1));
    }
    templatehisto->SetLineColor(kGreen+3);
    templatehisto->SetLineWidth(4);
    can2->cd(2);
    histoave->GetYaxis()->SetRangeUser(-.120,1.2);
    //templatehisto->Draw();
    histoave->SetStats(0);
    histoave->Draw();
	for (Int_t i=0;i<nevents;i++) {
    //for (Int_t i=0;i<10;i++) {
        if (vEvent[i]==513) continue; //event for which electronics die out for a bit halfway through
		TString histoname = "TempHisto_";
        histoname += i;
        TString histoname2 = "TempHisto2_";
        histoname2 += i;
        TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -40, 120, 1000, -120, 800); //nanoseconds
        TH2F* TempHisto2 = new TH2F (histoname2, "Temp Histo", nbins, -40, 120, 1000, -120, 800); //nanoseconds
        TString h1name = "h1001_";
        h1name += i;
        TString h1name2 = "h1002_";
        h1name2 += i;
        TString h1name2fft = "h1002fft_";
        h1name2 += i;
        TString h1name2ph = "h1002ph_";
        h1name2 += i;
        TString h1name3 = "h1003_";
        h1name3 += i;
        TString h1name4 = "h1004_";
        h1name4 += i;
        TH1F *h1001 = new TH1F(h1name,"Red = Unfiltered, Blue = Filtered", nbins, -40, 120);
        TH1F *h1002 = new TH1F(h1name2,"h1002", nbins, -40, 120);
        TH1F *h1002fft = new TH1F(h1name2fft,"h1002fft", nbins, 0, 5);
        TH1F *h1002ph = new TH1F(h1name2ph,"h1002ph", nbins, 0, 800);
        TH1F *h1003 = new TH1F(h1name3,"Filtered WF - Unfiltered WF", nbins, -40, 120);
        TH1F *h1004 = new TH1F(h1name4,"Percent Change in Filtered WF - Unfiltered WF", nbins, -40, 120);
        plot = "unfilteredwfval:(unfilteredwftime-unfilteredtimeref)>>";
        plot += histoname;
        plot2 = "filteredwfval:(filteredwftime-filteredtimeref)>>";
        plot2 += histoname2;
        cut = "abs(unfilteredbslope)<6 && unfilteredampfit>700 && unfilteredevent==";
        cut += vEvent[i];
        MyTree->Draw(plot, cut, "goff");
        TempHisto = (TH2F*) gDirectory->Get(histoname);
        h1001 = transform2Dto1D(TempHisto);
        MyTree->Draw(plot2, cut, "goff");
        TempHisto2 = (TH2F*) gDirectory->Get(histoname2);
        h1002 = transform2Dto1D(TempHisto2);
        if (h1002->GetBinCenter(h1002->GetMaximumBin()) < 30) continue;
        sprintf(name, "Good Events/Event%d", vEvent[i]);
        strcat(name, ".png");
        h1001->SetLineColor(kRed);
        h1002->SetLineColor(kBlue);
        h1001->SetStats(0);
        h1002->SetStats(0);
        //can1->cd(1);
        can1->cd();
        h1001->GetXaxis()->SetTitle("Time (ns)");
        h1001->GetYaxis()->SetTitle("Amplitude");
        h1001->Draw();
        h1002->Draw("same");
        //for (Int_t i=0;i<nbins;i++) {
        //    h1003->SetBinContent(i+1, (h1002->GetBinContent(i+1))-(h1001->GetBinContent(i+1)));
        //    if (h1001->GetBinContent(i+1) != 0) h1004->SetBinContent(i+1, ((h1002->GetBinContent(i+1))-(h1001->GetBinContent(i+1)))/(h1001->GetBinContent(i+1)));
        //    else h1004->SetBinContent(i+1, 0);
        //}
        //can1->cd(2);
        //h1003->Draw();
        //can1->cd(3);
        //h1004->Draw();
        //h1004->GetYaxis()->SetRangeUser(-0.1,0.1);
        //gPad->SetGrid();
        can1->SaveAs(name);
        //h1002->FFT(h1002fft, "MAG");
        //h1002->FFT(h1002ph, "PH");
        //histoavefft->Add(histoavefft, h1002fft);
        //histoaveph->Add(histoaveph, h1002ph);
        histoave->Add(histoave, h1002);
        can2->cd(2);
        h1002->Scale(1./(h1002->GetMaximum()));
        h1002->Draw("same");
        count++;
        delete TempHisto, TempHisto2, h1001, h1002, h1003, h1004, histoname, histoname2, h1name, h1name2, h1name3, h1name4;
        gDirectory->Clear();
	}
    can2->cd(1);
    //histoavefft->Scale(1./count);
    //histoaveph->Scale(1./count);
    histoave->Scale(1./count);
    
    //Double_t *re_full = new Double_t[nbins];
    //Double_t *im_full = new Double_t[nbins];
    //TH1 *Throwaway = 0;
    //TH1F *invhistoave = new TH1F(invhistoave, "Average Pulse", nbins, -40, 120);
    //TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
    //for (Int_t n=0; n<nbins; n++) {
    //    (re_full)[n]=(histoavefft->GetBinContent(n+1)*cos(histoaveph->GetBinContent(n+1)));
    //    (im_full)[n]=(histoavefft->GetBinContent(n+1)*sin(histoaveph->GetBinContent(n+1)));
    //}
    //invFFT->SetPointsComplex(re_full, im_full);
    //invFFT->Transform();
    //Throwaway = TH1::TransformHisto(invFFT, Throwaway, "Re");
    //for (Int_t p=0; p<nbins; p++) {
    //    histoave->SetBinContent(p+1, Throwaway->GetBinContent(p+1)/nbins);
    //}
    histoave->Scale(1./(histoave->GetMaximum()));
    histoave->SetLineColor(kRed);
    histoave->SetLineWidth(4);
    integral = templatehisto->Integral();
    templatehisto->Scale(1./integral);
    templatehisto->SetStats(0);
    templatehisto->Draw();
    integral = histoave->Integral();
    histoave->Scale(1./integral);
    histoave->DrawClone("same");
    histoave->Scale(integral);
    can2->cd(2);
    histoave->GetXaxis()->SetTitle("Time (ns)");
    histoave->GetYaxis()->SetTitle("Normalized Amplitude");
    histoave->DrawClone("same");
    can2->SaveAs("AmpSpread.png");
    can2->SaveAs("AmpSpreadRoot.root");
    TFile *output = new TFile("Alignment.root", "recreate");
    output->cd();
    originaltemplate->Write();
    histoave->Write();
    templatehisto->Write();
    output->Close();
    can3->cd(1);
    histoave->Scale(1./integral);
    for (int i=0;i<nbins;i++) {
        difference->SetBinContent(i+1, histoave->GetBinContent(i+1) - templatehisto->GetBinContent(i+1));
        if (templatehisto->GetBinContent(i+1) != 0) percentdifference->SetBinContent(i+1, (histoave->GetBinContent(i+1) - templatehisto->GetBinContent(i+1))/templatehisto->GetBinContent(i+1));
        else percentdifference->SetBinContent(i+1, 0);
    }
    difference->GetXaxis()->SetTitle("Time (ns)");
    difference->GetYaxis()->SetTitle("Average - Template");
    difference->SetStats(0);
    difference->Draw();
    can3->cd(2);
    percentdifference->GetXaxis()->SetTitle("Time (ns)");
    percentdifference->GetYaxis()->SetTitle("Percent Difference Average - Template");
    percentdifference->SetStats(0);
    percentdifference->Draw();
    percentdifference->GetYaxis()->SetRangeUser(-0.1,0.1);
    gPad->SetGrid();
    can3->SaveAs("Difference.png");
    can3->SaveAs("Difference.root");
}