#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TVirtualFFT.h"
#include <iostream>
#include <fstream>
using namespace std;

Float_t pedestalrms1D (TH1 *h1) {
    std::vector<Float_t> y;
    for (Int_t q=0;q<150;q++) {
            y.push_back(h1->GetBinContent(q+1));
    }
    Float_t newrms = TMath::RMS(y.begin(), y.end());
    y.clear();
    return newrms;
}

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
    TH1F *h1 = new TH1F (name, "1D Histogram", xbins, xmin, xmax);
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
    Int_t nspill = 1;
    Int_t nbins = 800;
    TFile f("/home/marko/Desktop/Data Files/analysis_4443.root");
    TFile f1("/home/marko/Music/analysis_4443.root");
    TTree* h4unfiltered = (TTree*) f.Get("h4");
    TTree* h4filtered = (TTree*) f1.Get("h4");
    h4unfiltered->SetEntryList(0);
    TString listcut = "WF_ch==2 && spill==1";
    h4unfiltered->Draw(">>myList", listcut, "entrylist");
    TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    h4unfiltered->SetEntryList(myList);
    Int_t nevents = myList->GetN();

    h4unfiltered->Draw("event", "spill==1", "goff");
    Double_t *vTemp = h4unfiltered->GetV1();
    Double_t *vEvent = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vEvent[iEntry] = vTemp[iEntry];
    }

    h4filtered->Draw("event", "spill==1", "goff");
    Double_t *vTemp1 = h4filtered->GetV1();
    Double_t *vfEvent = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfEvent[iEntry] = vTemp1[iEntry];
    }

    h4unfiltered->Draw("b_rms[APD2]", "spill==1", "goff");
    Double_t *vTemp_uBrms = h4unfiltered->GetV1();
    Double_t *vuBrms = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuBrms[iEntry] = vTemp_uBrms[iEntry];
    }

    h4filtered->Draw("b_rms[APD2]", "spill==1", "goff");
    Double_t *vTemp_fBrms = h4filtered->GetV1();
    Double_t *vfBrms = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfBrms[iEntry] = vTemp_fBrms[iEntry];
    }

    Float_t unfilteredrms, filteredrms, unfilteredbrms, filteredbrms;
    Int_t count = 0, counter;
    Float_t sum, sum2;
    //TFile out("FilterRMSComparison.root", "recreate");
    //TTree *MyTree = new TTree ("RMS", "RMS");
    TH1F *h1 = new TH1F("h1","h1", nbins, -0.1, 159.9);

    TString plot, plot2, cut;

    //for (Int_t j=0;j<nevents;j++) {
    for (Int_t j=0;j<1;j++) {
        TString histoname = "TempHisto_";
        histoname += j;
        TString histoname2 = "TempHisto2_";
        histoname2 += j;
        TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
        TH2F* TempHisto2 = new TH2F (histoname2, "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
        plot = "WF_val:WF_time>>";
        plot += histoname;
        plot2 = "WF_val:WF_time>>";
        plot2 += histoname2;
        cut = "WF_ch==2 && spill==";
        cut += nspill;
        cut += " && event==";
        cut += vEvent[j];
        sum = 0;
        sum2 = 0;
        counter = 0;
        h4unfiltered->Draw(plot, cut, "goff");
        TempHisto = (TH2F*) gDirectory->Get(histoname);
        TempHisto->GetXaxis()->SetRange(0,150);
        cout << "h2->GetRMS(2) = " << TempHisto->GetRMS(2) << endl;
        TempHisto->GetXaxis()->SetRange(0,nbins);
        h1 = transform2Dto1D(TempHisto);
        for(int i=0; i<150; i++) {
            sum += h1->GetBinContent(i+1);
            sum2 += (h1->GetBinContent(i+1))*(h1->GetBinContent(i+1));
            counter += 1;
        }
        unfilteredrms = sqrt(sum2/counter - pow(sum/counter, 2));
        cout << "Simone calculation RMS = " << unfilteredrms << endl;
        //TempHisto->GetXaxis()->SetRange(0,150); // sets the range to only view the pedestal
        //unfilteredrms = TempHisto->GetRMS(2);
        //TempHisto->GetXaxis()->SetRange(0,nbins);
        //sum = 0;
        //sum2 = 0;
        //counter = 0;
        //h4filtered->Draw(plot2, cut, "goff");
        //TempHisto2 = (TH2F*) gDirectory->Get(histoname2);
        //h1 = transform2Dto1D(TempHisto2);
        //for(int i=0; i<150; i++) {
        //    sum += h1->GetBinContent(i+1);
        //    sum2 += (h1->GetBinContent(i+1))*(h1->GetBinContent(i+1));
        //    counter += 1;
        //}
        //filteredrms = sqrt(sum2/counter - pow(sum/counter, 2));
        //TempHisto->GetXaxis()->SetRange(0,150); // sets the range to only view the pedestal
        //filteredrms = TempHisto->GetRMS(2);
        //TempHisto->GetXaxis()->SetRange(0,nbins);

        unfilteredbrms = vuBrms[j];

        cout << "b_rms = " << unfilteredbrms << endl;

        cout << "Event " << vEvent[j] << " out of " << nevents << endl;
        count += 1;

        //out.cd();
        //MyTree->Fill();
    }
    //MyTree->Write();
}