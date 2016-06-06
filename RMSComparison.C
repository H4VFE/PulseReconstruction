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
    Int_t *vEvent = new Int_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vEvent[iEntry] = vTemp[iEntry];
    }

    h4filtered->Draw("event", "spill==1", "goff");
    Double_t *vTemp1 = h4filtered->GetV1();
    Double_t *vfEvent = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfEvent[iEntry] = vTemp1[iEntry];
    }

    h4unfiltered->Draw("b_slope[APD2]", "spill==1", "goff");
    Double_t *vTemp_uBslope = h4unfiltered->GetV1();
    Double_t *vuBslope = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuBslope[iEntry] = vTemp_uBslope[iEntry];
    }

    h4filtered->Draw("b_slope[APD2]", "spill==1", "goff");
    Double_t *vTemp_fBslope = h4filtered->GetV1();
    Double_t *vfBslope = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfBslope[iEntry] = vTemp_fBslope[iEntry];
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

    h4unfiltered->Draw("maximum[APD2]", "spill==1", "goff");
    Double_t *vTemp_uMaximum = h4unfiltered->GetV1();
    Double_t *vuMaximum = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuMaximum[iEntry] = vTemp_uMaximum[iEntry];
    }

    h4filtered->Draw("maximum[APD2]", "spill==1", "goff");
    Double_t *vTemp_fMaximum = h4filtered->GetV1();
    Double_t *vfMaximum = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfMaximum[iEntry] = vTemp_fMaximum[iEntry];
    }

    h4unfiltered->Draw("amp_max[APD2]", "spill==1", "goff");
    Double_t *vTemp_uAmpmax = h4unfiltered->GetV1();
    Double_t *vuAmpmax = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuAmpmax[iEntry] = vTemp_uAmpmax[iEntry];
    }

    h4filtered->Draw("amp_max[APD2]", "spill==1", "goff");
    Double_t *vTemp_fAmpmax = h4filtered->GetV1();
    Double_t *vfAmpmax = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfAmpmax[iEntry] = vTemp_fAmpmax[iEntry];
    }

    h4unfiltered->Draw("time[MCP1]", "spill==1", "goff");
    Double_t *vTemp_uTimeref = h4unfiltered->GetV1();
    Double_t *vuTimeref = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuTimeref[iEntry] = vTemp_uTimeref[iEntry];
    }

    h4filtered->Draw("time[MCP1]", "spill==1", "goff");
    Double_t *vTemp_fTimeref = h4filtered->GetV1();
    Double_t *vfTimeref = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfTimeref[iEntry] = vTemp_fTimeref[iEntry];
    }

    h4unfiltered->Draw("fit_time[APD2]", "spill==1", "goff");
    Double_t *vTemp_uTimefit = h4unfiltered->GetV1();
    Double_t *vuTimefit = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuTimefit[iEntry] = vTemp_uTimefit[iEntry];
    }

    h4filtered->Draw("fit_time[APD2]", "spill==1", "goff");
    Double_t *vTemp_fTimefit = h4filtered->GetV1();
    Double_t *vfTimefit = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfTimefit[iEntry] = vTemp_fTimefit[iEntry];
    }

    h4unfiltered->Draw("fit_ampl[APD2]", "spill==1", "goff");
    Double_t *vTemp_uAmpfit = h4unfiltered->GetV1();
    Double_t *vuAmpfit = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vuAmpfit[iEntry] = vTemp_uAmpfit[iEntry];
    }

    h4filtered->Draw("fit_ampl[APD2]", "spill==1", "goff");
    Double_t *vTemp_fAmpfit = h4filtered->GetV1();
    Double_t *vfAmpfit = new Double_t[nevents];
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
        vfAmpfit[iEntry] = vTemp_fAmpfit[iEntry];
    }

    float* unfilteredwfval;
    float* filteredwfval;
    float* unfilteredwftime;
    float* filteredwftime;

    unfilteredwfval = new float[nbins];
    filteredwfval = new float[nbins];
    unfilteredwftime = new float[nbins];
    filteredwftime = new float[nbins];

    //h4unfiltered->Draw("WF_val", "WF_ch==APD2 && spill==1", "goff");
    //Double_t *vTemp_uWfval = h4unfiltered->GetV1();
    //Double_t *vuWfval = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vuWfval[iEntry] = vTemp_uWfval[iEntry];
    //}
//
    //h4filtered->Draw("WF_val", "WF_ch==APD2 && spill==1", "goff");
    //Double_t *vTemp_fWfval = h4filtered->GetV1();
    //Double_t *vfWfval = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vfWfval[iEntry] = vTemp_fWfval[iEntry];
    //}
//
    //h4unfiltered->Draw("WF_time", "WF_ch==APD2 && spill==1", "goff");
    //Double_t *vTemp_uWftime = h4unfiltered->GetV1();
    //Double_t *vuWftime = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vuWftime[iEntry] = vTemp_uWftime[iEntry];
    //}
//
    //h4filtered->Draw("WF_time", "WF_ch==APD2 && spill==1", "goff");
    //Double_t *vTemp_fWftime = h4filtered->GetV1();
    //Double_t *vfWftime = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vfWftime[iEntry] = vTemp_fWftime[iEntry];
    //}
//
    //h4unfiltered->Draw("charge_tot", "spill==1", "goff");
    //Double_t *vTemp_uChargetot = h4unfiltered->GetV1();
    //Double_t *vuChargetot = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vuChargetot[iEntry] = vTemp_uChargetot[iEntry];
    //}
//
    //h4filtered->Draw("charge_tot", "spill==1", "goff");
    //Double_t *vTemp_fChargetot = h4filtered->GetV1();
    //Double_t *vfChargetot = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vfChargetot[iEntry] = vTemp_fChargetot[iEntry];
    //}
//
    //h4unfiltered->Draw("charge_sig", "spill==1", "goff");
    //Double_t *vTemp_uChargesig = h4unfiltered->GetV1();
    //Double_t *vuChargesig = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vuChargesig[iEntry] = vTemp_uChargesig[iEntry];
    //}
//
    //h4filtered->Draw("charge_sig", "spill==1", "goff");
    //Double_t *vTemp_fChargesig = h4filtered->GetV1();
    //Double_t *vfChargesig = new Double_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vfChargesig[iEntry] = vTemp_fChargesig[iEntry];
    //}

    Float_t unfilteredevent, filteredevent, unfilteredrms, filteredrms, unfilteredbrms, filteredbrms, unfilteredbslope, filteredbslope, unfilteredmax, filteredmax, unfilteredampmax, filteredampmax, unfilteredtimefit, filteredtimefit, unfilteredtimeref, filteredtimeref, unfilteredampfit, filteredampfit;
    Int_t count = 0, counter;
    Float_t sum, sum2;
    TFile out("FilterRMSComparison.root", "recreate");
    TTree *MyTree = new TTree ("RMS", "RMS");

    MyTree->Branch("unfilteredevent", &unfilteredevent, "unfilteredevent/F");
    MyTree->Branch("filteredevent", &filteredevent, "filteredevent/F");

    MyTree->Branch("unfilteredrms", &unfilteredrms, "unfilteredrms/F");
    MyTree->Branch("filteredrms", &filteredrms, "filteredrms/F");

    MyTree->Branch("unfilteredbrms", &unfilteredbrms, "unfilteredbrms/F");
    MyTree->Branch("filteredbrms", &filteredbrms, "filteredbrms/F");

    MyTree->Branch("unfilteredbslope", &unfilteredbslope, "unfilteredbslope/F");
    MyTree->Branch("filteredbslope", &filteredbslope, "filteredbslope/F");

    MyTree->Branch("unfilteredmax", &unfilteredmax, "unfilteredmax/F");
    MyTree->Branch("filteredmax", &filteredmax, "filteredmax/F");

    MyTree->Branch("unfilteredampmax", &unfilteredampmax, "unfilteredampmax/F");
    MyTree->Branch("filteredampmax", &filteredampmax, "filteredampmax/F");

    MyTree->Branch("unfilteredtimeref", &unfilteredtimeref, "unfilteredtimeref/F");
    MyTree->Branch("filteredtimeref", &filteredtimeref, "filteredtimeref/F");

    MyTree->Branch("unfilteredtimefit", &unfilteredtimefit, "unfilteredtimefit/F");
    MyTree->Branch("filteredtimefit", &filteredtimefit, "filteredtimefit/F");

    MyTree->Branch("unfilteredampfit", &unfilteredampfit, "unfilteredampfit/F");
    MyTree->Branch("filteredampfit", &filteredampfit, "filteredampfit/F");

    //MyTree->Branch("unfilteredwfval", &unfilteredwfval, "unfilteredwfval/F");
    //MyTree->Branch("filteredwfval", &filteredwfval, "filteredwfval/F");
//
    //MyTree->Branch("unfilteredwftime", &unfilteredwftime, "unfilteredwftime/F");
    //MyTree->Branch("filteredwftime", &filteredwftime, "filteredwftime/F");

    MyTree->Branch("unfilteredwfval", unfilteredwfval, "unfilteredwfval[800]/F");
    MyTree->Branch("filteredwfval", filteredwfval, "filteredwfval[800]/F");

    MyTree->Branch("unfilteredwftime", unfilteredwftime, "unfilteredwftime[800]/F");
    MyTree->Branch("filteredwftime", filteredwftime, "filteredwftime[800]/F");
    
    //MyTree->Branch("unfilteredbcharge", &unfilteredbcharge, "unfilteredbcharge/F");
    //MyTree->Branch("filteredbcharge", &filteredbcharge, "filteredbcharge/F");

    //MyTree->Branch("unfilteredchargetot", &unfilteredchargetot, "unfilteredchargetot/F");
    //MyTree->Branch("filteredchargetot", &filteredchargetot, "filteredchargetot/F");

    //MyTree->Branch("unfilteredchargesig", &unfilteredchargesig, "unfilteredchargesig/F");
    //MyTree->Branch("filteredchargesig", &filteredchargesig, "filteredchargesig/F");

    TString plot, plot2, cut;
    char name[50];

    TCanvas *can1 = new TCanvas("can1", "canvas", 1200,800);

    for (Int_t j=0;j<nevents;j++) {
    //for (Int_t j=0;j<10;j++) {
        TString histoname = "TempHisto_";
        histoname += j;
        TString histoname2 = "TempHisto2_";
        histoname2 += j;
        TH2F* TempHisto = new TH2F (histoname, "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
        TH2F* TempHisto2 = new TH2F (histoname2, "Temp Histo", nbins, -0.1, 159.9, 1000, -120, 800); //nanoseconds
        TString h1name = "h1001_";
        h1name += j;
        TString h1name2 = "h1002_";
        h1name2 += j;
        TH1F *h1001 = new TH1F(h1name,"Red = Unfiltered, Blue = Filtered", nbins, -0.1, 159.9);
        TH1F *h1002 = new TH1F(h1name2,"h1002", nbins, -0.1, 159.9);
        plot = "WF_val:WF_time>>";
        plot += histoname;
        plot2 = "WF_val:WF_time>>";
        plot2 += histoname2;
        cut = "WF_ch==2 && spill==";
        cut += nspill;
        cut += " && event==";
        cut += vEvent[j];

        unfilteredevent = vEvent[j];
        filteredevent = vfEvent[j];

        sum = 0;
        sum2 = 0;
        counter = 0;
        h4unfiltered->Draw(plot, cut, "goff");
        TempHisto = (TH2F*) gDirectory->Get(histoname);
        h1001 = transform2Dto1D(TempHisto);
        for(int i=0; i<150; i++) {
            sum += h1001->GetBinContent(i+1);
            sum2 += (h1001->GetBinContent(i+1))*(h1001->GetBinContent(i+1));
            counter += 1;
        }
        unfilteredrms = sqrt(sum2/counter - pow(sum/counter, 2));
        //TempHisto->GetXaxis()->SetRange(0,150); // sets the range to only view the pedestal
        //unfilteredrms = TempHisto->GetRMS(2);
        //TempHisto->GetXaxis()->SetRange(0,nbins);
        sum = 0;
        sum2 = 0;
        counter = 0;
        h4filtered->Draw(plot2, cut, "goff");
        TempHisto2 = (TH2F*) gDirectory->Get(histoname2);
        h1002 = transform2Dto1D(TempHisto2);
        for(int i=0; i<150; i++) {
            sum += h1002->GetBinContent(i+1);
            sum2 += (h1002->GetBinContent(i+1))*(h1002->GetBinContent(i+1));
            counter += 1;
        }
        filteredrms = sqrt(sum2/counter - pow(sum/counter, 2));
        //TempHisto->GetXaxis()->SetRange(0,150); // sets the range to only view the pedestal
        //filteredrms = TempHisto->GetRMS(2);
        //TempHisto->GetXaxis()->SetRange(0,nbins);

        sprintf(name, "All Events/Event%d", vEvent[j]);
        strcat(name, ".png");
        h1001->SetLineColor(kRed);
        h1002->SetLineColor(kBlue);
        h1001->SetStats(0);
        h1002->SetStats(0);
        can1->cd();
        h1001->Draw();
        h1002->Draw("same");
        can1->SaveAs(name);

        unfilteredbslope = vuBslope[j];
        filteredbslope = vfBslope[j];

        unfilteredbrms = vuBrms[j];
        filteredbrms = vfBrms[j];

        unfilteredmax = vuMaximum[j];
        filteredmax = vfMaximum[j];

        unfilteredampmax = vuAmpmax[j];
        filteredampmax = vfAmpmax[j];

        unfilteredtimeref = vuTimeref[j];
        filteredtimeref = vfTimeref[j];

        unfilteredtimefit = vuTimefit[j];
        filteredtimefit = vfTimefit[j];

        unfilteredampfit = vuAmpfit[j];
        filteredampfit = vfAmpfit[j];

        for (Int_t i=0;i<nbins;i++) {
            unfilteredwfval[i] = h1001->GetBinContent(i+1);
            filteredwfval[i] = h1002->GetBinContent(i+1);

            unfilteredwftime[i] = h1001->GetBinCenter(i+1);
            filteredwftime[i] = h1002->GetBinCenter(i+1);
        }

        //unfilteredbcharge = vuBcharge[j];
        //filteredbcharge = vfBcharge[j];

        //unfilteredchargetot = vuChargetot[j];
        //filteredchargetot = vfChargetot[j];

        //unfilteredchargesig = vuChargesig[j];
        //filteredchargesig = vfChargesig[j];

        cout << "Event " << j << " out of " << nevents << endl;
        count += 1;

        out.cd();
        MyTree->Fill();

        delete TempHisto, TempHisto2, h1001, h1002;
    }
    MyTree->Write();
}