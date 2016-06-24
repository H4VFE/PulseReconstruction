#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEntryList.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TVirtualFFT.h"
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

void TimeResolution() {
	TCanvas* Time = new TCanvas("Time", "Time Resolution", 1600, 900);
	TCanvas* Energy = new TCanvas("Energy", "Energy", 1600, 900);
	TCanvas* Events = new TCanvas("Events", "Events", 1600, 900);
	Time->Divide(2,1);
	Energy->Divide(2,1);
	Events->Divide(2,1);
	TFile* f = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/New Template/40keventsanalysis_4443.root");
	TFile* f1 = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/analysis_4443original.root");
	TTree* h4filtered = (TTree*) f->Get("h4");
	TTree* h4unfiltered = (TTree*) f1->Get("h4");
	Int_t nbins = 3500;
	Float_t maxrange = 100;
	TH1F* TimeRes = new TH1F ("TimeRes", "Time Resolution", nbins, 0, maxrange);
	TH1F* TimeRes2 = new TH1F ("TimeRes2", "Time Resolution", nbins, 0, maxrange);
	TH1F* AmpMax = new TH1F ("AmpMax", "Energy", 250, 0, 1000);
	TH1F* AmpMax2 = new TH1F ("AmpMax2", "Energy", 250, 0, 1000);
	TH2F* TempHisto = new TH2F ("TempHisto", "Temp Histo", 800, -0.1, 159.9, 1000, -120, 800);
	TH1F* Unfiltered1D = new TH1F ("Unfiltered1D", "Unfiltered 1D", 800, -0.1, 159.9);
	TH1F* Filtered1D = new TH1F ("Filtered1D", "Filtered 1D", 800, -0.1, 159.9);
	Events->cd(1);
	h4unfiltered->Draw("WF_val:WF_time>>TempHisto", "WF_ch==APD2 && event==1 && spill==1");
	for (Int_t xbin=0; xbin<800; xbin++) {
        for (Int_t ybin=0; ybin<1000; ybin++) {
            if (TempHisto->GetBinContent(xbin+1, ybin) != 0) {
                Unfiltered1D->SetBinContent(xbin+1,ybin*(920./1000)-120);
            }
        }
    }
    Unfiltered1D->SetLineColor(kRed);
	Unfiltered1D->Draw();
	Events->cd(2);
	h4filtered->Draw("WF_val:WF_time>>TempHisto", "WF_ch==APD2 && event==1 && spill==1");
	for (Int_t xbin=0; xbin<800; xbin++) {
        for (Int_t ybin=0; ybin<1000; ybin++) {
            if (TempHisto->GetBinContent(xbin+1, ybin) != 0) {
                Filtered1D->SetBinContent(xbin+1,ybin*(920./1000)-120);
            }
        }
    }
    Filtered1D->SetLineColor(kBlue);
    Unfiltered1D->Draw();
	Filtered1D->Draw("same");
	Time->cd(1);
	h4unfiltered->Draw("fit_time[APD2]-time[MCP1]>>TimeRes","amp_max[MCP1]>100 && fit_ampl[APD2]>650");
	Int_t mean, maxbin, lowerlimitbin, upperlimitbin;
	Double_t prezoom_stddev, postzoom_stddev1, postzoom_stddev2, maxval, lowerlimit, upperlimit;
	mean = TimeRes->GetMean();
	maxbin = TimeRes->GetMaximumBin();
	maxval = TimeRes->GetBinContent(maxbin);
	prezoom_stddev = TimeRes->GetStdDev();
	TimeRes->GetXaxis()->SetRange(maxbin-2*prezoom_stddev/(maxrange/nbins), maxbin+2*prezoom_stddev/(maxrange/nbins));
	TimeRes2->GetXaxis()->SetRange(maxbin-2*prezoom_stddev/(maxrange/nbins), maxbin+2*prezoom_stddev/(maxrange/nbins));
	postzoom_stddev1 = TimeRes->GetStdDev();
	for (Int_t i=maxbin;i>0;i--) {
		lowerlimitbin = TimeRes->GetBinContent(i);
		if (lowerlimitbin <= maxval/2) {
			lowerlimit = TimeRes->GetBinCenter(i);
			break;
		}
	}
	for (Int_t i=maxbin;i<TimeRes->GetNbinsX();i++){
		upperlimitbin = TimeRes->GetBinContent(i);
		if (upperlimitbin <= maxval/2) {
			upperlimit = TimeRes->GetBinCenter(i);
			break;
		}
	}
	Time->cd(1);
	TimeRes->SetLineColor(kRed);
	TimeRes2->SetLineColor(kBlue);
	TimeRes->Scale(1./TimeRes->Integral());
	TimeRes->Draw();
	Energy->cd(1);
	h4unfiltered->Draw("fit_ampl[APD2]>>AmpMax");
	AmpMax->GetXaxis()->SetRange(125,225);
	AmpMax->SetLineColor(kRed);
	AmpMax->SetLineWidth(2);
	AmpMax->DrawClone();
	cout << "Unfiltered Full Width Half Max = " << abs(upperlimit-lowerlimit) << " ns" << endl;
	cout << "Unfiltered Standard Deviation = " << postzoom_stddev1 << " ns" << endl;


	Time->cd(2);
	h4filtered->Draw("fit_time[APD2]-time[MCP1]>>TimeRes2","amp_max[MCP1]>100 && fit_ampl[APD2]>650");
	mean = TimeRes->GetMean();
	maxbin = TimeRes->GetMaximumBin();
	maxval = TimeRes->GetBinContent(maxbin);
	postzoom_stddev2 = TimeRes2->GetStdDev();
	for (Int_t i=maxbin;i>0;i--) {
		lowerlimitbin = TimeRes->GetBinContent(i);
		if (lowerlimitbin <= maxval/2) {
			lowerlimit = TimeRes->GetBinCenter(i);
			break;
		}
	}
	for (Int_t i=maxbin;i<TimeRes->GetNbinsX();i++){
		upperlimitbin = TimeRes->GetBinContent(i);
		if (upperlimitbin <= maxval/2) {
			upperlimit = TimeRes->GetBinCenter(i);
			break;
		}
	}
	Time->cd(2);
	TimeRes2->Scale(1./TimeRes2->Integral());
	TimeRes2->Draw();
	TimeRes->DrawClone("same");
	Energy->cd(2);
	h4filtered->Draw("fit_ampl[APD2]>>AmpMax2");
	AmpMax2->GetXaxis()->SetRange(125,225);
	AmpMax2->SetLineColor(kBlue);
	AmpMax2->SetLineWidth(2);
	AmpMax->DrawClone();
	AmpMax2->DrawClone("same");
	cout << "Filtered Full Width Half Max = " << abs(upperlimit-lowerlimit) << " ns" << endl;
	cout << "Filtered Standard Deviation = " << postzoom_stddev2 << " ns" << endl;
	cout << "Std Dev ratio = " << 1.*postzoom_stddev2/postzoom_stddev1 << endl;
	cout << "before std dev = " << postzoom_stddev1 << endl;
	cout << "Bin Width = " << maxrange/nbins << endl;
}