#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TCanvas.h"
using namespace std;

void FilterUnfilter() {
	TCanvas *can1 = new TCanvas("can1", "canvas", 1600, 900);
	TCanvas *can2 = new TCanvas("can2", "canvas 2", 1200, 1200);
	can1->Divide(2,2);
	can2->Divide(2,2);
	TFile *input = new TFile("FilterRMSComparison.root");
	TTree *MyTree = (TTree*) input->Get("RMS");
	TH1F *utimeref = new TH1F("utimeref", "Unfiltered Time Ref", 200, 20, 60);
	TH1F *ftimeref = new TH1F("ftimeref", "Filtered Time Ref", 200, 20, 60);
	TH1F *utimefit = new TH1F("utimefit", "Unfiltered Fit Time", 200, 60, 110);
	TH1F *ftimefit = new TH1F("ftimefit", "Filtered Fit Time", 200, 60, 110);
	TH1F *uampfit = new TH1F("uampfit", "Unfiltered Fit Amp", 200, 400, 800);
	TH1F *fampfit = new TH1F("fampfit", "Filtered Fit Amp", 200, 400, 800);
	TH1F *udeltatime = new TH1F("udeltatime", "Unfiltered Delta Time", 200, 36, 42);
	TH1F *fdeltatime = new TH1F("fdeltatime", "Filtered Delta Time", 200, 36, 42);
	TH1F *ampfitpercent = new TH1F("ampfitpercent", "Amp Fit Percentage", 200, -1, 1);
	TH2F *timeref = new TH2F("timeref", "Time Ref", 200, 30, 50, 200, 30, 50);
	TH2F *timefit = new TH2F("timefit", "Time Fit", 200, 70, 90, 200, 70, 90);
	TH2F *ampfit = new TH2F("ampfit", "Amp Fit", 200, 500, 800, 200, 500, 800);
	TH2F *deltatime = new TH2F("deltatime", "Delta Time", 200, 36, 42, 200, 36, 42);
	utimeref->SetLineColor(kRed);
	ftimeref->SetLineColor(kBlue);
	utimefit->SetLineColor(kRed);
	ftimefit->SetLineColor(kBlue);
	uampfit->SetLineColor(kRed);
	fampfit->SetLineColor(kBlue);
	udeltatime->SetLineColor(kRed);
	fdeltatime->SetLineColor(kBlue);
	can1->cd(1);
	MyTree->Draw("unfilteredtimeref>>utimeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	MyTree->Draw("filteredtimeref>>ftimeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "same");
	can1->cd(2);
	MyTree->Draw("unfilteredtimefit>>utimefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	MyTree->Draw("filteredtimefit>>ftimefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "same");
	can1->cd(3);
	//MyTree->Draw("unfilteredampfit>>uampfit");
	//MyTree->Draw("filteredampfit>>fampfit", "", "same");
	MyTree->Draw("(filteredampfit-unfilteredampfit)/unfilteredampfit>>ampfitpercent", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	can1->cd(4);
	MyTree->Draw("unfilteredtimefit-unfilteredtimeref>>udeltatime", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	MyTree->Draw("filteredtimefit-filteredtimeref>>fdeltatime", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "same");
	can2->cd(1);
	MyTree->Draw("unfilteredtimeref:filteredtimeref>>timeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	can2->cd(2);
	MyTree->Draw("unfilteredtimefit:filteredtimefit>>timefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	can2->cd(3);
	MyTree->Draw("unfilteredampfit:filteredampfit>>ampfit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	can2->cd(4);
	MyTree->Draw("(unfilteredtimefit-unfilteredtimeref):(filteredtimefit-filteredtimeref)>>deltatime", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
}