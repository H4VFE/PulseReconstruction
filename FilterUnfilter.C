#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TCanvas.h"
using namespace std;

void FilterUnfilter() {
	TCanvas *can1 = new TCanvas("can1", "canvas", 1600, 900);
	TCanvas *can2 = new TCanvas("can2", "canvas 2", 1200, 1200);
	TCanvas *can3 = new TCanvas("can3", "canvas 3", 1200, 1200);
	can1->Divide(2,2);
	can2->Divide(2,2);
	can3->Divide(2,2);
	Float_t mean, stddev, int1, int2;
	TString drawtimeres;
	TFile *input = new TFile("FilterRMSComparison.root");
	TTree *MyTree = (TTree*) input->Get("RMS");
	TH1F *temp = new TH1F("temp", "temp", 200, 20, 60);

	TH1F *utimeref = new TH1F("utimeref", "Time Ref", 200, 20, 60);
	TH1F *ftimeref = new TH1F("ftimeref", "Filtered Time Ref", 200, 20, 60);

	TH1F *ubslope = new TH1F("ubslope", "Baseline Slope", 200, -60, 60);
	TH1F *fbslope = new TH1F("fbslope", "Filtered Baseline Slope", 200, -60, 60);

	TH1F *uampmaxmcp = new TH1F("uampmaxmcp", "MCP Amp Max", 200, 0, 2000);
	TH1F *fampmaxmcp = new TH1F("fampmaxmcp", "Filtered MCP Amp Max", 200, 0, 2000);

	TH1F *utimefit = new TH1F("utimefit", "Fit Time", 200, 60, 110);
	TH1F *ftimefit = new TH1F("ftimefit", "Filtered Fit Time", 200, 60, 110);

	TH1F *uampfit = new TH1F("uampfit", "Fit Amp", 200, 0, 800);
	TH1F *fampfit = new TH1F("fampfit", "Filtered Fit Amp", 200, 0, 800);

	TH1F *udeltatime = new TH1F("udeltatime", "Delta Time", 200, -3, 3);
	TH1F *fdeltatime = new TH1F("fdeltatime", "Filtered Delta Time", 200, -3, 3);

	TH1F *ampfitpercent = new TH1F("ampfitpercent", "Amp Fit Percentage", 200, -1, 1);

	TH2F *timeref = new TH2F("timeref", "Time Ref", 200, 30, 50, 200, 30, 50);

	TH2F *timefit = new TH2F("timefit", "Time Fit", 200, 70, 90, 200, 70, 90);

	TH2F *ampfit = new TH2F("ampfit", "Amp Fit", 200, 500, 800, 200, 500, 800);

	TH2F *deltatime = new TH2F("deltatime", "Delta Time", 200, 36, 42, 200, 36, 42);

	utimeref->SetLineColor(kRed);
	ftimeref->SetLineColor(kBlue);
	ubslope->SetLineColor(kRed);
	fbslope->SetLineColor(kBlue);
	uampmaxmcp->SetLineColor(kRed);
	fampmaxmcp->SetLineColor(kBlue);
	utimefit->SetLineColor(kRed);
	ftimefit->SetLineColor(kBlue);
	uampfit->SetLineColor(kRed);
	fampfit->SetLineColor(kBlue);
	udeltatime->SetLineColor(kRed);
	fdeltatime->SetLineColor(kBlue);
	can1->cd(1);
	//MyTree->Draw("unfilteredtimeref>>utimeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	//MyTree->Draw("filteredtimeref>>ftimeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "same");




	//MyTree->Draw("filteredbslope>>fbslope", "unfilteredampfit>700");
	//MyTree->Draw("unfilteredbslope>>ubslope", "filteredampfit>717", "same");
	MyTree->Draw("(filteredbslope/unfilteredbslope)>>temp2(1000,-10,10)", "abs(unfilteredbslope)>20");




	can1->cd(2);
	//MyTree->Draw("unfilteredtimefit>>utimefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400");
	//MyTree->Draw("filteredtimefit>>ftimefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "same");
	MyTree->Draw("unfilteredampmaxmcp>>uampmaxmcp");
	MyTree->Draw("filteredampmaxmcp>>fampmaxmcp", "", "same");
	can1->cd(3);
	MyTree->Draw("unfilteredampfit>>uampfit");
	MyTree->Draw("filteredampfit>>fampfit", "", "same");
	//MyTree->Draw("(filteredampfit-unfilteredampfit)/unfilteredampfit>>ampfitpercent", "abs(unfilteredbslope)<6 && unfilteredampfit>400");


	can1->cd(4);
	MyTree->Draw("unfilteredtimefit-unfilteredtimeref>>temp", "unfilteredampfit>500 && unfilteredampmaxmcp>100");
	//MyTree->Draw("unfilteredtimefit-unfilteredtimeref>>temp", "abs(unfilteredbslope)<6 && unfilteredampfit>700 && unfilteredampmaxmcp>100");
	mean = temp->GetMean();
	drawtimeres = Form("unfilteredtimefit-unfilteredtimeref-%f>>udeltatime", mean);
	MyTree->Draw(drawtimeres, "unfilteredampfit>500 && unfilteredampmaxmcp>100");
	//MyTree->Draw(drawtimeres, "abs(unfilteredbslope)<6 && unfilteredampfit>700 && unfilteredampmaxmcp>100");
	cout << "Pre-Filter Standard Deviation: " << udeltatime->GetStdDev()*1000 << " ps." << endl;

	MyTree->Draw("filteredtimefit-filteredtimeref>>temp", "filteredampfit>500 && filteredampmaxmcp>100");
	//MyTree->Draw("filteredtimefit-filteredtimeref>>temp", "abs(unfilteredbslope)<6 && unfilteredampfit>700 && unfilteredampmaxmcp>100");
	mean = temp->GetMean();
	drawtimeres = Form("filteredtimefit-filteredtimeref-%f>>fdeltatime", mean);
	MyTree->Draw(drawtimeres, "filteredampfit>500 && filteredampmaxmcp>100");
	//MyTree->Draw(drawtimeres, "abs(unfilteredbslope)<6 && unfilteredampfit>700 && unfilteredampmaxmcp>100", "same");

	int1 = udeltatime->Integral();
	int2 = fdeltatime->Integral();

	can2->cd(1);
	udeltatime->Scale(1./int1);
	fdeltatime->Scale(1./int2);
	fdeltatime->DrawClone();
	udeltatime->DrawClone("same");

	can2->cd(2);
	gPad->SetLogy();
	fdeltatime->DrawClone();
	udeltatime->DrawClone("same");

	can2->cd(3);
	MyTree->Draw("filteredbrms/unfilteredbrms");

	can1->cd(4);
	udeltatime->Scale(int1);
	fdeltatime->Scale(int2);
	mean = udeltatime->GetMean();
	stddev = udeltatime->GetStdDev();
	udeltatime->Fit("gaus", "", "", mean-stddev, mean+stddev);
	mean = fdeltatime->GetMean();
	stddev = fdeltatime->GetStdDev();
	fdeltatime->Fit("gaus", "", "", mean-stddev, mean+stddev);
	cout << "Post-Filter Standard Deviation: " << fdeltatime->GetStdDev()*1000 << " ps." << endl;
	//udeltatime->Scale(1./udeltatime->Integral());
	//fdeltatime->Scale(1./fdeltatime->Integral());
	fdeltatime->Draw();
	udeltatime->Draw("same");

	
	//can2->cd(1);
	//MyTree->Draw("unfilteredtimeref:filteredtimeref>>timeref", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	//can2->cd(2);
	//MyTree->Draw("unfilteredtimefit:filteredtimefit>>timefit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	can2->cd(3);
	MyTree->Draw("unfilteredampfit:filteredampfit>>ampfit", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");
	//can2->cd(4);
	//MyTree->Draw("(unfilteredtimefit-unfilteredtimeref):(filteredtimefit-filteredtimeref)>>deltatime", "abs(unfilteredbslope)<6 && unfilteredampfit>400", "colz");

	can3->cd(1);
	MyTree->Draw("abs(unfilteredbslope):unfilteredbrms>>h2001(50, 0, 45, 50, 0, 70)", "", "colz");
	can3->cd(3);
	MyTree->Draw("abs(filteredbslope):filteredbrms>>h2002(50, 0, 45, 50, 0, 70)", "", "colz");
}