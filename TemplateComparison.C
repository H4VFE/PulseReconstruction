#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
using namespace std;

void TemplateComparison() {
	TFile *input = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/analysis_4443original.root"); //original
	TFile *input2 = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/Old Template/analysis_4443.root"); //old template
	//TFile *input3 = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/New Template/analysis_4443.root"); //new template, lower limit 80 (default)
    TFile *input3 = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/New Template/analysis_4443lowerlimit80.root"); //new template, lower limit 80 (default) -- naming convention, same file as above line
    //TFile *input4 = new TFile("/media/marko/TOSHIBA EXT/CERN/TB Timing Res/All ntuples/4443/New Template/analysis_4443lowerlimit160.root"); //new template, lower limit 160
	TCanvas *c1 = new TCanvas("c1", "canvas", 1600, 900);
	TCanvas *c2 = new TCanvas("c2", "canvas", 1600, 900);
	TCanvas *c3 = new TCanvas("c3", "canvas", 1600, 900);
    //TCanvas *c4 = new TCanvas("c4", "canvas", 1600, 900);
	TTree *h4unf = (TTree*) input->Get("h4");
	TTree *h4old = (TTree*) input2->Get("h4");
	TTree *h4new = (TTree*) input3->Get("h4");
    //TTree *h4new2 = (TTree*) input4->Get("h4");
	TH1F *test = new TH1F ("test", "test", 200, 36, 42);
	TH1F *Unf = new TH1F ("Unf", "unfiltered", 600, -3, 3);
	TH1F *Old = new TH1F ("Old", "old template", 600, -3, 3);
	TH1F *New = new TH1F ("New", "new template", 600, -3, 3);
	TH1F *Unf2 = new TH1F ("Unf2", "unfiltered", 200, 0, 10);
	TH1F *Old2 = new TH1F ("Old2", "old template", 200, 0, 10);
	TH1F *New2 = new TH1F ("New2", "new template", 200, 0, 10);
	Float_t mean, stddev;
	TString drawtimeres, cut;
	Int_t nspill = 102;

    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200";
    //cut += nspill;
	h4unf->Draw("fit_time[APD2]-time[MCP1]>>test", cut);
	mean = test->GetMean();
	drawtimeres = Form("fit_time[APD2]-time[MCP1]-%f>>Unf", mean);
	cout << "unf mean = " << mean << endl;
	h4unf->Draw(drawtimeres, cut);
	
    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200";
    //cut += nspill;
	h4old->Draw("fit_time[APD2]-time[MCP1]>>test", cut);
	mean = test->GetMean();
	drawtimeres = Form("fit_time[APD2]-time[MCP1]-%f>>Old", mean);
	h4old->Draw(drawtimeres, cut);
	
	h4new->Draw("fit_time[APD2]-time[MCP1]>>test", "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200");
	mean = test->GetMean();
	drawtimeres = Form("fit_time[APD2]-time[MCP1]-%f>>New", mean);
	cout << "filt mean = " << mean << endl;
	h4new->Draw(drawtimeres, "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200");

    //cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200 && spill==";
    ////cut += nspill;
    //h4new2->Draw("fit_time[APD2]-time[MCP1]>>test", cut);
    //mean = test->GetMean();
    //drawtimeres = Form("fit_time[APD2]-time[MCP1]-%f>>New", mean);
    //cout << "filt mean = " << mean << endl;
    //h4new2->Draw(drawtimeres, cut);
	
	c1->cd();
    Unf->Scale(1./Unf->Integral());
    Old->Scale(1./Old->Integral());
    New->Scale(1./New->Integral());
	Unf->SetLineColor(kRed);
	Old->SetLineColor(kOrange);
	New->SetLineColor(kBlue);
	Unf->Draw();
	Old->Draw("same");
	New->Draw("same");
	mean = New->GetMean();
	stddev = New->GetStdDev();
	cout << "New Template Fit" << endl;
	New->Fit("gaus", "", "", mean-stddev, mean+stddev);
	mean = Old->GetMean();
	stddev = Old->GetStdDev();
	cout << "Old Template Fit" << endl;
	Old->Fit("gaus", "", "", mean-stddev, mean+stddev);
	mean = Unf->GetMean();
	stddev = Unf->GetStdDev();
	cout << "Unfiltered Fit" << endl;
	Unf->Fit("gaus", "", "", mean-stddev, mean+stddev);

	cout << "Mean = " << Unf->GetMean() << "   Unfiltered Std Dev = " << Unf->GetStdDev() << "+-" << Unf->GetStdDevError() << " & Fit Std Dev = " << Unf->GetFunction("gaus")->GetParameter(2) << "+-" << Unf->GetFunction("gaus")->GetParError(2) << endl;
	cout << "Mean = " << Old->GetMean() << "  Old Template Std Dev = " << Old->GetStdDev() << "+-" << Old->GetStdDevError() << " & Fit Std Dev = " << Old->GetFunction("gaus")->GetParameter(2) << "+-" << Old->GetFunction("gaus")->GetParError(2) << endl;
	cout << "Mean = " << New->GetMean() << " New Template Std Dev = " << New->GetStdDev() << "+-" << New->GetStdDevError() << " & Fit Std Dev = " << New->GetFunction("gaus")->GetParameter(2) << "+-" << New->GetFunction("gaus")->GetParError(2) << endl;
	cout << "percent change = " << (New->GetStdDev()-Unf->GetStdDev())/Unf->GetStdDev()*100 << endl;
	cout << "percent change fit = " << (New->GetFunction("gaus")->GetParameter(2)-Unf->GetFunction("gaus")->GetParameter(2))/Unf->GetFunction("gaus")->GetParameter(2)*100 << endl;



	c2->cd();
	h4unf->Draw("b_rms[APD2]>>Unf2", "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200");
	h4old->Draw("b_rms[APD2]>>Old2", "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200");
	h4new->Draw("b_rms[APD2]>>New2", "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200");
	Unf2->SetLineColor(kRed);
	Old2->SetLineColor(kOrange);
	New2->SetLineColor(kBlue);
	Unf2->Draw();
	Old2->Draw("same");
	New2->Draw("same");
	c3->cd();

	h4unf->SetEntryList(0);
    TString listcut = "fit_ampl[APD2]>700 && WF_ch==2";
    list//cut += nspill;
    h4unf->Draw(">>myList", listcut, "entrylist");
    TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    h4unf->SetEntryList(myList);
    Int_t nevents = myList->GetN();

	h4unf->SetEntryList(0);
    listcut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200 && WF_ch==2";
    list//cut += nspill;
    h4unf->Draw(">>myList", listcut, "entrylist");
    myList = (TEntryList*) gDirectory->Get("myList");
    h4unf->SetEntryList(myList);
    Int_t neventsunf = myList->GetN();

    h4new->SetEntryList(0);
    listcut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200 && WF_ch==2";
    list//cut += nspill;
    h4new->Draw(">>myList", listcut, "entrylist");
    myList = (TEntryList*) gDirectory->Get("myList");
    h4new->SetEntryList(myList);
    Int_t neventsnew = myList->GetN();

    //h4new2->SetEntryList(0);
    //listcut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200 && WF_ch==2 && spill==";
    //list//cut += nspill;
    //h4new2->Draw(">>myList", listcut, "entrylist");
    //myList = (TEntryList*) gDirectory->Get("myList");
    //h4new2->SetEntryList(myList);
    //Int_t neventsnew2 = myList->GetN();

    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200";
    //cut += nspill;
    h4unf->Draw("b_rms[APD2]", cut, "goff");
    Double_t *vTemp_uBrms = h4unf->GetV1();
    Double_t *vuBrms = new Double_t[neventsunf];
    for (int iEntry = 0; iEntry<neventsunf; iEntry++) {
        vuBrms[iEntry] = vTemp_uBrms[iEntry];
    }

    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200";
    //cut += nspill;
    h4new->Draw("b_rms[APD2]", cut, "goff");
    Double_t *vTemp_fBrms = h4new->GetV1();
    Double_t *vfBrms = new Double_t[neventsnew];
    for (int iEntry = 0; iEntry<neventsnew; iEntry++) {
        vfBrms[iEntry] = vTemp_fBrms[iEntry];
    }

    //cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200 && spill==";
    ////cut += nspill;
    //h4new2->Draw("b_rms[APD2]", cut, "goff");
    //Double_t *vTemp_fBrms2 = h4new2->GetV1();
    //Double_t *vfBrms2 = new Double_t[neventsnew2];
    //for (int iEntry = 0; iEntry<neventsnew2; iEntry++) {
    //    vfBrms2[iEntry] = vTemp_fBrms2[iEntry];
    //}

    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>700 && amp_max[MCP1]>200";
    //cut += nspill;
	h4unf->Draw("fit_ampl[APD2]", cut, "goff");
    Double_t *vTemp_uAmpfit = h4unf->GetV1();
    Double_t *vuAmpfit = new Double_t[neventsunf];
    for (int iEntry = 0; iEntry<neventsunf; iEntry++) {
        vuAmpfit[iEntry] = vTemp_uAmpfit[iEntry];
    }

    cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200";
    //cut += nspill;
    h4new->Draw("fit_ampl[APD2]", cut, "goff");
    Double_t *vTemp_fAmpfit = h4new->GetV1();
    Double_t *vfAmpfit = new Double_t[neventsnew];
    for (int iEntry = 0; iEntry<neventsnew; iEntry++) {
        vfAmpfit[iEntry] = vTemp_fAmpfit[iEntry];
    }

    //cut = "abs(b_slope[APD2])<6 && fit_ampl[APD2]>717 && amp_max[MCP1]>200";
    ////cut += nspill;
    //h4new2->Draw("fit_ampl[APD2]", cut, "goff");
    //Double_t *vTemp_fAmpfit2 = h4new->GetV1();
    //Double_t *vfAmpfit2 = new Double_t[neventsnew2];
    //for (int iEntry = 0; iEntry<neventsnew2; iEntry++) {
    //    vfAmpfit2[iEntry] = vTemp_fAmpfit2[iEntry];
    //}

    cut = "spill ==";
    //cut += nspill;
    h4unf->Draw("fit_time[APD2]", cut, "goff");
    Double_t *vTemp_uTime = h4unf->GetV1();
    Double_t *vuTime = new Double_t[neventsunf];
    for (int iEntry = 0; iEntry<neventsunf; iEntry++) {
        vuTime[iEntry] = vTemp_uTime[iEntry];
    }

    cut = "spill ==";
    //cut += nspill;
    h4new->Draw("fit_time[APD2]", cut, "goff");
    Double_t *vTemp_fTime = h4new->GetV1();
    Double_t *vfTime = new Double_t[neventsunf];
    for (int iEntry = 0; iEntry<neventsunf; iEntry++) {
        vfTime[iEntry] = vTemp_fTime[iEntry];
    }
	c3->Divide(2,2);
    //c3->cd(1);
    //TGraph *graph = new TGraph(nevents, vuAmpfit, vfAmpfit);
    //graph->Draw();
    //graph->Fit("pol1", "", "", 400, 900);
    //cout << "Filtered Cut Value = " << graph->GetFunction("pol1")->GetParameter(0) + graph->GetFunction("pol1")->GetParameter(1)*700 << endl;
    //c3->cd(2);
    //TGraph *graph2 = new TGraph(nevents, vuTime, vfTime);
    //graph2->Draw();


    c3->cd(3);
    TH1F *maxovernoiseratio = new TH1F("maxovernoiseratio", "maxovernoiseratio", 50, -20, 800);
    TH1F *maxovernoiseratio2 = new TH1F("maxovernoiseratio2", "maxovernoiseratio2", 50, -20, 800);
    for (int iEntry = 0; iEntry<neventsunf; iEntry++) {
    	maxovernoiseratio->Fill(vuAmpfit[iEntry]/vuBrms[iEntry]);
    }
    for (int iEntry = 0; iEntry<neventsnew; iEntry++) {
    	maxovernoiseratio2->Fill(vfAmpfit[iEntry]/vfBrms[iEntry]);
    }
    maxovernoiseratio->SetLineColor(kRed);
	maxovernoiseratio2->SetLineColor(kBlue);
    maxovernoiseratio->DrawClone();
    maxovernoiseratio2->DrawClone("same");

    c3->cd(4);
    TH1F *timeratio = new TH1F("timeratio", "timeratio", 1000, 0, 2);
    for (int iEntry = 0; iEntry<nevents; iEntry++) {
    	timeratio->Fill(vfTime[iEntry]/vuTime[iEntry]);
    }
    timeratio->Draw();

    c3->cd(1);
    TH1F *maxovernoiseratio160 = new TH1F("maxovernoiseratio160", "maxovernoiseratio160", 50, 0, 800);
    //for (int iEntry = 0; iEntry<neventsnew; iEntry++) {
    //    maxovernoiseratio160->Fill(vfAmpfit2[iEntry]/vfBrms2[iEntry]);
    //}
    maxovernoiseratio2->SetLineColor(kRed);
    //maxovernoiseratio160->SetLineColor(kBlue);
    maxovernoiseratio2->Draw();
    //maxovernoiseratio160->Draw("same");

    c3->cd(2);
    TH1F *max = new TH1F("max", "max", 100, 0, 2000);
    TH1F *max160 = new TH1F("max160", "max160", 100, 600, 1000);
    h4new->Draw("fit_ampl[APD2]>>max");
    //h4new2->Draw("fit_ampl[APD2]>>max160");
    max->SetLineColor(kRed);
    //max160->SetLineColor(kBlue);
    //max160->Draw();
    max->Draw("same");

    TString filename = "spill";
    filename += nspill;
    filename += ".root";
    TFile *outfile = new TFile(filename, "recreate");
    outfile->cd();
    Unf->Write();
    Old->Write();
    New->Write();
    Unf2->Write();
    Old2->Write();
    New2->Write();


	//h4unf->Draw("event", "spill==1", "goff");
    //Double_t *vTemp = h4unfiltered->GetV1();
    //Int_t *vEvent = new Int_t[nevents];
    //for (int iEntry = 0; iEntry<nevents; iEntry++) {
    //    vEvent[iEntry] = vTemp[iEntry];
    //}
    //h4unf->SetEntryList(0);
    //TString listcut = "WF_ch==2 && spill==1";
    //h4unfiltered->Draw(">>myList", listcut, "entrylist");
    //TEntryList *myList = (TEntryList*) gDirectory->Get("myList");
    //h4unfiltered->SetEntryList(myList);
    //Int_t nevents = myList->GetN();
	//for (int i=0; i<)
}