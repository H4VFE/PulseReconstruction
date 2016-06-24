#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;
void NormNoise2() {
    Int_t nbins = 800;
    Int_t count = 0;
    Float_t max;
    Float_t min;

    Float_t divergenceval10;
    Float_t divergenceval50;
    Float_t divergenceval100;
    Float_t divergenceval225;
    Float_t divergenceval450;
    Float_t divergenceval900;
    TString HistoName;
    
    Double_t mean;
    Double_t rms;

    std::vector<Double_t> v10Events1; 
    std::vector<Double_t> v50Events1;
    std::vector<Double_t> v100Events1;
    std::vector<Double_t> v225Events1;
    std::vector<Double_t> v450Events1;
    std::vector<Double_t> v900Events1;

    std::vector<Double_t> v10Events2; 
    std::vector<Double_t> v50Events2;
    std::vector<Double_t> v100Events2;
    std::vector<Double_t> v225Events2;
    std::vector<Double_t> v450Events2;
    std::vector<Double_t> v900Events2;

    std::vector<Double_t> v10EventsNorm1; 
    std::vector<Double_t> v50EventsNorm1;
    std::vector<Double_t> v100EventsNorm1;
    std::vector<Double_t> v225EventsNorm1;
    std::vector<Double_t> v450EventsNorm1;
    std::vector<Double_t> v900EventsNorm1;

    std::vector<Double_t> v10EventsNorm2; 
    std::vector<Double_t> v50EventsNorm2;
    std::vector<Double_t> v100EventsNorm2;
    std::vector<Double_t> v225EventsNorm2;
    std::vector<Double_t> v450EventsNorm2;
    std::vector<Double_t> v900EventsNorm2;

    std::vector<Double_t> v10EventsKolmog; 
    std::vector<Double_t> v50EventsKolmog;
    std::vector<Double_t> v100EventsKolmog;
    std::vector<Double_t> v225EventsKolmog;
    std::vector<Double_t> v450EventsKolmog;
    std::vector<Double_t> v900EventsKolmog;

    TCanvas* Divergence = new TCanvas("Divergence","Divergence",1800,1200);
    Divergence->Divide(2,2);

    TH1F* First10Events = new TH1F ("First10Events", "First 10 Events", nbins, 0, 5);
    TH1F* Second10Events = new TH1F ("Second10Events", "Second 10 Events", nbins, 0, 5);
    TH1F* First50Events = new TH1F ("First50Events", "First 50 Events", nbins, 0, 5);
    TH1F* Second50Events = new TH1F ("Second50Events", "Second 50 Events", nbins, 0, 5);
    TH1F* First100Events = new TH1F ("First100Events", "First 100 Events", nbins, 0, 5);
    TH1F* Second100Events = new TH1F ("Second100Events", "Second 100 Events", nbins, 0, 5);
    TH1F* First225Events = new TH1F ("First225Events", "First 225 Events", nbins, 0, 5);
    TH1F* Second225Events = new TH1F ("Second225Events", "Second 225 Events", nbins, 0, 5);
    TH1F* First450Events = new TH1F ("First450Events", "First 450 Events", nbins, 0, 5);
    TH1F* Second450Events = new TH1F ("Second450Events", "Second 450 Events", nbins, 0, 5);
    TH1F* First900Events = new TH1F ("First900Events", "First 900 Events", nbins, 0, 5);
    TH1F* Second900Events = new TH1F ("Second900Events", "Second 900 Events", nbins, 0, 5);

    TH1F* Subbed10Events = new TH1F ("Subbed10Events", "Subtracted 50 Events", nbins/2, 0, 2.5);
    TH1F* Subbed50Events = new TH1F ("Subbed50Events", "Subtracted 50 Events", nbins/2, 0, 2.5);
    TH1F* Subbed100Events = new TH1F ("Subbed100Events", "Subtracted 100 Events", nbins/2, 0, 2.5);
    TH1F* Subbed225Events = new TH1F ("Subbed225Events", "Subtracted 225 Events", nbins/2, 0, 2.5);
    TH1F* Subbed450Events = new TH1F ("Subbed450Events", "Subtracted 450 Events", nbins/2, 0, 2.5);
    TH1F* Subbed900Events = new TH1F ("Subbed900Events", "Subtracted 900 Events", nbins/2, 0, 2.5);

    TH1F* NormSubbed10Events = new TH1F ("NormSubbed10Events", "NormSubbed 50 Events", nbins/2, 0, 2.5);
    TH1F* NormSubbed50Events = new TH1F ("NormSubbed50Events", "NormSubbed 50 Events", nbins/2, 0, 2.5);
    TH1F* NormSubbed100Events = new TH1F ("NormSubbed100Events", "NormSubbed 100 Events", nbins/2, 0, 2.5);
    TH1F* NormSubbed225Events = new TH1F ("NormSubbed225Events", "NormSubbed 225 Events", nbins/2, 0, 2.5);
    TH1F* NormSubbed450Events = new TH1F ("NormSubbed450Events", "NormSubbed 450 Events", nbins/2, 0, 2.5);
    TH1F* NormSubbed900Events = new TH1F ("NormSubbed900Events", "NormSubbed 900 Events", nbins/2, 0, 2.5);

    TH1F* NormNoiseFFT100Events = new TH1F ("NormNoiseFFT100Events", "Normalized Noise 100 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT250Events = new TH1F ("NormNoiseFFT250Events", "Normalized Noise 250 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT500Events = new TH1F ("NormNoiseFFT500Events", "Normalized Noise 500 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT750Events = new TH1F ("NormNoiseFFT750Events", "Normalized Noise 750 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT1000Events = new TH1F ("NormNoiseFFT1000Events", "Normalized Noise 1000 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT1250Events = new TH1F ("NormNoiseFFT1250Events", "Normalized Noise 1250 Events", nbins, 0, 5);
    TH1F* NormNoiseFFT1500Events = new TH1F ("NormNoiseFFT1500Events", "Normalized Noise 1500 Events", nbins, 0, 5);

    TH1F* EventFFT;
    TFile* f1 = new TFile("/home/marko/Desktop/TB Timing Res/AllPedestalEventFFTs.root", "read");

    Int_t i;
    Int_t j;

    for (i=0;i<1800/20;i++) {
        for (j=(0+20*i);j<(10+20*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First10Events->Add(First10Events, EventFFT);
            count += 1;
        }
        First10Events->Scale(1./count);
        count = 0;
        for (j=(10+20*i);j<(20+20*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second10Events->Add(Second10Events, EventFFT);
            count += 1;
        }
        Second10Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed10Events->SetBinContent(k+1, First10Events->GetBinContent(k+1)-Second10Events->GetBinContent(k+1));
            NormSubbed10Events->SetBinContent(k+1, abs(2*(First10Events->GetBinContent(k+1)-Second10Events->GetBinContent(k+1))/(First10Events->GetBinContent(k+1)+Second10Events->GetBinContent(k+1))));
        }

        Subbed10Events->GetXaxis()->SetRange(0,100);
        max = Subbed10Events->GetMaximum();
        min = Subbed10Events->GetMinimum();
        if (max > abs(min)) divergenceval10 = max;
        if (abs(min) > max) divergenceval10 = abs(min);
        v10Events1.push_back(divergenceval10);

        Subbed10Events->GetXaxis()->SetRange(100,400);
        max = Subbed10Events->GetMaximum();
        min = Subbed10Events->GetMinimum();
        if (max > abs(min)) divergenceval10 = max;
        if (abs(min) > max) divergenceval10 = abs(min);
        v10Events2.push_back(divergenceval10);

        NormSubbed10Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed10Events->GetMaximum();
        min = NormSubbed10Events->GetMinimum();
        if (max > abs(min)) divergenceval10 = max;
        if (abs(min) > max) divergenceval10 = abs(min);
        v10EventsNorm1.push_back(divergenceval10);

        NormSubbed10Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed10Events->GetMaximum();
        min = NormSubbed10Events->GetMinimum();
        if (max > abs(min)) divergenceval10 = max;
        if (abs(min) > max) divergenceval10 = abs(min);
        v10EventsNorm2.push_back(divergenceval10);
    }

    for (i=0;i<1800/100;i++) {
        for (j=(0+100*i);j<(50+100*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First50Events->Add(First50Events, EventFFT);
            count += 1;
        }
        First50Events->Scale(1./count);
        count = 0;
        for (j=(50+100*i);j<(100+100*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second50Events->Add(Second50Events, EventFFT);
            count += 1;
        }
        Second50Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed50Events->SetBinContent(k+1, First50Events->GetBinContent(k+1)-Second50Events->GetBinContent(k+1));
            NormSubbed50Events->SetBinContent(k+1, abs(2*(First50Events->GetBinContent(k+1)-Second50Events->GetBinContent(k+1))/(First50Events->GetBinContent(k+1)+Second50Events->GetBinContent(k+1))));
        }

        Subbed50Events->GetXaxis()->SetRange(0,100);
        max = Subbed50Events->GetMaximum();
        min = Subbed50Events->GetMinimum();
        if (max > abs(min)) divergenceval50 = max;
        if (abs(min) > max) divergenceval50 = abs(min);
        v50Events1.push_back(divergenceval50);

        Subbed50Events->GetXaxis()->SetRange(100,400);
        max = Subbed50Events->GetMaximum();
        min = Subbed50Events->GetMinimum();
        if (max > abs(min)) divergenceval50 = max;
        if (abs(min) > max) divergenceval50 = abs(min);
        v50Events2.push_back(divergenceval50);

        NormSubbed50Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed50Events->GetMaximum();
        min = NormSubbed50Events->GetMinimum();
        if (max > abs(min)) divergenceval50 = max;
        if (abs(min) > max) divergenceval50 = abs(min);
        v50EventsNorm1.push_back(divergenceval50);

        NormSubbed50Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed50Events->GetMaximum();
        min = NormSubbed50Events->GetMinimum();
        if (max > abs(min)) divergenceval50 = max;
        if (abs(min) > max) divergenceval50 = abs(min);
        v50EventsNorm2.push_back(divergenceval50);
    }
    for (i=0;i<1800/200;i++) {
        for (j=(0+200*i);j<(100+200*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First100Events->Add(First100Events, EventFFT);
            count += 1;
        }
        First100Events->Scale(1./count);
        count = 0;
        for (j=(100+200*i);j<(200+200*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second100Events->Add(Second100Events, EventFFT);
            count += 1;
        }
        Second100Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed100Events->SetBinContent(k+1, First100Events->GetBinContent(k+1)-Second100Events->GetBinContent(k+1));
            NormSubbed100Events->SetBinContent(k+1, abs(2*(First100Events->GetBinContent(k+1)-Second100Events->GetBinContent(k+1))/(First100Events->GetBinContent(k+1)+Second100Events->GetBinContent(k+1))));
        }
        
        Subbed100Events->GetXaxis()->SetRange(0,100);
        max = Subbed100Events->GetMaximum();
        min = Subbed100Events->GetMinimum();
        if (max > abs(min)) divergenceval100 = max;
        if (abs(min) > max) divergenceval100 = abs(min);
        v100Events1.push_back(divergenceval100);

        Subbed100Events->GetXaxis()->SetRange(100,400);
        max = Subbed100Events->GetMaximum();
        min = Subbed100Events->GetMinimum();
        if (max > abs(min)) divergenceval100 = max;
        if (abs(min) > max) divergenceval100 = abs(min);
        v100Events2.push_back(divergenceval100);

        NormSubbed100Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed100Events->GetMaximum();
        min = NormSubbed100Events->GetMinimum();
        if (max > abs(min)) divergenceval100 = max;
        if (abs(min) > max) divergenceval100 = abs(min);
        v100EventsNorm1.push_back(divergenceval100);

        NormSubbed100Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed100Events->GetMaximum();
        min = NormSubbed100Events->GetMinimum();
        if (max > abs(min)) divergenceval100 = max;
        if (abs(min) > max) divergenceval100 = abs(min);
        v100EventsNorm2.push_back(divergenceval100);
    }
    for (i=0;i<1800/450;i++) {
        for (j=(0+450*i);j<(225+450*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First225Events->Add(First225Events, EventFFT);
            count += 1;
        }
        First225Events->Scale(1./count);
        count = 0;
        for (j=(225+450*i);j<(450+450*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second225Events->Add(Second225Events, EventFFT);
            count += 1;
        }
        Second225Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed225Events->SetBinContent(k+1, First225Events->GetBinContent(k+1)-Second225Events->GetBinContent(k+1));
            NormSubbed225Events->SetBinContent(k+1, abs(2*(First225Events->GetBinContent(k+1)-Second225Events->GetBinContent(k+1))/(First225Events->GetBinContent(k+1)+Second225Events->GetBinContent(k+1))));
        }
        
        Subbed225Events->GetXaxis()->SetRange(0,100);
        max = Subbed225Events->GetMaximum();
        min = Subbed225Events->GetMinimum();
        if (max > abs(min)) divergenceval225 = max;
        if (abs(min) > max) divergenceval225 = abs(min);
        v225Events1.push_back(divergenceval225);

        Subbed225Events->GetXaxis()->SetRange(100,400);
        max = Subbed225Events->GetMaximum();
        min = Subbed225Events->GetMinimum();
        if (max > abs(min)) divergenceval225 = max;
        if (abs(min) > max) divergenceval225 = abs(min);
        v225Events2.push_back(divergenceval225);

        NormSubbed225Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed225Events->GetMaximum();
        min = NormSubbed225Events->GetMinimum();
        if (max > abs(min)) divergenceval225 = max;
        if (abs(min) > max) divergenceval225 = abs(min);
        v225EventsNorm1.push_back(divergenceval225);

        NormSubbed225Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed225Events->GetMaximum();
        min = NormSubbed225Events->GetMinimum();
        if (max > abs(min)) divergenceval225 = max;
        if (abs(min) > max) divergenceval225 = abs(min);
        v225EventsNorm2.push_back(divergenceval225);
    }
    for (i=0;i<1800/900;i++) {
        for (j=(0+900*i);j<(450+900*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First450Events->Add(First450Events, EventFFT);
            count += 1;
        }
        First450Events->Scale(1./count);
        count = 0;
        for (j=(450+900*i);j<(900+900*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second450Events->Add(Second450Events, EventFFT);
            count += 1;
        }
        Second450Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed450Events->SetBinContent(k+1, First450Events->GetBinContent(k+1)-Second450Events->GetBinContent(k+1));
            NormSubbed450Events->SetBinContent(k+1, abs(2*(First450Events->GetBinContent(k+1)-Second450Events->GetBinContent(k+1))/(First450Events->GetBinContent(k+1)+Second450Events->GetBinContent(k+1))));
        }
        
        Subbed450Events->GetXaxis()->SetRange(0,100);
        max = Subbed450Events->GetMaximum();
        min = Subbed450Events->GetMinimum();
        if (max > abs(min)) divergenceval450 = max;
        if (abs(min) > max) divergenceval450 = abs(min);
        v450Events1.push_back(divergenceval450);

        Subbed450Events->GetXaxis()->SetRange(100,400);
        max = Subbed450Events->GetMaximum();
        min = Subbed450Events->GetMinimum();
        if (max > abs(min)) divergenceval450 = max;
        if (abs(min) > max) divergenceval450 = abs(min);
        v450Events2.push_back(divergenceval450);

        NormSubbed450Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed450Events->GetMaximum();
        min = NormSubbed450Events->GetMinimum();
        if (max > abs(min)) divergenceval450 = max;
        if (abs(min) > max) divergenceval450 = abs(min);
        v450EventsNorm1.push_back(divergenceval450);

        NormSubbed450Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed450Events->GetMaximum();
        min = NormSubbed450Events->GetMinimum();
        if (max > abs(min)) divergenceval450 = max;
        if (abs(min) > max) divergenceval450 = abs(min);
        v450EventsNorm2.push_back(divergenceval450);
    }
    for (i=0;i<1800/1800;i++) {
        for (j=(0+1800*i);j<(900+1800*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            First900Events->Add(First900Events, EventFFT);
            count += 1;
        }
        First900Events->Scale(1./count);
        count = 0;
        for (j=(900+1800*i);j<(1800+1800*i);j++) {
            HistoName = "NewHistoEventFFT";
            HistoName += j;
            EventFFT = (TH1F*) f1->Get(HistoName);
            Second900Events->Add(Second900Events, EventFFT);
            count += 1;
        }
        Second900Events->Scale(1./count);
        count = 0;
        for (Int_t k=0;k<nbins/2;k++) {
            Subbed900Events->SetBinContent(k+1, First900Events->GetBinContent(k+1)-Second900Events->GetBinContent(k+1));
            NormSubbed900Events->SetBinContent(k+1, abs(2*(First900Events->GetBinContent(k+1)-Second900Events->GetBinContent(k+1))/(First900Events->GetBinContent(k+1)+Second900Events->GetBinContent(k+1))));
        }
        
        Subbed900Events->GetXaxis()->SetRange(0,100);
        max = Subbed900Events->GetMaximum();
        min = Subbed900Events->GetMinimum();
        if (max > abs(min)) divergenceval900 = max;
        if (abs(min) > max) divergenceval900 = abs(min);
        v900Events1.push_back(divergenceval900);

        Subbed900Events->GetXaxis()->SetRange(100,400);
        max = Subbed900Events->GetMaximum();
        min = Subbed900Events->GetMinimum();
        if (max > abs(min)) divergenceval900 = max;
        if (abs(min) > max) divergenceval900 = abs(min);
        v900Events2.push_back(divergenceval900);

        NormSubbed900Events->GetXaxis()->SetRange(0,100);
        max = NormSubbed900Events->GetMaximum();
        min = NormSubbed900Events->GetMinimum();
        if (max > abs(min)) divergenceval900 = max;
        if (abs(min) > max) divergenceval900 = abs(min);
        v900EventsNorm1.push_back(divergenceval900);

        NormSubbed900Events->GetXaxis()->SetRange(100,400);
        max = NormSubbed900Events->GetMaximum();
        min = NormSubbed900Events->GetMinimum();
        if (max > abs(min)) divergenceval900 = max;
        if (abs(min) > max) divergenceval900 = abs(min);
        v900EventsNorm2.push_back(divergenceval900);
    }
    Float_t mean10a = TMath::Mean(v10Events1.begin(), v10Events1.end());
    Float_t rms10a = TMath::RMS(v10Events1.begin(), v10Events1.end());
 
    Float_t mean50a = TMath::Mean(v50Events1.begin(), v50Events1.end());
    Float_t rms50a = TMath::RMS(v50Events1.begin(), v50Events1.end());
 
    Float_t mean100a = TMath::Mean(v100Events1.begin(), v100Events1.end());
    Float_t rms100a = TMath::RMS(v100Events1.begin(), v100Events1.end());
 
    Float_t mean225a = TMath::Mean(v225Events1.begin(), v225Events1.end());
    Float_t rms225a = TMath::RMS(v225Events1.begin(), v225Events1.end());
 
    Float_t mean450a = TMath::Mean(v450Events1.begin(), v450Events1.end());
    Float_t rms450a = TMath::RMS(v450Events1.begin(), v450Events1.end());
 
    Float_t mean900a = TMath::Mean(v900Events1.begin(), v900Events1.end());
    Float_t rms900a = TMath::RMS(v900Events1.begin(), v900Events1.end());
 
    Float_t mean10b = TMath::Mean(v10Events2.begin(), v10Events2.end());
    Float_t rms10b = TMath::RMS(v10Events2.begin(), v10Events2.end());
 
    Float_t mean50b = TMath::Mean(v50Events2.begin(), v50Events2.end());
    Float_t rms50b = TMath::RMS(v50Events2.begin(), v50Events2.end());
 
    Float_t mean100b = TMath::Mean(v100Events2.begin(), v100Events2.end());
    Float_t rms100b = TMath::RMS(v100Events2.begin(), v100Events2.end());
 
    Float_t mean225b = TMath::Mean(v225Events2.begin(), v225Events2.end());
    Float_t rms225b = TMath::RMS(v225Events2.begin(), v225Events2.end());
 
    Float_t mean450b = TMath::Mean(v450Events2.begin(), v450Events2.end());
    Float_t rms450b = TMath::RMS(v450Events2.begin(), v450Events2.end());
 
    Float_t mean900b = TMath::Mean(v900Events2.begin(), v900Events2.end());
    Float_t rms900b = TMath::RMS(v900Events2.begin(), v900Events2.end());
 
    Float_t normmean10a = TMath::Mean(v10EventsNorm1.begin(), v10EventsNorm1.end());
    Float_t normrms10a = TMath::RMS(v10EventsNorm1.begin(), v10EventsNorm1.end());
 
    Float_t normmean50a = TMath::Mean(v50EventsNorm1.begin(), v50EventsNorm1.end());
    Float_t normrms50a = TMath::RMS(v50EventsNorm1.begin(), v50EventsNorm1.end());
 
    Float_t normmean100a = TMath::Mean(v100EventsNorm1.begin(), v100EventsNorm1.end());
    Float_t normrms100a = TMath::RMS(v100EventsNorm1.begin(), v100EventsNorm1.end());
 
    Float_t normmean225a = TMath::Mean(v225EventsNorm1.begin(), v225EventsNorm1.end());
    Float_t normrms225a = TMath::RMS(v225EventsNorm1.begin(), v225EventsNorm1.end());
 
    Float_t normmean450a = TMath::Mean(v450EventsNorm1.begin(), v450EventsNorm1.end());
    Float_t normrms450a = TMath::RMS(v450EventsNorm1.begin(), v450EventsNorm1.end());
 
    Float_t normmean900a = TMath::Mean(v900EventsNorm1.begin(), v900EventsNorm1.end());
    Float_t normrms900a = TMath::RMS(v900EventsNorm1.begin(), v900EventsNorm1.end());
 
    Float_t normmean10b = TMath::Mean(v10EventsNorm2.begin(), v10EventsNorm2.end());
    Float_t normrms10b = TMath::RMS(v10EventsNorm2.begin(), v10EventsNorm2.end());
 
    Float_t normmean50b = TMath::Mean(v50EventsNorm2.begin(), v50EventsNorm2.end());
    Float_t normrms50b = TMath::RMS(v50EventsNorm2.begin(), v50EventsNorm2.end());
 
    Float_t normmean100b = TMath::Mean(v100EventsNorm2.begin(), v100EventsNorm2.end());
    Float_t normrms100b = TMath::RMS(v100EventsNorm2.begin(), v100EventsNorm2.end());
 
    Float_t normmean225b = TMath::Mean(v225EventsNorm2.begin(), v225EventsNorm2.end());
    Float_t normrms225b = TMath::RMS(v225EventsNorm2.begin(), v225EventsNorm2.end());
 
    Float_t normmean450b = TMath::Mean(v450EventsNorm2.begin(), v450EventsNorm2.end());
    Float_t normrms450b = TMath::RMS(v450EventsNorm2.begin(), v450EventsNorm2.end());
 
    Float_t normmean900b = TMath::Mean(v900EventsNorm2.begin(), v900EventsNorm2.end());
    Float_t normrms900b = TMath::RMS(v900EventsNorm2.begin(), v900EventsNorm2.end());

    Double_t x[6] = {10, 50, 100, 225, 450, 900};
    Double_t y[6] = {mean10a, mean50a, mean100a, mean225a, mean450a, mean900a};
    Double_t ey[6] = {rms10a, rms50a, rms100a, rms225a, rms450a, rms900a};
    Double_t ex[6] = {0, 0, 0, 0, 0, 0,};

    Double_t y2[6] = {normmean10a, normmean50a, normmean100a, normmean225a, normmean450a, normmean900a};
    Double_t ey2[6] = {normrms10a, normrms50a, normrms100a, normrms225a, normrms450a, normrms900a};

    Double_t y3[6] = {mean10b, mean50b, mean100b, mean225b, mean450b, mean900b};
    Double_t ey3[6] = {rms10b, rms50b, rms100b, rms225b, rms450b, rms900b};

    Double_t y4[6] = {normmean10b, normmean50b, normmean100b, normmean225b, normmean450b, normmean900b};
    Double_t ey4[6] = {normrms10b, normrms50b, normrms100b, normrms225b, normrms450b, normrms900b};

    Divergence->cd(1);
    TGraphErrors *graph = new TGraphErrors(6, x, y, ex, ey);
    graph->SetTitle("Maximum Divergence from 0 as a Function of Number of Events (0-0.6 GHz)");
    graph->SetMarkerSize(2);
    graph->SetMarkerStyle(22);
    graph->GetXaxis()->SetTitle("Events");
    graph->GetYaxis()->SetTitle("Maximum Divergence");
    graph->SetMinimum(0);
    graph->SetMaximum(50);
    TF1 *myfunc = new TF1("myfunc", "[0]/x+[1]", 0, 1000);
    graph->Fit("myfunc", "", "");
    graph->Draw("AP");

    Divergence->cd(2);
    TGraphErrors *graph2 = new TGraphErrors(6, x, y2, ex, ey2);
    graph2->SetTitle("Normalized Maximum Divergence from 0 as a Function of Number of Events (0-0.6 GHz)");
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(22);
    graph2->GetXaxis()->SetTitle("Events");
    graph2->GetYaxis()->SetTitle("Normalized Maximum Divergence");
    graph2->SetMinimum(0);
    graph2->SetMaximum(1);
    TF1 *myfunc2 = new TF1("myfunc2", "[0]+[1]*exp(-x/[2])", 0, 500);
    myfunc2->SetParameter(2, 100);
    myfunc2->SetParameter(0, 0.1);
    graph2->Fit("myfunc2", "", "");
    graph2->Draw("AP");

    Divergence->cd(3);
    TGraphErrors *graph3 = new TGraphErrors(6, x, y3, ex, ey3);
    graph3->SetTitle("Maximum Divergence from 0 as a Function of Number of Events (0.6-2.5 GHz)");
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(22);
    graph3->GetXaxis()->SetTitle("Events");
    graph3->GetYaxis()->SetTitle("Maximum Divergence");
    graph3->SetMinimum(0);
    graph3->SetMaximum(50);
    graph3->Fit("myfunc", "", "");
    graph3->Draw("AP");

    Divergence->cd(4);
    TGraphErrors *graph4 = new TGraphErrors(6, x, y4, ex, ey4);
    graph4->SetTitle("Normalized Maximum Divergence from 0 as a Function of Number of Events (0.6-2.5 GHz)");
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(22);
    graph4->GetXaxis()->SetTitle("Events");
    graph4->GetYaxis()->SetTitle("Normalized Maximum Divergence");
    graph4->SetMinimum(0);
    graph4->SetMaximum(1);
    graph4->Fit("myfunc", "", "");
    graph4->Draw("AP");

    Divergence->SaveAs("DivergenceVsEvents.pdf");

 //   TFile* out = new TFile("NormNoiseFFTs.root", "recreate");
//
 //   for (j=0;j<100;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT100Events->Add(NormNoiseFFT100Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT100Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<250;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT250Events->Add(NormNoiseFFT250Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT250Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<500;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT500Events->Add(NormNoiseFFT500Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT500Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<750;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT750Events->Add(NormNoiseFFT750Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT750Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<1000;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT1000Events->Add(NormNoiseFFT1000Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT1000Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<1250;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT1250Events->Add(NormNoiseFFT1250Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT1250Events->Scale(1./count);
 //   count = 0;
 //   for (j=0;j<1500;j++) {
 //       HistoName = "NewHistoEventFFT";
 //       HistoName += j;
 //       EventFFT = (TH1F*) f1->Get(HistoName);
 //       NormNoiseFFT1500Events->Add(NormNoiseFFT1500Events, EventFFT);
 //       count += 1;
 //   }
 //   NormNoiseFFT1500Events->Scale(1./count);
 //   count = 0;
//
 //   NormNoiseFFT100Events->Write();
 //   NormNoiseFFT250Events->Write();
 //   NormNoiseFFT500Events->Write();
 //   NormNoiseFFT750Events->Write();
 //   NormNoiseFFT1000Events->Write();
 //   NormNoiseFFT1250Events->Write();
 //   NormNoiseFFT1500Events->Write();
}