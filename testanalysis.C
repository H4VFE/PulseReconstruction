#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;




//add phase comparison to canvas
void testanalysis() {
	TCanvas* canvas = new TCanvas("canvas","Time & Frequency",1800,1200);
    //TCanvas* canvas2 = new TCanvas("canvas2","Frequency",800,1200);
    //TCanvas* canvas3 = new TCanvas("canvas3","Frequency2", 800, 1200);
	canvas->Divide(2,3);
    Int_t nbins = 200*1;
    TF1 *func = new TF1("func", "(x<6)*((x>0.5)*800*(-exp(-3*(x-0.5)) + exp(-1.5*(x-0.5))))", 0, 8*1);
    canvas->cd(1);
    func->SetTitle("Shifted Function");
    func->SetMinimum(0);
    func->SetMaximum(210);
    func->DrawClone();
    TH1F* funchisto = new TH1F("funchisto", "Original Function", nbins, 0, 8*1);
    TH1F* funcFFT = new TH1F("funcFFT", "Func FFT", nbins, 0, 4);
    TH1F* funcFFTph = new TH1F("funcFFTph", "Func Phase", nbins, 0, 4);
    for (Int_t i=0; i<=nbins; i++) {
        //Double_t x = Double_t(i)/nbins*(4);
                Double_t x = Double_t(i)/nbins*(8*1);

        funchisto->SetBinContent(i+1, func->Eval(x));
    }
    funchisto->SetMinimum(0);
    funchisto->SetMaximum(210);
    funchisto->SetLineColor(kRed);
    funchisto->Draw();
    //canvas->cd(3);
    funchisto->FFT(funcFFT, "MAG");
    funchisto->FFT(funcFFTph, "PH");
    funcFFTph->SetLineColor(kRed);
    funcFFT->SetLineColor(kRed);
    //funcFFT->DrawClone();
    //gPad->SetLogy();
    TF1* noisy = new TF1("noisy", "(800*(-exp(-3*x) + exp(-1.5*x))) + (sin(3*x)/10 + x/10)", 0, 4);
    TF1* shiftedfunc = new TF1("shiftedfunc", "(x<6.01)*((x>0.51)*800*(-exp(-3*(x-0.51)) + exp(-1.5*(x-0.51))))", 0, 8*1);
    shiftedfunc->SetTitle("Shifted Function");
    canvas->cd(2);
    //noisy->SetTitle("Signal & Noise");
    //noisy->GetYaxis()->SetRange(0,210);
    //noisy->DrawClone();
    shiftedfunc->SetMinimum(0);
    shiftedfunc->SetMaximum(210);
    shiftedfunc->Draw();
    TH1F* noisyhisto = new TH1F("noisyhisto", "Signal & Noise ", nbins, 0, 4);
    TH1F* noisyFFT = new TH1F("noisyFFT", "Signal & Noise FFT S+N", nbins, 0, 5);
    TH1F* noisyPH = new TH1F("noisyPH", "Noisy Signal Phase", nbins, 0, 100);
    TH1F* shiftedhisto = new TH1F("shiftedhisto", "Shifted Function", nbins, 0, 8*1);
    TH1F* shiftedFFT = new TH1F("shiftedFFT", "Shifted Func FFT", nbins, 0, 4);
    TH1F* shiftedFFTph = new TH1F("shiftedFFTph", "Shifted Func Phase", nbins, 0, 4);
    for (Int_t j=0; j<=nbins; j++) {
        Double_t x = (Double_t(j)/nbins)*(4);
        noisyhisto->SetBinContent(j+1, noisy->Eval(x));
    }
    for (Int_t j=0; j<=nbins; j++) {
        Double_t x = (Double_t(j)/nbins)*(8*1);
        shiftedhisto->SetBinContent(j+1, shiftedfunc->Eval(x));
    }
    shiftedhisto->SetMinimum(0);
    shiftedhisto->SetMaximum(210);
    shiftedhisto->Draw();
    noisyhisto->FFT(noisyFFT, "MAG");
    noisyhisto->FFT(noisyPH, "PH");
    shiftedhisto->FFT(shiftedFFT, "MAG");
    shiftedhisto->FFT(shiftedFFTph, "PH");
    canvas->cd(3);
    //noisyFFT->DrawClone();
    shiftedFFT->Draw();
    funcFFT->Draw("same");
    gPad->SetLogy();
    canvas->cd(4);
    shiftedFFTph->Draw();
    funcFFTph->Draw("same");
    TF1 *noise = new TF1("noise", "sin(3*x)/10 + x/10", 0, 4);
    canvas->cd(5);
    TH1F *FFTsub = new TH1F("FFTsub", "FFT subtraction", nbins, 0, 4);
    TH1F *FFTsubnorm = new TH1F("FFTsubnorm", "Normalized FFT subtraction", nbins, 0, 4);
    for (Int_t i=0;i<nbins;i++) {
        FFTsub->SetBinContent(i+1, shiftedFFT->GetBinContent(i+1) - funcFFT->GetBinContent(i+1));
    }
    for (Int_t i=0;i<nbins;i++) {
        //if (shiftedFFT->GetBinContent(i+1) + funcFFT->GetBinContent(i+1)==0) FFTsubnorm->SetBinContent(i+1, 0);
        FFTsubnorm->SetBinContent(i+1, (shiftedFFT->GetBinContent(i+1) - funcFFT->GetBinContent(i+1))*2/(shiftedFFT->GetBinContent(i+1) + funcFFT->GetBinContent(i+1)));
    }
    FFTsub->SetLineColor(6);
    FFTsub->Draw();
    //noise->DrawClone();    
    TH1F* noisehisto = new TH1F("noisehisto", "Noise", nbins, 0, 4);
    TH1F* noiseFFT = new TH1F("noiseFFT", "Noise FFT", nbins, 0, 5);
    for (Int_t k=0; k<=nbins; k++) {
        Double_t x = Double_t(k)/nbins*(4);
        noisehisto->SetBinContent(k+1, noise->Eval(x));
    }
    noisehisto->FFT(noiseFFT, "MAG");
    canvas->cd(6);
    FFTsubnorm->SetLineColor(6);
    FFTsubnorm->Draw();
    //noiseFFT->DrawClone();
    gPad->SetLogy();
    //canvas2->Divide(1,3);
    //canvas2->cd(1);
    //noisyFFT->DrawClone();
    //gPad->SetLogy();
    //canvas2->cd(2);
    //noiseFFT->DrawClone();
    //gPad->SetLogy();
    //canvas2->cd(3); 
    //noisyFFT->Draw();
    //noiseFFT->Draw("same");
    //gPad->SetLogy();
    TH1F* normsignalFFT = new TH1F ("normsignalFFT", "Normalized Signal FFT S/(S+N)", nbins, 0, 5);
    for (Int_t l=0; l<nbins; l++) {
         normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));
    }
    TH1F* signalFFT = new TH1F ("signalFFT", "Signal FFT", nbins, 0, 5);
    for (Int_t q=0; q<nbins; q++) {
        signalFFT->SetBinContent(q+1, normsignalFFT->GetBinContent(q+1)*noisyFFT->GetBinContent(q+1));
    }
    //canvas3->Divide(1,3);
    //canvas3->cd(1);
    //normsignalFFT->Draw();
    //canvas3->cd(2);
    //signalFFT->Draw();
    //canvas3->cd(3);
    Double_t *re_full = new Double_t[nbins];
    Double_t *im_full = new Double_t[nbins];
    for (Int_t m=0; m<nbins; m++) {
        (re_full)[m]=(signalFFT->GetBinContent(m+1)*cos(noisyPH->GetBinContent(m+1)));
        (im_full)[m]=(signalFFT->GetBinContent(m+1)*sin(noisyPH->GetBinContent(m+1)));
    }
    TVirtualFFT *invFFT = TVirtualFFT::FFT(1, &nbins, "C2R M K");
    invFFT->SetPointsComplex(re_full, im_full);
    invFFT->Transform();
    TH1 *Signal = 0;
    Signal = TH1::TransformHisto(invFFT,Signal,"Re");
    Signal->SetTitle("Recovered Signal 'S'");
    TH1F* bettersignal = new TH1F ("bettersignal", "Recovered Signal", nbins, 0, 4);
    for (Int_t p=0; p<nbins; p++) {
        bettersignal->SetBinContent(p+1, Signal->GetBinContent(p+1)/nbins);
    }
    //func->Draw();
    //bettersignal->Draw("same");
    



	

    //Double_t noisyFFTIntegral = noisyFFT->Integral();
    //Double_t noiseFFTIntegral = noiseFFT->Integral();
    //noisyFFT->Scale(1/noisyFFTIntegral);
    //noiseFFT->Scale(1/noiseFFTIntegral);
    //normsignalFFT->SetBinContent(l+1, (noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1));
    //normsignalFFT->SetBinContent(l+1, pow((noisyFFT->GetBinContent(l+1)-noiseFFT->GetBinContent(l+1))/noisyFFT->GetBinContent(l+1),2));
    //signalFFT->SetBinContent(q+1, noisyFFT->GetBinContent(q+1)); //checking if output = input, without taking noise out


}