#ifndef __WFCLASS_H__
#define __WFCLASS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

#define MAX_INTERPOLATOR_POINTS 10000

struct WFBaseline
{
    float baseline;
    float rms;
    float slope;
    float k;
    float chi2;
};

struct WFFitResults
{
    double ampl;
    double time;
    double chi2;
};      

class WFClass : public TObject
{
public:
    //---ctors---
    WFClass() {};
    WFClass(int polarity, float tUnit);
    //---dtor---
    ~WFClass();

    //---getters---
    inline vector<float>* GetSamples() {return &samples_;};
    inline float          GetBaseline() {return baseline_;}
    float                 GetAmpMax(int min=-1, int max=-1);
    float                 GetInterpolatedAmpMax(int min=-1, int max=-1, int nFitSamples=7);
    pair<float, float>    GetTime(string method, vector<float>& params); 
    pair<float, float>    GetTimeCF(float frac, int nFitSamples=5, int min=-1, int max=-1);
    pair<float, float>    GetTimeLE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    float                 GetIntegral(int min=-1, int max=-1);
    float                 GetModIntegral(int min=-1, int max=-1);
    virtual float         GetSignalIntegral(int riseWin, int fallWin);
    //---setters---
    void                  SetSignalWindow(int min, int max);
    void                  SetBaselineWindow(int min, int max);
    void                  SetHisto(TString rootfilename, TString histoname);
    void                  SetTemplate(TH1* templateWF=NULL);
    //---utils---
    void                  Reset();
    void                  AddSample(float sample) {samples_.push_back(polarity_*sample);};
    WFBaseline            SubtractBaseline(int min=-1, int max=-1);
    WFFitResults          TemplateFit(float offset=0., int lW=0, int hW=0);
    void                  EmulatedWF(WFClass& wf, float rms, float amplitude, float time);
    void                  FFT(WFClass& wf, float tau, int cut);
    void                  FilterFFT();
    void                  Print();
    //---operators---
    WFClass&              operator=(const WFClass& origin);
    WFClass               operator-(const WFClass& sub);
    WFClass               operator+(const WFClass& add);
    WFClass&              operator-=(const WFClass& sub);
    WFClass&              operator+=(const WFClass& add);
    //---nasty trick to get code working---
    //TFile*                supportFile_;
    TFile*                inputFile_;
    TH1F*                 normNoiseFFT_;
    TH1F*                 h1_;
    TH1F*                 h1mag_;
    TH1F*                 h1phase_;
    TH1F*                 h1signalfft_;
protected:
    //---utils---
    float                 BaselineRMS();
    float                 LinearInterpolation(float& A, float& B, const int& min, const int& max);
    double                TemplateChi2(const double* par=NULL);
    
protected:
    vector<float> samples_;
    float         tUnit_;
    int           polarity_;
    int           sWinMin_;
    int           sWinMax_;
    int           bWinMin_;
    int           bWinMax_;
    int           maxSample_;
    float         fitAmpMax_;
    float         baseline_;
    float         bRMS_;
    int           cfSample_;
    float         cfFrac_;
    float         cfTime_;
    int           leSample_;
    float         leThr_;
    float         leTime_;
    float         chi2cf_;
    float         chi2le_;
    int           fWinMin_;
    int           fWinMax_;
    float         tempFitTime_;
    float         tempFitAmp_;
    
    int           nbinsFFT_;
    ROOT::Math::Interpolator* interpolator_;
};

#endif
