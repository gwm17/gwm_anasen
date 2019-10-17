#ifndef DALITZ_H
#define DALITZ_H

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "LookUp.h"

using namespace std;

class DalAnalyzer {

  public:
    DalAnalyzer(float min, float max, float step);
    ~DalAnalyzer();
    void run(char *dfilename, char *sfilename, char *output);

  private:
    int getPoints();
    void goCrossSection(Double_t minKE, Double_t maxKE);
    Double_t getCSError(Double_t cs, Double_t counts, Double_t e_perr, Double_t dist);
    void makePlots(char *output);
    void makeThetaPlot(char *output);
    Double_t interpBinContent(Int_t binx, Int_t biny, Int_t gbin);

    //Points for cross-section calc (efficiency, counts)
    vector<pair<Double_t,Int_t>> points;
    vector<Double_t> cross;
    vector<Double_t> Ecm;
    vector<Double_t> eff_perr;
    vector<Double_t> cstot_err;
    vector<Double_t> avg_eff;
    vector<Double_t> ae_Ecm;
    vector<Double_t> theta;
    float maxE, minE, stepE;
    int maxBin;
    Int_t depth;
    TFile *dataFile;
    TFile *simFile;
    TH2F *dataPlot, *simPlot, *freePlot, *detPlot;
    LookUp *beam_eloss;
    const Double_t maxBeamKE = 17.19;
    const Double_t density = 2.64e19;
    const Double_t nBParticles = 2.422e9;
    const Double_t cm2mb = 1e27;
    const Double_t m_beam = 6534.1836677282;
    const Double_t m_target = 1875.61291385342;
    const Double_t beam_perr = 0.09;
    Double_t dist_err;
    Int_t tot_counts;

};

#endif
