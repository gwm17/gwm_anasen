#ifndef SINGLECHAN_H
#define SINGLECHAN_H

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>
#include "LookUp.h"

using namespace std;

class singleChan {
  
  public:
    singleChan(float min, float max, float step, int chained);
    ~singleChan();
    void run(char *dfName, char *sfName, char *ofName);

  private:
    int getPoints(TH1F *dhisto, TH1F *shisto, TH1F *fhisto, TH1F *dethisto);
    void getCrossSection(Double_t minKE, Double_t maxKE);
    Double_t getCSError(Double_t cs,Double_t cnts,Double_t cntsf,Double_t cntsd,Double_t dist);
    void makeGraph(const char *gname);

    TFile *datafile, *simfile, *outfile;
    vector<pair<Double_t, Int_t>> points;
    vector<Double_t> Ecm, cross, avg_eff;
    vector<Double_t> cstot_err;
    vector<Double_t> detCounts;
    vector<Double_t> freeCounts;
    vector<Double_t> tot_counts;
    float maxE, minE, stepE;
    int maxBin, minBin, chainFlag;
    LookUp *beam_eloss;
    const Double_t maxBeamKE = 17.19;
    const Double_t density = 2.64e19;
    const Double_t nBParticles = 2.422e9;
    const Double_t cm2mb = 1e27;
    const Double_t m_beam = 6534.1836677282;
    const Double_t m_target = 1875.61291385342;
    const Double_t beam_perr = 0.09;
    Double_t dist_err;
};

#endif
