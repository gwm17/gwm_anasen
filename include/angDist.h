#ifndef ANGDIST_H
#define ANGDIST_H

#include <TROOT.h>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <vector>
#include <string>
#include <iostream>
#include "LookUp.h"

using namespace std;

class angDist {

  public:
    angDist(char *dfName, char *sfName, char *ofName, int chained);
    ~angDist();
    void effCorrection(TH1F *dhisto, TH1F *shisto);
    void calcErrors(TH1F *dhisto, TH1F *fhisto, TH1F *dethisto);
    void makeGraph(const char *gname);
    void run();

  private:
    TFile *datafile, *simfile, *outfile;
    vector<Double_t> cost;
    vector<Double_t> cor_Counts;
    vector<Double_t> err;
};

#endif
