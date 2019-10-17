#ifndef SFACTOR_H
#define SFACTOR_H

#include <TROOT.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

class SFactor {

  public:
    SFactor(int z1, int z2);
    ~SFactor();
    void Run(const char *filename);

  private:
    void GenerateData(TGraphErrors *graph);
    void MakeGraph(const char *gname);

    const Double_t alpha = 7.2973525693e-3;
    const Double_t m_be7 = 6534.1836677282;
    const Double_t m_d = 1875.61291385342;
    const Double_t lookup_frac_uncert = 2e-4;

    Double_t reduced_mass, charge_factor;
    vector<Double_t> Ecm;
    vector<Double_t> s;
    vector<Double_t> s_uncert;
    TFile *file;
};

#endif
