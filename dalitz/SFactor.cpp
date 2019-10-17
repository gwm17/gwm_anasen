#include "SFactor.h"
#include <TMath.h>
#include <TAxis.h>

using namespace std;

SFactor::SFactor(int z1, int z2) {
  reduced_mass = m_be7*m_d/(m_be7+m_d);
  charge_factor = alpha*z1*z2;
}

SFactor::~SFactor() {
  if(file->IsOpen()) file->Close();
}

void SFactor::GenerateData(TGraphErrors *graph) {
  vector<Double_t> ecm(graph->GetX(), graph->GetX()+graph->GetN());
  vector<Double_t> cs(graph->GetY(), graph->GetY()+graph->GetN());
  vector<Double_t> cs_uncert(graph->GetEY(), graph->GetEY()+graph->GetN());
  for(unsigned int i=0; i<cs.size(); i++) {
    Double_t eta = charge_factor*sqrt(reduced_mass/(2.0*ecm[i]));
    Double_t expon = TMath::Exp(2.0*TMath::Pi()*eta);
    Double_t sf = cs[i]*ecm[i]*expon*1e-3;
    Double_t sf_uncert = sf*sqrt((cs_uncert[i]/cs[i])*(cs_uncert[i]/cs[i])+
                                 lookup_frac_uncert*lookup_frac_uncert);
    Ecm.push_back(ecm[i]);
    s.push_back(sf);
    s_uncert.push_back(sf_uncert);
  }
}

void SFactor::MakeGraph(const char *gname) {
  TGraphErrors *graph = new TGraphErrors(s.size(), &(Ecm[0]), &(s[0]), 0, &(s_uncert[0]));
  graph->SetName(gname);
  graph->SetTitle(gname);
  graph->SetMarkerStyle(22);
  graph->SetMarkerColor(4);
  graph->GetXaxis()->SetLimits(0.18, 2.4);
  graph->GetYaxis()->SetRangeUser(0, 250);
  file->cd();
  graph->Write();
}

void SFactor::Run(const char *filename) {
  file = new TFile(filename, "UPDATE");
  TGraphErrors *graph = (TGraphErrors*) file->Get("3He CS");
  if(graph != NULL) {
    GenerateData(graph);
    string gname = "3He S-factor";
    MakeGraph(gname.c_str());
  } else {
    cout<<"Excitation function graph does not exist!! Check gname passed to SFactor::Run()"
        <<endl;
  }
  graph = (TGraphErrors*) file->Get("(d,p&a) Cross Section Ecm");
  if(graph != NULL) {
    Ecm.resize(0); s.resize(0); s_uncert.resize(0);
    GenerateData(graph);
    string gname = "(d,p&a) S-factor";
    MakeGraph(gname.c_str());
  } else {
    cout<<"Excitation function graph does not exist!! Check gname passed to SFactor::Run()"
        <<endl;
  }
  file->Close();
}
