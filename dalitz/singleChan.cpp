#include "singleChan.h"
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;

singleChan::singleChan(float min, float max, float step, int chained) {
  minE = min;
  maxE = max;
  stepE = step;
  maxBin = (maxE-minE)/stepE;
  minBin = 21;
  dist_err = 0.1;
  chainFlag = chained;
  beam_eloss = NULL;
}

singleChan::~singleChan() {
  delete beam_eloss;
}

int singleChan::getPoints(TH1F *dhisto, TH1F *shisto, TH1F *fhisto, TH1F *dethisto) {
  points.resize(0);
  if(shisto != NULL && dhisto != NULL) {
    int binsX = dhisto->GetNbinsX();
    if(binsX != shisto->GetNbinsX()) {
      cout<<"Non matching number of bins! Unable to compute cross section!"<<endl;
      return 0;
    }
    for(int i=0; i<binsX; i++) {
      Int_t bin = shisto->GetBin(i);
      Double_t eff = shisto->GetBinContent(bin);
      Int_t counts = dhisto->GetBinContent(bin);
      if(eff != 0) {
        points.push_back(make_pair(eff, counts));
        detCounts.push_back(dethisto->GetBinContent(bin));
        freeCounts.push_back(fhisto->GetBinContent(bin));
      }
    }
    return 1;
  } else {
    cout<<"Histograms are NULL! Unable to compute cross section!"<<endl;
    return 0;
  }
}

void singleChan::getCrossSection(Double_t minKE, Double_t maxKE) {
  Double_t distA = beam_eloss->GetDistance(maxBeamKE, minKE, dist_err);
  Double_t distB = beam_eloss->GetDistance(maxBeamKE, maxKE, dist_err);
  Double_t dist = distA-distB;
  Double_t nTParticles_per_A = dist*density;
  Double_t cs_total = 0;
  Double_t cst_err = 0;
  Double_t eff_bar = 0;
  Double_t tc = 0;
  float npoints = points.size();
  for(unsigned int i=0; i<points.size(); i++) {
    Double_t counts = points[i].second;
    Double_t eff = points[i].first;
    Double_t cs = cm2mb*counts/(eff*nTParticles_per_A*nBParticles);
    Double_t cs_err = getCSError(cs, counts, freeCounts[i], detCounts[i], dist);
    cst_err += cs_err*cs_err;
    cs_total += cs;
    eff_bar += eff;
    tc += counts;
  }
  cross.push_back(cs_total);
  cstot_err.push_back(sqrt(cst_err));
  avg_eff.push_back(eff_bar/npoints);
  tot_counts.push_back(tc);
  TLorentzVector beam_LV, target_LV;
  Double_t beamE = minKE+m_beam;
  Double_t beamP = sqrt(beamE*beamE-m_beam*m_beam);
  beam_LV.SetPxPyPzE(0.,0.,beamP,beamE);
  target_LV.SetPxPyPzE(0.,0.,0.,m_target);
  TLorentzVector parent_LV = beam_LV+target_LV;
  TVector3 boost = parent_LV.BoostVector();
  parent_LV.Boost(-boost);
  Ecm.push_back(parent_LV.E()-m_beam-m_target);
}

Double_t singleChan::getCSError(Double_t cs,Double_t cnts,Double_t cntsf,Double_t cntsd,
                                Double_t dist) {
  Double_t value = 0;
  if(cs != 0) {
    value =  cs*sqrt(1.0/cnts+1.0/cntsf+1.0/cntsd+beam_perr*beam_perr
                     +(dist_err/dist)*(dist_err/dist));
  }
  return value;
}

void singleChan::makeGraph(const char *gname) {
  TGraphErrors *graph = new TGraphErrors(cross.size(), &(Ecm[0]), &(cross[0]),
                                         0,&(cstot_err[0]));
  graph->SetName(gname);
  graph->SetTitle("3He Cross Section");
  graph->SetMarkerStyle(22);
  graph->SetMarkerColor(4);
  graph->GetXaxis()->SetLimits(0.8,2.4);
  outfile->cd();
  graph->Write();
  graph->SetTitle("");

  TGraph *egraph = new TGraph(avg_eff.size(), &(Ecm[0]), &(avg_eff[0]));
  egraph->SetTitle("");
  egraph->SetMarkerStyle(21);
  egraph->SetMarkerColor(2);
  egraph->GetXaxis()->SetLimits(0.8,2.4);

  TGraph *cgraph = new TGraph(tot_counts.size(),&(Ecm[0]),&(tot_counts[0]));
  cgraph->SetTitle("");
  cgraph->SetMarkerStyle(20);
  cgraph->SetMarkerColor(6);
  cgraph->GetXaxis()->SetLimits(0.8,2.4);
  


  TCanvas *c1 = new TCanvas("3He cs and eff","3He cs and eff",0,0,800,600);
  TPad *p1 = new TPad("pad1","",0,0,1,1);
  TPad *p2 = new TPad("pad2","",0,0,1,1);
  p2->SetFillStyle(4000);
  p2->SetFrameFillStyle(0);
  p1->Draw();
  p1->cd();
  graph->Draw("ALP");
  p2->Draw();
  p2->cd();
  egraph->Draw("ALPY+");
  TLegend *lc1 = new TLegend(0.1,0.7,0.48,0.9);
  lc1->AddEntry(graph,"Cross section");
  lc1->AddEntry(egraph,"Average Efficiency");
  lc1->Draw();
  outfile->cd();
  c1->Write();

  TCanvas *c2 = new TCanvas("3He counts and eff", "3He counts and eff",0,0,800,600);
  c2->cd();
  TPad *p3 = new TPad("pad3","",0,0,1,1);
  TPad *p4 = new TPad("pad4","",0,0,1,1);
  p4->SetFillStyle(4000);
  p4->SetFrameFillStyle(0);
  p4->SetLogy();
  p3->SetLogy();
  p3->Draw();
  p3->cd();
  cgraph->Draw("ALP");
  p4->Draw();
  p4->cd();
  egraph->Draw("ALPY+");
  TLegend *lc2 = new TLegend(0.1,0.7,0.48,0.9);
  lc2->AddEntry(cgraph,"Total Counts");
  lc2->AddEntry(egraph,"Average Efficiency");
  lc2->Draw();
  outfile->cd();
  c2->Write();

}

void singleChan::run(char *dfName, char *sfName, char *ofName) {
  datafile = new TFile(dfName, "READ");
  simfile = new TFile(sfName, "READ");
  if(chainFlag) {
    outfile = new TFile(ofName, "UPDATE");
  } else {
    outfile = new TFile(ofName, "RECREATE");
  }
  beam_eloss = new LookUp("./srim/7be_in_d2_290torr.eloss",m_beam);
  beam_eloss->InitializeLookupTables(30.0,200.0,0.01,0.04);
  TH1F *tempd, *temps, *tempf, *tempdet;
  for(int i=minBin; i<maxBin+1; i++) {
    string hname = "ExLi6_"+to_string(i)+"bin";
    tempd = (TH1F*) datafile->Get(hname.c_str());
    hname = "ExLi6e_"+to_string(i) + "bin";
    temps = (TH1F*) simfile->Get(hname.c_str());
    hname = "ExLi6f_"+to_string(i)+"bin";
    tempf = (TH1F*) simfile->Get(hname.c_str());
    hname = "ExLi6d_"+to_string(i)+"bin";
    tempdet = (TH1F*) simfile->Get(hname.c_str());
    getPoints(tempd, temps, tempf, tempdet);
    Double_t minKE = i*stepE; Double_t maxKE = minKE+stepE;
    getCrossSection(minKE, maxKE);
  }
  string gname("3He CS");
  makeGraph(gname.c_str());
  if(datafile->IsOpen()) datafile->Close();
  if(simfile->IsOpen()) simfile->Close();
  if(outfile->IsOpen()) outfile->Close();
}

