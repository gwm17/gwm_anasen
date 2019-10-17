#include "dalitz.h"
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TPad.h>

DalAnalyzer::DalAnalyzer(float min, float max, float step) {
  minE = min;
  maxE = max;
  stepE = step;
  maxBin = (maxE-minE)/stepE;
  dist_err = 0.1;
  beam_eloss = NULL;
  tot_counts = 0;
}

DalAnalyzer::~DalAnalyzer() {
  delete beam_eloss;
}

int DalAnalyzer::getPoints() {
  points.resize(0);
  eff_perr.resize(0);
  if(simPlot != NULL && dataPlot != NULL) {
    int binsX = simPlot->GetNbinsX();
    int binsY = simPlot->GetNbinsY();
    if(binsX != dataPlot->GetNbinsX() || binsY != dataPlot->GetNbinsY()) {
      cout<<"Data and simulation do not have matching number of bins!! Cannot"
          <<" calculate cross section!"<<endl;
      return 0;
    }
    for (int i=0; i<binsX; i++) {
      for(int j=0; j<binsY; j++) {
        Int_t bin = simPlot->GetBin(i,j);
        Double_t eff = simPlot->GetBinContent(bin);
        Int_t counts = dataPlot->GetBinContent(bin);
        tot_counts += counts;
        /*if(eff == 0 && counts != 0) {
          depth = 1;
          eff = interpBinContent(i,j,bin);
          if(eff != 0) 
          points.push_back(make_pair(eff, counts));
        } else*/ if(eff > 1e-2 && counts != 0) {
          points.push_back(make_pair(eff, counts));
          Double_t cntsd = detPlot->GetBinContent(bin);
          Double_t cntsf = freePlot->GetBinContent(bin);
          eff_perr.push_back(sqrt(1.0/cntsd+1.0/cntsf));
        }
      }
    }
    return 1;
  } else {
    cout<<"NULL histograms!!";
    if(simPlot == NULL) cout<<" simPlot!";
    if(dataPlot == NULL) cout<<" dataPlot!";
    cout<<endl;
    return 0;
  }
}

Double_t DalAnalyzer::interpBinContent(Int_t binx, Int_t biny, Int_t gbin) {
  Double_t sum = 0;
  Double_t count = 0;
  Double_t e_error = 0;
  for(int i = binx-depth; i<binx+depth+1; i++) {
    for(int j = biny-depth; j<biny+depth+1; j++) {
      if(i <= 0 || i> simPlot->GetNbinsX()) continue;
      if(j <= 0 || j> simPlot->GetNbinsY()) continue;
      if((i<binx+depth && i>binx-depth) && (j<biny+depth && j>biny-depth)) continue;
      Int_t this_bin = simPlot->GetBin(i,j);
      Double_t value = simPlot->GetBinContent(this_bin);
      if((value == 0)) {
        continue;
      } else {
        sum += value;
        count++;
        Double_t cntsf = freePlot->GetBinContent(this_bin);
        Double_t cntsd = detPlot->GetBinContent(this_bin);
        e_error += value*value*(1.0/cntsf + 1.0/cntsd);
      }
    }
  }
  if(depth >10.0) return 0;
  if(sum == 0 || (double)sum/count <1.0e-2) {
    depth++;
    return interpBinContent(binx, biny, gbin);
  } else {
    Double_t eff = sum/count;
    e_error = sqrt(e_error)/count;
    eff_perr.push_back(e_error/eff);
    return eff;
  }
}

void DalAnalyzer::goCrossSection(Double_t minKE, Double_t maxKE) {
  Double_t distA = beam_eloss->GetDistance(maxBeamKE, minKE, dist_err);
  Double_t distB = beam_eloss->GetDistance(maxBeamKE, maxKE, dist_err);
  Double_t dist = distA-distB;
  Double_t nTParticles_per_A = dist*density;
  Double_t cs_total = 0;
  Double_t cst_err = 0;
  Double_t eff_bar = 0;
  Double_t npoints = points.size();
  for(unsigned int i=0; i<points.size(); i++) {
    Int_t counts = points[i].second;
    Double_t eff = points[i].first;
    Double_t factor = cm2mb/(eff*nTParticles_per_A*nBParticles);
    Double_t cs = factor*counts;
    Double_t cs_err = getCSError(cs,counts,eff_perr[i],dist);
    cst_err += cs_err*cs_err;
    cs_total += cs;
    eff_bar += eff;
  }
  cross.push_back(cs_total);
  cstot_err.push_back(sqrt(cst_err));
  if(points.size() != 0)avg_eff.push_back(eff_bar/npoints);
  TLorentzVector beam_LV, target_LV;
  Double_t beamE = minKE + m_beam;
  Double_t beamP = sqrt(beamE*beamE-m_beam*m_beam);
  beam_LV.SetPxPyPzE(0.,0.,beamP,beamE);
  target_LV.SetPxPyPzE(0.,0.,0.,m_target);
  TLorentzVector parent_LV = beam_LV+target_LV;
  TVector3 boost = parent_LV.BoostVector();
  parent_LV.Boost(-boost);
  Double_t ecm = parent_LV.E()-m_beam-m_target;
  Ecm.push_back(ecm);
  if(points.size() != 0) ae_Ecm.push_back(ecm);
}

Double_t DalAnalyzer::getCSError(Double_t cs,Double_t cnts,Double_t e_perr,Double_t dist) {
  Double_t value = 0;
  if(cs != 0) {
    value =  cs*sqrt(1.0/cnts+e_perr*e_perr+beam_perr*beam_perr
                     +(dist_err/dist)*(dist_err/dist));
  }
  return value;
}

void DalAnalyzer::makePlots(char *output) {
  TFile *outputFile = new TFile(output, "RECREATE");
  TGraphErrors *cs_Ecm = new TGraphErrors(cross.size(), &(Ecm[0]), &(cross[0]),0, 
                                             &(cstot_err[0]));
  cs_Ecm->GetXaxis()->SetLimits(0.0,2.4);
  cs_Ecm->SetMarkerStyle(21);
  cs_Ecm->SetMarkerColor(2);
  cs_Ecm->SetName("(d,p&a) Cross Section Ecm");
  cs_Ecm->Write();

  TGraph *egraph = new TGraph(avg_eff.size(), &(ae_Ecm[0]), &(avg_eff[0]));
  egraph->GetXaxis()->SetLimits(0.0,2.4);
  egraph->SetMarkerStyle(22);
  egraph->SetMarkerColor(4);
  
  TCanvas *c1 = new TCanvas("7Be(d,ap)","7Be(d,ap)",0,0,800,600);
  TPad *p1 = new TPad("pad1","",0,0,1,1);
  TPad *p2 = new TPad("pad2","",0,0,1,1);
  p2->SetFillStyle(4000);
  p2->SetFrameFillStyle(0);
  p1->Draw();
  p1->cd();
  cs_Ecm->Draw("ALP");

  p2->Draw();
  p2->cd();
  egraph->Draw("ALPY+");

  outputFile->cd();
  c1->Write(); 
  outputFile->Close();
}

void DalAnalyzer::makeThetaPlot(char *output) {
  TFile *outputFile = new TFile(output, "UPDATE");
  TCanvas *c1 = new TCanvas("7Be(d,a)","7Be(d,a)",0,0,800,600);
  TGraph2D *graph = new TGraph2D(cross.size(), &(Ecm[0]), &(theta[0]), &(cross[0]));
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(2);
  graph->GetXaxis()->SetTitle("Ecm");
  graph->GetYaxis()->SetTitle("theta");
  graph->GetZaxis()->SetTitle("Cross Section");
  graph->Draw("surf1");
  outputFile->cd();
  c1->Write();
  outputFile->Close();
}

void DalAnalyzer::run(char *dfilename, char *sfilename, char *output) {
  dataFile = new TFile(dfilename, "READ");
  simFile = new TFile(sfilename, "READ");
  beam_eloss = new LookUp("./srim/7be_in_d2_290torr.eloss",m_beam);
  beam_eloss->InitializeLookupTables(30.0,200.0,0.01,0.04);
  for (int i=1; i<maxBin; i++) {
    string plotname = "DP_"+to_string(i)+"bin";
    dataPlot = (TH2F*) dataFile->Get(plotname.c_str());
    plotname = "DPeff_"+to_string(i)+"bin";
    simPlot = (TH2F*) simFile->Get(plotname.c_str());
    plotname = "DPdtot_"+to_string(i)+"bin";
    detPlot = (TH2F*) simFile->Get(plotname.c_str());
    plotname = "DPftot_"+to_string(i)+"bin";
    freePlot = (TH2F*) simFile->Get(plotname.c_str());
    getPoints();
    Double_t minKE = i*stepE; Double_t maxKE = minKE+stepE;
    goCrossSection(minKE, maxKE);
  }
  makePlots(output);

  cross.resize(0);
  Ecm.resize(0);
  cstot_err.resize(0);
  
  for (int i=1; i<maxBin; i++) {
    for(int j=0; j<10; j++) {
      string plotname = "DP_"+to_string(i)+"bbin_"+to_string(j)+"tbin";
      dataPlot = (TH2F*) dataFile->Get(plotname.c_str());
      plotname = "DPeff_"+to_string(i)+"bbin_"+to_string(j)+"tbin";
      simPlot = (TH2F*) simFile->Get(plotname.c_str());
      plotname = "DPdtot_"+to_string(i)+"bbin_"+to_string(j)+"tbin";
      detPlot = (TH2F*) simFile->Get(plotname.c_str());
      plotname = "DPftot_"+to_string(i)+"bbin_"+to_string(j)+"tbin";
      freePlot = (TH2F*) simFile->Get(plotname.c_str());
      if(dataPlot != NULL && simPlot != NULL) {
        getPoints();
      }
      Double_t minKE = i*stepE; Double_t maxKE = minKE+stepE;
      Double_t t = (j-5)*0.2;
      goCrossSection(minKE, maxKE);
      theta.push_back(t);
    }
  }
  makeThetaPlot(output);

  dataFile->Close();
  simFile->Close();
}



