#include "angDist.h"

using namespace std;

angDist::angDist(char *dfName, char *sfName, char *ofName, int chained) {
  datafile = new TFile(dfName, "READ");
  simfile = new TFile(sfName, "READ");
  if(chained) {
    outfile = new TFile(ofName, "UPDATE");
  } else {
    outfile = new TFile(ofName, "RECREATE");
  }
}

angDist::~angDist() {
}

void angDist::effCorrection(TH1F *dhisto, TH1F *shisto) {
  cost.resize(0);
  cor_Counts.resize(0);
  if(dhisto != NULL && shisto != NULL) {
    int binsX = dhisto->GetNbinsX();
    if(binsX != shisto->GetNbinsX()) {
      cout<<"Data and sim do not have the same number of bins!! Cannot perform correction!"
          <<endl;
      cout<<"BinsX: "<<binsX<<" other: "<<shisto->GetNbinsX()<<endl;
      return;
    }
    for(int i=0; i<binsX; i++) {
      Int_t bin = shisto->GetBin(i);
      Double_t eff = shisto->GetBinContent(bin);
      Int_t counts = dhisto->GetBinContent(bin);
      cost.push_back(shisto->GetXaxis()->GetBinUpEdge(bin));
      if(eff >5.0e-3) {
        cor_Counts.push_back((Double_t)counts/eff);
      } else {
        cor_Counts.push_back(counts);
      }
    }
  } else {
    cout<<"Histograms are NULL!"<<endl;
    if(shisto == NULL) cout<<"shisto"<<endl;
    if(dhisto == NULL) cout<<"dhisto"<<endl;
  }
}

void angDist::calcErrors(TH1F *dhisto, TH1F *fhisto, TH1F *dethisto) {
  err.resize(0);
  if(dhisto != NULL && fhisto != NULL && dethisto != NULL) {
    int binsX = dhisto->GetNbinsX();
    if(binsX != fhisto->GetNbinsX() || binsX != dethisto->GetNbinsX()) {
      cout<<"Non-matching number of bins!!"<<endl;
      return;
    }
    for(int i=0; i<binsX; i++) {
      Int_t bin = dhisto->GetBin(i);
      Double_t counts = dhisto->GetBinContent(bin);
      Double_t free_counts = fhisto->GetBinContent(bin);
      Double_t det_counts = dethisto->GetBinContent(bin);
      Double_t se = 0;
      if(counts != 0 && free_counts != 0 && det_counts != 0) {
        se = cor_Counts[i]*sqrt(1.0/counts+1.0/free_counts+1.0/det_counts);
      }
      err.push_back(se); 
    }
  } else {
    cout<<"Histograms are NULL!";
    if(dhisto == NULL) cout<<" dhisto!";
    if(fhisto == NULL) cout<<" fhisto!";
    if(dethisto == NULL) cout<<" dethisto!";
    cout<<endl;
  }
}

void angDist::makeGraph(const char *gname) {
  TGraphErrors *graph = new TGraphErrors(cost.size(),&(cost[0]),&(cor_Counts[0]),0,&(err[0]));
  graph->SetName(gname);
  outfile->cd();
  graph->Write();
}

void angDist::run() {

  TH1F *dhisto = (TH1F*) datafile->Get("AngDist_Li5_0.267_0.468");
  TH1F *shisto = (TH1F*) simfile->Get("AngDiste_5Li_0.267_0.468");
  effCorrection(dhisto, shisto);
  TH1F *fhisto = (TH1F*) simfile->Get("AngDistf_5Li_0.267_0.468");
  TH1F *dethisto = (TH1F*) simfile->Get("AngDistd_5Li_0.267_0.468");
  calcErrors(dhisto, fhisto, dethisto);
  string gname = "corAngDist_5Li_0.248_0.468";
  makeGraph(gname.c_str());

  dhisto = (TH1F*) datafile->Get("AngDist_Li5_1.11_1.2919");
  shisto = (TH1F*) simfile->Get("AngDiste_5Li_0.1.11_1.2919");
  effCorrection(dhisto,shisto);
  fhisto = (TH1F*) simfile->Get("AngDistf_5Li_1.11_1.2919");
  dethisto = (TH1F*) simfile->Get("AngDistd_5Li_1.11_1.2919");
  calcErrors(dhisto, fhisto, dethisto);
  gname = "corAngDist_5Li_1.11_1.2919";
  makeGraph(gname.c_str());

  dhisto = (TH1F*) datafile->Get("AngDist_Be8_1.07_1.244");
  shisto = (TH1F*) simfile->Get("AngDiste_8Be_1.07_1.244");
  effCorrection(dhisto,shisto);
  fhisto = (TH1F*) simfile->Get("AngDistf_8Be_1.07_1.244");
  dethisto = (TH1F*) simfile->Get("AngDistd_8Be_1.07_1.244");
  calcErrors(dhisto, fhisto, dethisto);
  gname = "corAngDist_8Be_1.07_1.244";
  makeGraph(gname.c_str()); 

  dhisto = (TH1F*) datafile->Get("AngDist_6Li");
  shisto = (TH1F*) simfile->Get("AngDiste_6Li");
  effCorrection(dhisto, shisto);
  fhisto = (TH1F*) simfile->Get("AngDistf_6Li");
  dethisto = (TH1F*) simfile->Get("AngDistd_6Li");
  calcErrors(dhisto, fhisto, dethisto);
  gname = "corAngDist_6Li";
  makeGraph(gname.c_str());

  if(datafile->IsOpen()) datafile->Close();
  if(simfile->IsOpen()) simfile->Close();
  if(outfile->IsOpen()) outfile->Close();
}

