/*gwmAnalyzer.h
 *Analyzer class for ANASEN Detector in Active Target mode. Contains all functions necessary to 
 *sort data and calculate physical values. Will do everything from making all of the tracking data
 *up to calling the experiment dependant reconstruction and sorting and storing all data. Current
 *asks user for the name of a data list file, an output file, and three cut files
 *
 *Gordon M. -- April 2019
 *Based on previous versions written by M. Anastasiou, N. Rijal, J. Parker, et al
 */
#ifndef GWM_ANALYZER_H
#define GWM_ANALYZER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <map>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom1.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "Track.h"
#include "SiHit.h"
#include "PCHit.h"
#include "CsIHit.h"
#include "LookUp.h"
#include "Reconstruct.h"

using namespace std;

class analyzer {

  public:
    analyzer();
    ~analyzer();
    void run();
    void SetFlag(string flagName);

  private:
    /**********Private Functions************/
    Int_t FindMaxPC(Double_t phi, PCHit& PC);
    Int_t FindGoodCsI(Double_t phi, CsIHit& CsI);
    Int_t FindMaxSi(Double_t phi, SiHit& Si);
    vector<Int_t> FindGoodCsI_Vector(Double_t phi, CsIHit& CsI);
    void getCut();
    void recoilReset();
    bool MCP_RF();
    void Track1();
    void Track2();
    void Track3();
    void PCPlotting();
    void EdEcor();
    void TrackCalc();
    void PCWireCalibration();
    void PCPlotting2();
    void CalculateResidE();
    void ReconstructMe();
    void AlphaScatter();
    void MyFill(string name, int binsX, double lowX, double highX, double valueX);
    void MyFill(string name, int binsX, double lowX, double highX, double valueX,
                int binsY, double lowY, double highY, double valueY);
    vector<Double_t> GetPCWireRadius();
    Double_t PhiDiff(Float_t phi1, Float_t phi2);
    
    
    /************Flags****************/
    int PCPlots, PCPlots_2, Beam_and_Eloss, FillTree, FillEdE_cor, CheckBasic, MCP_RF_Cut,
        ReadPCWire, PCWireCal, ResidualEnergy_Calc, ResidualEnergy_Calc_20Ne, 
        Reconstruction_Session, Elastic_Scat_Alpha, cutFlag;

    /***********Class Constants******/
    //Detector Class Parameters
    const int MaxSiHits = 500, MaxADCHits = 500, MaxTDCHits = 500, MaxTracks = 100, 
              MaxCsIHits = 500, DiffIP = 2;
    const static int NPCWires = 24;

    //Conversions, physical measurements
    const float rads2deg = 57.27272727, pcr = 3.846264509, ana_length = 55.0545;

    //Nuclear Masses (MeV)
    const float M_p = 938.27197, M_alpha = 3727.37892, M_160 = 14895.079, M_19F = 17692.29956,
                M_18Ne = 16767.09961, M_21Na = 19553.56837, M_14N = 13040.20242, 
                M_17O = 15830.50124, M_20Ne = 18617.72807;

    //Beam & Reaction Parameters (MeV)
    const float BeamE = 72.34, QValue = 2.63819, QValue_20Ne = 0.206519; 

    /**********Private Variables***********/
    //Storage
    TList *fhlist;
    map<string, TH1*> fhmap;
    TCutG *protonCut, *alphaCut, *needleCut;
    TObjArray *rootObj;
    vector<Double_t> WireRadii;
    
    //Detectors
    SiHit Si;
    PCHit PC;
    CsIHit CsI;

    //Energy loss lookups
    LookUp *Ne18_eloss, *Na21_eloss, *Ne20_eloss, *alpha_eloss, *proton_eloss;

    //Tree and Histogram Parameters
    Track tracks;
    Int_t RFTime, MCPTime, TDC2, input_RFTime, input_MCPTime, input_TDC2;
    
    
    Double_t IC_needle, R_IC_needle, ICne_E_diff, ICne_E_sum, ICne_T_diff, ICne_T_sum,
             input_ICne_E_diff, input_ICne_E_sum, input_ICne_T_diff, input_ICne_T_sum, 
             ICne_E_diff_cal;
    RecoilEvent recoil;

    Double_t PCGoodEnergy[NPCWires], PCGoodPCZ[NPCWires];
    vector<vector<Double_t>> SiEnergy_vec;
    Double_t SiGoodEnergy[28];

};


#endif
