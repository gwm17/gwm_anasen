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
    void getCut(string pcutfile, string acutfile, string he3cutfile, string dcutfile, string 
                jacutfile);
    void recoilReset(RecoilEvent &recoil);
    bool RecoverTrack1(TrackEvent track1, TrackEvent &track2);
    bool MCP_RF();
    void Track1();
    void Track2();
    void Track3();
    void PCPlotting();
    void EdEcor();
    void TrackCalc();
    void PCWireCalibration();
    void PCPlotting2();
    void CalcRecoilE();
    void CalcSmAngleAlphas();
    void ReconstructMe();
    void ElasticScatter();
    void MyFill(string name, int binsX, double lowX, double highX, double valueX);
    void MyFill(string name, int binsX, double lowX, double highX, double valueX,
                int binsY, double lowY, double highY, double valueY);
    vector<Double_t> GetPCWireRadius();
    Double_t PhiDiff(Float_t phi1, Float_t phi2);
    
    
    /************Flags****************/
    int PCPlots, PCPlots_2, Beam_and_Eloss, FillTree, FillEdE_cor, MCP_RF_Cut,
        ReadPCWire, PCWireCal, ResidualEnergy_Calc, ResidualEnergy_Calc_20Ne, 
        Reconstruction_Session, Elastic_Scat_Alpha, cutFlag;

    /***********Class Constants******/
    //Detector Class Parameters
    const int MaxSiHits = 500, MaxADCHits = 500, MaxTDCHits = 500, MaxTracks = 100, 
              MaxCsIHits = 500, DiffIP = 2;
    const static int NPCWires = 24;
    string be7eloss_name, he3eloss_name, he4eloss_name, peloss_name, deloss_name;

    //Conversions, physical measurements
    const float rads2deg = 57.27272727, pcr = 3.846264509, ana_length = 55.0545;

    //Nuclear Masses (MeV)
    const float m_p = 938.27206671856, m_alpha = 3727.37929745092, m_7be = 6534.1836677282,
                m_3he = 2808.3915032078, m_8be = 7454.85043438849, m_6li = 5601.518452737,
                m_7li = 6533.83277448969, m_d = 1875.61291385342, m_5he = 4667.67970996292,
                m_n = 1875.61291385342, m_5li = 4667.6163636366931;

    //Beam & Reaction Parameters (MeV)
    const float BeamE = 17.19, QValue_8Be = 16.674, QValue_6Li = -0.1133726; 

    /**********Private Variables***********/
    //Storage
    TList *fhlist;
    map<string, TH1*> fhmap;
    TCutG *protonCut, *alphaCut, *he3Cut, *deutCut, *joinedAlphaCut;
    TObjArray *rootObj;
    vector<Double_t> WireRadii;
    
    //Detectors
    SiHit Si;
    PCHit PC;
    CsIHit CsI;

    //Energy loss lookups
    LookUp *be7_eloss, *alpha_eloss, *proton_eloss, *deuteron_eloss, *he3_eloss;

    //Tree and Histogram Parameters
    Track tracks;
    Int_t RFTime, MCPTime, TDC2, input_RFTime, input_MCPTime, input_TDC2;
    
    RecoilEvent Li6, Be8_1p, Be8_1p_2a, Li6_qval, Be8_1p_qval, Be8_1p_2a_qval, Be8_1p_any,
                Be8_1p_any_qval, Li5, Li5_qval;

    Double_t PCGoodEnergy[NPCWires], PCGoodPCZ[NPCWires];
    vector<vector<Double_t>> SiEnergy_vec;
    Double_t SiGoodEnergy[28];

};


#endif
