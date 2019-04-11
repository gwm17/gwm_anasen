#ifndef __LOOKUP_H__
#define __LOOKUP_H__
#include "TROOT.h"
#include "TGraph.h"
#include <iostream>
#include <string>
#include <fstream>
/////////////////////////////////////////////////////////
//////////////Code: LookUp.h/////////////////////////////
////////////Date: June28_2016////////////////////////////
///////////Author: Nabin Rijal //////////////////////////
/////////////////////////////////////////////////////////
using namespace std;

class LookUp{

  public:

    LookUp();
    LookUp(string Eloss_file, Double_t InputMass);
    ~LookUp();
    Double_t GetEnergyLoss(Double_t initial_energy, Double_t distance);
    Double_t GetInitialEnergy(Double_t FinalEnergy, Double_t PathLength, Double_t StepSize);
    Double_t GetFinalEnergy(Double_t InitialEnergy, Double_t PathLength, Double_t StepSize);
    Double_t GetDistance(Double_t InitialE, Double_t FinalE, Double_t StepSize);
    Double_t GetPathLength(Float_t InitialEnergy, Float_t FinalEnergy, Float_t DeltaT);
    Double_t LoadRange(Float_t energy1);
    Double_t GetTimeOfFlight(Double_t InitialEnergy, Double_t PathLength, Double_t StepSize);
    void SetIonMass(Double_t IonMass);  
    void InitializeLookupTables(Double_t MaximumEnergy, Double_t MaximumDistance, Double_t DeltaE, Double_t DeltaD);
    void PrintLookupTables();
    Double_t GetLookupEnergy(Double_t InitialEnergy, Double_t distance);
    bool GoodELossFile;
    TGraph* EvD;

  private:

    Double_t c;
    Double_t IonMass;
    Double_t IonEnergy;
    Double_t dEdx_e;
    Double_t dEdx_n;
    Double_t EtoDtab;
    Double_t DtoEtab; 

    vector<Double_t> IonEnergy_v;
    vector<Double_t> dEdx_e_v;
    vector<Double_t> dEdx_n_v;

    Double_t MaximumEnergy;
    Double_t MaximumDistance;
    Double_t DeltaD;
    Double_t DeltaE;

    vector<Double_t> EtoDtab_v;
    vector<Double_t> DtoEtab_v;

    int points;
    int last_point;
    int points1;
    int last_point1;
    bool Energy_in_range;
};


#endif 
//////////////////////////////////////////////////////////////////////////
