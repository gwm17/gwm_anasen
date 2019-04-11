//Since to be included in dict, must be completely defined in header
#ifndef SIHIT_H
#define SIHIT_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct SiDetector {
  Int_t DetID, UpMult, DownMult, FrontMult, BackMult, HitType;
  vector<Int_t> UpChNum, DownChNum, FrontChNum, BackChNum;

  //Raw Energy (units of channels)
  vector<Double_t> EUp_Raw, EDown_Raw, EFront_Raw, EBack_Raw;

  //Energy after pulser cals (units of channels)
  vector<Double_t> EUp_Pulser, EDown_Pulser, EFront_Pulser, EBack_Pulser;

  //Energy after relative calibration (units of channels)
  vector<Double_t> EUp_Rel, EDown_Rel, EFront_Rel, EBack_Rel, SX3_ZUp, SX3_ZDown,
                   EUp_RelShift, EDown_RelShift, EBack_RelShift;
  //Energy after energy cal (units of MeV)
  vector<Double_t> EUp_Cal, EDown_Cal, EFront_Cal, EBack_Cal, TUp, TDown, TFront, TBack;
};

struct sievent {
  Int_t NHitsInDet, DetID, HitType;
  Int_t FrontChannel, BackChannel, UpChannel, DownChannel;
  Double_t EnergyBack, EnergyFront, EnergyUp, EnergyDown, Energy,
           Time, X, Y, Z, Z_linear, ZUp_Dummy, ZDown_Dummy, XW, YW, ZW, RW, PhiW;
  Int_t TrackType;
  Double_t RFTime, MCPTime;
};

class SiHit{

  public:

    Int_t NSiHits;
    struct SortByDetector:SiDetector {};
    struct SortByHit:sievent {};
    SortByDetector det_obj;
    SortByHit hit_obj;
    vector<SortByDetector> Detector;
    vector<SortByHit> Hit;
    vector<SortByDetector> *ReadDet;
    vector<SortByHit> *ReadHit; 

    SiHit() {};

    void ZeroSi_obj() {
      det_obj.DetID = -1;
      det_obj.UpMult = 0;
      det_obj.DownMult = 0;
      det_obj.BackMult = 0;
      det_obj.FrontMult = 0;
      det_obj.HitType = 0;

      det_obj.UpChNum.clear();
      det_obj.DownChNum.clear();
      det_obj.FrontChNum.clear();
      det_obj.BackChNum.clear();

      det_obj.EUp_Raw.clear();
      det_obj.EDown_Raw.clear();   
      det_obj.EFront_Raw.clear();
      det_obj.EBack_Raw.clear();

      det_obj.EUp_Pulser.clear();
      det_obj.EDown_Pulser.clear();    
      det_obj.EFront_Pulser.clear();
      det_obj.EBack_Pulser.clear();

      det_obj.EUp_Rel.clear();
      det_obj.EDown_Rel.clear();  
      det_obj.EFront_Rel.clear();
      det_obj.EBack_Rel.clear();

      det_obj.EUp_RelShift.clear();
      det_obj.EDown_RelShift.clear();  
      det_obj.EBack_RelShift.clear();

      det_obj.SX3_ZUp.clear();
      det_obj.SX3_ZDown.clear();

      det_obj.EUp_Cal.clear();
      det_obj.EDown_Cal.clear();    
      det_obj.EFront_Cal.clear();
      det_obj.EBack_Cal.clear();

      det_obj.TUp.clear();
      det_obj.TDown.clear();    
      det_obj.TFront.clear();
      det_obj.TBack.clear();

      hit_obj.NHitsInDet = 0;
      hit_obj.DetID = -1;
  
      hit_obj.HitType = 0;
      hit_obj.FrontChannel = -1;
      hit_obj.BackChannel = -1;
      hit_obj.UpChannel = -1;
      hit_obj.DownChannel = -1;
      hit_obj.EnergyBack = -1000;
      hit_obj.EnergyFront = -1000;
      hit_obj.EnergyUp = -1000;
      hit_obj.EnergyDown = -1000;
      hit_obj.Energy = -1000;
      hit_obj.Time = 0;
      hit_obj.X = 0;
      hit_obj.Y = 0;
      hit_obj.Z = -10;
      hit_obj.Z_linear = -10; // added 12/09/2016 for the ZPosCal_linear
      hit_obj.ZUp_Dummy = -10;
      hit_obj.ZDown_Dummy = -10;
      hit_obj.XW = 0;
      hit_obj.YW = 0;
      hit_obj.ZW = -10;
      hit_obj.RW = 0;
      hit_obj.PhiW = 0;
      hit_obj.TrackType = 0;
        //hit_obj.RFSubtract = 0;
    };

    void ZeroSiHit() {
      NSiHits = 0;
      Detector.clear();
      Hit.clear();  
      ZeroSi_obj();
    };
 
};

#endif
