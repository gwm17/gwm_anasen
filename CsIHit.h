//Since included in dict, must be completely defined in the header
#ifndef CSIHIT_H
#define CSIHIT_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct SortByCsI {
  Int_t CsI_ID, mADC_ID, mADC_Ch;
  Double_t Up, Down, QQQ, Up_Shift, Down_Shift, QQQ_Shift;
  Double_t WCsI_X, WCsI_Y, WCsI_Z, WCsI_R, WCsI_Phi;
  Double_t Energy;
};

class CsIHit {

  public:
    vector<Double_t> CsI_Energy;
    Int_t NCsIHits;
    SortByCsI csi_obj;
    vector<SortByCsI> Hit;
    vector<SortByCsI> *ReadHit;

    CsIHit() {};

    void ZeroCsI_obj() {
      csi_obj.CsI_ID = -1;
      csi_obj.mADC_Ch = -1;
      csi_obj.mADC_ID = -1;

      csi_obj.Up = -10;
      csi_obj.Down = -10;
      csi_obj.QQQ = -10;
  
  
      csi_obj.Up_Shift = -10;
      csi_obj.Down_Shift = -10;
      csi_obj.QQQ_Shift = -10;
  
      csi_obj.Energy = -10;
      csi_obj.WCsI_X = -1000;
      csi_obj.WCsI_Y = -1000;
      csi_obj.WCsI_R = -1000;
      csi_obj.WCsI_Z = -1000;
      csi_obj.WCsI_Phi = -1000;
    };

    void ZeroCsIHit(){
      NCsIHits = 0;
      Hit.clear();
      ZeroCsI_obj();
    };

};



#endif
