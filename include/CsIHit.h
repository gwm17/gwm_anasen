/*PCHit.h
 *Class and structure for generating ANASEN tracking dictionary 
 *in ROOT. This one is for the CsI, contains the larger CsIHit
 *class along with a SortByCsI structure for use in other codes
 *Some forced inheritance due to previous stages of analysis
 *
 *MUST be entirely defined in a header file due to the way that
 *ROOT compiles a dictionary
 *
 *Gordon M. -- April 2019
 *Based on previous version by M. Anastasiou, N. Rijal, J. Parker, et al
 */

#ifndef CSIHIT_H
#define CSIHIT_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct csievent {
  Int_t CsI_ID, mADC_ID, mADC_Ch;
  Double_t Up, Down, QQQ, Up_Shift, Down_Shift, QQQ_Shift;
  Double_t WCsI_X, WCsI_Y, WCsI_Z, WCsI_R, WCsI_Phi;
  Double_t Energy;
};

class CsIHit {

  public:
    vector<Double_t> CsI_Energy;
    Int_t NCsIHits;
    struct SortByCsI:csievent {};
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
