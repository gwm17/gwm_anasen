/*PCHit.h
 *Class and structure for generating ANASEN tracking dictionary 
 *in ROOT. This one is for the PC, contains the larger PCHit
 *class along with a pcevent structure for use in other codes
 *Some forced inheritance due to previous stages of analysis
 *
 *MUST be entirely defined in a header file due to the way that
 *ROOT compiles a dictionary
 *
 *Gordon M. -- April 2019
 *Based on previous version by M. Anastasiou, N. Rijal, J. Parker, et al
 */

#ifndef PCHIT_H
#define PCHIT_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct pcevent {
  Int_t WireID;
  Double_t Down, Up, DownVoltage, UpVoltage, Energy,
           Z, XW, YW, ZW, RW, PhiW;
  Int_t TrackType;
};

class PCHit {
  
  public: 
    Int_t NPCHits;
    struct SortByPC:pcevent {};
    SortByPC pc_obj;
    vector<PCHit::SortByPC> Hit;
    vector<SortByPC> *ReadHit;

    PCHit() {};

    void ZeroPC_obj() {
      pc_obj.WireID = -1;
      pc_obj.Down = -10.0;
      pc_obj.Up = -10.0;
      pc_obj.DownVoltage = -10.0;
      pc_obj.UpVoltage = -10.0;
      pc_obj.Energy = -10.0;
      pc_obj.Z = -100.0;
      pc_obj.XW = 0.0;
      pc_obj.YW = 0.0;
      pc_obj.ZW = -10.0;
      pc_obj.RW = 0.0;
      pc_obj.PhiW = 0.0;
      pc_obj.TrackType = 0;
    };
   
    void ZeroPCHit() {
      NPCHits = 0;
      Hit.clear();
      ZeroPC_obj();
    };
};



#endif
