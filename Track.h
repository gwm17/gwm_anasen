/*Track.h
 *Class and structure for generating ANASEN tracking dictionary 
 *in ROOT. This one is for the tracking, contains the larger Track
 *class along with TrackEvent and RecoilEvent structure for use in other codes
 *
 *MUST be entirely defined in a header file due to the way that
 *ROOT compiles a dictionary
 *
 *Gordon M. -- April 2019
 *Based on previous version by M. Anastasiou, N. Rijal, J. Parker, et al
 */

#ifndef TRACK_H
#define TRACK_H

#include <TROOT.h>
#include <vector>

using namespace std;

struct TrackEvent {

  Int_t TrackType;
  Int_t HitType;
  
  Double_t SiEnergy;
  Double_t SiTime;
  Double_t SiZ;
  Double_t SiX;
  Double_t SiY;
  Double_t SiR;
  Double_t SiPhi;
  Double_t SiBCh; //added by M.A. 02/23/2017
  Int_t DetID;    

  Double_t PCEnergy;    
  Double_t PCZ;
  Double_t PCX; 
  Double_t PCY;
  Double_t PCR;
  Double_t Z; //By M.A.01/16/2017
  Double_t PCPhi;
  Int_t WireID;

  Double_t Down, Up;
  Double_t DownVoltage, UpVoltage;
  

  Double_t IntPoint;
  Double_t IntPoint_X;
  Double_t IntPoint_Y;

  Double_t DiffIntPoint;
  Double_t DiffIntPoint_X;
  Double_t DiffIntPoint_Y;
  
  Double_t PathLength;

  Double_t Theta;
  Double_t Phi;

  Double_t BeamEnergy;
  Double_t EnergyLoss;
  Double_t LightParEnergy;
  Double_t BeamQvalue;
  Double_t ThetaQvalue;
  Double_t HeEnergyQvalue;
  Double_t ResidualEn;
  Double_t ResidualEn_20Ne;

  Double_t PCZ_Ref;
  Double_t pcz_ref;

  //-----CsI variables---------------//

  Double_t CsIEnergy; 
  Double_t CsI_Phi; 
  Double_t CsI_X; 
  Double_t CsI_Y; 
  Double_t CsI_R;    
  Int_t CsI_ID; 

  //-------for 2 proton events-----------//

  Double_t BeamE_p1;
  Double_t BeamE_p2;
  Double_t Ex_p1;
  Double_t Ex_p2;
  Double_t ResE_p1;
  Double_t ResE_p2;
  Double_t SiE_p1;
  Double_t SiE_p2;
  Double_t IntP_p1;
  Double_t IntP_p2;
  Double_t Theta_p1;
  Double_t Theta_p2;
  Double_t SiPhi_p1;
  Double_t SiPhi_p2;
  Double_t PCZ_p1;
  Double_t PCZ_p2;

};

///Not sure these structs are totally necessary; basically simplified versions of each
//class' sort by hit and already in TrackEvent...
struct Silicon_Event {
  Int_t TrackType;
  Int_t DetID;
  Double_t SiEnergy;
  Double_t SiTime;
  Double_t SiZ;
  Double_t SiR;
  Double_t SiPhi;
};

struct PropCounter_Event {
  Int_t TrackType;
  Int_t WireID;
  Double_t PCEnergy;
  Double_t PCZ;
  Double_t PCR;
  Double_t PCPhi;
  Double_t Down, Up;
  Double_t DownVoltage, UpVoltage;
};

struct CsI_Event {
  Int_t TrackType;
  Int_t CsI_ID;
  Double_t CsIEnergy; 
  Double_t CsI_Phi; 
  Double_t CsI_X; 
  Double_t CsI_Y; 
  Double_t CsI_R;    
};
//////////

struct RecoilEvent {

  Double_t IntPoint;
  Double_t BeamEnergy;
  Double_t BeamEnergy_20Ne;
  Double_t SiEnergy_tot;
  Double_t SiEnergy_calc;
  Double_t PCEnergy_tot;
  Double_t Energy_tot;
  Double_t Theta;
  Double_t Theta_20Ne;
  Double_t Theta_Qvalue_20Ne;
  Double_t Phi;
  Double_t Phi_20Ne;
  Double_t Phi_Qvalue_20Ne;
  Double_t KE;
  Double_t Ex;
  Double_t KE_20Ne;
  Double_t KE_Qvalue_20Ne;
  Double_t Ex_20Ne;
  Double_t Ex_Qvalue_20Ne;
  Double_t Ex_rec_small;
  Double_t Ex_rec_large;
  Double_t Ex_rec_qvalue_small;
  Double_t Ex_rec_qvalue_large;
  Double_t Delta_Theta;
  Double_t Delta_Phi;

  Double_t BeamWA_20Ne;
  Double_t BeamWA_Qvalue_20Ne;
  Double_t IntP_20Ne;
  Double_t IntP_Qvalue_20Ne;
  Double_t theta_p1_20Ne;
  Double_t theta_p2_20Ne;

  Double_t Ex_21Na;
  Double_t Beam_Qv_21Na;
    
};

class Track {
  
  public:

    Int_t NTracks, NTracks1, NTracks2, NTracks3, NTracks4,
          protonCounter, counterTrack1Cut;

    //stores number of tracks: max 3 (track1, track2, track3)
    vector<TrackEvent> TrEvent;
    vector<TrackEvent> *ReadTrEvent;
    vector<Silicon_Event> SiEvent;
    vector<Silicon_Event> *ReadSiEvent;
    vector<PropCounter_Event> PCEvent;
    vector<PropCounter_Event> *ReadPCEvent;
    vector<CsI_Event> CsIEvent;
    vector<CsI_Event> *ReadCsIEvent;

    Track() {};

    void ZeroTrack() {
      NTracks = 0;
      NTracks1 = 0;
      NTracks2 = 0;
      NTracks3 = 0;
      NTracks4 = 0;
  
      protonCounter = 0;
      counterTrack1Cut = 0;
  
      TrEvent.clear();
      SiEvent.clear();
      PCEvent.clear();
      CsIEvent.clear();
    };

    static bool Tr_Sisort_method(TrackEvent a, TrackEvent b) {
      if(a.SiEnergy > b.SiEnergy) return 1;
      else return 0;
    };

    static bool Tr_PCsort_method(TrackEvent a, TrackEvent b) {
      if(a.PCEnergy > b.PCEnergy) return 1;
      else return 0;
    };

    static bool Tr_CsIsort_method(TrackEvent a, TrackEvent b) {
      if(a.CsIEnergy > b.CsIEnergy) return 1;
      else return 0; 
    };
};

#endif
