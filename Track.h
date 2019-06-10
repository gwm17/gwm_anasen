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

/*TrackEvent
 *Stores tracking info including Si, PC, CsI,
 *interaction point, reaction geometery and energetics, 
 *and in the event of a two particle event info on both ejectiles from reconstruction
 */
struct TrackEvent {

  Int_t TrackType = -10, HitType = -10;    
 
  /****Si variables***/ 
  Double_t SiEnergy = -10, SiTime = -10, SiZ = -10, SiX = -10, SiY = -10, SiR = -10, SiPhi = -10,
           SiBCh = -10; //added by M.A. 02/23/2017
  Int_t DetID = -10;

  /*****PC variables****/
  //Z added by M.A.01/16/2017
  Double_t PCEnergy = -10, PCZ = -10, PCX = -10, PCY = -10, PCR = -10, Z, PCPhi = -10;
  Int_t WireID = -10;
  Double_t Down = -10, Up = -10, DownVoltage = -10, UpVoltage = -10;
  Double_t pcz_ref = -10; //For calibrations

  /****Tracking Variables*****/
  //Interaction point coords (IntPoint is Z coord (beam-axis))
  Double_t IntPoint = -10, IntPoint_X = -10, IntPoint_Y = -10;
  Double_t DiffIntPoint = -10, DiffIntPoint_X = -10, DiffIntPoint_Y = -10;
 
  //Geometery 
  Double_t PathLength = -10, Theta = -10, Phi = -10;

  //Energetics
  Double_t BeamEnergy = -10, EnergyLoss = -10, LightParEnergy = -10, BeamQvalue = -10, ThetaQvalue = -10, HeEnergyQvalue = -10,
           ResidualEn = -10, ResidualEn_20Ne = -10;

  /****CsI variables****/
  Double_t CsIEnergy = -10, CsI_Phi = -10, CsI_X = -10, CsI_Y = -10, CsI_R = -10;    
  Int_t CsI_ID = -10; 

  /****for 2 ejectile (particle) events*****/
  Double_t BeamE_p1 = -10, BeamE_p2 = -10, Ex_p1 = -10, Ex_p2 = -10, ResE_p1 = -10, ResE_p2 = -10, SiE_p1 = -10, SiE_p2 = -10, IntP_p1 = -10, IntP_p2 = -10,
           Theta_p1 = -10, Theta_p2 = -10, SiPhi_p1 = -10, SiPhi_p2 = -10, PCZ_p1 = -10, PCZ_p2 = -10;

};

///Not sure these structs are totally necessary; basically simplified versions of each
//class' sort by hit and already in TrackEvent... Turns out it is necessary for backwards 
//compatibility with dictionaries
struct Silicon_Event {
  Int_t TrackType, DetID; 
  Double_t SiEnergy, SiTime, SiZ, SiR, SiPhi;
};

struct PropCounter_Event {
  Int_t TrackType, WireID;
  Double_t PCEnergy, PCZ, PCR, PCPhi, Down, Up, DownVoltage, UpVoltage;
};

struct CsI_Event {
  Int_t TrackType, CsI_ID;
  Double_t CsIEnergy,  CsI_Phi, CsI_X, CsI_Y, CsI_R;    
};
//////////

/*RecoilEvent
 *Storage for reconstructed recoil information; note that there are several ways that the 
 *reconstruction is performed (with tracking, without tracking, ect.)
 */
struct RecoilEvent {

  Double_t IntPoint, SiEnergy_tot, SiEnergy_calc, PCEnergy_tot,
           Energy_eject_tot, Theta, Phi, KE_recoil, Ex_recoil, Ex_temp_small, Ex_temp_large,
           Ex_temp_qvalue_small, Ex_temp_qvalue_large, Delta_Theta, Delta_Phi, Ex_recoil_2a, Ex_recoil_a1p, Ex_recoil_a2p;
  
  Double_t BeamKE, BeamKE_WAIntP, BeamKE_eject, theta_eject1, theta_eject2;
  Double_t BeamPz, a1_ip, a2_ip, p_ip, recoil_mass_sq;
    
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
