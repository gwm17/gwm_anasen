#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H

#include <TROOT.h>
#include <vector> 
#include <TList.h>
#include "LookUp.h"
#include "Track.h"
#include "CsIHit.h"
#include "SiHit.h"
#include "PCHit.h"

using namespace std;

class Reconstruct {

  public:
    /*******Class wide variables********/
    Double_t m_beam, m_target, m_ejectile, m_recoil;
    LookUp* ELoss_light;
    LookUp* ELoss_beam;
    const int MaxPCHits = 24, MaxSiHits = 500, MaxTracks = 100, NPCWires = 24;

    /********Construction Area********/
    Reconstruct();
    Reconstruct(Double_t beam, Double_t target, Double_t light, Double_t heavy);
    ~Reconstruct();

    /*********Functions*************/
    void SetELossFile(const char* name1, const char* name2);
    Double_t ReconstructElastic(Double_t &EnergyIncoming,Double_t &LightTheta,Track Tr,Int_t c);
    void ReconstructHeavy(Track Tr, Int_t proton, RecoilEvent &recoil);
    Double_t DeltaZ(Double_t Zpc, Double_t Zsi, Double_t Rpc, Double_t Rsi);
    Double_t WIntP(Track Tr, vector<Int_t> IsProton);
    void ThetaPath(Double_t &Theta, Double_t &Path, Double_t Zsi, Double_t Rsi, Track Tr,
                   vector<Int_t> IsProton);
    void Reconstruct20Ne(Track Tr, vector<Int_t> IsProton, Double_t Qvalue, 
                                bool IntPnt_flag, RecoilEvent &recoil);
    void ReconstructHeavy_Qvalue(Double_t Qvalue, Track Tr, Int_t proton, RecoilEvent &recoil);

};


#endif
