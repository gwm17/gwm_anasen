/*Reconstruct.h
 *Class for implementing kinematic reconstruction of reaction events in ANASEN detector
 *Should be customized for each unique experiment
 *For use in the analyzer class 
 *
 * Gordon M. -- April 2019
 * Based on previous version by M. Anastasiou
 */

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
    Double_t m_beam, m_target, m_ejectile, m_recoil, m_ejectile2;
    LookUp* ELoss_eject;
    LookUp* ELoss_eject2;
    LookUp* ELoss_beam;
    const int MaxPCHits = 24, MaxSiHits = 500, MaxTracks = 100, NPCWires = 24;
    const float BeamE = 17.19;

    /********Construction Area********/
    Reconstruct();
    Reconstruct(Double_t beam, Double_t target, Double_t ejectile, Double_t recoil);
    Reconstruct(Double_t beam, Double_t target, Double_t eject, Double_t eject2, Double_t recoil);
    ~Reconstruct();

    /*********Functions*************/
    void SetELossFile(const char* name1, const char* name2, const char* name3);
    Double_t CalcElastic(Double_t &EnergyProj,Double_t &EjectTheta,Track Tr,Int_t c);
    void CalcRecoil(Track Tr, Int_t proton, RecoilEvent &recoil);
    Double_t DeltaZ(Double_t Zpc, Double_t Zsi, Double_t Rpc, Double_t Rsi);
    Double_t WIntP(Track Tr, vector<Int_t> IsEject);
    void ThetaPath(Double_t &Theta, Double_t &Path, Double_t Zsi, Double_t Rsi, Track Tr,
                   vector<Int_t> IsEject);
    void CalcRecoil_MultiParticle(Track Tr, vector<Int_t> IsEject, RecoilEvent &recoil_qval,
                                  RecoilEvent &recoil);
    void CalcRecoil_from_Qvalue(Double_t Qvalue, Track Tr, Int_t proton,RecoilEvent &recoil_qval);
    void CalcE_ex(Track Tr, vector<int> IsEject, RecoilEvent &recoil, RecoilEvent &recoil_qval);
    void CalcRecoil_InvMassSelect(Track Tr, vector<int> IsEject, RecoilEvent &recoil, 
                                  RecoilEvent &recoil_qval);
};


#endif
