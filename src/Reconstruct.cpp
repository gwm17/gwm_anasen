/*Reconstruct.cpp
 *Class for implementing kinematic reconstruction of reaction events in ANASEN detector
 *Should be customized for each unique experiment
 *For use in the analyzer class 
 *
 * Gordon M. -- April 2019
 * Based on previous version by M. Anastasiou
 */

#include "Reconstruct.h"
#include "TLorentzVector.h"

Reconstruct::Reconstruct() {
  m_beam = 0.0;
  m_target = 0.0;
  m_ejectile = 0.0;
  m_recoil = 0.0;
  ELoss_light = NULL;
  ELoss_beam = NULL;
}

Reconstruct::Reconstruct(Double_t beam,Double_t target,Double_t light,Double_t heavy) {
  m_beam = beam;
  m_target = target;
  m_ejectile = light;
  m_recoil = heavy;
  ELoss_light = NULL;
  ELoss_beam = NULL;
}

Reconstruct::~Reconstruct() {
  delete ELoss_light;
  delete ELoss_beam;
}

void Reconstruct::SetELossFile(const char* name1, const char* name2) {
  ELoss_light = new LookUp(name1, m_ejectile);
  ELoss_beam = new LookUp(name2, m_beam);
}

Double_t Reconstruct::ReconstructElastic(Double_t &EnergyIncoming, Double_t &LightTheta,
                                         Track Tr, Int_t c) {
  Double_t Elight_si = Tr.TrEvent[c].SiEnergy;
  Double_t d_light = Tr.TrEvent[c].PathLength;

  //Get the inital energy of the light projectile
  Double_t Elight_rxn = ELoss_light->GetInitialEnergy(Elight_si, d_light, 0.1);

  //Calculate the energy of the beam initially, and the angle of the scattered beam
  EnergyIncoming = ((m_beam+m_ejectile)*(m_beam+m_ejectile)*Elight_rxn)/
                   ((4.0*Tr.TrEvent[c].Theta*Tr.TrEvent[c].Theta*m_beam*m_ejectile));
  LightTheta = asin((sin(Tr.TrEvent[c].Theta)*sqrt(m_ejectile/m_beam))/
                    (sqrt(EnergyIncoming/Elight_rxn-1.0)));
  return Elight_rxn;
}
                            
void Reconstruct::ReconstructHeavy(Track Tr, Int_t proton, RecoilEvent &recoil) {

  TrackEvent ejectile = Tr.TrEvent[proton];

  TLorentzVector Light_LV(0., 0., 0., 0.), Beam_LV(0., 0., 0., 0.), Target_LV(0., 0., 0., 0.),
                 Parent_LV(0., 0., 0., 0.), Heavy_LV(0., 0., 0., 0.);

  Float_t Elight_si = ejectile.SiEnergy;
  Float_t d_light = ejectile.PathLength;
  Float_t theta_light = ejectile.Theta;
  Float_t phi_light = ejectile.SiPhi;
  Float_t Elight_rxn = ELoss_light->GetLookupEnergy(Elight_si, -d_light);
  Float_t Elight_tot = Elight_rxn + m_ejectile;
  Float_t Plight = sqrt(Elight_tot*Elight_tot - m_ejectile*m_ejectile);

  Light_LV.SetPxPyPzE(Plight*sin(theta_light)*cos(phi_light), 
                      Plight*sin(theta_light)*sin(phi_light),
                      Plight*cos(theta_light),
                      Elight_tot);
  
  Float_t BeamE_tot = ejectile.BeamEnergy + m_beam;
  Float_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - m_beam*m_beam);
  Float_t TargetE = m_target;

  Beam_LV.SetPxPyPzE(0.0, 0.0, BeamP_z, BeamE_tot);
  Target_LV.SetPxPyPzE(0.0,0.0,0.0, TargetE);

  Parent_LV = Beam_LV + Target_LV;

  Heavy_LV = Parent_LV - Light_LV;

  recoil.IntPoint = ejectile.IntPoint;
  recoil.BeamEnergy = ejectile.BeamEnergy;
  recoil.SiEnergy_tot = ejectile.SiEnergy;
  recoil.PCEnergy_tot = ejectile.PCEnergy*ejectile.Theta; //Why Theta here???
  recoil.Energy_tot = Elight_rxn;
  recoil.Theta = Heavy_LV.Theta()*180.0/TMath::Pi();
  recoil.Phi = Heavy_LV.Phi()*180.0/TMath::Pi();
  recoil.KE = Heavy_LV.E()-Heavy_LV.M();
  recoil.Ex = Heavy_LV.M() - m_recoil;
  
}

Double_t Reconstruct::DeltaZ(Double_t Zpc, Double_t Zsi, Double_t Rpc, Double_t Rsi) {
  
  Double_t dZsi = -(Rpc/Rsi-Rpc);//Check to make sure that denom is correct for all of these
  Double_t dZpc = 1+dZsi; 
  Double_t dRsi = (Rpc*(Zsi-Zpc))/((Rsi-Rpc)*(Rsi-Rpc));
  Double_t dRpc = -(Rsi*(Zsi-Zpc))/((Rsi-Rpc)*(Rsi-Rpc));

  Double_t dz = sqrt(dZpc*dZpc + dRpc*dRpc + dZsi*dZsi*0.1*0.1 +  dRsi*dRsi*0.1*0.1);
  return dz;

}

Double_t Reconstruct::WIntP(Track Tr, vector<Int_t> IsProton) {

  TrackEvent ejectile1 = Tr.TrEvent[IsProton[0]];
  TrackEvent ejectile2 = Tr.TrEvent[IsProton[1]];

  Double_t IntPoint_e1 = ejectile1.IntPoint;
  Double_t IntPoint_e2 = ejectile2.IntPoint;

  Double_t dz1 = DeltaZ(ejectile1.PCZ, ejectile1.SiZ, ejectile1.PCR, ejectile1.SiR);
  Double_t dz2 = DeltaZ(ejectile2.PCZ, ejectile2.SiZ, ejectile2.PCR, ejectile2.SiR);

  Double_t IntPoint = (IntPoint_e1/dz1 + IntPoint_e2/dz2)/(1.0/dz1+1.0/dz2);
  return IntPoint;

}

//Used to ref track, but never updated track... Was it supposed to?
void Reconstruct::ThetaPath(Double_t &Theta, Double_t &Path, Double_t Zsi, Double_t Rsi, Track Tr,
                            vector<Int_t> IsProton) {

  Double_t IntPoint = WIntP(Tr, IsProton);

  if ((IntPoint-Zsi) > 0.0) {
    Theta = atan(Rsi/(IntPoint - Zsi));
    Path = Rsi/sin(Theta);
  } else if ((IntPoint-Zsi) < 0.0) {
    Theta = TMath::Pi() + atan(Rsi/(IntPoint - Zsi));
    Path = Rsi/sin(Theta);
  } else {
    Theta = TMath::Pi()/2.0;
    Path = Rsi;
  }

}

void Reconstruct::Reconstruct20Ne(Track Tr, vector<Int_t> IsProton, Double_t Qvalue,  
                                         bool IntPnt_flag, RecoilEvent &recoil) {

  TLorentzVector Beam_LV(0.,0.,0.,0.), Target_LV(0.,0.,0.,0.), Parent_LV(0.,0.,0.,0.),
                 Recoil_LV(0.,0.,0.,0.), eject1_LV(0.,0.,0.,0.), eject2_LV(0.,0.,0.,0.),
                 Sum_eject_LV(0.,0.,0.,0.), Recoil_LV_beam_rcnstrct(0.,0.,0.,0.),
                 Beam_LV_rcnstrct(0.,0.,0.,0.), Parent_LV_beam_rcnstrct(0.,0.,0.,0.);

  TVector3 eject1_V(0.,0.,0.), eject2_V(0.,0.,0.), Sum_eject_V(0.,0.,0.);

  Target_LV.SetPxPyPzE(0.0,0.0,0.0, m_target);
  Double_t IntPoint = 0.0, BeamE_tot = 0.0, BeamE_tot_rcnstrct = 0.0,BeamWA = 0.0, 
           E_eject1_rxn = 0.0, E_eject2_rxn = 0.0, E_eject_tot = 0.0;
  
  if (IntPnt_flag) {
    IntPoint = WIntP(Tr, IsProton);
  }

  for (unsigned int i = 0; i<IsProton.size(); i++) {
    TrackEvent ejectile = Tr.TrEvent[IsProton[i]];
    Double_t E_eject_si = ejectile.SiEnergy;
    Double_t theta_eject, d_eject;

    if (IntPnt_flag) {
      ThetaPath(theta_eject, d_eject, ejectile.SiZ, ejectile.SiR, Tr, IsProton);
    } else {
      theta_eject = ejectile.Theta;
      d_eject = ejectile.PathLength;
    }
    
    if (i == 0) {
      E_eject1_rxn = ELoss_light->GetLookupEnergy(E_eject_si, -d_eject);
      E_eject_tot = E_eject1_rxn + m_ejectile;
    } else if (i == 1) {
      E_eject2_rxn = ELoss_light->GetLookupEnergy(E_eject_si, -d_eject);
      E_eject_tot = E_eject2_rxn + m_ejectile;
    }
    Double_t P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile*m_ejectile);

    Double_t p_x = P_eject*cos(ejectile.SiPhi)*sin(theta_eject);
    Double_t p_y = P_eject*sin(ejectile.SiPhi)*sin(theta_eject);
    Double_t p_z = P_eject*cos(theta_eject);

    if (i == 0) {
      eject1_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
      eject1_V.SetXYZ(p_x, p_y, p_z);
    } else if (i == 1) {
      eject2_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
      eject2_V.SetXYZ(p_x, p_y, p_z);
    }
  }

  Sum_eject_LV = eject1_LV + eject2_LV;
  Sum_eject_V = eject1_V + eject2_V;
 
  Double_t a = (m_beam/m_recoil-1.0);
  //Check that theta used bellow is correct theta
  Double_t b = -(sqrt(2.0*m_beam)/m_recoil)*Sum_eject_V.Mag()*cos(Sum_eject_V.Theta());
  Double_t c = E_eject1_rxn + E_eject2_rxn + Sum_eject_V.Mag()*Sum_eject_V.Mag()/(2.0*m_ejectile) 
               - Qvalue;

  Double_t s1 = (-b + sqrt(b*b-4*a*c))/(2.0*a);
  Double_t s2 = (-b - sqrt(b*b-4*a*c))/(2.0*a);

  if(s1>s2) {
    cout <<"BB s1: "<<s1<<" s2: "<<s2<<" Qvalue: "<<Qvalue<<endl;
    swap(s1, s2);
  }
  Double_t Beam_KE = s2*s2;
  if (IntPnt_flag) {
  }
  else recoil.BeamEnergy_20Ne = Beam_KE;
   
  if (IntPnt_flag) {
    BeamWA = ELoss_beam->GetLookupEnergy(72.34, (55.0545-IntPoint));
    BeamE_tot = BeamWA + m_beam;
    recoil.BeamWA_Qvalue_20Ne = Beam_KE;
    BeamE_tot_rcnstrct = Beam_KE;
  } else { // in orig, was hard coded proton 0; necessary?
    BeamE_tot = Tr.TrEvent[IsProton[0]].BeamEnergy + m_beam; 
    recoil.BeamEnergy_20Ne = Beam_KE;
  }

  Double_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - m_beam*m_beam);
  Double_t BeamP_z_rcnstrct = sqrt(BeamE_tot_rcnstrct*BeamE_tot_rcnstrct - m_beam*m_beam);
  Beam_LV.SetPxPyPzE(0.0,0.0, BeamP_z, BeamE_tot);
  Beam_LV_rcnstrct.SetPxPyPzE(0.0,0.0, BeamP_z_rcnstrct, BeamE_tot_rcnstrct);
  Parent_LV = Target_LV + Beam_LV;
  Parent_LV_beam_rcnstrct = Target_LV + Beam_LV_rcnstrct;
  Recoil_LV = Parent_LV - Sum_eject_LV;
  Recoil_LV_beam_rcnstrct = Parent_LV_beam_rcnstrct - Sum_eject_LV;

  ////So many repeats from Beam calc and 20Ne stored in recoil, are they all necessary? 
  if (IntPnt_flag) {
    recoil.IntP_20Ne = IntPoint;
    recoil.IntP_Qvalue_20Ne = IntPoint;
    recoil.BeamWA_20Ne = BeamWA;
    recoil.Theta_Qvalue_20Ne = Recoil_LV_beam_rcnstrct.Theta()*180.0/TMath::Pi();
    recoil.Phi_Qvalue_20Ne = Recoil_LV_beam_rcnstrct.Phi()*180.0/TMath::Pi();
    recoil.KE_Qvalue_20Ne = Recoil_LV_beam_rcnstrct.E()-Recoil_LV_beam_rcnstrct.M();
    recoil.Ex_Qvalue_20Ne = Recoil_LV_beam_rcnstrct.M() - m_recoil;
  } else {
    recoil.BeamWA_20Ne = BeamE_tot - m_beam;
  }

  recoil.theta_p1_20Ne = eject1_LV.Theta()*180.0/TMath::Pi();
  recoil.theta_p2_20Ne = eject2_LV.Theta()*180.0/TMath::Pi();
  recoil.Theta_20Ne = Recoil_LV.Theta()*180.0/TMath::Pi();
  recoil.Phi_20Ne = Recoil_LV.Phi()*180.0/TMath::Pi();
  recoil.KE_20Ne = Recoil_LV.E() - Recoil_LV.M();
  recoil.Ex_20Ne = Recoil_LV.M() - m_recoil;

  recoil.Delta_Phi = abs(eject1_LV.Phi()-eject2_LV.Phi());
  if(recoil.Delta_Phi > TMath::Pi() && recoil.Delta_Phi <= 2*TMath::Pi()) {
    recoil.Delta_Phi = 2*TMath::Pi() - recoil.Delta_Phi;
  } 
  recoil.Delta_Theta = abs(eject1_LV.Theta() - eject2_LV.Theta())*180.0/TMath::Pi();

  Double_t m_21Na = 19553.56837;
  TLorentzVector Na21_p1_LV = Recoil_LV + eject1_LV, Na21_p2_LV = Recoil_LV + eject2_LV;
  Double_t Ex_p1 = Na21_p1_LV.M() - m_21Na, Ex_p2 = Na21_p2_LV.M() - m_21Na;

  if (Ex_p1 > Ex_p2){
    recoil.Ex_rec_small = Ex_p2;
    recoil.Ex_rec_qvalue_small = Ex_p2;
    recoil.Ex_rec_large = Ex_p1;
    recoil.Ex_rec_qvalue_large = Ex_p1;
  } else {
    recoil.Ex_rec_small = Ex_p1;
    recoil.Ex_rec_qvalue_small = Ex_p1;
    recoil.Ex_rec_large = Ex_p2;
    recoil.Ex_rec_qvalue_large = Ex_p2;
  }

}

void Reconstruct::ReconstructHeavy_Qvalue(Double_t Qvalue, Track Tr, Int_t proton, 
                                          RecoilEvent &recoil) {

  TrackEvent ejectile = Tr.TrEvent[proton];

  Double_t E_eject_si = ejectile.SiEnergy;
  Double_t d_eject = ejectile.PathLength;
  Double_t E_eject_rxn = ELoss_light->GetLookupEnergy(E_eject_si, -d_eject);
  Double_t E_eject_tot = E_eject_rxn + m_ejectile;
  Double_t P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile*m_ejectile);
  Double_t px = P_eject*sin(ejectile.Theta)*cos(ejectile.SiPhi),
           py = P_eject*sin(ejectile.Theta)*sin(ejectile.SiPhi),
           pz = P_eject*cos(ejectile.Theta);

  TLorentzVector Eject_LV(0.,0.,0.,0.), Recoil_LV(0.,0.,0.,0.), Beam_LV(0.,0.,0.,0.),
                 Target_LV(0.,0.,0.,0.), Parent_LV(0.,0.,0.,0.);

  Double_t a = -(1.0-m_beam/m_recoil);
  Double_t b = -2.0*sqrt(m_beam*m_ejectile*E_eject_rxn)*cos(ejectile.Theta)/m_recoil;
  Double_t c = -Qvalue + E_eject_rxn*(1+m_ejectile/m_recoil);
  Double_t s1 = (-b+sqrt(b*b-4*a*c))/(2.0*a);
  Double_t s2 = (-b-sqrt(b*b-4*a*c))/(2.0*a);

  if(s1>s2) swap(s1,s2);
  Double_t Beam_KE = s2*s2;
  //Double_t E_Recoil = (m_target+m_beam+Beam_KE)-m_ejectile-E_eject_rxn-m_recoil; Unused?
  
  Eject_LV.SetPxPyPzE(px, py, pz, E_eject_tot);
  Double_t E_beam_tot = Beam_KE + m_beam;
  Double_t pz_beam = sqrt(E_beam_tot*E_beam_tot-m_beam*m_beam);
  Beam_LV.SetPxPyPzE(0.,0.,pz_beam, E_beam_tot);
  Target_LV.SetPxPyPzE(0.,0.,0.,m_target);
  Parent_LV = Beam_LV+Target_LV;
  Recoil_LV = Parent_LV-Eject_LV;

  recoil.Beam_Qv_21Na = Beam_KE;
  recoil.Ex_21Na = Recoil_LV.M()-m_recoil;
  
}





