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
#include "TMath.h"

//Constructor
Reconstruct::Reconstruct() {
  m_beam = 0.0;
  m_target = 0.0;
  m_ejectile = 0.0;
  m_recoil = 0.0;
  ELoss_eject = NULL;
  ELoss_beam = NULL;
}

//Constructor overload
Reconstruct::Reconstruct(Double_t beam,Double_t target,Double_t ejectile,Double_t recoil) {
  m_beam = beam;
  m_target = target;
  m_ejectile = ejectile;
  m_ejectile2 = ejectile;
  m_recoil = recoil;
  ELoss_eject = NULL;
  ELoss_eject2 = NULL;
  ELoss_beam = NULL;
}
Reconstruct::Reconstruct(Double_t beam,Double_t target,Double_t eject,Double_t eject2, 
                         Double_t recoil) {
  m_beam = beam;
  m_target = target;
  m_ejectile = eject;
  m_ejectile2 = eject2;
  m_recoil = recoil;
  ELoss_eject = NULL;
  ELoss_eject2 = NULL;
  ELoss_beam = NULL;
}


//Destructor
Reconstruct::~Reconstruct() {
  delete ELoss_eject;
  delete ELoss_eject2;
  delete ELoss_beam;
}

//SetELossFile() Takes in names for eloss files and and sets class vars
void Reconstruct::SetELossFile(const char* name1, const char* name2, const char* name3) {
  ELoss_eject = new LookUp(name1, m_ejectile);
  ELoss_eject2 = new LookUp(name3, m_ejectile2);
  ELoss_beam = new LookUp(name2, m_beam);
}

/*Elastic()
 *Reconstructs an elastic scattering event. 
 *Returns the energy of the ejectile at the reaction point, the energy of the incoming 
 *particle, and the angle of the ejectile path.
 */
Double_t Reconstruct::CalcElastic(Double_t &EnergyProj, Double_t &EjectTheta, Track Tr, Int_t c) {

  Double_t E_eject_si = Tr.TrEvent[c].SiEnergy;
  Double_t d_eject = Tr.TrEvent[c].PathLength;
  Double_t theta = Tr.TrEvent[c].Theta;

  //Get the inital energy of the ejectile
  Double_t E_eject_rxn = ELoss_eject->GetInitialEnergy(E_eject_si, d_eject, 0.1);

  //Calculate the energy of the proj initially, and the angle of the scattered proj
  EnergyProj = ((m_beam+m_ejectile)*(m_beam+m_ejectile)*E_eject_rxn)/
                   ((4.0*theta*theta*m_beam*m_ejectile));
  EjectTheta = asin((sin(theta)*sqrt(m_ejectile/m_beam))/(sqrt(EnergyProj/E_eject_rxn-1.0)));
  return E_eject_rxn;
}

/*CalcRecoil()
 *Reconstructs a recoil nucleus using the tracking information of a single ejectile, and assumed
 *knowledge of the beam and target. Should be compared to CaclRecoil_from_Qvalue, which does not
 *rely as much on the tracking to calculate the same value.
 */                            
void Reconstruct::CalcRecoil(Track Tr, Int_t eject, RecoilEvent &recoil) {

  TrackEvent ejectile = Tr.TrEvent[eject];

  TLorentzVector Eject_LV(0., 0., 0., 0.), Beam_LV(0., 0., 0., 0.), Target_LV(0., 0., 0., 0.),
                 Parent_LV(0., 0., 0., 0.), Recoil_LV(0., 0., 0., 0.);

  Float_t E_eject_si = ejectile.SiEnergy;
  Float_t d_eject = ejectile.PathLength;
  Float_t theta_eject = ejectile.Theta;
  Float_t phi_eject = ejectile.SiPhi;
  Float_t E_eject_rxn = ELoss_eject->GetLookupEnergy(E_eject_si, -d_eject);
  Float_t E_eject_tot = E_eject_rxn + m_ejectile;
  Float_t P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile*m_ejectile);

  Eject_LV.SetPxPyPzE(P_eject*sin(theta_eject)*cos(phi_eject), 
                      P_eject*sin(theta_eject)*sin(phi_eject),
                      P_eject*cos(theta_eject),
                      E_eject_tot);
   
  Float_t BeamE_tot = ejectile.BeamEnergy + m_beam;
  Float_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - m_beam*m_beam);
  Float_t TargetE = m_target;

  Beam_LV.SetPxPyPzE(0.0, 0.0, BeamP_z, BeamE_tot);
  Target_LV.SetPxPyPzE(0.0,0.0,0.0, TargetE);
  
  Parent_LV = Beam_LV + Target_LV;

  Recoil_LV = Parent_LV - Eject_LV;

  recoil.IntPoint = ejectile.IntPoint;
  recoil.BeamKE = ejectile.BeamEnergy;
  recoil.BeamKE_eject = ejectile.BeamEnergy;
  recoil.SiEnergy_tot = ejectile.SiEnergy;
  recoil.PCEnergy_tot = ejectile.PCEnergy*ejectile.Theta; //Why Theta here???
  recoil.Energy_eject_tot = E_eject_rxn;
  recoil.Theta = Recoil_LV.Theta()*180.0/TMath::Pi();
  recoil.Phi = Recoil_LV.Phi()*180.0/TMath::Pi();
  recoil.KE_recoil = Recoil_LV.E()-Recoil_LV.M();
  recoil.Ex_recoil = Recoil_LV.M() - m_recoil;
}

/*DeltaZ()
 *Caclulates uncertainties for interaction point average. In principle not really necessary for 
 *7Be+d since alphas & 3He are pretty easy, but good to keep just in case
 */
Double_t Reconstruct::DeltaZ(Double_t Zpc, Double_t Zsi, Double_t Rpc, Double_t Rsi) {
  
  Double_t dZsi = -(Rpc/(Rsi-Rpc));
  Double_t dZpc = 1-dZsi; 
  Double_t dRsi = (Rpc*(Zsi-Zpc))/((Rsi-Rpc)*(Rsi-Rpc));
  Double_t dRpc = -(Rsi*(Zsi-Zpc))/((Rsi-Rpc)*(Rsi-Rpc));

  Double_t dz = sqrt(dZpc*dZpc + dRpc*dRpc + dZsi*dZsi*0.1*0.1 +  dRsi*dRsi*0.1*0.1);
  return dz;

}

/*WIntP()
 *Caclulates weighted average interaction point
 *Weights are from the uncertainty calculation DeltaZ
 */
Double_t Reconstruct::WIntP(Track Tr, vector<Int_t> IsEject) {
  vector<Double_t> dz;
  vector<Double_t> intp;
  Double_t denom = 0;
  Double_t IntPoint = 0;
  for(unsigned int i = 0; i<IsEject.size(); i++) {
    TrackEvent ejectile = Tr.TrEvent[IsEject[i]];
    intp.push_back(ejectile.IntPoint);
    Double_t dz_i = DeltaZ(ejectile.PCZ, ejectile.SiZ, ejectile.PCR, ejectile.SiR);
    dz.push_back(dz_i);
    denom += 1.0/dz_i;
  }
  for(unsigned int i = 0; i<IsEject.size(); i++) {
    IntPoint += (intp[i]/dz[i])/denom;
  }
  return IntPoint;

}

/*ThetaPath()
 *Calculates angle Theta and the path length
 */
void Reconstruct::ThetaPath(Double_t &Theta, Double_t &Path, Double_t Zsi, Double_t Rsi, Track Tr,
                            vector<Int_t> IsEject) {

  Double_t IntPoint = WIntP(Tr, IsEject);

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
  float test = Path;
  if(isnan(test)) {
    cout<<"Gotcha"<<endl;
    Theta = TMath::Pi()/2.0;
    Path = Rsi;
  }

}

/*CalcRecoil_MultiParticle()
 *This is the big boy... Reconstructs an event from 3 outgoing particles, in this specific ver,
 *1 p and 2 alpha. Then calculates the various properties of the residual nucleus from the 
 *detected particles. 
 */
void Reconstruct::CalcRecoil_MultiParticle(Track Tr, vector<Int_t> IsEject, 
                                           RecoilEvent &recoil_qval, RecoilEvent &recoil) 
{
  if (IsEject.size() == 3) {
    TLorentzVector Beam_LV(0.,0.,0.,0.), Target_LV(0.,0.,0.,0.), Parent_LV(0.,0.,0.,0.),
                   Recoil_LV(0.,0.,0.,0.), eject1_LV(0.,0.,0.,0.), eject2_LV(0.,0.,0.,0.),
                   eject3_LV(0.,0.,0.,0.), Sum_eject_LV(0.,0.,0.,0.),
                   Recoil_LV_beam_qval(0.,0.,0.,0.), Beam_LV_qval(0.,0.,0.,0.), 
                   Parent_LV_beam_qval(0.,0.,0.,0.);
 
    TVector3 eject1_V(0.,0.,0.), eject2_V(0.,0.,0.), eject3_V(0.,0.,0.), Sum_eject_V(0.,0.,0.);
 
    Target_LV.SetPxPyPzE(0.0,0.0,0.0, m_target);
    Double_t BeamE_tot = 0.0, BeamE_tot_qval = 0.0, 
             E_eject1_rxn = 0.0, E_eject2_rxn = 0.0, E_eject3_rxn = 0.0, E_eject_tot = 0.0;
    
    Double_t BeamKE_avg = 0.0;
    for (int i = 0; i<3; i++) {
      TrackEvent ejectile = Tr.TrEvent[IsEject[i]];
      Double_t E_eject_si = ejectile.SiEnergy;
      Double_t theta_eject, d_eject;
 
      theta_eject = ejectile.Theta;
      d_eject = ejectile.PathLength;
      
      Double_t P_eject = 0.0;
      if (i == 0) {
        E_eject1_rxn = ELoss_eject->GetLookupEnergy(E_eject_si, -d_eject);
        E_eject_tot = E_eject1_rxn + m_ejectile;
        P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile*m_ejectile);
        Double_t p_x = P_eject*cos(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_y = P_eject*sin(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_z = P_eject*cos(theta_eject);
        eject1_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
        eject1_V.SetXYZ(p_x, p_y, p_z);
      } else if (i == 1) {
        E_eject2_rxn = ELoss_eject2->GetLookupEnergy(E_eject_si, -d_eject);
        E_eject_tot = E_eject2_rxn + m_ejectile2;
        P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile2*m_ejectile2);
        Double_t p_x = P_eject*cos(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_y = P_eject*sin(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_z = P_eject*cos(theta_eject);
        eject2_LV.SetPxPyPzE(p_x,p_y,p_z, E_eject_tot);
        eject2_V.SetXYZ(p_x, p_y, p_z);
        BeamKE_avg += ejectile.BeamEnergy;
      } else if (i == 2) {
        E_eject3_rxn = ELoss_eject2->GetLookupEnergy(E_eject_si, -d_eject);
        E_eject_tot = E_eject3_rxn+m_ejectile2;
        P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile2*m_ejectile2);
        Double_t p_x = P_eject*cos(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_y = P_eject*sin(ejectile.SiPhi)*sin(theta_eject);
        Double_t p_z = P_eject*cos(theta_eject);
        eject3_LV.SetPxPyPzE(p_x, p_y, p_z, E_eject_tot);
        eject3_V.SetXYZ(p_x, p_y, p_z);
        BeamKE_avg += ejectile.BeamEnergy;
      }
    }
 
    Sum_eject_LV = eject1_LV + eject2_LV + eject3_LV;
    Sum_eject_V = eject1_V + eject2_V + eject3_V;
  
    Double_t Beam_KE = Sum_eject_LV.E()-m_target-m_beam;//Qvalue included in this calc implicitly
 
    BeamKE_avg = BeamKE_avg/2.0;
    BeamE_tot = BeamKE_avg + m_beam; 
    recoil.BeamKE = BeamKE_avg;
    Double_t BeamP_z = sqrt(BeamE_tot*BeamE_tot - m_beam*m_beam);
    BeamE_tot_qval = Beam_KE + m_beam;
 
    Double_t BeamP_z_qval = sqrt(BeamE_tot_qval*BeamE_tot_qval - m_beam*m_beam);
    Beam_LV_qval.SetPxPyPzE(0.0,0.0, BeamP_z_qval, BeamE_tot_qval);
    Parent_LV_beam_qval = Target_LV + Beam_LV_qval;
    Recoil_LV_beam_qval = Parent_LV_beam_qval - eject1_LV;
 
    Beam_LV.SetPxPyPzE(0.0,0.0, BeamP_z, BeamE_tot);
    Parent_LV = Target_LV + Beam_LV;
    Recoil_LV = Parent_LV - Sum_eject_LV;
    TLorentzVector Recoil_LV_2a = eject2_LV+eject3_LV;
 
    recoil.BeamPz = BeamP_z;
    recoil_qval.recoil_mass_sq = Recoil_LV_2a.M()*Recoil_LV_2a.M()/1.0e6;//convert to GeV^2
    recoil.recoil_mass_sq = Recoil_LV_2a.M()*Recoil_LV_2a.M()/1.0e6;
    recoil.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil.a2_ip = Tr.TrEvent[IsEject[2]].IntPoint;
    recoil.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    recoil_qval.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil_qval.a2_ip = Tr.TrEvent[IsEject[2]].IntPoint;
    recoil_qval.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    recoil_qval.BeamPz = BeamP_z;
    recoil_qval.BeamKE = Beam_KE;
    recoil_qval.BeamKE_eject = BeamKE_avg;
    recoil_qval.Theta = Recoil_LV_beam_qval.Theta()*180.0/TMath::Pi();
    recoil_qval.Phi = Recoil_LV_beam_qval.Phi()*180.0/TMath::Pi();
    recoil_qval.KE_recoil = Recoil_LV_beam_qval.E()-Recoil_LV_beam_qval.M();
    recoil_qval.Ex_recoil = Recoil_LV_beam_qval.M() - m_recoil;
    recoil_qval.Ex_recoil_2a = Recoil_LV_2a.M() - m_recoil;
    recoil_qval.theta_eject1 = eject2_LV.Theta()*180.0/TMath::Pi();
    recoil_qval.theta_eject2 = eject3_LV.Theta()*180.0/TMath::Pi();
    recoil.theta_eject1 = eject1_LV.Theta()*180.0/TMath::Pi();
    recoil.theta_eject2 = eject2_LV.Theta()*180.0/TMath::Pi();
    recoil.Theta = Recoil_LV.Theta()*180.0/TMath::Pi();
    recoil.Phi = Recoil_LV.Phi()*180.0/TMath::Pi();
    recoil.KE_recoil = Recoil_LV.E() - Recoil_LV.M();
    recoil.Ex_recoil = Recoil_LV.M() - m_recoil;
    recoil.Ex_recoil_2a = Recoil_LV_2a.M()-m_recoil;
    recoil.Ex_recoil_2a = Recoil_LV_2a.M() - m_recoil;
 
    recoil.Delta_Phi = abs(eject1_LV.Phi()-eject2_LV.Phi());
    if(recoil.Delta_Phi > TMath::Pi() && recoil.Delta_Phi <= 2*TMath::Pi()) {
      recoil.Delta_Phi = 2*TMath::Pi() - recoil.Delta_Phi;
    } 
    recoil.Delta_Theta = abs(eject1_LV.Theta() - eject2_LV.Theta())*180.0/TMath::Pi();

  } else {
    cout<<"Error CalcRecoil_MultiParticle! Has to have 3 ejectiles"<<endl;
  }
}

/*CalcRecoil_from_Qvalue()
 *Reconstructs a recoil nucleus from a single ejectile assuming knowledge of the beam,
 *reaction Qvalue, and target. This is mainly meant to be a comparison to CalcRecoil() to 
 *ensure that the tracking isn't super screwed up. 
 */
void Reconstruct::CalcRecoil_from_Qvalue(Double_t Qvalue, Track Tr, Int_t proton, 
                                         RecoilEvent &recoil_qval) {

  TrackEvent ejectile = Tr.TrEvent[proton];

  Double_t E_eject_si = ejectile.SiEnergy;
  Double_t d_eject = ejectile.PathLength;
  Double_t E_eject_rxn = ELoss_eject->GetLookupEnergy(E_eject_si, -d_eject);
  Double_t E_eject_tot = E_eject_rxn + m_ejectile;
  Double_t P_eject = sqrt(E_eject_tot*E_eject_tot - m_ejectile*m_ejectile);
  Double_t px = P_eject*sin(ejectile.Theta)*cos(ejectile.SiPhi),
           py = P_eject*sin(ejectile.Theta)*sin(ejectile.SiPhi),
           pz = P_eject*cos(ejectile.Theta);

  TLorentzVector Eject_LV(0.,0.,0.,0.), Recoil_LV(0.,0.,0.,0.), Beam_LV(0.,0.,0.,0.),
                 Target_LV(0.,0.,0.,0.), Parent_LV(0.,0.,0.,0.);
  Double_t a = (m_beam-m_recoil);
  Double_t b = (m_recoil+m_ejectile);
  Double_t c = -Qvalue*m_recoil;
  Double_t d = 2.0*sqrt(m_beam*m_ejectile);
  Double_t cos2 = cos(ejectile.Theta)*cos(ejectile.Theta);
  Double_t a1 = a*a;
  Double_t a2 = (2.0*a*(b*E_eject_rxn+c)-d*d*E_eject_rxn*cos2);
  Double_t a3 = (b*E_eject_rxn+c)*(b*E_eject_rxn+c);
  Double_t s1 = (-a2+sqrt(a2*a2-4.0*a1*a3))/(2.0*a1);
  Double_t s2 = (-a2-sqrt(a2*a2-4.0*a1*a3))/(2.0*a1);

  if(s1<s2) swap(s1,s2);
  Double_t Beam_KE = s2;
  if(TMath::IsNaN(s2)) {
    return;
  }
  
  Eject_LV.SetPxPyPzE(px, py, pz, E_eject_tot);
  Double_t E_beam_tot = Beam_KE + m_beam;
  Double_t pz_beam = sqrt(E_beam_tot*E_beam_tot-m_beam*m_beam);
  Beam_LV.SetPxPyPzE(0.,0.,pz_beam, E_beam_tot);
  Target_LV.SetPxPyPzE(0.,0.,0.,m_target);
  Parent_LV = Beam_LV+Target_LV;
  Recoil_LV = Parent_LV-Eject_LV;

  recoil_qval.BeamKE = Beam_KE;
  recoil_qval.BeamKE_eject = ejectile.BeamEnergy;
  recoil_qval.Ex_recoil = Recoil_LV.M()-m_recoil;
}

/*CalcE_ex()
 *Function for calculating excitation energy of a recoil 
 *Assumption is that the ejectile AND a double hit are tracked
 *Current implementation is for 7Be(d,p)2alpha where the two alphas are so close together
 *they are counted as one track; energy at si is divided equally between them
 */
void Reconstruct::CalcE_ex(Track Tr, vector<int> IsEject, RecoilEvent &recoil, 
                           RecoilEvent &recoil_qval) {
  if(IsEject.size() == 2) {
    Double_t BeamKE_avg = 0.0;
    TLorentzVector eject1_LV(0.,0.,0.,0.), eject2_LV(0.,0.,0.,0.); 
    for(int i =0; i<2; i++) {
      TrackEvent ejectile = Tr.TrEvent[IsEject[i]];
      Double_t pl = ejectile.PathLength;
      Double_t sie = ejectile.SiEnergy;
      Double_t theta = ejectile.Theta;
      Double_t phi = ejectile.SiPhi;
      Double_t be = ejectile.BeamEnergy;
      if(i == 0) {
        Double_t E_eject_rxn = ELoss_eject->GetLookupEnergy(sie, -pl);
        Double_t E_eject_tot = E_eject_rxn + m_ejectile;
        Double_t P_eject_tot = sqrt(E_eject_tot*E_eject_tot-m_ejectile*m_ejectile);
        eject1_LV.SetPxPyPzE(P_eject_tot*sin(theta)*cos(phi),
                             P_eject_tot*sin(theta)*sin(phi),
                             P_eject_tot*cos(theta),
                             E_eject_tot);
        BeamKE_avg += be;
      } else {
        Double_t E_eject_rxn = ELoss_eject2->GetLookupEnergy(sie/2.0, -pl);
        Double_t E_eject_tot = E_eject_rxn + m_ejectile2;
        Double_t P_eject_tot = sqrt(E_eject_tot*E_eject_tot-m_ejectile2*m_ejectile2);
        eject2_LV.SetPxPyPzE(P_eject_tot*sin(theta)*cos(phi),
                             P_eject_tot*sin(theta)*sin(phi),
                             P_eject_tot*cos(theta),
                             E_eject_tot);
        BeamKE_avg += be;
      }
    }
    Double_t BeamKE = eject1_LV.E()+2.0*eject2_LV.E()-m_target-m_beam;
    BeamKE_avg = BeamKE_avg/2.0;
    Double_t E_beam_tot_qval = BeamKE+m_beam;
    Double_t E_beam_tot = BeamKE_avg+m_beam;
    Double_t Pz_beam_qval = sqrt(E_beam_tot_qval*E_beam_tot_qval-m_beam*m_beam);
    Double_t Pz_beam = sqrt(E_beam_tot*E_beam_tot-m_beam*m_beam);
    TLorentzVector recoil_2a_LV = eject2_LV + eject2_LV;
    TLorentzVector beam_LV(0.,0.,0.,0.), beam_qval_LV(0.,0.,0.,0.), targ_LV(0.,0.,0.,0.);
    //cout<<"M: "<<recoil_LV.M()<<" eject M: "<<eject2_LV.M()<<endl;
    beam_LV.SetPxPyPzE(0.0,0.0,Pz_beam,E_beam_tot);
    beam_qval_LV.SetPxPyPzE(0.0,0.0,Pz_beam_qval,E_beam_tot_qval);
    targ_LV.SetPxPyPzE(0.0,0.0,0.0,m_target);
    TLorentzVector parent_LV = beam_LV+targ_LV;
    TLorentzVector parent_qval_LV = beam_qval_LV+targ_LV;
    TLorentzVector recoil_LV = parent_LV-eject1_LV;
    TLorentzVector recoil_qval_LV = parent_qval_LV-eject1_LV;

    recoil.BeamKE = BeamKE_avg;
    recoil.BeamPz = Pz_beam;
    recoil.Ex_recoil = recoil_LV.M() - m_recoil;
    recoil_qval.Ex_recoil = recoil_qval_LV.M()-m_recoil;
    //for this ex, need to add on qval of 8Be breakup. Without individ tracks, this info is lost
    recoil.Ex_recoil_2a = recoil_2a_LV.M() - m_recoil+0.091;
    recoil.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    recoil_qval.BeamKE = BeamKE;
    recoil_qval.BeamPz = Pz_beam_qval;
    recoil_qval.Ex_recoil_2a = recoil_2a_LV.M() - m_recoil + 0.091;
    recoil_qval.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil_qval.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    return;
  } else {
    cout<<"Incorrect Reconstruction in CalcE_ex! Should only have two ejectiles"<<endl;
    return;
  }
}

/*CalcRecoil_InvMassSelect()
 *Reconstruction where one of the identical ejectiles comes from the inital reaction while the 
 *other comes from the breakup. Choses the ejectile which gives the lowest ex energy of the final
 *recoil. 
 */
void Reconstruct::CalcRecoil_InvMassSelect(Track Tr, vector<int> IsEject, RecoilEvent &recoil,
                                           RecoilEvent &recoil_qval) {
  if(IsEject.size() == 3) {
    TLorentzVector eject1_LV(0,0,0,0), eject2_LV(0,0,0,0), eject3_LV(0,0,0,0);
    Double_t BeamKE_avg = 0.0;
    for (int i=0; i<3; i++) {
      TrackEvent ejectile = Tr.TrEvent[IsEject[i]];
      Double_t pl = ejectile.PathLength;
      Double_t sie = ejectile.SiEnergy;
      Double_t theta = ejectile.Theta;
      Double_t phi = ejectile.SiPhi;
      Double_t be = ejectile.BeamEnergy;
      if(i==0) {
        Double_t E_eject_rxn = ELoss_eject->GetLookupEnergy(sie, -pl);
        Double_t E_eject_tot = E_eject_rxn+m_ejectile;
        Double_t P_eject_tot = sqrt(E_eject_tot*E_eject_tot-m_ejectile*m_ejectile);
        eject1_LV.SetPxPyPzE(P_eject_tot*sin(theta)*cos(phi),
                             P_eject_tot*sin(theta)*sin(phi),
                             P_eject_tot*cos(theta),
                             E_eject_tot);
        BeamKE_avg += be;
      }else if(i == 1) {
        Double_t E_eject_rxn = ELoss_eject2->GetLookupEnergy(sie, -pl);
        Double_t E_eject_tot = E_eject_rxn+m_ejectile2;
        Double_t P_eject_tot = sqrt(E_eject_tot*E_eject_tot-m_ejectile2*m_ejectile2);
        eject2_LV.SetPxPyPzE(P_eject_tot*sin(theta)*cos(phi),
                             P_eject_tot*sin(theta)*sin(phi),
                             P_eject_tot*cos(theta),
                             E_eject_tot);
        BeamKE_avg += be;
      }else {
        Double_t E_eject_rxn = ELoss_eject2->GetLookupEnergy(sie, -pl);
        Double_t E_eject_tot = E_eject_rxn+m_ejectile2;
        Double_t P_eject_tot = sqrt(E_eject_tot*E_eject_tot-m_ejectile2*m_ejectile2);
        eject3_LV.SetPxPyPzE(P_eject_tot*sin(theta)*cos(phi),
                             P_eject_tot*sin(theta)*sin(phi),
                             P_eject_tot*cos(theta),
                             E_eject_tot);
        BeamKE_avg += be;
      }
    }
    TLorentzVector Sum_eject_LV = eject1_LV+eject2_LV+eject3_LV;
    Double_t BeamKE = Sum_eject_LV.E()-m_target-m_beam;
    BeamKE_avg = BeamKE_avg/3.0;
    Double_t E_beam_tot = BeamKE_avg +m_beam;
    Double_t E_beam_tot_qval = BeamKE +m_beam;
    Double_t Pz_beam = sqrt(E_beam_tot*E_beam_tot - m_beam*m_beam);
    Double_t Pz_beam_qval = sqrt(E_beam_tot_qval*E_beam_tot_qval - m_beam*m_beam);
    TLorentzVector beam_LV(0,0,0,0), beam_qval_LV(0,0,0,0), targ_LV(0,0,0,0);
    beam_LV.SetPxPyPzE(0,0,Pz_beam,E_beam_tot);
    beam_qval_LV.SetPxPyPzE(0,0,Pz_beam_qval,E_beam_tot_qval);
    targ_LV.SetPxPyPzE(0,0,0,m_target);
    TLorentzVector parent_LV = targ_LV+beam_LV;
    TLorentzVector parent_qval_LV = targ_LV+beam_qval_LV;

    TLorentzVector recoil_LV;
    TLorentzVector rec1_LV = eject1_LV+eject2_LV;
    TLorentzVector rec2_LV = eject1_LV+eject3_LV;
    Double_t ex1 = rec1_LV.M()-m_recoil;
    Double_t ex2 = rec2_LV.M()-m_recoil;
    if(ex1<ex2) {
      recoil_LV = rec1_LV;
      recoil.Ex_recoil = ex1;
      recoil_qval.Ex_recoil = ex1;
    } else {
      recoil_LV = rec2_LV;
      recoil.Ex_recoil = ex2;
      recoil_qval.Ex_recoil = ex2;
    }

    recoil.recoil_mass_sq = recoil_LV.M()*recoil_LV.M()/1.0e6;
    recoil_qval.recoil_mass_sq = recoil_LV.M()*recoil_LV.M()/1.0e6;//convert to GeV^2
    recoil.BeamPz = Pz_beam;
    recoil.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil.a2_ip = Tr.TrEvent[IsEject[2]].IntPoint;
    recoil.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    recoil_qval.a1_ip = Tr.TrEvent[IsEject[1]].IntPoint;
    recoil_qval.a2_ip = Tr.TrEvent[IsEject[2]].IntPoint;
    recoil_qval.p_ip = Tr.TrEvent[IsEject[0]].IntPoint;
    recoil_qval.BeamPz = Pz_beam_qval;
    recoil_qval.BeamKE = BeamKE;
    recoil_qval.BeamKE_eject = BeamKE_avg;
    recoil_qval.Theta = recoil_LV.Theta()*180.0/TMath::Pi();
    recoil_qval.Phi = recoil_LV.Phi()*180.0/TMath::Pi();
    recoil_qval.KE_recoil = recoil_LV.E()-recoil_LV.M();
    recoil_qval.theta_eject1 = eject2_LV.Theta()*180.0/TMath::Pi();
    recoil_qval.theta_eject2 = eject3_LV.Theta()*180.0/TMath::Pi();
    recoil.theta_eject1 = eject2_LV.Theta()*180.0/TMath::Pi();
    recoil.theta_eject2 = eject3_LV.Theta()*180.0/TMath::Pi();
    recoil.Theta = recoil_LV.Theta()*180.0/TMath::Pi();
    recoil.Phi = recoil_LV.Phi()*180.0/TMath::Pi();
    recoil.KE_recoil = recoil_LV.E() - recoil_LV.M();
  } else {
    cout<<"Incorrect Reconstruction in InvMassSelect! Should have 3 ejectiles"<<endl;
  }
}

