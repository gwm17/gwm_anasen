/*gwmAnalyzer.cpp
 *Analyzer class for ANASEN Detector in Active Target mode. Contains all functions necessary to 
 *sort data and calculate physical values. Will do everything from making all of the tracking data
 *up to calling the experiment dependant reconstruction and sorting and storing all data. Current
 *asks user for the name of a data list file, an output file, and three cut files
 *
 *Gordon M. -- April 2019
 *Based on previous versions written by M. Anastasiou, N. Rijal, J. Parker, et al
 */
#include "gwmAnalyzer.h"

using namespace std;

//Constructor
analyzer::analyzer() {

  //Set all of the flags
  PCPlots = 0;
  PCPlots_2 = 0;
  Beam_and_Eloss = 0;
  FillTree = 0;
  FillEdE_cor = 0;
  CheckBasic = 0;
  MCP_RF_Cut = 0;
  ReadPCWire = 0;
  PCWireCal = 0;
  ResidualEnergy_Calc = 0;
  ResidualEnergy_Calc_20Ne = 0;
  Reconstruction_Session = 0;
  Elastic_Scat_Alpha = 0;
  cutFlag = 0;

}

//Destructor
analyzer::~analyzer() {
  //free dynamic memory
  delete Ne18_eloss;
  delete Ne20_eloss;
  delete Na21_eloss;
  delete proton_eloss;
  delete alpha_eloss;
}

/*SetFlag()
 *To compartmentalize the analysis, uses many different functions
 *To make the program still customizable, use flags in the main run progress to pick
 *and choose which analysis is run; Flags should be set by feeding an all LOWERCASE name string
 */
void analyzer::SetFlag(string flagName) {
  
  if (flagName == "pcplots") PCPlots = 1;
  else if (flagName == "pcplots2") PCPlots_2 = 1;
  else if (flagName == "beamandeloss") Beam_and_Eloss = 1;
  else if (flagName == "filltree") FillTree = 1;
  else if (flagName == "filledecor") FillEdE_cor = 1;
  else if (flagName == "checkbasic") CheckBasic = 1;
  else if (flagName == "mcprfcut") MCP_RF_Cut = 1;
  else if (flagName == "readpcwire") ReadPCWire = 1;
  else if (flagName == "pcwirecal") PCWireCal = 1;
  else if (flagName == "residualenergycalc") ResidualEnergy_Calc = 1;
  else if (flagName == "residualenergycalc20ne") ResidualEnergy_Calc_20Ne = 1;
  else if (flagName == "reconstructionsession") Reconstruction_Session = 1;
  else if (flagName == "elasticscatalpha") Elastic_Scat_Alpha = 1;
  else if (flagName == "cutflag") cutFlag = 1;
  else {
    cout<<"Invalid flag given as argument!"<<endl;
  }

}

/*getCut()
 *Grabs the different cuts for the program
 *Cuts should be stored in a file and feed through this function
 *Cuts can all be in one file or in many files, just make sure names match
 */
void analyzer::getCut(string pcutfile, string acutfile, string ncutfile) {
  cout<<"Proton cut file: "<<pcutfile<<endl;
  cout<<"Alpha cut file: "<<acutfile<<endl;
  cout<<"Needle cut file: "<<ncutfile<<endl;
  char protonname[pcutfile.length()], alphaname[acutfile.length()], needlename[ncutfile.length()];
  strcpy(protonname, pcutfile.c_str());
  strcpy(alphaname, acutfile.c_str());
  strcpy(needlename, ncutfile.c_str());

  TFile *protonfile = new TFile(protonname, "READ");
  protonCut = (TCutG*) protonfile->Get("protons_edetheta_q3r1_run2007_2044_18Necut_09132018");
  protonfile->Close();

  TFile *alphafile = new TFile(alphaname, "READ");
  alphaCut = (TCutG*) alphafile->Get("alphas_edetheta_q3r1_run2007_2044_18Necut_09132018");
  alphafile->Close();

  TFile *needlefile = new TFile(needlename, "READ");
  needleCut = (TCutG*) needlefile->Get("resEn_vs_neEn_348Torr_PCZfactor_09202019");
  needlefile->Close();
}

/*recoilReset()
 *With RecoilEvent being promoted to a global structure there needs to be
 *a method for setting the class instance of RecoilEvent back to zero
 *Was originally contained in the Track class
 */
void analyzer::recoilReset() {
  recoil.IntPoint = -10;
  recoil.BeamEnergy = -10;
  recoil.BeamEnergy_20Ne = -10;
  recoil.SiEnergy_tot = -10;
  recoil.SiEnergy_calc = -10;
  recoil.PCEnergy_tot = -10;
  recoil.Energy_tot = -10;
  recoil.Theta = -10;
  recoil.Theta_20Ne = -10;
  recoil.Theta_Qvalue_20Ne = -10;
  recoil.Phi = -10;
  recoil.Phi_20Ne = -10;
  recoil.Phi_Qvalue_20Ne = -10;
  recoil.KE = -10;
  recoil.Ex = -10;
  recoil.KE_20Ne = -10;
  recoil.KE_Qvalue_20Ne = -10;
  recoil.Ex_20Ne = -10;
  recoil.Ex_Qvalue_20Ne = -10;
  recoil.Ex_rec_small = -10;
  recoil.Ex_rec_large = -10;
  recoil.Ex_rec_qvalue_small = -10;
  recoil.Ex_rec_qvalue_large = -10;
  recoil.Delta_Theta = -10;
  recoil.Delta_Phi = -10;
  recoil.BeamWA_20Ne = -10;
  recoil.BeamWA_Qvalue_20Ne = -10;
  recoil.IntP_20Ne = -10;
  recoil.IntP_Qvalue_20Ne = -10;
  recoil.theta_p1_20Ne = -10;
  recoil.theta_p2_20Ne = -10;
  recoil.Ex_21Na = -10;
  recoil.Beam_Qv_21Na = -10;
}

/*GetPCWireRadius()
 *Reads the pc wire radius file, and stores in a vector
 *for use in the program
 *PCWireRadius file is a txt file not a root file
 */
vector<Double_t> analyzer::GetPCWireRadius() {
  ifstream pcrfile;
  string pcrName;
  cout<<"Enter full path name to PC wire radius file: ";
  cin>>pcrName;
  char pcrname[200];
  strcpy(pcrname, pcrName.c_str());
  pcrfile.open(pcrname);
  if(pcrfile.is_open()) {
    cout<<"Reading PC wire radius file: "<<pcrName<<" ..."<<endl;
    string line1;
    getline(pcrfile, line1);
    cout<<line1<<endl;

    Int_t wireNum;
    Double_t rad;
    vector<Double_t> WireRad;
    WireRad.resize(NPCWires);
    while(pcrfile >> wireNum) {
      pcrfile >> rad;
      WireRad[wireNum] = rad;
    }
    return WireRad;
  } else {
    cout<<"PC wire radius file "<<pcrName<<" does not exist!"<<endl;
    exit(EXIT_FAILURE);
  }
}

/*MyFill()
 *A method for creating, storing, and filling ROOT histograms
 *Looks for histogram stored in a map, if doesn't exist creates a new histogram which is then
 *stored in the map and in the ROOT object array
 *IMPROVEMENT: Look and see if TObjArray has any built in method for checking if there is an
 * object; would cut down on amount of dynamic memory allocated
 */
void analyzer::MyFill(string name, int binsX, double lowX, double highX, double valueX) {
  try{
    fhmap.at(name)->Fill(valueX);
  } catch(out_of_range e) {
    TH1F* histo = new TH1F(name.c_str(), name.c_str(), binsX, lowX, highX);
    histo->Fill(valueX);
    rootObj->Add(histo);
    fhmap[name] = histo;
  }
}

void analyzer::MyFill(string name, int binsX, double lowX, double highX, double valueX,
                                   int binsY, double lowY, double highY, double valueY) {
  try{ 
    fhmap.at(name)->Fill(valueX, valueY);
  } catch(out_of_range e) {
    TH2F* histo = new TH2F(name.c_str(), name.c_str(), binsX, lowX, highX, binsY, lowY, highY);
    histo->Fill(valueX, valueY);
    rootObj->Add(histo);
    fhmap[name] = histo;
  }
}

//PhiDiff() Calculates angular difference from 0 to pi
Double_t analyzer::PhiDiff(Float_t phi1, Float_t phi2) {
  Double_t phiDiff = abs(phi1-phi2);
  if(phiDiff>TMath::Pi() && phiDiff<= 2*TMath::Pi()) {
    phiDiff = 2*TMath::Pi() - phiDiff;
  }
  return phiDiff;
}

/*FindMaxPC()
 *Looks for the PC wire hit in the event with the maximum PC Energy
 */
Int_t analyzer::FindMaxPC(Double_t phi, PCHit &PC) {//Pretty sure the PCHit doesn't need ref
  
  Int_t MaxPCindex = -1;
  Double_t MaxPC = -10;

  //phi window
  Double_t phiMin = 30.0*TMath::Pi()/180.0;

  //loop over hits
  for (int k=0; k<PC.NPCHits; k++) {
    PC.pc_obj = (*PC.ReadHit)[k];
    if(PC.pc_obj.Energy>=0.0 && PhiDiff(PC.pc_obj.PhiW, phi) < phiMin) {//toss empties & check win
      if( PC.pc_obj.Energy >= MaxPC) { 
        MaxPC = PC.pc_obj.Energy;
        MaxPCindex = k;
      } 
    }
  }
  return MaxPCindex;

}

/*MCP_RF()
 *First of the analysis methods; Makes and fills histograms for the MCP and RF
 *timing. Also performs the major cut on the data, the MCP_RF cut, where the correct 
 *beam for the experiment is selected
 */
bool analyzer::MCP_RF() {

  ICne_E_sum = input_ICne_E_sum;
  ICne_E_diff = input_ICne_E_diff;
  ICne_T_diff = input_ICne_T_diff;
  ICne_T_sum = input_ICne_T_sum;
  ICne_E_diff_cal = 0.0116*ICne_E_diff;
  MyFill("ICne_En_sum",1024,0,4096,ICne_E_sum);
  MyFill("ICne_En_diff",1024,0,4096,ICne_E_diff);
  MyFill("ICne_En_diff_cal",200,0,70,ICne_E_diff_cal);
  MyFill("ICne_Ti_diff",1024,0,4096,ICne_T_diff);
  MyFill("ICne_Ti_sum",1024,0,4096,ICne_T_sum);
  MyFill("Sum_vs_Diff_Time_ICne",256,0,4096,ICne_T_diff,256,0,4096,ICne_T_sum);
  MyFill("Energy_vs_Time_Diff_ICne",256,0,4096,ICne_T_diff,256,0,4096,ICne_E_diff);

  Double_t correct = 1.004009623;
  MCPTime = input_MCPTime;
  RFTime = input_RFTime;
  TDC2 = input_TDC2;
  MyFill("MCPTime",1024,1,4096,MCPTime);
  MyFill("RFTime",1024,1,4096,RFTime);
  MyFill("TDC2",1024,1,4096,TDC2);
  MyFill("Timing",600,0,600,fmod((MCPTime*correct-RFTime),545));
  if (MCP_RF_Cut) {
    if (MCPTime>0 && RFTime>0) {
      float wrap_val = fmod((MCPTime*correct-RFTime), 545);
      if((wrap_val<40 || wrap_val>100) && (wrap_val<310 || wrap_val>370)){
        return false;
      } 
      if (MCPTime<2870 || MCPTime>3050) {
        return false;
      }
    } else {
        return false;
    }
  }
  MyFill("MCPTime_cut",1000,1,5000,MCPTime);
  MyFill("RFTime_cut",1000,1,5000,RFTime);
  MyFill("TDC2_cut",1000,1,5000,TDC2);
  MyFill("Timing_cut",600,0,600,fmod((MCPTime*correct-RFTime),545));
  
  MyFill("Timing_cut_ICne_En_sum",1024,0,4096,ICne_E_sum);
  MyFill("Timing_cut_ICne_En_diff",1024,0,4096,ICne_E_diff);
  MyFill("Timing_cut_ICne_En_diff_cal",200,0,70,ICne_E_diff_cal);
  MyFill("Timing_cut_ICne_Ti_diff",1024,0,4096,ICne_T_diff);
  MyFill("Timing_cut_ICne_Ti_sum",1024,0,4096,ICne_T_sum);
  MyFill("Timing_cut_Sum_vs_Diff_Time_ICne",256,0,4096,ICne_T_diff,256,0,4096,ICne_T_sum);
  MyFill("Timing_cut_Energy_vs_Time_Diff_ICne",256,0,4096,ICne_T_diff,256,0,4096,ICne_E_diff);
  return true;
}

/*Track1()
 *Creates and stores what are refered to as track1 type events. This is where there is both
 *a good Si event and a good PC event. Actually is a bit tricky to do. Basically use FindMaxPC()
 *to pick out if there is a good PC in range of the Si event; else move on and keep looking
 *For later reference the good PC event information is stored *in PCGoodEnergy & PCGoodPCZ and 
 *Si energy in SiEnergy_vec
 *Tracks are then sorted by the method specified in the Track class, and used events are marked 
 *by setting the Si energy and PC energy in the data to an unphysical value  
 */
void analyzer::Track1() {

  Int_t GoodPC = -1;

  for (int j = 0; j<24; j++) {
    PCGoodEnergy[j] = 0;
    PCGoodPCZ[j] = 0;
  }

  MyFill("Si_ReadHit_size",500,0,50,(*Si.ReadHit).size());
  
  for(unsigned int i=0; i<(*Si.ReadHit).size(); i++) {
    sievent sihit = (*Si.ReadHit)[i];
    if (sihit.Energy <= 0 ) {
      continue;
    }
    GoodPC = FindMaxPC(sihit.PhiW, PC);
    if(GoodPC > -1) {
      pcevent pchit = (*PC.ReadHit)[GoodPC];
      TrackEvent trackhit;
      trackhit.TrackType = 1;
      trackhit.SiEnergy = sihit.Energy;
      trackhit.SiTime = sihit.Time;
      trackhit.SiPhi = sihit.PhiW;
      trackhit.SiZ = sihit.ZW;
      trackhit.SiR = sihit.RW;
      trackhit.DetID = sihit.DetID;
      trackhit.SiBCh = sihit.BackChannel;
      trackhit.HitType = sihit.HitType;

      trackhit.PCEnergy = pchit.Energy;
      if (trackhit.PCEnergy < 0.0) {
        cout << "PCEnergy is too low: "<<trackhit.PCEnergy<<endl;
        exit(EXIT_FAILURE);
      }
      trackhit.PCZ = pchit.ZW;
      trackhit.PCR = pchit.RW;
      trackhit.WireID = pchit.WireID;
      trackhit.Down = pchit.Down;
      trackhit.Up = pchit.Up;
      trackhit.DownVoltage = pchit.DownVoltage;
      trackhit.UpVoltage = pchit.UpVoltage;
      if(trackhit.WireID==6 || trackhit.WireID==12 || trackhit.WireID==17) {
        trackhit.PCZ = -10.0;
        trackhit.PCEnergy = -10.0;
        continue;
      }
      
      PCGoodEnergy[trackhit.WireID] = trackhit.PCEnergy;
      PCGoodPCZ[trackhit.WireID] = trackhit.PCZ;
      SiEnergy_vec[trackhit.DetID].push_back(trackhit.SiEnergy);

      (*Si.ReadHit)[i].Energy = -1000.0;
      (*PC.ReadHit)[GoodPC].Energy = -10.0;
      
      tracks.NTracks++;
      tracks.NTracks1++;
      tracks.TrEvent.push_back(trackhit);
    } 
  }
  sort(tracks.TrEvent.begin(), tracks.TrEvent.begin()+tracks.NTracks1, tracks.Tr_Sisort_method);

}

/*Track2()
 *Do for Track2 as was done to Track1, where track2 type events are just Si hits
 *sans the PC info
 *Again, sorted with method specified in Track
 */
void analyzer::Track2() {
  for (unsigned int i = 0; i<(*Si.ReadHit).size(); i++) {
    sievent sihit  = (*Si.ReadHit)[i];
    if (sihit.Energy == -1000) {
      continue;
    } else if (sihit.Energy <= 0) {//Why two conditions? They overlap
      continue;
    } else {
      TrackEvent trackhit;
      trackhit.TrackType = 2;
      trackhit.SiEnergy = sihit.Energy;
      trackhit.SiTime = sihit.Time;
      trackhit.SiPhi = sihit.PhiW;
      trackhit.SiZ = sihit.ZW;
      trackhit.SiR = sihit.RW;
      trackhit.DetID = sihit.DetID;
      trackhit.SiBCh = sihit.BackChannel;
      trackhit.HitType = sihit.HitType;
      SiEnergy_vec[trackhit.DetID].push_back(trackhit.SiEnergy);
      tracks.NTracks2++;
      tracks.NTracks++;
      tracks.TrEvent.push_back(trackhit);
    }
  }
  sort(tracks.TrEvent.begin()+tracks.NTracks1, 
       tracks.TrEvent.begin()+tracks.NTracks1+tracks.NTracks2, 
       tracks.Tr_Sisort_method);
}

/*Track3()
 *And lo there do be track3 type events as well. Recipie is mostly the same, where track3's 
 *are defined as PC only events
 *Again, sort as specified in Track
 */
void analyzer::Track3() {

  MyFill("PC_ReadHit_size",500,0,50,(*PC.ReadHit).size());
  for (unsigned int i=0; i<(*PC.ReadHit).size(); i++) {
    pcevent pchit = (*PC.ReadHit)[i];
    if(pchit.Energy == -10) {
      continue;
    } else if (pchit.Energy <= 0) {//Again, not exclusive conditions. Necessary?
      continue;
    } else {
      TrackEvent trackhit;
      trackhit.TrackType = 3;
      trackhit.PCEnergy = pchit.Energy;
      if(trackhit.PCEnergy < 0.0) {
        cout<<" PC Energy is too low! "<<trackhit.PCEnergy<<endl;
        exit(EXIT_FAILURE);
      }
      trackhit.WireID = pchit.WireID;
      trackhit.PCZ = pchit.ZW;
      trackhit.PCR = pchit.RW;
      trackhit.PCPhi = pchit.PhiW;
      trackhit.Down = pchit.Down;
      trackhit.Up = pchit.Up;
      trackhit.DownVoltage = pchit.DownVoltage;
      trackhit.UpVoltage = pchit.UpVoltage;

      if (trackhit.WireID==6 || trackhit.WireID==12 || trackhit.WireID==17) {
        trackhit.PCZ = -10.0;
        trackhit.PCEnergy = -10.0;
        continue;
      }
      
      tracks.NTracks3++;
      tracks.NTracks++;
      if(abs(PCGoodEnergy[trackhit.WireID])>1.e-7) {//is this bad?
        cout<<"wire double fired"<<endl;
      }
      PCGoodEnergy[trackhit.WireID] = trackhit.PCEnergy;
      PCGoodPCZ[trackhit.WireID] = trackhit.PCZ;
      tracks.TrEvent.push_back(trackhit);
    }
  }
  sort(tracks.TrEvent.begin()+tracks.NTracks1+tracks.NTracks2,
       tracks.TrEvent.end(), tracks.Tr_PCsort_method);
 
}

/*PCPlotting()
 *Makes plots for basic PC properties
 *Really nothing too crazy
 */
void analyzer::PCPlotting() {

  for(int i=0; i<tracks.NTracks; i++) {
    MyFill("WireID_vs_PCEnergy",24,0,24,tracks.TrEvent[i].WireID,500,0,1.5,
           tracks.TrEvent[i].PCEnergy);
    for(int j=0; j<tracks.NTracks1; j++) {
      Int_t value = (tracks.TrEvent[i].WireID-tracks.TrEvent[j].WireID)%24;
      MyFill("WireID_mod1_vs_PCEnergy_0",24,0,24,value,500,0,1.5,tracks.TrEvent[i].PCEnergy);
    } 
    if(tracks.TrEvent[i].WireID>-1 && tracks.TrEvent[i].WireID<23) {
      int id1 = tracks.TrEvent[i].WireID;
      int id2 = id1+1;
      MyFill(Form("PCEnergy_wire%i_wire%i_ntracksAllPlus",id2,id1), 200,0,1.5,PCGoodEnergy[id1],
                  200,0,1.5,PCGoodEnergy[id2]);
    }
  }

  if(tracks.NTracks1 == 1) {
    for(int wire=0; wire<24; wire++) {
      MyFill(Form("PCEnergy_vs_PCWireID_%i",tracks.TrEvent[0].WireID),24,0,24,wire,200,0,1.5,
             PCGoodEnergy[wire]/PCGoodEnergy[tracks.TrEvent[0].WireID]);
    }
    if(tracks.TrEvent[0].WireID>-1 && tracks.TrEvent[0].WireID<23) {
      int id1 = tracks.TrEvent[0].WireID;
      int id2 = id1+1;
      MyFill(Form("PCEnergyPlus_wire%i_wire%i_ntracks1_1",id2,id1), 200,0,1.5,PCGoodEnergy[id1],
                  200,0,1.5,PCGoodEnergy[id2]);
    }
  }  
}

/*TrackCalc()
 *Method of checks! Goal is to run and make sure that all of the tracking info
 *is good up to this point (which should be MCPRF, Tracks 1,2,3). Calculates and stores
 *the interaction point and beam energy for each event 
 *IMPROVEMENT: CheckBasic should probably be removed, these tests should be done everytime,
 *and removing the extra if statements would help boost the speed (probably not major, but every
 *little bit helps)
 *OPTIMIZATION: Look at how some of the later methods work and make this like those. Probably is
 *more readable/clean to locally store all of these track values at the top of the loop instead
 *of constantly reaccesing the same TrackEvent (Probably no performance gain tho :( )
 */
void analyzer::TrackCalc() {

  for(int i=0; i<tracks.NTracks1; i++) {
    if(CheckBasic) {
      if(tracks.TrEvent[i].SiZ < 0.0 || tracks.TrEvent[i].PCZ < 0.0) {
        tracks.TrEvent[i].PCZ = -100;
        tracks.TrEvent[i].SiZ = -100;
        continue;
      }
      if(tracks.TrEvent[i].PCR>4.0 || tracks.TrEvent[i].PCR<3.5) {
        tracks.TrEvent[i].PCR = -100;
        continue;
      }
      if(tracks.TrEvent[i].SiR>11.0 || tracks.TrEvent[i].SiR<4.0) {
        tracks.TrEvent[i].SiR = -100;
        continue;
      }
    }
    
    Double_t m = (tracks.TrEvent[i].PCR-tracks.TrEvent[i].SiR)/
                 (tracks.TrEvent[i].PCZ-tracks.TrEvent[i].SiZ);
    
    if(CheckBasic && (m==0 || tracks.TrEvent[i].PCZ<=0 || tracks.TrEvent[i].SiZ<0 || 
      tracks.TrEvent[i].SiR<=0 || (tracks.TrEvent[i].PCZ-tracks.TrEvent[i].SiZ) ==0)) {
      cout<<"m for Interaction Point is zero "<<m<<endl;
      continue;
    }

    Double_t b = (tracks.TrEvent[i].PCR - m*tracks.TrEvent[i].PCZ);
    tracks.TrEvent[i].IntPoint = -b/m;

    if(CheckBasic && (tracks.TrEvent[i].IntPoint<0.0 || tracks.TrEvent[i].IntPoint>55.0)) {
      tracks.TrEvent[i].IntPoint = -100;
      continue;
    }
    
    if((tracks.TrEvent[i].IntPoint - tracks.TrEvent[i].SiZ)>0) {
      tracks.TrEvent[i].Theta = atan(tracks.TrEvent[i].SiR/
                                     (tracks.TrEvent[i].IntPoint-tracks.TrEvent[i].SiZ));
      tracks.TrEvent[i].PathLength = tracks.TrEvent[i].SiR/sin(tracks.TrEvent[i].Theta);
    } else if ((tracks.TrEvent[i].IntPoint-tracks.TrEvent[i].SiZ) < 0) {
      tracks.TrEvent[i].Theta = TMath::Pi() + atan(tracks.TrEvent[i].SiR/
                                (tracks.TrEvent[i].IntPoint-tracks.TrEvent[i].SiZ));
      tracks.TrEvent[i].PathLength = tracks.TrEvent[i].SiR/sin(tracks.TrEvent[i].Theta);
    } else { 
      tracks.TrEvent[i].Theta = TMath::Pi()/2.0;
      tracks.TrEvent[i].PathLength = tracks.TrEvent[i].SiR;
    }

    if (Beam_and_Eloss) {
      float length_check = ana_length - tracks.TrEvent[i].IntPoint;
      if(length_check>0.0 && length_check<ana_length) {   
        tracks.TrEvent[i].BeamEnergy = Ne18_eloss->GetLookupEnergy(BeamE, length_check);
        if(CheckBasic && (tracks.TrEvent[i].BeamEnergy<0.0 || 
                          tracks.TrEvent[i].BeamEnergy>BeamE)){
          tracks.TrEvent[i].BeamEnergy = -100.0;
          continue;
        }
      }
    }
  }
}

/*PCWireCalibration()
 *Option for calibrating pc wires; 
 *Makes a very large number of plots, so only run if necessary
 *Requires gold_pos, which has to be entered manually
 */
void analyzer::PCWireCalibration() {

  Double_t mpc, bpc;
  Double_t qqq_r1_Emin = 11.3;
  Double_t r2_Emin = 10.3;
  Double_t r2_Emax = 11.2;

  Double_t gold_pos = 0.0; //set by experiment

  for (int i=0; i<tracks.NTracks1; i++) {
    if((tracks.TrEvent[i].DetID>-1 && tracks.TrEvent[i].DetID<16 && 
      tracks.TrEvent[i].SiEnergy>qqq_r1_Emin) || (tracks.TrEvent[i].DetID>15 && 
      tracks.TrEvent[i].DetID<28 && tracks.TrEvent[i].SiEnergy>r2_Emin && 
      tracks.TrEvent[i].SiEnergy<r2_Emax)) {
        
      mpc = tracks.TrEvent[i].SiR/(tracks.TrEvent[i].SiZ-gold_pos);
      bpc = tracks.TrEvent[i].SiR - mpc*tracks.TrEvent[i].SiZ;
      tracks.TrEvent[i].pcz_ref = (WireRadii[tracks.TrEvent[i].WireID]-bpc)/mpc;

      Int_t wireid = tracks.TrEvent[i].WireID;
      Int_t detid = tracks.TrEvent[i].DetID;
      Double_t siz = tracks.TrEvent[i].SiZ;
      Double_t pcz = tracks.TrEvent[i].PCZ;
      Double_t sir = tracks.TrEvent[i].SiR;
      Double_t sie = tracks.TrEvent[i].SiEnergy;
      Double_t pce = tracks.TrEvent[i].PCEnergy;
      Double_t pczref = tracks.TrEvent[i].pcz_ref;
      Double_t pczdiff = pcz-pczref;
      Double_t theta = tracks.TrEvent[i].Theta;
      Double_t siphi = tracks.TrEvent[i].SiPhi;
      Double_t pcphi = tracks.TrEvent[i].PCPhi;
      Double_t intp = tracks.TrEvent[i].IntPoint;

      MyFill(Form("SiZ_vs_SiEnergy_WireID_%i",wireid),600,-1,30,sie,600,-1,25,siz);
      if(detid<28 && detid>3)
	MyFill("SiZ_vs_SiEnergy_R1_R2",600,-1,30,sie,600,-1,25,siz);
      if(detid<4 && detid>-1)
	  MyFill("SiR_vs_SiEnergy_QQQ3",600,-1,30,sie,600,3,12,sir);

      MyFill(Form("PCZ_vs_Z_afterCal_%i",wireid),400,-10.0,30.0,pcz,400,-10.0,30.0,pczref); 
      MyFill(Form("PCZoffset_vs_PCZ_afterCal_%i",wireid),600,-1,30,pcz,600,-20,20,(pczdiff));
      MyFill(Form("PCZoffset_vs_PCEnergy_%i",wireid),600,-0.1,0.5,pce,600,-20,20,(pczdiff));
      MyFill(Form("PCZoffset_vs_Theta_%i",wireid),600,0,190,theta*rads2deg,600,-20,20,(pczdiff));

      if(detid<4 && detid>-1){ // qqq
	MyFill("PCZ_vs_Z_afterCal_q3",400,-10.0,30.0,pcz,400,-10.0,30.0,pczref);
	MyFill(Form("PCZ_vs_Z_afterCal_q3_%i",wireid),400,-10.0,30.0,pcz,400,-10.0,30.0,pczref); 
	MyFill(Form("PCZoffset_vs_PCZ_afterCal_q3_%i",wireid),600,-1,30,pcz,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_vs_PCEnergy_q3_%i",wireid),600,-0.1,0.5,pce,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_Theta_q3_%i",wireid),600,0,190,theta*rads2deg,600,-20,20,pczdiff);
      }

      if(detid<16 && detid>3){ // r1
	MyFill(Form("PCZ_vs_Z_afterCal_r1_%i",wireid),400,-10.0,30.0,pcz,400,-10.0,30.0,pczref); 
	MyFill("PCZ_vs_Z_afterCal_r1",400,-10.0,30.0,pcz,400,-10.0,30.0,pczref);
	MyFill(Form("PCZoffset_vs_PCZ_afterCal_r1_%i",wireid),600,-1,30,pcz,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_vs_PCEnergy_r1_%i",wireid),600,-0.1,0.5,pce,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_Theta_r1_%i",wireid),600,0,190,theta*rads2deg,600,-20,20,pczdiff);
      }

      if(detid<28 && detid>15){ // r2
	MyFill(Form("PCZ_vs_Z_afterCal_r2_%i",wireid),400,-10.0,30.0,pcz,400,-10.0,30.0,pczref); 
	MyFill("PCZ_vs_Z_afterCal_r2",400,-10.0,30.0,pcz,400,-10.0,30.0,pczref);
	MyFill(Form("PCZoffset_vs_PCZ_afterCal_r2_%i",wireid),600,-1,30,pcz,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_vs_PCEnergy_r2_%i",wireid),600,-0.1,0.5,pce,600,-20,20,(pczdiff));
	MyFill(Form("PCZoffset_Theta_r2_%i",wireid),600,0,190,theta*rads2deg,600,-20,20,pczdiff);
      }

      if(detid<16 && detid>-1){
	MyFill("PCEnergy_vs_PCZ_q3_r1",600,-10,60,pcz,600,-0.1,0.5,pce);
	MyFill("PCZoffset_vs_PCEnergy_q3_r1",600,-0.1,0.5,pce,600,-20,20,(pczdiff));
	MyFill("PCZoffset_vs_PCZ_q3_r1",600,-1,30,pcz,600,-20,20,(pczdiff));
	MyFill("PCZoffset_vs_Theta_q3_r1",600,0,190,theta*rads2deg,600,-20,20,(pczdiff));
      }

      MyFill("PCZ_vs_Z_afterCal_All",400,-10.0,30.0,pcz,400,-10.0,30.0,pczref);
      MyFill("PCEnergy_vs_PCZ_All",600,-10,60,pcz,600,-0.1,0.5,pce);
      MyFill("PCZoffset_vs_PCEnergy_All",600,-0.1,0.5,pce,600,-20,20,(pczdiff));
      MyFill("PCZoffset_vs_PCZ_All",600,-1,30,pcz,600,-20,20,(pczdiff));
      MyFill("PCZoffset_vs_Theta_All",600,0,190,theta*180/TMath::Pi(),600,-20,20,(pczdiff));

      // for all the wires together the Phi_s make more sense
      MyFill("PCZoffset_vs_PCPhi",600,0,360,pcphi*180/TMath::Pi(),600,-20,20,(pczdiff));
      MyFill("PCZoffset_vs_SiPhi",600,0,360,siphi*180/TMath::Pi(),600,-20,20,(pczdiff));

      //E_de && E_Theta && IntPoint
      MyFill("IntPoint_cut_ALL",600,-10,60,intp);
      MyFill(Form("IntPoint_cut_ALL_WireID_%i",wireid),600,-10,60,intp);  
      MyFill("E_de_corrected_cut_ALL",600,-1,35,sie,600,-0.01,0.35,pce *sin(theta));
      MyFill("E_si_vs_Theta_cut_ALL",500,0,200,theta*rads2deg,500,0,35,sie);
      MyFill("E_si_vs_IntPoint_cut_ALL",500,0,60,intp,500,0,35,sie);

      if(detid<4 && detid>-1){
	MyFill("IntPoint_cut_Q3",600,-10,60,intp);
	MyFill(Form("IntPoint_cut_Q3_WireID_%i",wireid),600,-10,60,intp);
	MyFill("E_de_corrected_cut_Q3",600,-1,35,sie,600,-0.01,0.35,pce *sin(theta));
	MyFill("E_si_vs_Theta_cut_Q3",500,0,200,theta*rads2deg,500,0,35,sie);
	MyFill("E_si_vs_IntPoint_cut_Q3",500,0,60,intp,500,0,35,sie);
      }
      if(detid<16 && detid>3){
	MyFill("IntPoint_cut_SX3_R1",600,-10,60,intp);
	MyFill(Form("IntPoint_cut_R1_WireID_%i",wireid),600,-10,60,intp);
	MyFill("E_de_corrected_cut_SX3_R1",600,-1,35,sie,600,-0.01,0.35,pce *sin(theta));
	MyFill("E_si_vs_Theta_cut_SX3_R1",500,0,200,theta*rads2deg,500,0,35,sie);
	MyFill("E_si_vs_IntPoint_cut_SX3_R1",500,0,60,intp,500,0,35,sie);
      }

      if(detid<28 && detid>15){
	MyFill("IntPoint_cut_SX3_R2",600,-10,60,intp);
	MyFill(Form("IntPoint_cut_R2_WireID_%i",wireid),600,-10,60,intp);
	MyFill("E_de_corrected_cut_SX3_R2",600,-1,35,sie,600,-0.01,0.35,pce *sin(theta));
	MyFill("E_si_vs_Theta_cut_SX3_R2",500,0,200,theta*rads2deg,500,0,35,sie);
	MyFill("E_si_vs_IntPoint_cut_SX3_R2",500,0,60,intp,500,0,35,sie);
      }
    }
  }
}

/*EdEcor()
 *Makes corrected EdE plots
 *Slices by detector type
 */
void analyzer::EdEcor() {
  for (int i=0; i<tracks.NTracks1; i++) {
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t theta = tracks.TrEvent[i].Theta;
    Int_t detid = tracks.TrEvent[i].DetID;

    MyFill("BeamEnergy_vs_IntPoint",100,-20,80,intp,100,-10,90,be);
    MyFill("BeamEnergy",100,-10,90,be);
    MyFill("E_de",200,-1,35,sie,200,-0.01,1.5,pce);
    MyFill("E_de_corrected_ALL",200,-1,35,sie,200,-0.01,1.5,pce *sin(theta));
    MyFill("InteractionPoint",300,-10,56,intp);	
    MyFill("E_si_vs_Theta",500,0,200,theta*rads2deg,500,0,35,sie);
    if(detid<4 && detid>-1){
      MyFill("E_de_Q3",200,-1,35,sie,200,-0.01,1.5,pce);
      MyFill("E_de_corrected_Q3",200,-1,35,sie,200,-0.01,1.5,pce *sin(theta));
      MyFill("InteractionPoint_Q3",300,-10,56,intp);
      MyFill("E_si_vs_Theta_Q3",500,0,200,theta*rads2deg,500,0,35,sie);
    }
    if(detid<16 && detid>3){
      MyFill("E_de_SX3_1",200,-1,35,sie,200,-0.01,1.5,pce);
      MyFill("E_de_corrected_SX3_1",200,-1,35,sie,200,-0.01,1.5,pce *sin(theta));
      MyFill("InteractionPoint_SX3_1",300,-10,56,intp);
      MyFill("E_si_vs_Theta_SX3_1",500,0,200,theta*rads2deg,500,0,35,sie);
    }
    if(detid<28 && detid>15){
      MyFill("E_de_SX3_2",200,-1,35,sie,200,-0.01,1.5,pce);
      MyFill("E_de_corrected_SX3_2",200,-1,35,sie,200,-0.01,1.5,pce *sin(theta));
      MyFill("InteractionPoint_SX3_2",300,-10,56,intp);
      MyFill("E_si_vs_Theta_SX3_2",500,0,200,theta*rads2deg,500,0,35,sie);
    }
  }
}

/*PCPlottting2()
 *PC plots again, but this one is  for the later stages of analysis
 *Pretty intensive, so again only run when necessary
 */
void analyzer::PCPlotting2() {
  Double_t PCMaxE = -1.0;
  Int_t PCMaxIndex = -1;

  for (int i=0; i<tracks.NTracks; i++) {
    if(tracks.TrEvent[i].PCEnergy>0.005 && tracks.TrEvent[i].PCEnergy<1.0) {
      PCMaxE = tracks.TrEvent[i].PCEnergy;
      PCMaxIndex = i;
    }
    for(int j=0; j<24; j++) {
      if(PCGoodEnergy[j]>0.005 && PCGoodEnergy[j]<1.0) {
        for(int k=0; k<24; k++) {
          if(PCGoodEnergy[k]>0.005 && PCGoodEnergy[k]<1.0) {
            MyFill("Wire_vs_wire_allTracks",24,0,24,j,24,0,24,k);
          }
        }
      }
    }
  }
  if(PCMaxIndex != -1) {
    MyFill("PCEnergy_vs_WireID",24,0,24,tracks.TrEvent[PCMaxIndex].WireID,500,0,1.0,PCMaxE);
    for(int j=0; j<24; j++) {
      if(PCGoodEnergy[j]>0.005 && PCGoodEnergy[j]<1.0) {
        MyFill("Wire_Vs_WireMax",24,0,24,j,24,0,24,tracks.TrEvent[PCMaxIndex].WireID);
      }
    }
  }

  for(int i=0; i<tracks.NTracks1; i++) {
    int wireplus = (tracks.TrEvent[i].WireID+1)%24;
    int wireminus = (tracks.TrEvent[i].WireID-1+24)%24;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t theta = tracks.TrEvent[i].Theta;
    Int_t wireid = tracks.TrEvent[i].WireID;
    if(PCGoodEnergy[wireplus]<0 || PCGoodEnergy[wireminus]<0) {
      MyFill("EdE_corrected_all_flag_lowerthresh_0",500,-1,35,sie,600,-0.01,0.8,pce*sin(theta));
    }
    int wire_fired = 0;
    for(int wire=0; wire<24; wire++) {
      if((wire != wireid) && (PCGoodEnergy[wire]>0.005)) {
        wire_fired =1;
      }
    }
    if(!wire_fired) {
      MyFill("Ede_corrected_ALL_flag_thres0.005_exludedWire",600,-1,35,sie,
             600,-0.01,0.8,pce*sin(theta));
    }
    if(protonCut->IsInside(sie, pce*sin(theta))) {
      for(int sid=0; sid<28; sid++) {
        for(unsigned int svid; svid<SiEnergy_vec[sid].size(); svid++) {
          Double_t sie_v = SiEnergy_vec[sid][svid];
          if((sid != tracks.TrEvent[i].DetID && abs(sie-sie_v)>1.0e-10) || sie_v>0 ) {
            MyFill("SiEnergy_track1_vs_SiEnergy_All",300,-1,35,sie_v,300,-1,35,sie);
          }
        } 
      }
      if(sie>0) {
        for(int t2=tracks.NTracks1; t2<(tracks.NTracks1+tracks.NTracks2); t2++) {
          Double_t sie2 = tracks.TrEvent[t2].SiEnergy;
          if(sie2>0) {
            MyFill("SiEnergy_track2_vs_SiEnergy_track1_inPCUT",300,0,35,sie,300,0,35,sie2);
          }
        }
      }
    }
  }
}

/*CalculateResidE()
 *First reconstruction method
 *Uses ReconstructHeavy from Reconstruction class
 *Focuses on Energy of the residual
 */
void analyzer::CalculateResidE() {
  Reconstruct Elastic2(M_18Ne, M_alpha, M_p, M_21Na);
  Elastic2.ELoss_light = proton_eloss;
  Elastic2.ELoss_beam = Ne18_eloss;

  for(int i=0; i<tracks.NTracks1; i++) {
    Int_t detid = tracks.TrEvent[i].DetID;
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t theta = tracks.TrEvent[i].Theta;
   
    if(detid>-1 && detid<28 && pcz>0.0 && siz>=0.0 && intp>0.0 && intp<ana_length) {
      if(protonCut->IsInside(sie, pce*sin(theta))) {
        Elastic2.ReconstructHeavy(tracks, i, recoil);
        Double_t etot = recoil.Energy_tot;
        Double_t ke = recoil.KE;
        Double_t rtheta = recoil.Theta;

	MyFill("SiEnergy_tot_vs_Needle_4vec",256,0,4096,ICne_E_diff,300,-0.1,30,etot);
	MyFill("SiEnergy_tot_vs_NeedleCal_4vec",256,0,70,ICne_E_diff_cal,300,-0.1,30,etot);
	MyFill("SiEnergy_tot_vs_RecoilEn_4vec",300,-0.1,80,ke,300,-0.1,30,etot);
	MyFill("LightTheta_vs_RecoilTheta",300,-0.1,180,rtheta,300,-0.1,180,theta*rads2deg);
	MyFill("RecoilEnergy_vs_Theta_4vec",300,-0.1,180,rtheta,300,-0.1,80,ke);
	MyFill("RecoilEnergy_vs_Needle_4vec",256,0,4096,ICne_E_diff,300,-0.1,80,ke);
	MyFill("RecoilEnergy_vs_NeedleCal_4vec",256,0,70,ICne_E_diff_cal,300,-0.1,80,ke);
	MyFill("RecoilEnergy_vs_IntP_4vec",300,-0.1,60,intp,300,-0.1,80,ke);
        
        Float_t r_fullpath, r_IntPtoNeedle, ResidE_recoil;
        if(rtheta != 0.0 && rtheta != 90.0 && rtheta != 180.0) { 
          r_fullpath = 2.2/sin(rtheta/rads2deg); //2.2? radius of ic
          if(intp > 17.8) {
            r_IntPtoNeedle  = (intp-17.8)/cos(rtheta/rads2deg);
          }
        }
        if(intp<=17.8) {
          Float_t E_fullpath = Na21_eloss->GetLookupEnergy(ke, r_fullpath);
          Float_t BeamE_atneedle = Ne18_eloss->GetLookupEnergy(BeamE, (ana_length-17.8));
          Float_t deltaBeamE = BeamE_atneedle - be;

          if(E_fullpath <= 0.0) {
            ResidE_recoil = ke+deltaBeamE;
            if(detid>-1 && detid<4) {
              MyFill("ResidualE_vs_Needle_4vec_Q3_insideVol_FullLoss",256,0,4096,ICne_E_diff,
                     300,-11,50,ResidE_recoil);
            }
            if(detid>3 && detid<16) {
              MyFill("ResidualE_vs_Needle_4vec_R1_insideVol_FullLoss",256,0,4096,ICne_E_diff,
                     300,-11,50,ResidE_recoil);
            }
          } else {
            ResidE_recoil = ke-E_fullpath+deltaBeamE;
            if(detid>-1 && detid<4) {
              MyFill("ResidualE_vs_Needle_4vec_Q3_insideVol_someLoss",256,0,4096,ICne_E_diff,
                     300,-11,50,ResidE_recoil);
            }
            if(detid>3 && detid<16) {
              MyFill("ResidualE_vs_Needle_4vec_R1_insideVol_someLoss",256,0,4096,ICne_E_diff,
                     300,-11,50,ResidE_recoil);
            }
          }
          
          MyFill("ReisdualE_Recoil_vs_Needle_NeedleVolume",256,0,4096,ICne_E_diff,
                  300,-11,50,ResidE_recoil); 
        } else{
          Float_t E_needlestart = Na21_eloss->GetLookupEnergy(ke, r_IntPtoNeedle);
          if(E_needlestart <= 0.0) {
            ResidE_recoil = -10.0;
          } else if (r_fullpath<r_IntPtoNeedle) {
            ResidE_recoil = -5.0;
          } else {
            Float_t Efinal = Na21_eloss->GetLookupEnergy(E_needlestart,
                                                         (r_fullpath-r_IntPtoNeedle));
            if(Efinal <= 0.0) {
              ResidE_recoil = E_needlestart;
            } else{
              ResidE_recoil = E_needlestart-Efinal;
            }
          }
          MyFill("ResidualE_recoil_vs_Needle_outside_NeedleVol",256,0,4096,ICne_E_diff,
                                                                300,-11,50,ResidE_recoil);
        } 
        tracks.TrEvent[i].ResidualEn = ResidE_recoil;
        
	MyFill("ResidualE_vs_Path_4vec",300,-0.1,60,r_fullpath,300,-11,50,ResidE_recoil);
	MyFill("ResidualE_vs_IntP_4vec",300,-0.1,60,intp,300,-11,50,ResidE_recoil);
	MyFill("ResidualEn_IntP_vs_Path_4vec",300,-0.1,60,r_fullpath,300,-0.1,60,intp);
	MyFill("ResidualEn_vs_Needle_4vec",256,0,4096,ICne_E_diff,300,-11,50,ResidE_recoil);
	MyFill("ResidualEn_NeedleCal4vec",256,0,70,ICne_E_diff_cal,300,-11,50,ResidE_recoil);	    
	MyFill("ResidualEn_vs_RecoilE_4vec",300,-0.1,80,ke,300,-11,50,ResidE_recoil);
	MyFill("ResidualEn_vs_Theta_4vec",300,-0.1,180,rtheta,300,-11,50,ResidE_recoil);

	if(detid>-1 && detid<4){
	  MyFill("LightParEn_4VEC_vs_BeamEnergy_Tracked_4vec_Q3",300,-0.1,80,be,300,-0.1,30,etot);
	  MyFill("RecoilEn_4VEC_vs_BeamEnergy_Tracked_4vec_Q3",300,-0.1,80,be,300,-0.1,30,ke);
	  MyFill("SiEnergy_tot_vs_Needle_4vec_Q3",256,0,4096,ICne_E_diff,300,-0.1,30,etot);
	  MyFill("SiEnergy_tot_vs_NeedleCal_4vec_Q3",256,0,70,ICne_E_diff_cal,300,-0.1,30,etot);
	  MyFill("SiEnergy_tot_vs_RecoilEn_4vec_Q3",300,-0.1,80,ke,300,-0.1,30,etot);
	  MyFill("RecoilEnergy_vs_Theta_4vec_Q3",300,-0.1,180,rtheta,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_Path_4vec_Q3",300,-0.1,60,r_fullpath,300,-0.1,80,ke);		
	  MyFill("RecoilEnergy_vs_Needle_4vec_Q3",256,0,4096,ICne_E_diff,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_NeedleCal_4vec_Q3",256,0,70,ICne_E_diff_cal,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_IntP_4vec_Q3",300,-0.1,60,intp,300,-0.1,80,ke);
	  MyFill("ResidualEn_vs_Needle_4vec_Q3",256,0,4096,ICne_E_diff,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_NeedleCal_4vec_Q3",256,0,70,ICne_E_diff_cal,300,-11,50,
                                                                            ResidE_recoil);
	  MyFill("ResidualEn_vs_RecoilEn_4vec_Q3",300,-0.1,80,ke,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_Theta_4vec_Q3",300,-0.1,180,rtheta,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_Path_4vec_Q3",300,-0.1,60,r_fullpath,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_IntP_4vec_Q3",300,-0.1,60,intp,300,-11,50,ResidE_recoil);
	} else if (detid>3 && detid<16) {
	  MyFill("LightParEn_4VEC_vs_BeamEnergy_Tracked_4vec_R1",300,-0.1,80,be,300,-0.1,30,etot);
	  MyFill("RecoilEn_4VEC_vs_BeamEnergy_Tracked_4vec_R1",300,-0.1,80,be,300,-0.1,80,ke);
	  MyFill("SiEnergy_tot_vs_Needle_4vec_R1",256,0,4096,ICne_E_diff,300,-0.1,30,etot);
	  MyFill("SiEnergy_tot_vs_NeedleCal_4vec_R1",256,0,70,ICne_E_diff_cal,300,-0.1,30,etot);
	  MyFill("SiEnergy_tot_vs_RecoilEn_4vec_R1",300,-0.1,80,ke,300,-0.1,30,etot);
	  MyFill("RecoilEnergy_vs_Theta_4vec_R1",300,-0.1,180,rtheta,300,-0.1,30,ke);
	  MyFill("RecoilEnergy_vs_Path_4vec_R1",300,-0.1,60,r_fullpath,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_Needle_4vec_R1",256,0,4096,ICne_E_diff,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_NeedleCal_4vec_R1",256,0,70,ICne_E_diff_cal,300,-0.1,80,ke);
	  MyFill("RecoilEnergy_vs_IntP_4vec_R1",300,-0.1,60,intp,300,-0.1,80,ke);
	  MyFill("ResidualEn_vs_Needle_4vec_R1",256,0,4096,ICne_E_diff,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_NeedleCal_4vec_R1",256,0,70,ICne_E_diff_cal,300,-11,50,
                                                                            ResidE_recoil);
	  MyFill("ResidualEn_vs_RecoilEn_4vec_R1",300,-0.1,80,ke,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_Theta_4vec_R1",300,-0.1,180,rtheta,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_Path_4vec_R1",300,-0.1,60,r_fullpath,300,-11,50,ResidE_recoil);
	  MyFill("ResidualEn_vs_IntP_4vec_R1",300,-0.1,60,intp,300,-11,50,ResidE_recoil);
        }
      }
    }
  }
  Elastic2.ELoss_light = NULL;
  Elastic2.ELoss_beam = NULL;
}

/*ReconstructMe()
 *This is the main reconstruction method; reconstructs properties of the beam,
 *interaction point, and recoil using both kinematics and tracking information
 *Utilizes ReconstructHeavy, Reconstruct20Ne, etc. This is an absolute monster of a code
 *and will take the most time of the bunch to run.
 *Have fun
 */
void analyzer::ReconstructMe() {
  Reconstruct Elastic1(M_18Ne,M_alpha,M_p,M_21Na);
  Elastic1.ELoss_light = proton_eloss;
  Elastic1.ELoss_beam = Ne18_eloss;
  vector<Int_t> IsProton;
  vector<Int_t> IsAlpha;
  vector<Double_t> RecoilE;

  for(int i=0; i<tracks.NTracks1; i++) {
    if(tracks.TrEvent[i].TrackType !=1) {
      cout<<"gwmAnalyzer::ReconstructMe() error: Track is not type 1!"<<endl;
    }
    Int_t detid = tracks.TrEvent[i].DetID;
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t sir = tracks.TrEvent[i].SiR;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t pcphi = tracks.TrEvent[i].PCPhi;
    Double_t siphi = tracks.TrEvent[i].SiPhi;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t theta = tracks.TrEvent[i].Theta;
    Double_t resE = tracks.TrEvent[i].ResidualEn;
    if(detid>-1 && detid<28 && pcz>0.0 && siz>=0.0 && intp>0.0 && intp<ana_length && 
       be>0.0 && be<BeamE && sir>4.0 && sir<11.0) {
      if(protonCut->IsInside(sie, pce*sin(theta))) {
        IsProton.push_back(i);
        Elastic1.ReconstructHeavy(tracks, i, recoil);
        RecoilE.push_back(recoil.Ex);
        Double_t rex = recoil.Ex;
        Double_t retot = recoil.Energy_tot;
        Double_t rbe = recoil.BeamEnergy;

        if(rex>-1.0 && rex<1.0 && needleCut->IsInside(ICne_E_diff*0.0116, resE)) {
          Elastic1.ReconstructHeavy_Qvalue(QValue, tracks, i, recoil);
	  MyFill("Ex_21Na_Qvalue_Rec",600,-10,40,recoil.Ex_21Na);		   
	  MyFill("Ex_21Na_vs_Beam_Qvalue_Rec",600,-10,80,recoil.Beam_Qv_21Na,
                                              400,-10,20,recoil.Ex_21Na);
	  MyFill("BeamEnergy_Tracked_vs_Beam_Qv_21Na_Qvalue_Rec",600,-10,80,recoil.Beam_Qv_21Na,
                                                                 600,-10,80,be);
	  MyFill("ResidualEn_vs_NCal_Qvalue_Rec",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_Track_vs_NCal_Qvalue_Rec",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
	  MyFill("Beam_Qv_21Na_vs_NCal_Qvalue_Rec",600,0,40,ICne_E_diff*0.0116,
                                                   600,-10,80,recoil.Beam_Qv_21Na);
        } else if (rex<0.0) {
	  MyFill("PCZ_Ex<0",400,-10,40,pcz);
	  MyFill("IntPoint_Ex<0",400,-10,80,intp);
	  MyFill("BeamEnergy_Ex<0",400,-10,80,be);
	  MyFill("ProtonEnergy_vs_ThetaEx<0",400,0,190,theta*rads2deg,400,0,15,retot);
	  MyFill("BeamEnergy_vs_IntPoint_Ex<0",400,-10,80,intp,400,-10,80,be);
	  MyFill("ProtonEnergy_vs_BeamEnergy_Ex<0",400,-10,80,be,400,0,15,retot);
	  MyFill("ProtonEnergy_vs_SiR_Ex<0",400,-1,11,sir,400,0,15,retot);
	  MyFill("ProtonEnergy_vs_IntPoint_Ex<0",400,-10,80,intp,400,0,15,retot);
	  MyFill("ProtonEnergy_vs_NeedleEn_Ex<0",600,0,40,ICne_E_diff*0.0116,400,0,15,retot);
	  MyFill("IntPoint_vs_DetID_Ex<0",28,0,28,detid,400,-10,80,intp);
	  MyFill("IntPoint_vs_PCZ_Ex<0",400,-10,40,pcz,400,-10,80,intp);
	  MyFill("SiPhi_vs_PCPhi_Ex<0",400,-0,360,pcphi*rads2deg,400,-0,360,siphi*rads2deg);
	  MyFill("ResidualEn_vs_Needle_Ex<0",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_vs_Needle_Ex<0",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          if (needleCut->IsInside(ICne_E_diff*0.0116, resE)) {
	    MyFill("BeamEnergy_for_cross_section_Ex<0_Needlecut",300,0,55,rbe);
	    MyFill("BeamEnergy_atCM_cross_section_Ex<0Needlecut",300,0,10,rbe*4/22);
	    MyFill("ResidualEn_vs_Needle_Ex<0NeedleC",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	    MyFill("BeamEnergy_vs_Needle_Ex<0NeedleC",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          }
        }
        if(needleCut->IsInside(ICne_E_diff*0.0116, resE)) {
	  MyFill("BeamEnergy_for_cross_section_NeedleCutAllDet",300,0,55,rbe);
	  MyFill("BeamEnergy_atCM_cross_sectionNeedleCutAllDet",300,0,10,rbe*4/22);
	  MyFill("ExEnergy_vs_BeamEnergy_NeedleCutAllDet",600,-10,80,rbe,400,-10,40,rex);
	  MyFill("ResidualEn_vs_NCal_NeedleCutAllDet",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_vs_NCal_NeedleCutAllDet",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          if(detid>-1 && detid<4) {
	    MyFill("BeamEnergy_for_cross_section_NeedleCutQ3",300,0,55,rbe);
	    MyFill("BeamEnergy_atCM_cross_section_NeedleCutQ3",300,0,10,rbe*4/22);
	    MyFill("ExEnergy_vs_BeamEnergy_NeedleCutQ3",600,-10,80,rbe,400,-10,40,rex);
	    MyFill("ResidualEn_vs_NCal_NeedleCutQ3",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	    MyFill("BeamEnergy_vs_NCal_NeedleCutQ3",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          }else if (detid>3 && detid<16) {
	    MyFill("BeamEnergy_for_cross_section_NeedleCutR1",300,0,55,rbe);
	    MyFill("BeamEnergy_atCM_cross_sectionNeedleCutR1",300,0,10,rbe*4/22);
	    MyFill("ExEnergy_vs_BeamEnergy_NeedleCutR1",600,-10,80,rbe,400,-10,40,rex);
 	    MyFill("ResidualEn_vs_NCal_NeedleCutR1",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	    MyFill("BeamEnergy_vs_NCal_NeedleCutR1",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          }else if (detid>15 && detid<28) {
	    MyFill("BeamEnergy_for_cross_section_NeedleCutR2",300,0,55,rbe);
	    MyFill("BeamEnergy_atCM_cross_sectionNeedleCutR2",300,0,10,rbe*4/22);
	    MyFill("ExEnergy_vs_BeamEnergy_NeedleCutR2",600,-10,80,rbe,400,-10,40,rex);
	    MyFill("ResidualEn_vs_NCal_NeedleCutR2",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	    MyFill("BeamEnergy_vs_NCal_NeedleCutR2",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
          }
        }
	MyFill("Ex_of_HeavyRecoil",600,-10,40,rex);
	MyFill("BeamEnergy_apRec",600,-10,80,rbe);
	MyFill("BeamEnergy_apRec_cm",600,-1,15,rbe*4/22);
	MyFill("BeamEnergy_apRec_Tracked",600,-10,80,be);
	MyFill("ProtonEnergy_Tracked_protons",600,-10,20,retot);
	MyFill("ExEnergy_vs_BeamEnergy",600,-10,80,rbe,400,-10,20,rex);
	MyFill("ProtonEnergy_vs_ExEnergy",600,-10,60,rex,400,-10,20,retot);
	MyFill("ProtonEnergy_vs_Theta",600,0,180,theta*rads2deg,600,0,20,retot);
	MyFill("ProtonEnergy_vs_BeamEnergy",600,-10,80,rbe,600,0,20,retot);
	MyFill("BeamEnergy_Tr_vs_BeamQval_apRecHeavy",600,-10,80,rbe,600,-10,80,be);
	MyFill("ProtonEnergy_vs_NCalAllDet",600,0,40,ICne_E_diff*0.0116,600,-10,20,retot);
	MyFill("ResidualEn_vs_NeedleCalAllDet",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	MyFill("BeamEnergy_vs_NeedleCalAllDet",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
        if(detid>-1 && detid<4) {
	  MyFill("Ex_of_HeavyRecoilQ3",600,-10,40,rex);
	  MyFill("BeamEnergy_apRecQ3",600,-10,80,rbe);
	  MyFill("BeamEnergy_apRecQ3_cm",600,-1,15,rbe*4/22);
	  MyFill("BeamEnergy_apRec_TrackedQ3",600,-10,80,be);
	  MyFill("ProtonEnergy_Tracked_protonsQ3",600,-10,20,retot);
	  MyFill("ExEnergy_vs_BeamEnergyQ3",600,-10,80,rbe,400,-10,20,rex);
	  MyFill("ProtonEnergy_vs_ExEnergyQ3",600,-10,60,rex,400,-10,20,retot);
	  MyFill("ProtonEnergy_vs_ThetaQ3",600,0,180,theta*rads2deg,600,0,20,retot);
	  MyFill("ProtonEnergy_vs_BeamEnergyQ3",600,-10,80,rbe,600,0,20,retot);
	  MyFill("BeamEnergy_Track_vs_BeamQval_apRecHeavyQ3",600,-10,80,rbe,600,-10,80,be);
	  MyFill("ProtonEnergy_vs_NCalQ3",600,0,40,ICne_E_diff*0.0116,600,-10,20,retot);
	  MyFill("ResidualEn_vs_NCalQ3",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_vs_NCalQ3",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
        } else if (detid>3 && detid <16) {
	  MyFill("Ex_of_HeavyRecoilR1",600,-10,40,rex);
	  MyFill("BeamEnergy_apRecR1",600,-10,80,rbe);
	  MyFill("BeamEnergy_apRecR1_cm",600,-1,15,rbe*4/22);
	  MyFill("BeamEnergy_apRec_TrackR1",600,-10,80,rbe);
	  MyFill("ProtonEnergy_Track_protonsR1",600,-10,20,retot);
	  MyFill("ExEnergy_vs_BeamEnergyR1",600,-10,80,rbe,400,-10,20,rex);
	  MyFill("ProtonEnergy_vs_ExEnergyR1",600,-10,60,rex,400,-10,20,retot);
	  MyFill("ProtonEnergy_vs_ThetaR1",600,0,180,theta*rads2deg,600,0,20,retot);
	  MyFill("ProtonEnergy_vs_BeamEnergyR1",600,-10,80,rbe,600,0,20,retot);
	  MyFill("BeamEnergy_Track_vs_BeamQval_apRecHeavyR1",600,-10,80,rbe,600,-10,80,be);
	  MyFill("ProtonEnergy_vs_NCalR1",600,0,40,ICne_E_diff*0.0116,600,-10,20,retot);
	  MyFill("ResidualEn_vs_NCalR1",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_vs_NCalR1",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
        } else if (detid>15 && detid<28) {
	  MyFill("Ex_of_HeavyRecoilR2",600,-10,40,rex);
	  MyFill("BeamEnergy_apRecR2",600,-10,80,rbe);
	  MyFill("BeamEnergy_apRecR2_cm",600,-1,15,rbe*4/22);
	  MyFill("BeamEnergy_apRec_TrackR2",600,-10,80,be);
	  MyFill("ProtonEnergy_Track_protonsR2",600,-10,20,retot);
	  MyFill("ExEnergy_vs_BeamEnergyR2",600,-10,80,rbe,400,-10,20,rex);
	  MyFill("ProtonEnergy_vs_ExEnergyR2",600,-10,60,rex,400,-10,20,retot);
	  MyFill("ProtonEnergy_vs_ThetaR2",600,0,180,theta*rads2deg,600,0,20,retot);
	  MyFill("ProtonEnergy_vs_BeamEnergyR2",600,-10,80,rbe,600,0,20,retot);
	  MyFill("BeamEnergy_Track_vs_BeamQval_apRecHeavyR2",600,-10,80,rbe,600,-10,80,be);
	  MyFill("ProtonEnergy_vs_NCalR2",600,0,40,ICne_E_diff*0.0116,600,-10,20,retot);
	  MyFill("ResidualEn_vs_NCalR2",600,0,40,ICne_E_diff*0.0116,600,0,40,resE);
	  MyFill("BeamEnergy_vs_NCalR2",600,0,40,ICne_E_diff*0.0116,600,-10,80,be);
        }
      }//if protonCut
      if (IsProton.size() == 1) { //This is a check case; REMOVE POST TESTING
	MyFill("SiProtonE_de_p0",400,-1,35,tracks.TrEvent[IsProton[0]].SiEnergy,
               400,-0.01,0.8,tracks.TrEvent[IsProton[0]].PCEnergy*sin(tracks.TrEvent[IsProton[0]].Theta));
      }
      if(IsProton.size()==1 && tracks.NTracks1==1) {//only single proton events
        Double_t p1theta = tracks.TrEvent[IsProton[0]].Theta;
        Double_t p1sie = tracks.TrEvent[IsProton[0]].SiEnergy;
        Double_t p1be = tracks.TrEvent[IsProton[0]].BeamEnergy;
        Double_t p1resE = tracks.TrEvent[IsProton[0]].ResidualEn;
        Double_t p1intp = tracks.TrEvent[IsProton[0]].IntPoint;
        Int_t p1detid = tracks.TrEvent[IsProton[0]].DetID;
        Double_t rex = RecoilE[0];

	MyFill("SiProtonE_theta_p0_ntracks1_1",400,-1,190,p1theta*rads2deg,400,-1,15,p1sie);
	MyFill("SiProtonE_vs_IntP_p0_ntracks1_1",400,-1,70,p1intp,400,-1,15,p1sie);
	MyFill("SiProtonE_vs_BeamEnergy_p0_ntracks1_1",400,-1,80,p1be,400,-1,15,p1sie);
	MyFill("ExEnergy_vs_BeamEnergy_p0_ntracks1_1",600,-10,80,p1be,400,-10,20,rex);
	MyFill("ExEnergy_vs_BeamE_cm_p0_ntracks1_1",600,-10,20,p1be*4/22,400,-10,20,rex);
	MyFill("ExcitationEnergy_p0_ntracks1_1",400,-10,20,rex);
	MyFill("SiProtonEnergy_vs_N_p0_ntracks1_1",600,0,40,ICne_E_diff*0.0116,400,-10,20,p1sie);
	MyFill("ResidualEn_vs_N_p0_ntracks1_1",600,0,40,ICne_E_diff*0.0116,600,0,40,p1resE);
	MyFill("BeamEnergy_vs_N_p0_ntrack1_1",600,0,40,ICne_E_diff*0.0116,600,-10,80,p1be);
	MyFill("BeamE_cm_p0_ntracks1_1",180,-1,15,be*4/22);
        if(needleCut->IsInside(ICne_E_diff*0.0116, p1resE)) {//what is this 0.0116 factor?
	  MyFill("BeamEnergy_for_cross_section_p0_ntracks1_1_Ncut",300,0,55,p1be);
	  MyFill("BeamEnergy_at_CM_for_cross_section_p0_ntracks1_1_Ncut",300,0,10,p1be*4/22);
	  MyFill("ExEnergy_vs_BeamEnergy_p0_ntracks1_1_Ncut",600,-10,80,p1be,400,-10,20,rex);
          MyFill("ExEnergy_vs_BeamE_cm_p0_ntracks1_1_Ncut",600,-10,20,p1be*4/22,400,-10,20,rex);
	  MyFill("ExcitationEnergy_p0_ntracks1_1_Ncut",400,-10,20,rex);
          if(p1detid>-1 && p1detid<4) {
	    MyFill("BeamEnergy_for_cross_section_p0_ntracks1_1_NcutQ3",300,0,55,p1be);
	    MyFill("BeamEnergy_atCM_cross_section_p0_ntracks1_1_NcutQ3",300,0,10,p1be*4/22); 
	    MyFill("ExEnergy_vs_BeamEnergy_p0_ntracks1_1_NcutQ3",600,-10,80,p1be,400,-10,20,rex);
	    MyFill("ExEnergy_vs_BeamEcm_p0ntracks1_1_NcutQ3",600,-10,20,p1be*4/22,400,-10,20,rex);
	    MyFill("ExEnergy_p0_ntracks1_1_NcutQ3",400,-10,20,rex);
          } else if (p1detid>3 && p1detid<16) {
	    MyFill("BeamEnergy_for_cross_section_p0_ntracks1_1_NcutR1",300,0,55,p1be);
	    MyFill("BeamEnergy_at_CM_for_cross_section_p0_ntracks1_1_NcutR1",300,0,10,p1be*4/22); 
	    MyFill("ExEnergy_vs_BeamEnergy_p0_ntracks1_1_NcutR1",600,-10,80,p1be,400,-10,20,rex);
	    MyFill("ExEnergy_vs_BeamEcm_p0ntracks1_1_NcutR1",600,-10,20,p1be*4/22,400,-10,20,rex);
	    MyFill("ExEnergy_p0_ntracks1_1_NcutR1",400,-10,20,rex);
          } else if (p1detid>15 && p1detid<27) {
	    MyFill("BeamEnergy_for_cross_section_p0_ntracks1_1_NcutR2",300,0,55,p1be);
	    MyFill("BeamEnergy_at_CM_for_cross_section_p0_ntracks1_1_NcutR2",300,0,10,p1be*4/22); 
	    MyFill("EEnergy_vs_BeamEnergy_p0_ntracks1_1_NcutR2",600,-10,80,p1be,400,-10,20,rex);
	    MyFill("ExEnergy_vs_BeamEcm_p0ntracks1_1_NcutR2",600,-10,20,p1be*4/22,400,-10,20,rex);
	    MyFill("ExcitationEnergy_p0_ntracks1_1_NcutR2",400,-10,20,rex);
          }
        }
      } else if (IsProton.size()==2 && tracks.NTracks1==2) {//only 2 proton events
        Double_t p1theta = tracks.TrEvent[IsProton[0]].Theta;
        Double_t p1siphi = tracks.TrEvent[IsProton[0]].SiPhi;
        Double_t p1pcphi = tracks.TrEvent[IsProton[0]].PCPhi;
        Double_t p1pcz = tracks.TrEvent[IsProton[0]].PCZ;
        Double_t p1sie = tracks.TrEvent[IsProton[0]].SiEnergy;
        Double_t p1pce = tracks.TrEvent[IsProton[0]].PCEnergy;
        Double_t p1be = tracks.TrEvent[IsProton[0]].BeamEnergy;
        Double_t p1resE = tracks.TrEvent[IsProton[0]].ResidualEn;
        Double_t p1intp = tracks.TrEvent[IsProton[0]].IntPoint;
        Int_t p1detid = tracks.TrEvent[IsProton[0]].DetID;
        Double_t p2theta = tracks.TrEvent[IsProton[1]].Theta;
        Double_t p2siphi = tracks.TrEvent[IsProton[1]].SiPhi;
        Double_t p2pcphi = tracks.TrEvent[IsProton[1]].PCPhi;
        Double_t p2pcz = tracks.TrEvent[IsProton[1]].PCZ;
        Double_t p2sie = tracks.TrEvent[IsProton[1]].SiEnergy;
        Double_t p2pce = tracks.TrEvent[IsProton[1]].PCEnergy;
        Double_t p2be = tracks.TrEvent[IsProton[1]].BeamEnergy;
        Double_t p2resE = tracks.TrEvent[IsProton[1]].ResidualEn;
        Double_t p2intp = tracks.TrEvent[IsProton[1]].IntPoint;
        Int_t p2detid = tracks.TrEvent[IsProton[1]].DetID;
        Double_t rex1 = RecoilE[0];
        Double_t rex2 = RecoilE[1];

	MyFill("IntPoint_p1_vs_IntPoint_p2",400,-10,70,p1intp,400,-10,70,p2intp);
	MyFill("SiProtonE_de_p1",400,-1,35,p1sie,400,-0.01,0.8,p1pce*sin(p1theta));
	MyFill("SiProtonE_theta_p1",400,-1,190,p1theta*rads2deg,400,-1,15,p1sie);
	MyFill("SiProtonE_de_p2",400,-1,35,p2sie,400,-0.01,0.8,p2pce*sin(p2theta));
	MyFill("SiProtonE_theta_p2",400,-1,190,p2theta*rads2deg,400,-1,15,p2sie);
	MyFill("Theta_p1_vs_Theta_p2",400,-1,190,p1theta*rads2deg,400,-1,190,p2theta*rads2deg);
	MyFill("SiPhi_p1_vs_SiPhi_p2",400,-1,360,p1siphi*rads2deg,400,-1,360,p2siphi*rads2deg);
	MyFill("SiProtonE_p1_vs_SiProtonE_p2",400,-1,20,p1sie,400,-1,20,p2sie);
	MyFill("BeamE_p1_vs_BeamE_p2",400,-1,80,p1be,400,-1,80,p2be);
	MyFill("PCZ_p1_vs_PCZ_p2",600,-10,60,p1pcz,600,-1,60,p2pcz);
	MyFill("Ex_p1_vs_Ex_p2",400,-10,20,rex1,400,-10,20,rex2);
        if(needleCut->IsInside(ICne_E_diff*0.0116,p1resE) || 
           needleCut->IsInside(ICne_E_diff*0.0116,p2resE)) {
	  MyFill("IntPoint_p1_vs_IntPoint_p2_NeedleCut",400,-10,70,p1intp,400,-10,70,p2intp);
	  MyFill("Theta1_vs_Theta2_NCut",400,-1,190,p1theta*rads2deg,400,-1,190,p2theta*rads2deg);
	  MyFill("SiPhi1_vs_SiPhi2_NCut",400,-1,360,p1siphi*rads2deg,400,-1,360,p2siphi*rads2deg);
	  MyFill("SiProtonE_p1_vs_SiProtonE_p2_NCut",400,-1,20,p1sie,400,-1,20,p2sie);
	  MyFill("BeamE_p1_vs_BeamE_p2_NeedleCut",400,-1,80,p1be,400,-1,80,p2be);
	  MyFill("PCZ_p1_vs_PCZ_p2_NeedleCut",600,-10,60,p1pcz,600,-1,60,p2pcz);
	  MyFill("Ex_p1_vs_Ex_p2_NeedleCut",400,-10,20,rex1,400,-10,20,rex2);
        }
	MyFill("ExEnergy_vs_BeamEnergy_p1_p2_forP0",600,-10,80,p1be,400,-10,20,rex1);
	MyFill("ExEnergy_vs_BeamEnergy_p1_p2_forP1",600,-10,80,p2be,400,-10,20,rex2);
	MyFill("ExEnergy_vs_BeamEnergy_p1_p2",600,-10,80,p1be,400,-10,20,rex1);//dupe; RM POST
	MyFill("ExEnergy_vs_BeamEnergy_p1_p2",600,-10,80,p2be,400,-10,20,rex2);//dupe; RM POST
	MyFill("ExEnergy_p1_p2_forP0",400,-10,20,rex1);
	MyFill("ExEnergy_p1_p2_forP1",400,-10,20,rex2);
	MyFill("ExEnergy_p1_p2",400,-10,20,rex1);//dupe; REMOVE POST TEST
	MyFill("ExEnergy_p1_p2",400,-10,20,rex2);//dupe; REMOVE POST TEST
	MyFill("SiProtonEnergy_vs_N_p1_p2_forP0",600,0,40,ICne_E_diff*0.0116,400,-10,20,p1sie);
	MyFill("SiProtonEnergy_vs_N_p1_p2_forP1",600,0,40,ICne_E_diff*0.0116,400,-10,20,p2sie);
	MyFill("SiProtonEnergy_vs_N_p1_p2",600,0,40,ICne_E_diff*0.0116,400,-10,20,p1sie);//dupe
	MyFill("SiProtonEnergy_vs_N_p1_p2",600,0,40,ICne_E_diff*0.0116,400,-10,20,p2sie);//dupe
	MyFill("BeamE_cm_p1_p2_forP0",180,-1,15,p1be*4/22);
	MyFill("BeamE_cm_p1_p2_forP1",180,-1,15,p2be*4/22);
	MyFill("BeamE_cm_p1_p2",180,-1,15,p1be*4/22);//dupe
	MyFill("BeamE_cm_p1_p2",180,-1,15,p2be*4/22);//dupe; Factor of 4/22?
	MyFill("ResidualE_vs_N_p1_p2_forP0",600,0,40,ICne_E_diff*0.0116,600,0,40,p1resE);
	MyFill("ResidualE_vs_N_p1_p2_forP1",600,0,40,ICne_E_diff*0.0116,600,0,40,p2resE);
	MyFill("ResidualE_vs_N_p1_p2",600,0,40,ICne_E_diff*0.0116,600,0,40,p1resE);
	MyFill("ResidualE_vs_N_p1_p2",600,0,40,ICne_E_diff*0.0116,600,0,40,p2resE);
	MyFill("ResidualE_vs_N_p1_p2_largelim",600,-11,100,ICne_E_diff*0.0116,600,-11,60,p1resE);
	MyFill("ResidualE_vs_N_p1_p2_largelim",600,-11,100,ICne_E_diff*0.0116,600,-11,60,p2resE);

        //Why do these exist? Never used ever again; terminal plotting tools
        tracks.TrEvent[IsProton[0]].BeamE_p1 = tracks.TrEvent[IsProton[0]].BeamEnergy;
        tracks.TrEvent[IsProton[1]].BeamE_p2 = tracks.TrEvent[IsProton[1]].BeamEnergy;
        tracks.TrEvent[IsProton[0]].SiE_p1 = tracks.TrEvent[IsProton[0]].SiEnergy;
        tracks.TrEvent[IsProton[1]].SiE_p2 = tracks.TrEvent[IsProton[1]].SiEnergy;
        tracks.TrEvent[IsProton[0]].ResE_p1 = tracks.TrEvent[IsProton[0]].ResidualEn;
        tracks.TrEvent[IsProton[1]].ResE_p2 = tracks.TrEvent[IsProton[1]].ResidualEn;
        tracks.TrEvent[IsProton[0]].Ex_p1 = RecoilE[0];
        tracks.TrEvent[IsProton[1]].Ex_p2 = RecoilE[1];
        tracks.TrEvent[IsProton[0]].IntP_p1 = tracks.TrEvent[IsProton[0]].IntPoint;
        tracks.TrEvent[IsProton[1]].IntP_p2 = tracks.TrEvent[IsProton[1]].IntPoint;
        tracks.TrEvent[IsProton[0]].Theta_p1 = tracks.TrEvent[IsProton[0]].Theta;
        tracks.TrEvent[IsProton[1]].Theta_p2 = tracks.TrEvent[IsProton[1]].Theta;
        tracks.TrEvent[IsProton[0]].SiPhi_p1 = tracks.TrEvent[IsProton[0]].SiPhi;
        tracks.TrEvent[IsProton[1]].SiPhi_p2 = tracks.TrEvent[IsProton[1]].SiPhi;
        tracks.TrEvent[IsProton[0]].PCZ_p1 = tracks.TrEvent[IsProton[0]].PCZ;
        tracks.TrEvent[IsProton[1]].PCZ_p2 = tracks.TrEvent[IsProton[1]].PCZ;
        Double_t deltaPhi = abs(p1siphi - p2siphi);
        if(deltaPhi>TMath::Pi() && deltaPhi<=2*TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
        Double_t deltaTheta = abs(p1theta-p2theta);

        Reconstruct Elastic3(M_18Ne, M_alpha, M_p, M_20Ne);
        Elastic3.ELoss_light = proton_eloss;
        Elastic3.ELoss_beam = Ne18_eloss;
        Elastic3.Reconstruct20Ne(tracks, IsProton, QValue_20Ne, 0, recoil);//no intp reconstruct
        Double_t rbe20Ne = recoil.BeamEnergy_20Ne;
        //dupes throughout
	MyFill("BeamEnergy_Track_vs_BeamQval_20Ne_RHeavy_P0",600,-10,80,rbe20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_Track_vs_BeamQval_20Ne_RHeavy_P1",600,-10,80,rbe20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_Track_vs_BeamQval_20Ne_RHeavy_all",600,-10,80,rbe20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_Track_vs_BeamQval_20Ne_RHeavy_all",600,-10,80,rbe20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_20Ne_RHeavy",600,-10,80,rbe20Ne);
	MyFill("BeamEnergy_20Ne_RHeavy_cm",180,-1,15,rbe20Ne*4/22);
	MyFill("BeamEnergy_20Ne_RHeavy_vs_NEnergy",128,0,4096,ICne_E_diff,600,-10,80,rbe20Ne);
	MyFill("ExEnergy_vs_BeamEnergy_20Ne_RHeavy_forP0",600,-10,80,rbe20Ne,400,-10,20,rex1);
	MyFill("ExEnergy_vs_BeamEnergy_20Ne_RHeavy_forP1",600,-10,80,rbe20Ne,400,-10,20,rex2);
	MyFill("ExEnergy_vs_BeamEnergy_20Ne_RHeavy_p1_p2",600,-10,80,rbe20Ne,400,-10,20,rex1);
	MyFill("ExEnergy_vs_BeamEnergy_20Ne_RHeavy_p1_p2",600,-10,80,rbe20Ne,400,-10,20,rex2);
        Elastic3.Reconstruct20Ne(tracks, IsProton, QValue_20Ne, 1, recoil);//intp reconstruct
        Double_t exQval20Ne = recoil.Ex_Qvalue_20Ne;
        Double_t exrQvalsm = recoil.Ex_rec_qvalue_small;
        Double_t exrQvallg = recoil.Ex_rec_qvalue_large;
        Double_t bmWAQval20Ne = recoil.BeamWA_Qvalue_20Ne;
        Double_t ex20Ne = recoil.Ex_20Ne;
        Double_t exrlg = recoil.Ex_rec_large;
        Double_t exrsm = recoil.Ex_rec_small;
        Double_t bmWA20Ne = recoil.BeamWA_20Ne;

	MyFill("ExEnergy_20Ne_QvalRec",600,-20,50,exQval20Ne);
	MyFill("21Na_Ex_lg_vs_21Na_Ex_sm_QvalRec",600,-10,20,exrQvalsm,600,-10,20,exrQvallg);
	MyFill("21Na_Ex_lg_vs_BeamEnergy_QvalRec",600,-10,80,bmWAQval20Ne,600,-10,20,exrQvallg);
	MyFill("21Na_Ex_sm_vs_BeamEnergy_QvalRec",600,-10,80,bmWAQval20Ne,600,-10,20,exrQvalsm);
	MyFill("21Na_Ex_large_QvalueRec",600,-10,20,exrQvallg);
	MyFill("21Na_Ex_small_QvalueRec",600,-10,20,exrQvalsm);
	MyFill("ExEnergy_20Ne_vs_BeamWA_QvalRec",600,-10,80,bmWAQval20Ne,600,-20,50,exQval20Ne);
	MyFill("BeamEnergy_Track_vs_BeamWA_Qval_20Ne_P0",600,-10,80,bmWAQval20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_Track_vs_BeamWA_Qval_20Ne_P1",600,-10,80,bmWAQval20Ne,600,-10,80,p2be);
	MyFill("BeamEnergy_Track_vs_BeamWA_Qval_20Ne",600,-10,80,bmWAQval20Ne,600,-10,80,p1be);
	MyFill("BeamEnergy_Track_vs_BeamWA_Qval_20Ne",600,-10,80,bmWAQval20Ne,600,-10,80,p2be);
	MyFill("BeamWA_QvalueRec_20Ne",600,-10,80,bmWAQval20Ne);
	MyFill("BeamWA_QvalueRec_20Ne_cm",180,-1,15,bmWAQval20Ne*4/22);
    
        //test spectra; REMOVE POST TEST? Yes
	if((p1detid>-1 && p1detid<4) && (p2detid>-1 && p2detid<4)) {
	  MyFill("ExEnergy_20Ne_Q3_testAnd",600,-20,50,ex20Ne);
        }
	if((p1detid>-1 || p1detid<4) || (p2detid>-1 || p2detid<4)) {
	  MyFill("ExEnergy_20Ne_Q3_testOr",600,-20,50,ex20Ne);
        }
	if(p1detid>-1 && p1detid<4) {
	  MyFill("ExEnergy_20Ne_Q3_testP0",600,-20,50,ex20Ne);
        }
	if(p2detid>-1 && p2detid<4) {
	  MyFill("ExEnergy_20Ne_Q3_testP1",600,-20,50,ex20Ne);
        }
	if(needleCut->IsInside(ICne_E_diff*0.0116,p1resE)) {
	   MyFill("ExEnergy_20Ne_NCut_cutP0",600,-20,50,ex20Ne);
        }
	if(needleCut->IsInside(ICne_E_diff*0.0116,p2resE)) {
	   MyFill("ExEnergy_20Ne_NCut_cutP1",600,-20,50,ex20Ne);
        }
	if(needleCut->IsInside(ICne_E_diff*0.0116,p1resE)||
           needleCut->IsInside(ICne_E_diff*0.0116,p2resE)) {
	  MyFill("ExEnergy_20Ne_NCut_cutOR",600,-20,50,ex20Ne);
        }
	if(needleCut->IsInside(ICne_E_diff*0.0116,p1resE)&&
           needleCut->IsInside(ICne_E_diff*0.0116,p2resE)) {
	  MyFill("ExEnergy_20Ne_NCut_cutAND",600,-20,50,ex20Ne);
        }
	MyFill("21Na_Ex_lg_vs_21Na_Ex_sm",600,-10,20,exrsm,600,-10,20,exrlg);
	MyFill("21Na_Ex_lg_vs_BeamEnergy",600,-10,80,bmWA20Ne,600,-10,20,exrlg);
	MyFill("21Na_Ex_sm_vs_BeamEnergy",600,-10,80,bmWA20Ne,600,-10,20,exrsm);
	MyFill("21Na_Ex_large",600,-10,20,exrlg);
	MyFill("21Na_Ex_small",600,-10,20,exrsm);
	MyFill("ExEnergy_20Ne_vs_BeamWA",600,-10,80,bmWA20Ne,600,-20,50,ex20Ne);
	MyFill("ExEnergy_20Ne_vs_BeamEnergy_P0",600,-10,80,p1be,600,-20,50,ex20Ne);
	MyFill("ExEnergy_20Ne_vs_BeamEnergy_P1",600,-10,80,p2be,600,-20,50,ex20Ne);
	MyFill("ExEnergy_20Ne_vs_BeamEnergy_all",600,-10,80,p1be,600,-20,50,ex20Ne);//dupe
	MyFill("ExEnergy_20Ne_vs_BeamEnergy_all",600,-10,80,p2be,600,-20,50,ex20Ne); //dupe
	MyFill("ExEnergy_20Ne_vs_DeltaSiPhi",400,-10,360,deltaPhi*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_DeltaTheta",400,-10,190,deltaTheta*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_PCZ_P0_all",400,-10,30,p1pcz,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_PCZ_P1_all",400,-10,30,p2pcz,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_PCZ_all",400,-10,30,p1pcz,400,-10,30,ex20Ne);//dupe 
	MyFill("ExEnergy_20Ne_vs_PCZ_all",400,-10,30,p2pcz,400,-10,30,ex20Ne); //dupe
	MyFill("ExEnergy_20Ne_vs_SiPhi_P0",400,-10,360,p1siphi*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_SiPhi_P1",400,-10,360,p2siphi*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_SiPhi_all",400,-10,360,p1siphi*rads2deg,400,-10,30,ex20Ne);//d 
	MyFill("ExEnergy_20Ne_vs_SiPhi_all",400,-10,360,p2siphi*rads2deg,400,-10,30,ex20Ne);//d 
	MyFill("ExEnergy_20Ne_vs_PCPhi_P0",400,-10,360,p1pcphi*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_PCPhi_P1",400,-10,360,p2pcphi*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_PCPhi_all",400,-10,360,p1pcphi*rads2deg,400,-10,30,ex20Ne);//d 
	MyFill("ExEnergy_20Ne_vs_PCPhi_all",400,-10,360,p2pcphi*rads2deg,400,-10,30,ex20Ne);//d 
	MyFill("ExEnergy_20Ne_vs_Theta_P0",400,-10,180,p1theta*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_Theta_P1",400,-10,180,p2theta*rads2deg,400,-10,30,ex20Ne); 
	MyFill("ExEnergy_20Ne_vs_Theta_all",400,-10,180,p1theta*rads2deg,400,-10,30,ex20Ne);//d 
	MyFill("ExEnergy_20Ne_vs_Theta_all",400,-10,180,p2theta*rads2deg,400,-10,30,ex20Ne);//d
        if(ResidualEnergy_Calc_20Ne) { 
          Float_t recoilFullpath = 0.0;
          Float_t recoilIP2Npath = 0.0;
          Float_t beam_needle = 0.0;
          Float_t deltaEbeam = 0.0; 
          Float_t Efinal = 0.0;
          Float_t E_fullpath = 0.0;
          Float_t E_needleStart = 0.0;
          Float_t resErecoil20Ne = 0.0; 
          Double_t ke20Ne = recoil.KE_20Ne;
          Double_t theta20Ne = recoil.Theta_20Ne;
          Double_t intp20Ne = recoil.IntP_20Ne;
	  MyFill("RecoilEnergy_vs_Theta_20Ne",300,-0.1,180,theta20Ne,300,-0.1,80,ke20Ne);
	  MyFill("RecoilEnergy_vs_NCal_20Ne",256,0,70,ICne_E_diff_cal,300,-0.1,80,ke20Ne);
	  if(theta20Ne != 0.0 && theta20Ne != 90.0 && theta20Ne != 180.0){
	    recoilFullpath = 2.2/sin(theta20Ne/rads2deg);
	    if(intp20Ne > 17.8) recoilIP2Npath = (intp20Ne - 17.8)/cos(theta20Ne/rads2deg);
	  }
	  if(intp20Ne <= 17.8) {
	    E_fullpath = Ne20_eloss->GetLookupEnergy(ke20Ne,recoilFullpath);
	    beam_needle = Ne18_eloss->GetLookupEnergy(BeamE,(ana_length-17.8));
	    deltaEbeam = beam_needle - Ne18_eloss->GetLookupEnergy(BeamE,(ana_length-intp20Ne));
	    if(E_fullpath <= 0.0){ 
	      resErecoil20Ne = ke20Ne + deltaEbeam;
	      if(p1detid>-1 && p1detid<4) {
	        MyFill("ResE_vs_N_20Ne_Q3insideVolFullLoss",256,0,4096,ICne_E_diff,
                                                            300,-11,50,resErecoil20Ne);
              } else if(p1detid>3 && p1detid<16) {
	        MyFill("ResEn_vs_N_20Ne_R1insideVolFullLoss",256,0,4096,ICne_E_diff,
                                                             300,-11,50,resErecoil20Ne);
              }
	    } else{
	      resErecoil20Ne = (ke20Ne - E_fullpath) + deltaEbeam;
	      if(p1detid>-1 && p1detid<4) {
	        MyFill("ResE_vs_N_20Ne_Q3insideVolSomeLoss",256,0,4096,ICne_E_diff,
                                                            300,-11,50,resErecoil20Ne);
              } else if(p1detid>3 && p1detid<16) {
	        MyFill("ResE_vs_N_20Ne_R1insideVolSomeLoss",256,0,4096,ICne_E_diff,
                                                            300,-11,50,resErecoil20Ne);
	      }
	    }
	    MyFill("ResE_Recoil_vs_N_NVol_20Ne",256,0,4096,ICne_E_diff,300,-11,50,resErecoil20Ne);
	  } else{ 
	    E_needleStart = Ne20_eloss->GetLookupEnergy(ke20Ne,recoilIP2Npath);
	    if(E_needleStart <= 0.0) { 
	      resErecoil20Ne = -10.0;
	    } else if(recoilFullpath < recoilIP2Npath) { 
	      resErecoil20Ne = -5.0;
	    } else { 
	      Efinal=Ne20_eloss->GetLookupEnergy(E_needleStart,(recoilFullpath-recoilIP2Npath));
	      if(Efinal <= 0.0){ 
	        resErecoil20Ne = E_needleStart;
	      } else{ 
	        resErecoil20Ne = E_needleStart - Efinal;
	      }
	    } 
	    MyFill("ResE_Recoil_vs_N_OUTSIDENeedleVol20Ne",256,0,4096,ICne_E_diff,
                                                           300,-11,50,resErecoil20Ne);
	  } 
	  tracks.TrEvent[IsProton[0]].ResidualEn_20Ne = resErecoil20Ne;	     	     
	  MyFill("ResE_vs_N_20Ne",256,0,4096,ICne_E_diff,300,-11,50,resErecoil20Ne);
	  MyFill("ResE_vs_NCal_20Ne",600,0,40,ICne_E_diff_cal,600,0,40,resErecoil20Ne);
	  MyFill("ResE_vs_RecoilE_20Ne",300,-0.1,80,ke20Ne,300,-11,50,resErecoil20Ne);
	  MyFill("ResE_vs_Theta_20Ne",300,-0.1,180,theta20Ne,300,-11,50,resErecoil20Ne);
        }
        Elastic3.ELoss_light = NULL;
        Elastic3.ELoss_beam = NULL;
      }
    }
  }
  Elastic1.ELoss_light = NULL;
  Elastic1.ELoss_beam = NULL;
}

/*AlphaScatter()
 *Method for elastic scattering reconstruction
 *Need to check, but my guess is its for a calibration
 *or as a way to normalize for cross sections
 */
void analyzer::AlphaScatter() {

  for(int i=0; i<tracks.NTracks1; i++) {
    Int_t detid = tracks.TrEvent[i].DetID;
    Int_t wireid = tracks.TrEvent[i].DetID;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t theta = tracks.TrEvent[i].Theta;
    Double_t pathL = tracks.TrEvent[i].PathLength;
    Int_t wirePlus3 = (wireid+2)%24;
    Int_t wireMinus3 = (wireid-2)%24;
    if(PCGoodEnergy[wirePlus3]>0 || PCGoodEnergy[wireMinus3]>0) continue;
    if(alphaCut->IsInside(sie, pce*sin(theta))) {
      MyFill("4He_E_si",600,0,35,sie);
      MyFill("4He_E_si_vs_Theta",600,0,600,theta*rads2deg,600,0,35,sie);
      MyFill("4He_E_si_vs_IntPoint",600,0,60,intp,600,0,35,sie);
      if(detid<4 && detid>-1){
        MyFill("4He_E_si_Q3",600,0,35,sie);
        MyFill("4He_E_si_vs_Theta_Q3",600,0,600,theta*rads2deg,600,0,35,sie);
        MyFill("4He_E_si_vs_IntPoint_Q3",600,0,50,intp,600,0,35,sie);
      }else if(detid<16 && detid>3){
        MyFill("4He_E_si_SX3_1",600,0,35,sie);
        MyFill("4He_E_si_vs_Theta_SX3_1",600,0,600,theta*rads2deg,600,0,35,sie);
        MyFill("4He_E_si_vs_IntPoint_SX3_1",600,0,60,intp,600,0,35,sie);
      }else if(detid<28 && detid>15){
        MyFill("4He_E_si_SX3_2",600,0,35,sie);
        MyFill("4He_E_si_vs_Theta_SX3_2",600,0,600,theta*rads2deg,600,0,35,sie);
        MyFill("4He_E_si_vs_IntPoint_SX3_2",600,0,60,intp,600,0,35,sie);
      }
      if (sie>0.0 && sie<40.0 && pathL>0.0 && pathL>100.0 && intp>0.0 && intp<ana_length && 
          be>0.0 && be<BeamE && pcz>0.0 && pce>0.0) {
        Float_t E_4Herxn = alpha_eloss->GetLookupEnergy(sie, -pathL);
        tracks.TrEvent[i].LightParEnergy = E_4Herxn;
	MyFill("4He_E_rxn",600,0,35,E_4Herxn);
	MyFill("4He_E_rxn_vs_Theta",600,0,600,theta*rads2deg,600,0,30,E_4Herxn);
	MyFill("4He_E_rxn_vs_IntPoint",600,0,50,intp,600,0,30,E_4Herxn);
	if(detid<4 && detid>-1){
	  MyFill("4He_E_rxn_Q3",600,0,30,E_4Herxn);
	  MyFill("4He_E_rxn_vs_Theta_Q3",600,0,80,theta*rads2deg,600,0,30,E_4Herxn);
	  MyFill("4He_E_rxn_vs_IntPoint_Q3",600,0,50,intp,600,0,30,E_4Herxn);
	}else if(detid<16 && detid>3){
	  MyFill("4He_E_rxn_SX3_1",600,0,30,E_4Herxn);
	  MyFill("4He_E_rxn_vs_Theta_SX3_1",600,0,80,theta*rads2deg,600,0,30,E_4Herxn);
	  MyFill("4He_E_rxn_vs_IntPoint_SX3_1",600,0,50,intp,600,0,30,E_4Herxn);
	}else if(detid<28 && detid>15){
	  MyFill("4He_E_rxn_SX3_2",600,0,30,E_4Herxn);
	  MyFill("4He_E_rxn_vs_Theta_SX3_2",600,0,600,theta*rads2deg,600,0,30,E_4Herxn);	
	  MyFill("4He_E_rxn_vs_IntPoint_SX3_2",600,0,50,intp,600,0,30,E_4Herxn);
	}

	//Beam Energy Calculation from the Elastic scattering: 
	Float_t EBeam4He = ((M_18Ne+M_alpha)*(M_18Ne+M_alpha)*E_4Herxn)/
                          (4.0*M_18Ne*M_alpha*cos(theta)*cos(theta));
	tracks.TrEvent[i].BeamQvalue = EBeam4He;
	//Beam Theta Calculation from the Elastic scattering:
	Float_t ThetaB4He = asin((sin(theta)*sqrt(M_alpha/M_18Ne))/
                                sqrt((EBeam4He/E_4Herxn)-1.0))*rads2deg;
	tracks.TrEvent[i].ThetaQvalue = ThetaB4He;
	//4He energy calculated from elastic scattering assuming known SRIM beam energy
	Float_t E4He_elas = (4.0*M_alpha*M_18Ne*be*cos(theta)*cos(theta))/
                                     ((M_18Ne+M_alpha)*(M_18Ne+M_alpha));
	tracks.TrEvent[i].HeEnergyQvalue =  E4He_elas;
	// Recoil Kinematics
	Float_t ERecoil = EBeam4He - E_4Herxn;
	MyFill("4He_Elastic",600,-10,40,E4He_elas);
	if(detid<4 && detid>-1) MyFill("4He_Elastic_q3",600,-10,40,E4He_elas);
	else if(detid<16 && detid>3) MyFill("4He_Elastic_r1",600,-10,40,E4He_elas);
	else if(detid<28 && detid>15) MyFill("4He_Elastic_r2",600,-10,40,E4He_elas);

	MyFill("LightParrxn_vs_Elastic_4He_E",600,0,30,E4He_elas,600,0,30,E_4Herxn);
	if(detid<4 && detid>-1){
	  MyFill("LightParrxn_vs_Elastic4He_EQ3",600,0,30,E4He_elas,600,0,30,E_4Herxn);
	  MyFill("RecE_v_(BeamQval-LightParE)Q3",600,-5,30,(EBeam4He-E_4Herxn),600,-5,30,ERecoil);
	  MyFill("RecE_vs_(BeamTrack - LightParE)Q3",600,-5,30,(be - E_4Herxn),600,-5,30,ERecoil);
	  MyFill("RecE_vs_BeamTrack_ElasticScatQ3",600,-5,80,be,600,-5,30,ERecoil);
	  MyFill("RecE_vs_BeamQval_ElasticScatQ3",600,-5,80,EBeam4He,600,-5,30,ERecoil);
	  MyFill("RecE_vs_LightParE_ElasticScatQ3",600,-5,30,E_4Herxn,600,-5,30,ERecoil);
	  MyFill("RecE_vs_LightParQval_ElasticScatQ3",600,-5,30,E4He_elas,600,-5,30,ERecoil);
	  MyFill("RecE_vs_N_ElasticQ3",256,0,4096,ICne_E_diff,600,-0.1,20,ERecoil);
	  MyFill("SiE_tot_vs_N_ElasticQ3",256,0,4096,ICne_E_diff,600,-0.1,30,E_4Herxn);
	  MyFill("SiE_tot_vs_NCal_ElasticQ3",256,0,70,ICne_E_diff_cal,600,-0.1,30,E_4Herxn);
	}else if(detid<16 && detid>3){
	  MyFill("LightParrxn_vs_Elastic4He_ESX31",600,0,30,E4He_elas,600,0,30,E_4Herxn);
	  MyFill("RecE_vs_N_Elastic_SX3_1",256,0,4096,ICne_E_diff,600,-0.1,30,ERecoil);
	  MyFill("SiE_tot_vs_N_Elastic_SX3_1",256,0,4096,ICne_E_diff,600,-0.1,30,E_4Herxn);
	  MyFill("SiE_tot_vs_NCal_Elastic_SX3_1",256,0,70,ICne_E_diff_cal,600,-0.1,30,E_4Herxn);
	}else if(detid<28 && detid>15){
	  MyFill("LightParrxn_vs_Elastic4He_ESX32",600,0,30,E4He_elas,600,0,30,E_4Herxn);
	}if(detid<16 && detid>-1){
	  MyFill("LightParrxn_vs_Elastic4He_E_Q3&&SX31",600,0,30,E4He_elas,600,0,30,E_4Herxn);
	  MyFill("RecE_vs_N_Elastic_Q3&&SX31",256,0,4096,ICne_E_diff,600,-5,30,ERecoil);
	  MyFill("SiE_tot_vs_N_Elastic_Q3&&SX31",256,0,4096,ICne_E_diff,600,-5,30,E_4Herxn);
	  MyFill("SiE_tot_vs_NCal_Elastic_Q3&&SX31",256,0,70,ICne_E_diff_cal,600,-5,30,E_4Herxn);
	}
	MyFill("Beam_4He_Energy",600,-1,80,EBeam4He);
	MyFill("BeamTrack_Energy_ElasticScat",600,-1,80,be);
	MyFill("BeamEnergy_vs_Beam_4He_Energy",600,-1,80,EBeam4He,600,-1,80,be);
	if(detid<4 && detid>-1){
	  MyFill("Beam_4He_Energy_Q3",600,-1,80,EBeam4He);
	  MyFill("BeamTrack_Energy_ElasticScat_Q3",600,-1,80,be);
	  MyFill("BeamEnergy_vs_Beam_4He_Energy_Q3",600,-1,80,EBeam4He,600,-1,80,be);
	}else if(detid<16 && detid>3){
	  MyFill("Beam_4He_Energy_SX3_1",600,-1,80,EBeam4He);
	  MyFill("BeamTrack_Energy_ElasticScat_SX3_1",600,-1,80,be);
	  MyFill("BeamEnergy_vs_Beam_4He_Energy_SX3_1",600,-1,80,EBeam4He,600,-1,80,be);
	}else if(detid<28 && detid>15){
	  MyFill("Beam_4He_Energy_SX3_2",600,-1,80,EBeam4He);
	  MyFill("BeamTrack_Energy_ElasticScat_SX3_2",600,-1,80,be);
	  MyFill("BeamEnergy_vs_Beam_4He_Energy_SX3_2",600,-1,80,EBeam4He,600,-1,80,be);
	}if(detid<16 && detid>-1){
	  MyFill("Beam_4He_Energy_Q3_&&_SX3_1",600,-1,80,EBeam4He);
	  MyFill("BeamTrack_Energy_ElasticScat_Q3_&&_SX3_1",600,-1,80,be);
	  MyFill("BeamEnergy_vs_Beam_4He_Energy_Q3_&&_SX3_1",600,-1,80,EBeam4He,600,-1,80,be);
	}
      }
    } 
  }
}

/*run()
 *Method that runs the code
 *Takes in a name for a .txt data list and a single output file for histograms,
 *tracking data, and recoil data
 *This is really the only function outside of SetFlag() that should ever be called by the 
 *main.
 */
void analyzer::run() {


  Si.ReadDet = 0;
  Si.ReadHit = 0;
  PC.ReadHit = 0;
  CsI.ReadHit = 0;


  if(ReadPCWire) {
    WireRadii = GetPCWireRadius();
  } else {
    for(int wire=0; wire<NPCWires; wire++) {
      WireRadii.push_back(pcr);
    }
  }

  Ne18_eloss = new LookUp("../anasen_software/srim_files/18Ne_in_HeCO2_348Torr_18Nerun_09122018.eloss", M_18Ne);  
  Ne18_eloss->InitializeLookupTables(90.0,760.0,0.02,0.04);  
  Na21_eloss = new LookUp("../anasen_software/srim_files/21Na_in_HeCO2_348Torr_18Nerun_02192019.eloss", M_21Na);  
  Na21_eloss->InitializeLookupTables(90.0,760.0,0.02,0.04); 
  Ne20_eloss = new LookUp("../anasen_software/srim_files/20Ne_in_HeCO2_348Torr_18Nerun_02192019.eloss", M_20Ne);  
  Ne20_eloss->InitializeLookupTables(90.0,760.0,0.02,0.04); 
  alpha_eloss = new LookUp("../anasen_software/srim_files/He_in_HeCO2_379Torr_allruns_09122018.eloss", M_alpha);
  alpha_eloss->InitializeLookupTables(60.0,2400.0,0.02,0.04); 
  proton_eloss = new LookUp("../anasen_software/srim_files/H_in_HeCO2_379Torr_allruns_09122018.eloss",M_p);
  proton_eloss->InitializeLookupTables(40.0,15000.0,0.02,0.04);

  string inputlistName;
  string outputName;
  cout<<"Enter inputlist file name: ";
  cin>>inputlistName;
  cout<<"Enter output ROOTFile name: ";
  cin>>outputName;
  char inputlistname[100], outputname[100];
  strcpy(inputlistname, inputlistName.c_str());
  strcpy(outputname, outputName.c_str());

  TFile *outFile = new TFile(outputname, "RECREATE");
  rootObj = new TObjArray();
  TTree *outTree = new TTree("TrackTree", "TrackTree");
  
  outTree->Branch("NTracks1", &tracks.NTracks1, "NTracks1/I");
  outTree->Branch("NTracks", &tracks.NTracks, "NTracks/I");
  outTree->Branch("NTracks2", &tracks.NTracks2, "NTracks2/I");
  outTree->Branch("NTracks3", &tracks.NTracks3, "NTracks3/I");
  outTree->Branch("NTracks4", &tracks.NTracks4, "NTracks4/I");
  outTree->Branch("TrackEvents", &tracks.TrEvent);
  outTree->Branch("recoil", &recoil);
  outTree->Branch("RFTime", &RFTime, "RFTime/I");
  outTree->Branch("MCPTime", &MCPTime, "MCPTime/I");
  outTree->Branch("TDC2", &TDC2, "TDC2/I");
  outTree->Branch("ICne_En_diff", &ICne_E_diff, "ICne_En_diff/D");
  outTree->Branch("ICne_E_sum", &ICne_E_sum, "ICne_E_sum/D");
  outTree->Branch("ICne_T_diff", &ICne_T_diff, "ICne_T_diff/D");
  outTree->Branch("ICne_T_sum", &ICne_T_sum, "ICne_T_sum/D");

  rootObj->Add(protonCut); 
  rootObj->Add(alphaCut); rootObj->Add(needleCut);

  ifstream inputList;
  inputList.open(inputlistname);
  if(!inputList.is_open()){
    cout<<"List of rootfiles could not be opened!"<<endl;
    exit(EXIT_FAILURE);
  }
  if(cutFlag) {
    string protonCutfile, alphaCutfile, needleCutfile;
    getline(inputList, protonCutfile);
    getline(inputList, alphaCutfile);
    getline(inputList, needleCutfile);
    getCut(protonCutfile, alphaCutfile, needleCutfile);
  }
  string rootName;
  char rootname[100];

  while(getline(inputList, rootName)) {
  
    if(rootName.empty()) {
      cout<<"Root file list is not formated correctly!"<<endl;
      exit(EXIT_FAILURE);
    }
    strcpy(rootname, rootName.c_str());
    TFile *inputFile = new TFile(rootname, "READ");
    if (!inputFile->IsOpen()) {
      cout << "ROOT file: "<<rootName<<" does not exist!"<<endl;
      exit(EXIT_FAILURE);
    }
    cout<<"Processsing "<<rootName<<" ..."<<endl;
    cout<<"Beam Energy: "<<BeamE<<endl;

    TTree *inputTree = (TTree*) inputFile->Get("MainTree");

    inputTree->SetBranchAddress("Si.NSiHits", &Si.NSiHits);
    inputTree->SetBranchAddress("Si.Detector", &Si.ReadDet);
    inputTree->SetBranchAddress("Si.Hit", &Si.ReadHit);
    inputTree->SetBranchAddress("PC.NPCHits", &PC.NPCHits);
    inputTree->SetBranchAddress("PC.Hit", &PC.ReadHit);
    //inputTree->SetBranchAddress("CsI.NCsIHits", &CsI.NCsIHits);
    //inputTree->SetBranchAddress("CsI.Hit", &CsI.ReadHit);
    inputTree->SetBranchAddress("RFTime", &input_RFTime);
    inputTree->SetBranchAddress("MCPTime", &input_MCPTime);
    //inputTree->SetBranchAddress("ICne_En_diff", &input_ICne_E_diff);
    //inputTree->SetBranchAddress("ICne_En_sum", &input_ICne_E_sum);
    //inputTree->SetBranchAddress("ICne_Ti_diff", &input_ICne_T_diff);
    //inputTree->SetBranchAddress("ICne_Ti_sum", &input_ICne_T_sum);

    Int_t nentries = inputTree->GetEntries();
    float blentries = nentries;
    cout<<"Number of entries: "<<nentries<<endl;

    for (int entry = 0; entry<nentries; entry++) {
      cout<<"\rPercent of file completed: "<<entry/blentries*100.0<<"%  "<<flush;
      inputTree->GetEvent(entry);

      SiEnergy_vec = vector<vector<Double_t>>(28, vector<Double_t>(0));
    
      bool keep_event = MCP_RF();
      if(!keep_event) continue;
       
      tracks.ZeroTrack();
      recoilReset();
      Track1();
      Track2();
      Track3();
      
      if(PCPlots) PCPlotting();
      
      TrackCalc();

      if(PCWireCal) PCWireCalibration();
      if(FillEdE_cor) EdEcor();
      if(PCPlots_2) PCPlotting2();
      if(ResidualEnergy_Calc) CalculateResidE();
      if(Reconstruction_Session) ReconstructMe();
      if(Elastic_Scat_Alpha) AlphaScatter();

      if(FillTree) outTree->Fill();
    }//end event loop
  }//end file list loop
  inputList.close();
  outFile->cd();
  outTree->Write();
  rootObj->Write();
  outFile->Close();
  cout<<"All Root objects are written to: "<<outputName<<endl;
  cout<<outputName<<" is closed"<<endl;
}
