/*gwmAnalyzer.cpp
 *Analyzer class for ANASEN Detector in Active Target mode. Contains all functions necessary to 
 *sort data and calculate physical values. Will do everything from making all of the tracking data
 *up to calling the experiment dependant reconstruction and sorting and storing all data. Current
 *asks user for the name of a data list file and an output file
 *
 *Gordon M. -- April 2019
 *Based on previous versions written by M. Anastasiou, N. Rijal, J. Parker, et al
 */
#include "gwmAnalyzer.h"

using namespace std;

//Constructor
analyzer::analyzer() {

  be7eloss_name = "./srim/7be_in_d2_290torr.eloss";
  he4eloss_name = "./srim/4he_in_d2_400torr.eloss";
  he3eloss_name = "./srim/3he_in_d2_400torr.eloss";
  peloss_name = "./srim/p_in_d2_400torr.eloss";
  deloss_name = "./srim/d_in_d2_400torr.eloss";

  //Set all of the flags
  PCPlots = 0;
  PCPlots_2 = 0;
  Beam_and_Eloss = 0;
  FillTree = 0;
  FillEdE_cor = 0;
  MCP_RF_Cut = 0;
  ReadPCWire = 0;
  PCWireCal = 0;
  ResidualEnergy_Calc = 0;
  Reconstruction_Session = 0;
  Elastic_Scat_Alpha = 0;
  cutFlag = 0;

}

//Destructor
analyzer::~analyzer() {
  //free dynamic memory
  delete be7_eloss;
  delete deuteron_eloss;
  delete he3_eloss;
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
void analyzer::getCut(string pcutfile, string acutfile, string he3cutfile, string dcutfile, string                      jacutfile) {
  cout<<"Proton cut file: "<<pcutfile<<endl;
  cout<<"Alpha cut file: "<<acutfile<<endl;
  cout<<"He3 cut file: "<<he3cutfile<<endl;
  cout<<"Deuteron cut file: "<<dcutfile<<endl;
  cout<<"Joined Alpha cut file: "<<jacutfile<<endl;
  char protonname[pcutfile.length()], alphaname[acutfile.length()], he3name[he3cutfile.length()];
  char deuteronname[dcutfile.length()], janame[jacutfile.length()];

  strcpy(protonname, pcutfile.c_str());
  strcpy(alphaname, acutfile.c_str());
  strcpy(he3name, he3cutfile.c_str());
  strcpy(deuteronname, dcutfile.c_str());
  strcpy(janame, jacutfile.c_str());

  TFile *protonfile = new TFile(protonname, "READ");
  protonCut = (TCutG*) protonfile->Get("CUTG");
  protonCut->SetName("protonCut");

  TFile *alphafile = new TFile(alphaname, "READ");
  alphaCut = (TCutG*) alphafile->Get("CUTG");
  alphaCut->SetName("alphaCut");

  TFile *he3file = new TFile(he3name, "READ");
  he3Cut = (TCutG*) he3file->Get("CUTG");
  he3Cut->SetName("he3Cut");

  TFile *deutfile = new TFile(deuteronname, "READ");
  deutCut = (TCutG*) deutfile->Get("CUTG");
  deutCut->SetName("deutCut");

  TFile *jafile = new TFile(janame, "READ");
  joinedAlphaCut = (TCutG*) jafile->Get("CUTG");
  joinedAlphaCut->SetName("joinedAlphaCut");
}

/*recoilReset()
 *With RecoilEvent being promoted to a global structure there needs to be
 *a method for setting the class instance of RecoilEvent back to zero
 *Was originally contained in the Track class
 */
void analyzer::recoilReset(RecoilEvent &recoil) {
  recoil.IntPoint = -10;
  recoil.SiEnergy_tot = -10;
  recoil.SiEnergy_calc = -10;
  recoil.PCEnergy_tot = -10;
  recoil.Energy_eject_tot = -10;
  recoil.Theta = -10;
  recoil.Phi = -10;
  recoil.KE_recoil = -10;
  recoil.Ex_recoil = -10;
  recoil.Ex_recoil_2a = -10;
  recoil.Ex_recoil_a1p = -10;
  recoil.Ex_recoil_a2p = -10;
  recoil.Ex_temp_small = -10;
  recoil.Ex_temp_large = -10;
  recoil.Ex_temp_qvalue_small = -10;
  recoil.Ex_temp_qvalue_large = -10;
  recoil.Delta_Theta = -10;
  recoil.Delta_Phi = -10;
  recoil.BeamKE_WAIntP = -10;
  recoil.theta_eject1 = -10;
  recoil.theta_eject2 = -10;
  recoil.BeamKE = -10;
  recoil.BeamKE_eject = -10;
  recoil.recoil_mass_sq = -10;
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
  Int_t NexttoMaxPCindex = -1;
  Double_t NexttoMaxPC = -10;
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
      } else if (PC.pc_obj.Energy >= NexttoMaxPC) {
        NexttoMaxPCindex = k;
        NexttoMaxPC = PC.pc_obj.Energy;
      } 
    }
  }
  if (NexttoMaxPCindex>0 && (PhiDiff((*PC.ReadHit)[MaxPCindex].PhiW, phi) >
                             PhiDiff((*PC.ReadHit)[NexttoMaxPCindex].PhiW,phi))) {
    return NexttoMaxPCindex;
  }
  return MaxPCindex;

}

/*RecoverTrack1()
 *This method is for recovering a Track of type 1 from a Track of type 2 when there are two 
 *tracks that are so close together that the PC cannot distinguish them. When the tracks are this
 *close the sorting will assign one of the Si hits as a track1 with all of the PC info, along with
 *the individual Si info, while the other Si hit will be left as a track2 with only Si info
 */
bool analyzer::RecoverTrack1(TrackEvent track1, TrackEvent &track2) {
  if(track1.TrackType == 1 && track2.TrackType == 2) {
    Double_t pce1 = track1.PCEnergy;
    Double_t pcz1 = track1.PCZ;
    Double_t pcr1 = track1.PCR;
    Double_t wid1 = track1.WireID;
    Double_t down1 = track1.Down;
    Double_t up1 = track1.Up;
    Double_t downV1 = track1.DownVoltage;
    Double_t upV1 = track1.UpVoltage;
    Double_t siz2 = track2.SiZ;
    Double_t sir2 = track2.SiR;
   
    Double_t m = (pcr1-sir2)/(pcz1-siz2);
    if(m == 0.0) {
      return false;
    }
    Double_t b =  pcr1-m*pcz1;
    Double_t intp = -b/m;
    
    Double_t theta, pl;      
    if((intp - siz2)>0) {
      theta = atan(sir2/(intp-siz2));
      pl = sir2/sin(theta);
    } else if ((intp-siz2) < 0) {
      theta = TMath::Pi() + atan(sir2/(intp-siz2));
      pl = sir2/sin(theta);
    } else { 
      theta = TMath::Pi()/2.0;
      pl = sir2;
    }
    float length_check = ana_length - intp;
    Double_t be;
    if(length_check>0.0 && length_check<ana_length) {   
      be = be7_eloss->GetLookupEnergy(BeamE, length_check);
    } else {
      return false;
    }
    
    track2.PCEnergy = pce1;
    track2.PCZ = pcz1;
    track2.PCR = pcr1;
    track2.WireID = wid1;
    track2.Down = down1;
    track2.Up = up1;
    track2.DownVoltage = downV1;
    track2.UpVoltage = upV1;
    track2.IntPoint = intp;
    track2.Theta = theta;
    track2.PathLength = pl;
    track2.BeamEnergy = be;
    track2.TrackType = 1;
    return true;
  } else {
    cout<<"Error in RecoverTrack1! track1 must be type 1 and track2 must be type 2"<<endl;
    return false;
  }
}

/*MCP_RF()
 *First of the analysis methods; Makes and fills histograms for the MCP and RF
 *timing. Also performs the major cut on the data, the MCP_RF cut, where the correct 
 *beam for the experiment is selected
 */
bool analyzer::MCP_RF() {


  Double_t correct = 1.004009623;
  MCPTime = input_MCPTime;
  RFTime = input_RFTime;
  TDC2 = input_TDC2;
  MyFill("MCPTime",1024,1,4096,MCPTime);
  MyFill("RFTime",1024,1,4096,RFTime);
  MyFill("TDC2",1024,1,4096,TDC2);
  MyFill("Timing",600,0,600,fmod((MCPTime*correct-RFTime),545));
  if (MCP_RF_Cut) {
    if ((MCPTime>2800&&MCPTime<3200) && RFTime>0) {
      float wrap_val = fmod((MCPTime*correct-RFTime), 545);
      if((wrap_val<47 || wrap_val>118) && (wrap_val<320 || wrap_val>384)){
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
      if(trackhit.WireID==6) {
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
 *the interaction point and beam energy for each event calculated purely through tracking
 */
void analyzer::TrackCalc() {

  for(int i=0; i<tracks.NTracks1; i++) {
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t pcr = tracks.TrEvent[i].PCR;
    Double_t sir = tracks.TrEvent[i].SiR;
    
    if(siz < 0.0 || pcz < 0.0) {
      tracks.TrEvent[i].PCZ = -100;
      tracks.TrEvent[i].SiZ = -100;
      continue;
    }
    if(pcr>4.0 || pcr<3.5) {
      tracks.TrEvent[i].PCR = -100;
      continue;
    }
    if(sir>11.0 || sir<4.0) {
      tracks.TrEvent[i].SiR = -100;
      continue;
    }
    
    Double_t m = (pcr-sir)/(pcz-siz);
    if(m==0 || pcz<=0 || siz<0 || sir<=0 || (pcz-siz) ==0) {
      cout<<"m for Interaction Point is zero "<<m<<endl;
      continue;
    }

    Double_t b = (tracks.TrEvent[i].PCR - m*tracks.TrEvent[i].PCZ);
    tracks.TrEvent[i].IntPoint = -b/m;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    if(intp<0.0 || intp>55.0) {
      tracks.TrEvent[i].IntPoint = -100;
      continue;
    }
    
    if((intp - siz)>0) {
      tracks.TrEvent[i].Theta = atan(sir/(intp-siz));
      Double_t theta = tracks.TrEvent[i].Theta;
      tracks.TrEvent[i].PathLength = sir/sin(theta);
    } else if ((intp-siz) < 0) {
      tracks.TrEvent[i].Theta = TMath::Pi() + atan(sir/(intp-siz));
      Double_t theta = tracks.TrEvent[i].Theta;
      tracks.TrEvent[i].PathLength = sir/sin(theta);
    } else { 
      tracks.TrEvent[i].Theta = TMath::Pi()/2.0;
      tracks.TrEvent[i].PathLength = sir;
    }

    if (Beam_and_Eloss) {
      float length_check = ana_length - intp;
      if(length_check>0.0 && length_check<ana_length) {   
        tracks.TrEvent[i].BeamEnergy = be7_eloss->GetLookupEnergy(BeamE, length_check);
        if(tracks.TrEvent[i].BeamEnergy<0.0 || tracks.TrEvent[i].BeamEnergy>BeamE){
          tracks.TrEvent[i].BeamEnergy = -100.0;
          continue;
        }
      }
    }
  }
}

/*PCWireCalibration()
 *Option for looking at pc wires calibration; 
 *Makes a very large number of plots, so only run if necessary
 *Requires gold_pos, which has to be entered manually
 */
void analyzer::PCWireCalibration() {

  Double_t mpc, bpc;
  Double_t qqq_r1_Emin = 11.3;
  Double_t r2_Emin = 10.3;
  Double_t r2_Emax = 11.2;

  Double_t gold_pos = 28.9; //set by experiment

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

    MyFill("BeamEnergy_vs_IntPoint",200,-10,60,intp,200,-10,30,be);
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
 *Focuses on Energy of the recoil; single ejectile reconstruction
 */
void analyzer::CalcRecoilE() {
  //Make a be8 calculator
  Reconstruct be8_calc(m_7be, m_d, m_p, m_8be);
  be8_calc.ELoss_eject = proton_eloss;
  be8_calc.ELoss_beam = be7_eloss;
  
  //Make a li6 calculator
  Reconstruct li6_calc(m_7be, m_d, m_3he, m_6li);
  li6_calc.ELoss_eject = he3_eloss;
  li6_calc.ELoss_beam = be7_eloss; 

  for(int i=0; i<tracks.NTracks1; i++) {
    Int_t detid = tracks.TrEvent[i].DetID;
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t siphi = tracks.TrEvent[i].SiPhi;
    Double_t theta = tracks.TrEvent[i].Theta;
    Double_t pcphi = tracks.TrEvent[i].PCPhi;
    Double_t sir = tracks.TrEvent[i].SiR;

    //All histos in this section tagged with _cre (CalcRecoilE) to avoid confusion
    //in rootfile
    if(detid>-1 && detid<28 && pcz>0.0 && siz>=0.0 && intp>0.0 && intp<ana_length) {
      if(protonCut->IsInside(sie, pce*sin(theta)) && tracks.NTracks1 == 3) {
        be8_calc.CalcRecoil(tracks, i, Be8_1p);
        Double_t eetot = Be8_1p.Energy_eject_tot;
        Double_t ker = Be8_1p.KE_recoil;
        Double_t rex = Be8_1p.Ex_recoil;
        Double_t rbe = Be8_1p.BeamKE;
        Double_t rtheta = Be8_1p.Theta;

	if(detid>-1 && detid<4){
	  MyFill("ProtonE_vs_BeamE_Tracked_Q3_cre",300,-0.1,80,be,300,-0.1,30,eetot);
	  MyFill("8BeKE_vs_BeamE_Tracked_Q3_cre",300,-0.1,80,be,300,-0.1,30,ker);
	  MyFill("ProtonE_vs_8BeKE_Q3_cre",300,-0.1,80,ker,300,-0.1,30,eetot);
	  MyFill("8BeKE_vs_8BeTheta_Q3_cre",300,-0.1,180,rtheta,300,-0.1,80,ker);
	  MyFill("8BeKE_vs_IntP_Tracked_Q3_cre",300,-0.1,60,intp,300,-0.1,80,ker);
	} else if (detid>3 && detid<16) {
	  MyFill("ProtonE_vs_BeamE_Tracked_R1_cre",300,-0.1,80,be,300,-0.1,30,eetot);
	  MyFill("8BeKE_vs_BeamE_Tracked_R1_cre",300,-0.1,80,be,300,-0.1,80,ker);
	  MyFill("ProtonE_vs_8BeKE_R1_cre",300,-0.1,80,ker,300,-0.1,30,eetot);
	  MyFill("8BeKE_vs_8BeTheta_R1_cre",300,-0.1,180,rtheta,300,-0.1,30,ker);
	  MyFill("8BeKE_vs_IntP_Tracked_R1_cre",300,-0.1,60,intp,300,-0.1,80,ker);
        }
        be8_calc.CalcRecoil_from_Qvalue(QValue_8Be, tracks, i, Be8_1p_qval);
        MyFill("Ex_8Be_Qvalue_Rec_cre",600,-10,40,Be8_1p_qval.Ex_recoil);  
        MyFill("Ex_8Be_vs_BeamKE_fromQv_cre",600,-10,80,Be8_1p_qval.BeamKE,
                                            400,-10,20,Be8_1p_qval.Ex_recoil);
        MyFill("8Be_BeamEnergy_Tracked_vs_BeamKE_fromQv_cre",600,-10,20,Be8_1p_qval.BeamKE,
                                                             600,-10,20,be);
        if (rex<0.0) {
          MyFill("PCZ_Ex8Be<0_cre",400,-10,40,pcz);
          MyFill("IntPoint_Ex8Be<0_cre",400,-10,80,intp);
          MyFill("BeamEnergy_Ex8Be<0_cre",400,-10,80,be);
          MyFill("ProtonEnergy_vs_ThetaEx8Be<0_cre",400,0,190,theta*rads2deg,400,0,15,eetot);
          MyFill("BeamEnergy_vs_IntPoint_Ex8Be<0_cre",400,-10,80,intp,400,-10,80,be);
          MyFill("ProtonEnergy_vs_BeamEnergy_Ex8Be<0_cre",400,-10,80,be,400,0,15,eetot);
          MyFill("ProtonEnergy_vs_SiR_Ex8Be<0_cre",400,-1,11,sir,400,0,15,eetot);
          MyFill("ProtonEnergy_vs_IntPoint_Ex8Be<0_cre",400,-10,80,intp,400,0,15,eetot);
          MyFill("IntPoint_vs_DetID_Ex8Be<0_cre",28,0,28,detid,400,-10,80,intp);
          MyFill("IntPoint_vs_PCZ_Ex8Be<0_cre",400,-10,40,pcz,400,-10,80,intp);
          MyFill("SiPhi_vs_PCPhi_Ex8Be<0_cre",400,0,360,pcphi*rads2deg,400,-0,360,siphi*rads2deg);
        }
        MyFill("Ex_of_8BeRecoil_cre",600,-10,40,rex);
        MyFill("BeamKE_8Be_cre",600,-10,80,rbe);
        MyFill("BeamKE_8Be_cm_cre",600,-1,15,rbe*4/22);
        MyFill("BeamEnergy_8Be_Tracked_cre",600,-10,80,be);
        MyFill("ProtonE_Tracked_8BeProtons_cre",600,-10,20,eetot);
        MyFill("Ex8Be_vs_BeamKE_cre",600,-10,80,rbe,400,-10,20,rex);
        MyFill("ProtonE_vs_Ex8Be_cre",600,-10,60,rex,400,-10,20,eetot);
        MyFill("ProtonE_vs_ProtonTheta_8Be_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
        MyFill("ProtonE_vs_BeamKE_8Be_cre",600,-10,80,rbe,600,0,20,eetot);
        MyFill("BeamKE_Tr_vs_BeamKE_8Be_cre",600,-10,80,rbe,600,-10,80,be);
        if(detid>-1 && detid<4) {
          MyFill("Ex8Be_Q3_cre",600,-10,40,rex);
          MyFill("BeamKE_8Be_Q3_cre",600,-10,80,rbe);
          MyFill("BeamKE_8Be_Q3_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_Q3_cre",600,-10,80,be);
          MyFill("ProtonE_Tracked_8BeProtons_Q3_cre",600,-10,20,eetot);
          MyFill("Ex8Be_vs_BeamKE_Q3_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex8Be_Q3_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_8Be_Q3_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_8Be_Q3_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_8Be_Q3_cre",600,-10,80,rbe,600,-10,80,be);
        } else if (detid>3 && detid <16) {
          MyFill("Ex8Be_R1_cre",600,-10,40,rex);
          MyFill("BeamKE_8Be_R1_cre",600,-10,80,rbe);
          MyFill("BeamKE_8Be_R1_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_R1_cre",600,-10,80,be);
          MyFill("ProtonE_Track_8BeProtons_R1_cre",600,-10,20,eetot);
          MyFill("Ex8Be_vs_BeaKE_R1_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex8Be_R1_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_8Be_R1_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_8Be_R1_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_8Be_R1_cre",600,-10,80,rbe,600,-10,80,be);
        } else if (detid>15 && detid<28) {
          MyFill("Ex8Be_R2_cre",600,-10,40,rex);
          MyFill("BeamKE_8Be_R2_cre",600,-10,80,rbe);
          MyFill("BeamKE_8Be_R2_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_R2_cre",600,-10,80,be);
          MyFill("ProtonE_Track_8BeProtons_R2_cre",600,-10,20,eetot);
          MyFill("Ex8Be_vs_BeaKE_R2_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex8Be_R2_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_8Be_R2_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_8Be_R2_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_8Be_R2_cre",600,-10,80,rbe,600,-10,80,be);
        }
        if(tracks.NTracks1==1) {
          MyFill("ExcitationEnergy_p0_ntracks1_1_cre",400,-10,20,rex);
          MyFill("ExEnergy_vs_BeamEnergy_p0_ntracks1_1_cre",600,-10,80,be,400,-10,20,rex);
          MyFill("ExEnergy_vs_BeamE_cm_p0_ntracks1_1_cre",600,-10,20,be*4/22,400,-10,20,rex);
        }
      } else if (he3Cut->IsInside(sie, pce*sin(theta))&&tracks.NTracks1==1) {
        li6_calc.CalcRecoil(tracks, i, Li6);
        Double_t eetot = Li6.Energy_eject_tot;
        Double_t ker = Li6.KE_recoil;
        Double_t rtheta = Li6.Theta;
        Double_t rex = Li6.Ex_recoil;
        Double_t rbe = Li6.BeamKE;

        MyFill("6LiExEnergy",100,0,15,rex);
        MyFill("BeamE_Tracked_vs_BeamKE_6Li_cre", 300,-0.1,20,be,300,-0.1,20,rbe);
        MyFill("Ex6Li_vs_BeamKE_cre",300, -0.1,20,rbe,300,0,15,rex);
	if(detid>-1 && detid<4){
	  MyFill("3HeE_vs_BeamE_Tracked_Q3_cre",300,-0.1,80,be,300,-0.1,30,eetot);
	  MyFill("6LiKE_vs_BeamE_Tracked_Q3_cre",300,-0.1,80,be,300,-0.1,30,ker);
	  MyFill("3HeE_vs_6LiKE_Q3_cre",300,-0.1,80,ker,300,-0.1,30,eetot);
	  MyFill("6LiKE_vs_6LiTheta_Q3",300,-0.1,180,rtheta,300,-0.1,80,ker);
	  MyFill("6LiKE_vs_IntP_Tracked_Q3_cre",300,-0.1,60,intp,300,-0.1,80,ker);
	} else if (detid>3 && detid<16) {
	  MyFill("3HeE_vs_BeamE_Tracked_R1_cre",300,-0.1,80,be,300,-0.1,30,eetot);
	  MyFill("6LiKE_vs_BeamE_Tracked_R1_cre",300,-0.1,80,be,300,-0.1,80,ker);
	  MyFill("3HeE_vs_6LiKE_R1_cre",300,-0.1,80,ker,300,-0.1,30,eetot);
	  MyFill("6LiKE_vs_6LiTheta_R1_cre",300,-0.1,180,rtheta,300,-0.1,30,ker);
	  MyFill("6LiKE_vs_IntP_Tracked_R1_cre",300,-0.1,60,intp,300,-0.1,80,ker);
        }
        li6_calc.CalcRecoil_from_Qvalue(QValue_6Li, tracks, i, Li6_qval);
        //cout<<"Li6 BeamKE: "<<Li6_qval.BeamKE<<endl;
        MyFill("Ex_6Li_Qvalue_Rec_cre",600,-10,40,Li6_qval.Ex_recoil);  
        MyFill("Ex_6Li_vs_BeamKE_fromQv_cre",600,-10,80,Li6_qval.BeamKE,
                                             400,-10,20,Li6_qval.Ex_recoil);
        MyFill("6Li_BeamEnergy_Tracked_vs_BeamKE_fromQv_cre",600,-10,20,Li6_qval.BeamKE,
                                                             600,-10,20,be);
        if (rex<0.0) {
          MyFill("PCZ_Ex6Li<0_cre",400,-10,40,pcz);
          MyFill("IntPoint_Ex6Li<0_cre",400,-10,80,intp);
          MyFill("BeamEnergy_Ex6Li<0_cre",400,-10,80,be);
          MyFill("ProtonEnergy_vs_ThetaEx6Li<0_cre",400,0,190,theta*rads2deg,400,0,15,eetot);
          MyFill("BeamEnergy_vs_IntPoint_Ex6Li<0_cre",400,-10,80,intp,400,-10,80,be);
          MyFill("ProtonEnergy_vs_BeamEnergy_Ex6Li<0_cre",400,-10,80,be,400,0,15,eetot);
          MyFill("ProtonEnergy_vs_SiR_Ex6Li<0_cre",400,-1,11,sir,400,0,15,eetot);
          MyFill("ProtonEnergy_vs_IntPoint_Ex6Li<0_cre",400,-10,80,intp,400,0,15,eetot);
          MyFill("IntPoint_vs_DetID_Ex6Li<0_cre",28,0,28,detid,400,-10,80,intp);
          MyFill("IntPoint_vs_PCZ_Ex6Li<0_cre",400,-10,40,pcz,400,-10,80,intp);
          MyFill("SiPhi_vs_PCPhi_Ex6Li<0_cre",400,0,360,pcphi*rads2deg,400,-0,360,siphi*rads2deg);
        }
        MyFill("Ex_of_6LiRecoil_cre",600,-10,40,rex);
        MyFill("BeamKE_6Li_cre",600,-10,80,rbe);
        MyFill("BeamKE_6Li_cm_cre",600,-1,15,rbe*4/22);
        MyFill("BeamEnergy_6Li_Tracked_cre",600,-10,80,be);
        MyFill("ProtonE_Tracked_6LiProtons_cre",600,-10,20,eetot);
        MyFill("Ex6Li_vs_BeamKE_cre",600,-10,80,rbe,400,-10,20,rex);
        MyFill("ProtonE_vs_Ex6Li_cre",600,-10,60,rex,400,-10,20,eetot);
        MyFill("ProtonE_vs_ProtonTheta_6Li_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
        MyFill("ProtonE_vs_BeamKE_6Li_cre",600,-10,80,rbe,600,0,20,eetot);
        MyFill("BeamKE_Tr_vs_BeamKE_6Li_cre",600,-10,80,rbe,600,-10,80,be);
        if(detid>-1 && detid<4) {
          MyFill("Ex6Li_Q3_cre",600,-10,40,rex);
          MyFill("BeamKE_6Li_Q3_cre",600,-10,80,rbe);
          MyFill("BeamKE_6Li_Q3_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_Q3_cre",600,-10,80,be);
          MyFill("ProtonE_Tracked_6LiProtons_Q3_cre",600,-10,20,eetot);
          MyFill("Ex6Li_vs_BeamKE_Q3_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex6Li_Q3_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_6Li_Q3_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_6Li_Q3_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_6Li_Q3_cre",600,-10,80,rbe,600,-10,80,be);
        } else if (detid>3 && detid <16) {
          MyFill("Ex6Li_R1_cre",600,-10,40,rex);
          MyFill("BeamKE_6Li_R1_cre",600,-10,80,rbe);
          MyFill("BeamKE_6Li_R1_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_R1_cre",600,-10,80,be);
          MyFill("ProtonE_Track_6LiProtons_R1_cre",600,-10,20,eetot);
          MyFill("Ex6Li_vs_BeaKE_R1_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex6Li_R1_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_6Li_R1_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_6Li_R1_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_6Li_R1_cre",600,-10,80,rbe,600,-10,80,be);
        } else if (detid>15 && detid<28) {
          MyFill("Ex6Li_R2_cre",600,-10,40,rex);
          MyFill("BeamKE_6Li_R2_cre",600,-10,80,rbe);
          MyFill("BeamKE_6Li_R2_cm_cre",600,-1,15,rbe*4/22);
          MyFill("BeamEnergy_Tracked_R2_cre",600,-10,80,be);
          MyFill("ProtonE_Track_6LiProtons_R2_cre",600,-10,20,eetot);
          MyFill("Ex6Li_vs_BeaKE_R2_cre",600,-10,80,rbe,400,-10,20,rex);
          MyFill("ProtonE_vs_Ex6Li_R2_cre",600,-10,60,rex,400,-10,20,eetot);
          MyFill("ProtonE_vs_ProtonTheta_6Li_R2_cre",600,0,180,theta*rads2deg,600,0,20,eetot);
          MyFill("ProtonE_vs_BeamKE_6Li_R2_cre",600,-10,80,rbe,600,0,20,eetot);
          MyFill("BeamEnergy_Track_vs_BeamKE_6Li_R2_cre",600,-10,80,rbe,600,-10,80,be);
        }
        if(tracks.NTracks1==1) {
          MyFill("ExcitationEnergy_3he_ntracks1_1_cre",400,-10,20,rex);
          MyFill("ExEnergy_vs_BeamEnergy_3he_ntracks1_1_cre",600,-10,80,be,400,-10,20,rex);
          MyFill("ExEnergy_vs_BeamE_cm_3he_ntracks1_1_cre",600,-10,20,be*4/22,400,-10,20,rex);
        }
      }
    }
  }
  be8_calc.ELoss_eject = NULL;
  be8_calc.ELoss_beam = NULL;
  li6_calc.ELoss_eject = NULL;
  li6_calc.ELoss_beam = NULL;
}

void analyzer::CalcSmAngleAlphas() {
  Reconstruct be8_calc(m_7be, m_d, m_p, m_alpha, m_8be);
  be8_calc.ELoss_eject = proton_eloss;
  be8_calc.ELoss_eject2 = alpha_eloss;
  be8_calc.ELoss_beam = be7_eloss;
  
  vector<int> IsProton;
  vector<int> IsT1alpha;
  vector<int> IsT2alpha;
  vector<int> IsEject;

  Double_t alpha_pce = 0.0;
  Double_t alpha_t = 0.0;
  Double_t alpha_be = 0.0;
  Double_t alpha_pl = 0.0;
  Double_t alpha_phi = 0.0;

  for(int i=0; i<tracks.NTracks1; i++) {
    Int_t detid = tracks.TrEvent[i].DetID;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t sir = tracks.TrEvent[i].SiR;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t theta = tracks.TrEvent[i].Theta;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t pl = tracks.TrEvent[i].PathLength;
    Double_t phi = tracks.TrEvent[i].SiPhi;
    if(detid>-1 && detid<28 && pcz>0.0 && siz>=0.0 && intp>0.0 && intp<ana_length && 
       be>0.0 && be<BeamE && sir>4.0 && sir<11.0 && tracks.NTracks1 == 2) {
      MyFill("Ede_8be_glob_p",300,0,30,sie/2.0,300,0,0.3,pce/2.0*sin(theta));
      if(protonCut->IsInside(sie, pce*sin(theta))) {
        IsProton.push_back(i);
      } else if(alphaCut->IsInside(sie/2.0, (pce/2.0)*sin(theta))){
        tracks.TrEvent[i].PCEnergy = pce/2.0;
        IsT1alpha.push_back(i);
        alpha_pce = pce/2.0;
        alpha_t = theta;
        alpha_be = be;
        alpha_pl = pl;
        alpha_phi = phi;
      }
    }
  }
  /*if(IsProton.size() == 1 && IsT1alpha.size() == 1) {
    for (int i=tracks.NTracks1; i<(tracks.NTracks1+tracks.NTracks2); i++) {
      bool valid = RecoverTrack1(tracks.TrEvent[IsT1alpha[0]], tracks.TrEvent[i]);
      if(valid) {
        TrackEvent hope = tracks.TrEvent[i];
        Double_t sie2 = hope.SiEnergy;
        Double_t theta2 = hope.Theta;
        Double_t pce2 = hope.PCEnergy;
        Double_t intp1 = tracks.TrEvent[IsT1alpha[0]].IntPoint;
        Double_t intp2 = tracks.TrEvent[i].IntPoint;
        Double_t theta1 = tracks.TrEvent[IsT1alpha[0]].Theta;
        //MyFill("Ede_Recov",300,0,30,sie2,300,0,0.3,alpha_pce*sin(alpha_t));
        MyFill("Ede_Recov",300,0,30,sie2,300,0,0.3,pce2*sin(theta2));
        //if (alphaCut->IsInside(sie2, alpha_pce*sin(alpha_t))){
        if (joinedAlphaCut->IsInside(sie2, pce2*sin(theta2))){
          MyFill("RecovAlphaIP_vs_T1AlphaIP",550,0,55,intp1,550,0,55,intp2);
          MyFill("RecovAlphaTheta_vs_T1AlphaTheta",300,0,2*TMath::Pi(),theta1,300,0,2*TMath::Pi(),theta2);
          IsT2alpha.push_back(i);
          tracks.TrEvent[i].PathLength = alpha_pl;
          tracks.TrEvent[i].BeamEnergy = alpha_be;
          tracks.TrEvent[i].Theta = alpha_t;
          tracks.TrEvent[i].PCEnergy = alpha_pce; //like a flag for later
        }
      }
    }
  }*/
  if(IsT1alpha.size() == 1 && IsProton.size() == 1) {
    IsEject.push_back(IsProton[0]);
    IsEject.push_back(IsT1alpha[0]);
    be8_calc.CalcE_ex(tracks, IsEject, Be8_1p_any, Be8_1p_any_qval);
  }

  be8_calc.ELoss_eject = NULL;
  be8_calc.ELoss_eject2 = NULL;
  be8_calc.ELoss_beam = NULL;
}

/*ReconstructMe()
 *Reconstruction method for multiparticle reconstruction; in the case of this version, used for 
 *reconstruct of 7Be+d->p+2*alpha
 */
void analyzer::ReconstructMe() {
  Reconstruct be8_calc(m_7be,m_d,m_p,m_alpha,m_8be);
  Reconstruct li5_calc(m_7be,m_d,m_p,m_alpha,m_5li);
  be8_calc.ELoss_eject = proton_eloss;
  be8_calc.ELoss_eject2 = alpha_eloss;
  be8_calc.ELoss_beam = be7_eloss;
  li5_calc.ELoss_eject = proton_eloss;
  li5_calc.ELoss_eject2 = alpha_eloss;
  li5_calc.ELoss_beam = be7_eloss;
  vector<Int_t> IsEject;
  vector<Int_t> IsProton;
  vector<Int_t> IsAlpha;
  //Grab all of the protons and alphas 
  for(int i=0; i<tracks.NTracks1; i++) {
    if(tracks.TrEvent[i].TrackType !=1) {
      cout<<"gwmAnalyzer::ReconstructMe() error: Track is not type 1!"<<endl;
    }
    Int_t detid = tracks.TrEvent[i].DetID;
    Double_t pcz = tracks.TrEvent[i].PCZ;
    Double_t intp = tracks.TrEvent[i].IntPoint;
    Double_t siz = tracks.TrEvent[i].SiZ;
    Double_t sir = tracks.TrEvent[i].SiR;
    Double_t be = tracks.TrEvent[i].BeamEnergy;
    Double_t theta = tracks.TrEvent[i].Theta;
    Double_t pce = tracks.TrEvent[i].PCEnergy;
    Double_t sie = tracks.TrEvent[i].SiEnergy;
    if(detid>-1 && detid<28 && pcz>0.0 && siz>=0.0 && intp>0.0 && intp<ana_length && 
       be>0.0 && be<BeamE && sir>4.0 && sir<11.0 /*&& tracks.NTracks1 == 3*/) {
      if(alphaCut->IsInside(sie, pce*sin(theta))) {
        IsAlpha.push_back(i);
      } else if (protonCut->IsInside(sie, pce*sin(theta))) {
        IsProton.push_back(i);
      }
    }
  }
  //all histos in this section tagged with _rm (ReconstructMe) to avoid confusion in 
  //rootfile
  if (IsProton.size() == 1) {
    TrackEvent proton = tracks.TrEvent[IsProton[0]];
    IsEject.push_back(IsProton[0]);
    Double_t sie = proton.SiEnergy;
    Double_t be = proton.BeamEnergy;
    Double_t intp = proton.IntPoint;
    Double_t theta = proton.Theta;

    if(tracks.NTracks1==1) {//only single proton events
      MyFill("SiProtonE_theta_p0_ntracks1_1_rm",400,-1,190,theta*rads2deg,400,-1,15,sie);
      MyFill("SiProtonE_vs_IntP_p0_ntracks1_1_rm",400,-1,70,intp,400,-1,15,sie);
      MyFill("SiProtonE_vs_BeamEnergy_p0_ntracks1_1_rm",400,-1,80,be,400,-1,15,sie);
      MyFill("BeamE_cm_p0_ntracks1_1_rm",180,-1,15,be*4/22);
    } else if (IsAlpha.size()==2) {//only 1 proton, 2 alpha events
      IsEject.push_back(IsAlpha[0]);
      IsEject.push_back(IsAlpha[1]);
      TrackEvent alpha1 = tracks.TrEvent[IsAlpha[0]];
      TrackEvent alpha2 = tracks.TrEvent[IsAlpha[1]];
      Double_t a1theta = alpha1.Theta;
      Double_t a1siphi = alpha1.SiPhi;
      Double_t a1pcphi = alpha1.PCPhi;
      Double_t a1pcz = alpha1.PCZ;
      Double_t a1sie = alpha1.SiEnergy;
      Double_t a1pce = alpha1.PCEnergy;
      Double_t a1be = alpha1.BeamEnergy;
      Double_t a1intp = alpha1.IntPoint;
      Int_t a1detid = alpha1.DetID;
      Double_t a2theta = alpha2.Theta;
      Double_t a2siphi = alpha2.SiPhi;
      Double_t a2pcphi = alpha2.PCPhi;
      Double_t a2pcz = alpha2.PCZ;
      Double_t a2sie = alpha2.SiEnergy;
      Double_t a2pce = alpha2.PCEnergy;
      Double_t a2be = alpha2.BeamEnergy;
      Double_t a2intp = alpha2.IntPoint;
      Int_t a2detid = alpha2.DetID;

      MyFill("IntPoint_a1_vs_IntPoint_a2_rm",400,-10,70,a1intp,400,-10,70,a2intp);
      MyFill("SiProtonE_de_a1_rm",400,-1,35,a1sie,400,-0.01,0.8,a1pce*sin(a1theta));
      MyFill("SiProtonE_theta_a1_rm",400,-1,190,a1theta*rads2deg,400,-1,15,a1sie);
      MyFill("SiProtonE_de_a2_rm",400,-1,35,a2sie,400,-0.01,0.8,a2pce*sin(a2theta));
      MyFill("SiProtonE_theta_a2_rm",400,-1,190,a2theta*rads2deg,400,-1,15,a2sie);
      MyFill("Theta_a1_vs_Theta_a2_rm",400,-1,190,a1theta*rads2deg,400,-1,190,a2theta*rads2deg);
      MyFill("SiPhi_a1_vs_SiPhi_a2_rm",400,-1,360,a1siphi*rads2deg,400,-1,360,a2siphi*rads2deg);
      MyFill("SiProtonE_a1_vs_SiProtonE_a2_rm",400,-1,20,a1sie,400,-1,20,a2sie);
      MyFill("BeamE_a1_vs_BeamE_a2_rm",400,-1,80,a1be,400,-1,80,a2be);
      MyFill("PCZ_a1_vs_PCZ_a2_rm",600,-10,60,a1pcz,600,-1,60,a2pcz);
      MyFill("BeamE_cm_a1_rm",180,-1,15,a1be*4/22);
      MyFill("BeamE_cm_a1_rm",180,-1,15,a2be*4/22);

      //Why do these exist? Never used ever again; terminal plotting tools
      tracks.TrEvent[IsAlpha[0]].BeamE_p1 = alpha1.BeamEnergy;
      tracks.TrEvent[IsAlpha[1]].BeamE_p2 = alpha2.BeamEnergy;
      tracks.TrEvent[IsAlpha[0]].SiE_p1 = alpha1.SiEnergy;
      tracks.TrEvent[IsAlpha[1]].SiE_p2 = alpha2.SiEnergy;
      tracks.TrEvent[IsAlpha[0]].IntP_p1 = alpha1.IntPoint;
      tracks.TrEvent[IsAlpha[1]].IntP_p2 = alpha2.IntPoint;
      tracks.TrEvent[IsAlpha[0]].Theta_p1 = alpha1.Theta;
      tracks.TrEvent[IsAlpha[1]].Theta_p2 = alpha2.Theta;
      tracks.TrEvent[IsAlpha[0]].SiPhi_p1 = alpha1.SiPhi;
      tracks.TrEvent[IsAlpha[1]].SiPhi_p2 = alpha2.SiPhi;
      tracks.TrEvent[IsAlpha[0]].PCZ_p1 = alpha1.PCZ;
      tracks.TrEvent[IsAlpha[1]].PCZ_p2 = alpha2.PCZ;
      Double_t deltaPhi = abs(a1siphi - a2siphi);
      if(deltaPhi>TMath::Pi() && deltaPhi<=2*TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
      Double_t deltaTheta = abs(a1theta-a2theta);

      be8_calc.CalcRecoil_MultiParticle(tracks, IsEject, Be8_1p_2a_qval, Be8_1p_2a);
      li5_calc.CalcRecoil_InvMassSelect(tracks, IsEject, Li5, Li5_qval);
      Double_t exQval = Be8_1p_2a_qval.Ex_recoil;
      Double_t bmWAQval = Be8_1p_2a_qval.BeamKE_WAIntP;
      Double_t ex = Be8_1p_2a.Ex_recoil;
      Double_t bmWA = Be8_1p_2a.BeamKE_WAIntP;

      MyFill("Ex8Be_Qval_rm",600,-20,50,exQval);
      MyFill("Ex8Be_vs_BeamKE_Qval_rm",600,-10,80,bmWAQval,600,-20,50,exQval);
      MyFill("BeamKE_Track_vs_BeamKE_Qval_a1_rm",600,-10,80,bmWAQval,600,-10,80,a1be);
      MyFill("BeamKE_Track_vs_BeamKE_Qval_a2_rm",600,-10,80,bmWAQval,600,-10,80,a2be);
      MyFill("BeamKE_Qvalue_8Be_rm",600,-10,80,bmWAQval);
      MyFill("BeamKE_Qvalue_8Be_cm_rm",180,-1,15,bmWAQval*4/22);
   
      //test spectra; REMOVE POST TEST? Yes
      if((a1detid>-1 && a1detid<4) && (a2detid>-1 && a2detid<4)) {
        MyFill("Ex8Be_Q3_testAnd_rm",600,-20,50,ex);
      }
      if((a1detid>-1 || a1detid<4) || (a2detid>-1 || a2detid<4)) {
        MyFill("Ex8Be_Q3_testOr_rm",600,-20,50,ex);
      }
      if(a1detid>-1 && a1detid<4) {
        MyFill("Ex8Be_Q3_testa1_rm",600,-20,50,ex);
      }
      if(a2detid>-1 && a2detid<4) {
        MyFill("Ex8Be_Q3_testa2_rm",600,-20,50,ex);
      }
      MyFill("BeamKE_Track_vs_BeamKE_a1_rm",600,-10,20,bmWA,600,-10,20,a1be);
      MyFill("BeamKE_Track_vs_BeamKE_a2_rm",600,-10,20,bmWA,600,-10,20,a2be);
      MyFill("Ex8Be_vs_BeamKE_rm",600,-10,80,bmWA,600,-20,50,ex);
      MyFill("Ex8Be_vs_BeamKE_tracked_a1_rm",600,-10,80,a1be,600,-20,50,ex);
      MyFill("Ex8Be_vs_BeamKE_tracked_a1_rm",600,-10,80,a2be,600,-20,50,ex);
      MyFill("Ex8Be_vs_DeltaSiPhi_rm",400,-10,360,deltaPhi*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_DeltaTheta_rm",400,-10,190,deltaTheta*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_PCZ_a1_all_rm",400,-10,30,a1pcz,400,-10,30,ex); 
      MyFill("Ex8Be_vs_PCZ_a2_all_rm",400,-10,30,a2pcz,400,-10,30,ex); 
      MyFill("Ex8Be_vs_SiPhi_a1_rm",400,-10,360,a1siphi*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_SiPhi_a2_rm",400,-10,360,a2siphi*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_PCPhi_a1_rm",400,-10,360,a1pcphi*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_PCPhi_a2_rm",400,-10,360,a2pcphi*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_Theta_a1_rm",400,-10,180,a1theta*rads2deg,400,-10,30,ex); 
      MyFill("Ex8Be_vs_Theta_a2_rm",400,-10,180,a2theta*rads2deg,400,-10,30,ex); 

      Double_t ker = Be8_1p_2a.KE_recoil;
      Double_t theta = Be8_1p_2a.Theta;
      MyFill("8BeKE_vs_Theta_8Be_rm",300,-0.1,180,theta,300,-0.1,80,ker);
    }
  }
  be8_calc.ELoss_eject = NULL;
  be8_calc.ELoss_eject2 = NULL;
  be8_calc.ELoss_beam = NULL;
  li5_calc.ELoss_eject = NULL;
  li5_calc.ELoss_eject2 = NULL;
  li5_calc.ELoss_beam = NULL;
}

/*AlphaScatter()
 *Method for elastic scattering reconstruction
 *Need to check, but my guess is its for a calibration
 *or as a way to normalize for cross sections
 */
void analyzer::ElasticScatter() {

  for(int i=0; i<tracks.NTracks1; i++) {
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
    if(deutCut->IsInside(sie, pce*sin(theta))) {
      MyFill("Elastic_E_si",600,0,35,sie);
      MyFill("Elastic_E_si_vs_Theta",600,0,600,theta*rads2deg,600,0,35,sie);
      MyFill("Elastic_E_si_vs_IntPoint",600,0,60,intp,600,0,35,sie);
      if (sie>0.0 && sie<40.0 && pathL>0.0 && pathL<100.0 && intp>0.0 && intp<ana_length && 
          be>0.0 && be<BeamE && pcz>0.0 && pce>0.0) {
        Float_t E_erxn = deuteron_eloss->GetLookupEnergy(sie, -pathL);
	MyFill("Elastic_E_rxn",600,0,35,E_erxn);
	MyFill("Elastic_E_rxn_vs_Theta",600,0,600,theta*rads2deg,600,0,30,E_erxn);
	MyFill("4He_E_rxn_vs_IntPoint",600,0,50,intp,600,0,30,E_erxn);

	Float_t EBeame = ((m_7be+m_d)*(m_7be+m_d)*E_erxn)/
                          (4.0*m_7be*m_d*cos(theta)*cos(theta));
	// Recoil Kinematics
	Float_t E7be_elas = EBeame - E_erxn;
	MyFill("Elastic_7Be",600,-10,40,E7be_elas);

	MyFill("EjectE_vs_RecoilE_Elastic",600,0,30,E7be_elas,600,0,30,E_erxn);
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

  //for init:  Large values that are within the SRIM file range and the ion has to go far enough 
  //that it stops. Step size should be small
  be7_eloss = new LookUp(be7eloss_name, m_7be);  
  be7_eloss->InitializeLookupTables(30.0,200.0,0.01,0.04);  
  alpha_eloss = new LookUp(he4eloss_name, m_alpha);
  alpha_eloss->InitializeLookupTables(30.0,900.0,0.01,0.04); 
  he3_eloss = new LookUp(he3eloss_name, m_3he);
  he3_eloss->InitializeLookupTables(30.0,1200.0,0.02,0.04);
  proton_eloss = new LookUp(peloss_name, m_p);
  proton_eloss->InitializeLookupTables(30.0, 11900.0, 0.05, 0.01);
  deuteron_eloss = new LookUp(deloss_name, m_d);
  deuteron_eloss->InitializeLookupTables(18.0, 6600.0, 0.02, 0.04);
  
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
  outTree->Branch("Be8_1p", &Be8_1p);
  outTree->Branch("Be8_1p_qval", &Be8_1p_qval);
  outTree->Branch("Be8_1p_2a", &Be8_1p_2a);
  outTree->Branch("Be8_1p_2a_qval", &Be8_1p_2a_qval);
  outTree->Branch("Li6", &Li6);
  outTree->Branch("Li6_qval", &Li6_qval);
  outTree->Branch("Li5", &Li5);
  outTree->Branch("Li5_qval", &Li5_qval);
  outTree->Branch("Be8_1p_any", &Be8_1p_any);
  outTree->Branch("Be8_1p_any_qval", &Be8_1p_any_qval);
  outTree->Branch("RFTime", &RFTime, "RFTime/I");
  outTree->Branch("MCPTime", &MCPTime, "MCPTime/I");
  outTree->Branch("TDC2", &TDC2, "TDC2/I");


  ifstream inputList;
  inputList.open(inputlistname);
  if(!inputList.is_open()){
    cout<<"List of rootfiles could not be opened!"<<endl;
    exit(EXIT_FAILURE);
  }

  if(cutFlag) {
    string protonCutfile, alphaCutfile, he3Cutfile, deutCutfile, joinedAlphaCutfile;
    getline(inputList, protonCutfile);
    getline(inputList, alphaCutfile);
    getline(inputList, he3Cutfile);
    getline(inputList, deutCutfile);
    getline(inputList, joinedAlphaCutfile);
    getCut(protonCutfile, alphaCutfile, he3Cutfile, deutCutfile, joinedAlphaCutfile);
  }
  rootObj->Add(protonCut); rootObj->Add(alphaCut); rootObj->Add(he3Cut);
  rootObj->Add(joinedAlphaCut);

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
    inputTree->SetBranchAddress("RFTime", &input_RFTime);
    inputTree->SetBranchAddress("MCPTime", &input_MCPTime);

    Int_t nentries = inputTree->GetEntries();
    float blentries = nentries;
    cout<<"Number of entries: "<<nentries<<endl;

    for (int entry = 0; entry<nentries; entry++) {
      cout<<"\rPercent of file completed: "<<entry/blentries*100.0<<"%  "<<flush;
      inputTree->GetEvent(entry);

      SiEnergy_vec = vector<vector<Double_t>>(28, vector<Double_t>(0));
      tracks.ZeroTrack();
      recoilReset(Li6);
      recoilReset(Li6_qval);
      recoilReset(Li5);
      recoilReset(Li5_qval);
      recoilReset(Be8_1p);
      recoilReset(Be8_1p_qval);
      recoilReset(Be8_1p_2a);
      recoilReset(Be8_1p_2a_qval);
      recoilReset(Be8_1p_any);
      recoilReset(Be8_1p_any_qval);
    
      if(MCP_RF()){        
        Track1();
        Track2();
        Track3();
      
        if(PCPlots) PCPlotting();
      
        TrackCalc();

        if(PCWireCal) PCWireCalibration();
        if(FillEdE_cor) EdEcor();
        if(PCPlots_2) PCPlotting2();
        if(ResidualEnergy_Calc) {
          CalcRecoilE();
          //CalcSmAngleAlphas();
        }
        if(Reconstruction_Session) ReconstructMe();
        if(Elastic_Scat_Alpha) ElasticScatter();

        if(FillTree) outTree->Fill();
      }
    }//end event loop
    cout<<endl;
  }//end file list loop
  inputList.close();
  outFile->cd();
  outTree->Write(outTree->GetName(), TObject::kOverwrite);
  rootObj->Write();
  outFile->Close();
  cout<<"All Root objects are written to: "<<outputName<<endl;
  cout<<outputName<<" is closed"<<endl;
}
