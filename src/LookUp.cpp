/*LookUp.cpp
 *Class designed to take in srim files and calculate the energy loss of a particle 
 *moving through gas a certain distance. Intended for use in ANASEN analysis. 
 *
 *N. Rijal -- June 2016
 *Facelift by G.M. April 2019
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <TGraph.h> 

#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>

#include <string>
#include <math.h>   
#include <TRandom3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TCanvas.h>
#include <TMath.h> 
#include <TROOT.h> 
#include <vector> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>
#include "LookUp.h"

using namespace std;

LookUp::LookUp() {
  c = 29.9792458;           // Speed of light in cm/ns.
  dEdx_e = 0;
  dEdx_n = 0;
  Energy_in_range = 1;
  EvD = new TGraph();
  GoodELossFile = 0;
  IonEnergy = 0;
  IonMass = 0;
  last_point = 0;
  points = 0;
  last_point1 = 0;
  points1 = 0;
  EtoDtab = 0;
  DtoEtab = 0;
  dEdx_e = 0;
  dEdx_n = 0;
}
LookUp::LookUp(string Eloss_file, Double_t InputMass) {

  string aux;
  ifstream Read(Eloss_file.c_str());
  last_point = 0;

  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: File " << Eloss_file << " was not found." << endl;
    GoodELossFile = 0;
  } else {
    GoodELossFile = 1;        
    do{
      Read >> IonEnergy >> dEdx_e >> dEdx_n ;
      IonEnergy_v.push_back(IonEnergy);
      dEdx_e_v.push_back(dEdx_e);
      dEdx_n_v.push_back(dEdx_n);
    }
    while(!Read.eof());
    Read.close();

    points = IonEnergy_v.size();
    Energy_in_range = 1;
    IonMass = InputMass;          // In MeV/c^2
    c = 29.9792458;                   // Speed of light in cm/ns.
    EvD = new TGraph();
  }
}
LookUp::~LookUp(){
}
/////////////////////////////////// SPLINE INTERPOLATION /////////////////////////////////////////
double LookUp::GetEnergyLoss(double energy /*MeV*/, double distance /*cm*/)
{
  Float_t a11=0.0, a12=0.0, a21=0.0, a22=0.0, a23=0.0, a32=0.0, a33=0.0;
  Float_t b11=0.0, b22=0.0, b33=0.0;
  Float_t a1=0.0,  b1=0.0;
  Float_t K0=0.0, K1=0.0; 
  Float_t N1=0.0, N2=0.0, N3=0.0;
  Float_t T1=0.0,q1=0.0;

  int i = -1;
  if(energy < 0.01) return(0);
  
  
  // Look for two points for which the initial energy lies in between.  
  //This for-loop should find the points 
  //unless there was a big jump from the energy used in the 
  //last point and the energy used now.
  
  for(int p=0; p<points-1; p++){
    if(energy>=IonEnergy_v[p]  && energy<IonEnergy_v[p+1]){
      i = p+1;
      last_point = p;
      break;
    }
  }
  // If after this two loop i is still -1 it means the energy was out of range.
  
  if(i==-1){
    cout << "*** EnergyLoss Error: energy not within range: " << energy << 
            " " << IonEnergy_v[0] << endl;
    Energy_in_range = 0;
    return 0;
  } 
  
  //Ion Energy
  double x0=IonEnergy_v[i-1];
  double x1=IonEnergy_v[i];
  double x2=IonEnergy_v[i+1];
  
  //Total Energy Loss (electric + nuclear) for one step
  double y0=dEdx_e_v[i-1]+dEdx_n_v[i-1];
  double y1=dEdx_e_v[i]+dEdx_n_v[i];
  double y2=dEdx_e_v[i+1]+dEdx_n_v[i+1];
  
  a11=2/(x1-x0);
  a12=1/(x1-x0);
  a21=1/(x1-x0);
  a22=2*((1/(x1-x0))+(1/(x2-x1)));
  a23=1/(x2-x1);
  a32=1/(x2-x1);
  a33=2/(x2-x1);
  
  b11=3*((y1-y0)/((x1-x0)*(x1-x0)));
  b22=3*(((y1-y0)/((x1-x0)*(x1-x0)))+((y2-y1)/((x2-x1)*(x2-x1))));
  b33=3*((y2-y1)/((x2-x1)*(x2-x1)));  
  
  //mathematical terms to calculate curvatures.
  N1=(a21*a33*a12-a11*(a22*a33-a23*a32))/(a33*a12);
  N2=(b22*a33-a23*b33)/a33;
  N3=b11*(a22*a33-a23*a32)/(a33*a12);
  
  //curvatures
  K0=(N2-N3)/N1;
  K1=(b11-a11*K0)/a12;
  a1=K0*(x1-x0)-(y1-y0);  
  b1=-K1*(x1-x0)+(y1-y0);  
  T1=(energy-x0)/(x1-x0);
  
  //polynomials //which gives the value of energy loss for given energy.
  q1=(1-T1)*y0+T1*y1+T1*(1-T1)*(a1*(1-T1)+b1*T1);
  
  return (q1*10*distance);
}

Double_t LookUp::GetInitialEnergy(Double_t FinalEnergy /*MeV*/, Double_t PathLength /*cm*/,
                                  Double_t StepSize/*cm*/)
{
  Double_t Energy = FinalEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  last_point = 0;

  // The function starts by assuming FinalEnergy is within the energy range, 
  //but this could be changed in the GetEnergyLoss() function.

  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {
    Energy = Energy + GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range){
      break;
    } 
  }

  Energy = Energy + GetEnergyLoss(Energy,PathLength-Steps*StepSize);

  if (!Energy_in_range)
    Energy = -1000.0; // Return an unrealistic value. 

  return Energy;
}

Double_t LookUp::GetFinalEnergy(Double_t InitialEnergy /*MeV*/, Double_t PathLength /*cm*/,
                                Double_t StepSize/*cm*/)
{

  Double_t Energy = InitialEnergy;
  int Steps = (int)floor(PathLength/StepSize);

  // The function starts by assuming InitialEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.

  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {
    Energy = Energy - GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range){
      break;
    }
  }  

  Energy = Energy - GetEnergyLoss(Energy,PathLength-Steps*StepSize);

  if (!Energy_in_range){
    Energy = -1000.0;
  }
  
  return Energy;

}

Double_t LookUp::GetDistance(Double_t InitialE, Double_t FinalE, Double_t StepSize)
{
  
  Double_t dist = 0;
  Double_t E = 0, Elast=0;

  E = InitialE;
  
  while(E>FinalE){
    dist += StepSize;
    Elast=E;
    E = E - GetEnergyLoss(E,StepSize);
  }

  return ((dist-StepSize)-(StepSize*(Elast-FinalE)/(E-Elast)));
}

// Calulates the ion's path length in cm.
Double_t LookUp::GetPathLength(Float_t InitialEnergy /*MeV*/, Float_t FinalEnergy /*MeV*/, Float_t DeltaT /*ns*/)
{
  Double_t L = 0, DeltaX = 0;
  Double_t Kn = InitialEnergy;
  Double_t Kn1 = InitialEnergy;
  Int_t n=0;


  if (IonMass==0)
    cout << "*** EnergyLoss Error: Path length cannot be calculated for IonMass = 0." << endl;
  else {

    // The path length (L) is proportional to sqrt(Kn). 
    //After the sum, L will be multiplied by the proportionality factor.

    while (Kn > FinalEnergy && Kn1 > FinalEnergy && n < (int)pow(10.0,6)){
      L += sqrt((Kn+Kn1)/2)*sqrt(2/IonMass)*DeltaT*c;  // DeltaL going from point n to n+1.
      DeltaX = sqrt((Kn+Kn1)/IonMass)*DeltaT*c;
      Kn1 = Kn;
      Kn -= GetEnergyLoss((Kn+Kn1)/2, DeltaX);// After L is incremented the kinetic energy at n+1 is calculated.
      n++;    
    }
    if (n>=(int)pow(10.0,6)) {
      cout << "*** EnergyLoss Warning: Full path length wasn't reached after 10^6 iterations." 
           << endl;
      L = 0;
    } else L *= 1.0;
  }
  return L;
}

// Calulates the ion's time of flight in ns.
Double_t LookUp::GetTimeOfFlight(Double_t InitialEnergy, Double_t PathLength, Double_t StepSize)
{
  Double_t TOF = 0;
  Double_t Kn = InitialEnergy;
  int Steps = (int)(PathLength/StepSize);
  
  if (IonMass==0) {
    cout << "Error: Time of flight cannot be calculated because mass is zero." << endl;
  }

  else {

    for (int n=0; n<Steps; n++) {

      TOF += sqrt(IonMass/(2*Kn))*StepSize/c; // DeltaT going from point n to n+1.
      Kn -= GetEnergyLoss(Kn, StepSize); // After the TOF is added the K.E. at point n+1 is calc 
    }
    return TOF;
  }
  return(0);
}

void LookUp::SetIonMass(Double_t InputMass)
{
  IonMass = InputMass;
  return;
}

void LookUp::InitializeLookupTables(Double_t MaxE, Double_t MaxDist, 
				    Double_t DelE, Double_t DelD){

  MaximumEnergy = MaxE;
  MaximumDistance = MaxDist;
  DeltaD = DelD;
  DeltaE = DelE; 

  int noE = (int)ceil(MaximumEnergy / DeltaE );
  int noD = (int)ceil(MaximumDistance / DeltaD );
  
  EtoDtab_v.resize(noE); DtoEtab_v.resize(noD);

  DtoEtab_v[0] = MaximumEnergy;
  cout << " Number of distance entries " << noD << endl;
  
  for (int i=1; i<noD; i++){
    DtoEtab_v[i] = GetFinalEnergy(DtoEtab_v[i-1],DeltaD,0.05*DeltaD);
  }
  
  cout << " Number of Energy entries " << noE << endl;
  EtoDtab_v[0] = 0.;
  for (int j=1;j<noE;j++){

    EtoDtab_v[j] = EtoDtab_v[j-1] + GetDistance((MaximumEnergy-(j-1)*DeltaE),
                                                (MaximumEnergy-(j)*DeltaE), (0.05*DeltaD)); 
  }
}

void LookUp::PrintLookupTables(){
  int noE = (int)ceil(MaximumEnergy / DeltaE );
  int noD = (int)ceil(MaximumDistance / DeltaD );
  cout << "Maximum Energy = " << MaximumEnergy << " " << DeltaE << noE << endl;
  for (int i=0;i<noE;i++){
    cout << "E,D= "<< MaximumEnergy - DeltaE*i << "," <<EtoDtab_v[i] <<endl; 
  }
  cout << "Maximum Distance = "<< MaximumDistance << " " << DeltaD << noD << endl;
  for (int i=0;i<noD;i++){
    cout << "D,E= "<< DeltaD*i << "," <<DtoEtab_v[i]<<endl; 
  }
}

Double_t LookUp::GetLookupEnergy(Double_t InitialEnergy, Double_t distance){

  Double_t D1,D2,D;
  Double_t E1,E2,E;
  int index;

  if (InitialEnergy<0 || InitialEnergy>MaximumEnergy){
    return(-1.);
  }

  // Find the distance for which the initial energy is matched, interpolating 
  index = (int)floor((MaximumEnergy-InitialEnergy) / DeltaE) ;
  D1=EtoDtab_v[index]; D2=EtoDtab_v[index+1];
  E1=MaximumEnergy-index*DeltaE; E2=MaximumEnergy-(index+1)*DeltaE; 

  D = ( InitialEnergy - E1 ) / ( E2 - E1 ) * (D2-D1) + D1 ;

  //Still in the table ?
  if(((D+distance <=0)) || ((D+distance)> MaximumDistance)){
    cout<<"i m here"<<(distance)<<" "<<D<<" Energy: "<<InitialEnergy<<endl;
    return 0.;
  }

  // Lookup what energy is reached for (D + distance) and interpolate
  index = (int)floor((D + distance) / DeltaD );

  E1=DtoEtab_v[index]; E2=DtoEtab_v[index+1];
  D1=index*DeltaD;D2=(index+1)*DeltaD;

  E = ( (D+distance) - D1 ) / ( D2 - D1 ) * (E2-E1) + E1 ;

  return(E);
}

