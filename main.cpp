#include "gwmAnalyzer.h"
#include <TROOT.h>
#include <TApplication.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  TApplication *app = new TApplication("app",&argc, argv);
  analyzer pleasework;
  /*set the analyzer flags here... options are pcplots, pcplots2, beamandeloss, filltree,
   *filledecor, checkbasic, mcprfcut, readpcwire, pcwirecal, residualenergycalc, 
   *residualenergycalc20ne, reconstructionsession, and elasticalpha. See readme for details 
   *on each flag
   */
  pleasework.SetFlag("filledecor");
  pleasework.SetFlag("mcprfcut");
  pleasework.SetFlag("checkbasic");
  pleasework.SetFlag("pcplots");
  //pleasework.SetFlag("pcplots2");
  pleasework.SetFlag("filltree");
  pleasework.run();
}
