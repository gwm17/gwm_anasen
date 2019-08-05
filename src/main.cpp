/*main.cpp
 *
 *Runs the ANASEN analysis
 *This is where the user should set all of the option flags
 *that they need for their analysis
 *
 * Gordon M. -- April 2019
 */

#include "gwmAnalyzer.h"
#include <TROOT.h>
#include <TApplication.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  TApplication *app = new TApplication("app",&argc, argv);
  analyzer pleasework;
  /*set the analyzer flags here... options are pcplots, pcplots2, beamandeloss, filltree,
   *filledecor, mcprfcut, readpcwire, pcwirecal, residualenergycalc, 
   *reconstructionsession, and elasticscatalpha. See readme for details 
   *on each flag
   */
  pleasework.SetFlag("filledecor");
  pleasework.SetFlag("mcprfcut");
  //pleasework.SetFlag("pcplots");
  pleasework.SetFlag("beamandeloss");
  //pleasework.SetFlag("pcplots2");
  pleasework.SetFlag("cutflag");
  pleasework.SetFlag("residualenergycalc");
  pleasework.SetFlag("reconstructionsession");
  //pleasework.SetFlag("pcwirecal");
  pleasework.SetFlag("filltree");
  //pleasework.SetFlag("elasticscatalpha");
  pleasework.run();
}
