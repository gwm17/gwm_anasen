#ifndef BAD_DETECTOR_H
#define BAD_DETECTOR_H

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "SiHit.h"
#include "CsIHit.h"
#include "PCHit.h"
#include "Track.h"

using namespace std;

class SiSearch {
  
  public:
    SiSearch();
    void run();

  private:
    void writeBadList();
    void chanFound(Int_t detID, Int_t chanNum);
    int isGood(Int_t detID, Int_t chanNum);
    void getDataList();
    void makeBadList();
    unordered_map<string, int> det_map;
    vector<string> data_list;
    SiHit Si;
};


#endif
