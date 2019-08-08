#include "BadDetector.h"

SiSearch::SiSearch() {
  for (int i=0; i<4; i++) {
    for(int j=0; j<32; j++) {
      string q3 = to_string(i)+"\t"+to_string(j);
      pair<string, int> entry(q3, 0);
      det_map.insert(entry);
    }
  }
  for(int i=6; i<30; i++) {
    for(int j=0; j<8; j++) {
      string sx3 = to_string(i) +"\t"+ to_string(j);
      pair<string, int> entry(sx3, 0);
      det_map.insert(entry);
    }
  }
}

int SiSearch::isGood(Int_t detID, Int_t chanNum) {
  string det = to_string(detID)+"\t"+to_string(chanNum);
  auto entry = det_map.find(det);
  if(entry->second) return 1;
  else return 0;
}

void SiSearch::chanFound(Int_t detID, Int_t chanNum) {
  string det = to_string(detID)+"\t"+to_string(chanNum);
  det_map[det] = 1;
}

void SiSearch::getDataList() {
  ifstream input("DataList.txt");
  if(input.is_open()) {
    string junk;
    string file;
    input>>junk>>junk>>junk>>junk>>junk;
    while(input>>file) {
      data_list.push_back(file);
    } 
  } else {
    cout<<"Unable to open data list!"<<endl;
  }
  input.close();
}

void SiSearch::makeBadList() {
  ofstream output("badDet/badDetectors.txt");
  output<<"DetID"<<"\t"<<"ChanNum"<<endl;
  for (auto entry:det_map) {
    if(entry.second == 0) {
      output<<entry.first<<endl;
    }
  }
  output.close();
}

void SiSearch::run() {
  getDataList();
  for(unsigned int i=0; i<data_list.size(); i++) {
    string fname = data_list[i];
    char filename[fname.size()+1];
    strcpy(filename, fname.c_str());
    cout<<"filename: "<<filename<<endl;
    TFile *file = new TFile(filename, "READ");
    if(!file->IsOpen()) {
      cout<<"blerg"<<endl;
    }
    TTree *tree = (TTree*) file->Get("MainTree");
    
    tree->SetBranchAddress("Si.Hit", &Si.ReadHit);
    tree->SetBranchAddress("Si.Detector", &Si.ReadDet);
    int nentries = tree->GetEntries();
    for (int e=0; e<nentries; e++) {
      cout<<"\rPercent processed: "<<e/((float)nentries)*100<<"% "<<flush;
      tree->GetEvent(e);
      for (unsigned int hit=0; hit<Si.ReadHit->size(); hit++) {
        sievent sihit = (*Si.ReadHit)[hit];
        Int_t detid = sihit.DetID;
        if(detid<4) {
          chanFound(detid, sihit.FrontChannel);
          chanFound(detid, sihit.BackChannel+16);
        } else {
          chanFound(detid+2, sihit.FrontChannel);
          chanFound(detid+2, sihit.BackChannel+4);
        }
      }
    }
    cout<<endl;
    file->Close();
  }
  makeBadList();
}

int main() {
  SiSearch generator;
  generator.run();
}
