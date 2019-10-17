#include "dalitz.h"
#include "angDist.h"
#include "singleChan.h"
#include "SFactor.h"
#include  <string>
#include <iostream>

using namespace std;

struct options {
  int dalitz = 0;
  int angcor = 0;
  int he3 = 0;
  int sf = 0;
  int counter = 0;
} options;

static const char *optString = "adhs";

int main(int argc, char **argv) {
  if(argc >= 5) {
    int opt = 0;
    while( (opt = getopt(argc,argv,optString)) != -1) {
      switch(opt) {
        case 'a':
          options.angcor = 1;
          options.counter++;
          break;
        case 'd':
          options.dalitz = 1;
          options.counter++;
          break;
        case 'h':
          options.he3 = 1;
          options.counter++;
          break;
        case 's':
          options.sf = 1;
          options.counter++;
          break;
      }
    }
    if(options.counter == 4) {
      cout<<"Performing data analysis..."<<endl;
      DalAnalyzer dal_it(0.0,10.0,0.2);
      dal_it.run(argv[2],argv[3],argv[4]);
      cout<<"Dalitz analysis completed."<<endl;
      angDist fix_it(argv[2],argv[3],argv[4],1);
      fix_it.run();
      cout<<"Angular distribution analysis completed."<<endl;
      singleChan he_it(0.2,10.2,0.2,1);
      he_it.run(argv[2],argv[3],argv[4]);
      cout<<"Single channel analysis completed."<<endl;
      SFactor s_it(1,4);
      s_it.Run(argv[4]); 
      cout<<"S-factor analysis completed."<<endl;
    } else if(options.angcor && options.counter == 1) {
      cout<<"Perfoming angular distribution correction..."<<endl;
      angDist fix_it(argv[2],argv[3],argv[4],0);
      fix_it.run();
      cout<<"Completed"<<endl;
    } else if(options.dalitz && options.counter == 1) {
      cout<<"Peforming Dalitz plot cross section..."<<endl;
      DalAnalyzer dal_it(0.0,10.2,0.2);
      dal_it.run(argv[2],argv[3],argv[4]);
      cout<<"Completed"<<endl;
    } else if(options.he3 && options.counter == 1) {
      cout<<"Peforming 3He cross section..."<<endl;
      singleChan he_it(0.2,10.2,0.2,0);
      he_it.run(argv[2],argv[3],argv[4]);
      cout<<"Completed"<<endl;
    } else if(options.sf && options.counter == 1) {
      cout<<"Performing S-factor calculation..."<<endl;
      SFactor s_it(1,4);
      s_it.Run(argv[4]);
      cout<<"Completed"<<endl;
    } else {
      cout<<"No method given... Terminating."<<endl;
    }
  } else {
    cout<<"Incorrect number of arguments!! Expects 4"<<endl;
  }
}
