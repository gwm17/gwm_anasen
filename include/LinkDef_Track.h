#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

#pragma link C++ struct TrackEvent+;
#pragma link C++ struct RecoilEvent+;
#pragma link C++ struct Silicon_Event+;
#pragma link C++ struct PropCounter_Event+;
#pragma link C++ struct CsI_Event+;
#pragma link C++ struct SiDetector+;
#pragma link C++ struct sievent+;
#pragma link C++ struct pcevent+;
#pragma link C++ struct csievent;
#pragma link C++ struct CsIHit::SortByCsI+;
#pragma link C++ struct SiHit::SortByDetector+;
#pragma link C++ struct SiHit::SortByHit+;
#pragma link C++ struct PCHit::SortByPC+;

#pragma link C++ class std::vector<TrackEvent>+;
#pragma link C++ class std::vector<Silicon_Event>+;
#pragma link C++ class std::vector<PropCounter_Event>+;
#pragma link C++ class std::vector<CsI_Event>+;
#pragma link C++ class std::vector<SiHit::SortByHit>+;
#pragma link C++ class std::vector<SiHit::SortByDetector>+;
#pragma link C++ class std::vector<PCHit::SortByPC>+;
#pragma link C++ class std::vector<CsIHit::SortByCsI>+;

#pragma link C++ defined_in "Track.h";
#pragma link C++ defined_in "PCHit.h";
#pragma link C++ defined_in "SiHit.h";
#pragma link C++ defined_in "CsIHit.h";

#endif
