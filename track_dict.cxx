// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME track_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "lib/SiHit.h"
#include "lib/PCHit.h"
#include "lib/Track.h"
#include "lib/CsIHit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TrackEvent_Dictionary();
   static void TrackEvent_TClassManip(TClass*);
   static void *new_TrackEvent(void *p = 0);
   static void *newArray_TrackEvent(Long_t size, void *p);
   static void delete_TrackEvent(void *p);
   static void deleteArray_TrackEvent(void *p);
   static void destruct_TrackEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TrackEvent*)
   {
      ::TrackEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TrackEvent));
      static ::ROOT::TGenericClassInfo 
         instance("TrackEvent", "lib/Track.h", 21,
                  typeid(::TrackEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TrackEvent_Dictionary, isa_proxy, 4,
                  sizeof(::TrackEvent) );
      instance.SetNew(&new_TrackEvent);
      instance.SetNewArray(&newArray_TrackEvent);
      instance.SetDelete(&delete_TrackEvent);
      instance.SetDeleteArray(&deleteArray_TrackEvent);
      instance.SetDestructor(&destruct_TrackEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TrackEvent*)
   {
      return GenerateInitInstanceLocal((::TrackEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TrackEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TrackEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TrackEvent*)0x0)->GetClass();
      TrackEvent_TClassManip(theClass);
   return theClass;
   }

   static void TrackEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RecoilEvent_Dictionary();
   static void RecoilEvent_TClassManip(TClass*);
   static void *new_RecoilEvent(void *p = 0);
   static void *newArray_RecoilEvent(Long_t size, void *p);
   static void delete_RecoilEvent(void *p);
   static void deleteArray_RecoilEvent(void *p);
   static void destruct_RecoilEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RecoilEvent*)
   {
      ::RecoilEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RecoilEvent));
      static ::ROOT::TGenericClassInfo 
         instance("RecoilEvent", "lib/Track.h", 138,
                  typeid(::RecoilEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RecoilEvent_Dictionary, isa_proxy, 4,
                  sizeof(::RecoilEvent) );
      instance.SetNew(&new_RecoilEvent);
      instance.SetNewArray(&newArray_RecoilEvent);
      instance.SetDelete(&delete_RecoilEvent);
      instance.SetDeleteArray(&deleteArray_RecoilEvent);
      instance.SetDestructor(&destruct_RecoilEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RecoilEvent*)
   {
      return GenerateInitInstanceLocal((::RecoilEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RecoilEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RecoilEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RecoilEvent*)0x0)->GetClass();
      RecoilEvent_TClassManip(theClass);
   return theClass;
   }

   static void RecoilEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_TrackEvent(void *p) {
      return  p ? new(p) ::TrackEvent : new ::TrackEvent;
   }
   static void *newArray_TrackEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TrackEvent[nElements] : new ::TrackEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TrackEvent(void *p) {
      delete ((::TrackEvent*)p);
   }
   static void deleteArray_TrackEvent(void *p) {
      delete [] ((::TrackEvent*)p);
   }
   static void destruct_TrackEvent(void *p) {
      typedef ::TrackEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TrackEvent

namespace ROOT {
   // Wrappers around operator new
   static void *new_RecoilEvent(void *p) {
      return  p ? new(p) ::RecoilEvent : new ::RecoilEvent;
   }
   static void *newArray_RecoilEvent(Long_t nElements, void *p) {
      return p ? new(p) ::RecoilEvent[nElements] : new ::RecoilEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_RecoilEvent(void *p) {
      delete ((::RecoilEvent*)p);
   }
   static void deleteArray_RecoilEvent(void *p) {
      delete [] ((::RecoilEvent*)p);
   }
   static void destruct_RecoilEvent(void *p) {
      typedef ::RecoilEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RecoilEvent

namespace {
  void TriggerDictionaryInitialization_track_dict_Impl() {
    static const char* headers[] = {
"lib/SiHit.h",
"lib/PCHit.h",
"lib/Track.h",
"lib/CsIHit.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/data1/gwm17/7Be_d/new/gwm_analysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "track_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$lib/Track.h")))  TrackEvent;
struct __attribute__((annotate("$clingAutoload$lib/Track.h")))  RecoilEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "track_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "lib/SiHit.h"
#include "lib/PCHit.h"
#include "lib/Track.h"
#include "lib/CsIHit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RecoilEvent", payloadCode, "@",
"TrackEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("track_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_track_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_track_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_track_dict() {
  TriggerDictionaryInitialization_track_dict_Impl();
}
