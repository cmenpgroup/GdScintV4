//
/// \file GdScint/src/EventAction.hh
/// \brief Definition of the GdScint::EventAction class
#ifndef GdScintEventAction_H
#define GdScintEventAction_H 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


#include "G4THitsMap.hh"


//2017: follwing have been changed:
#ifdef ANALYSIS_USE
 #include "TROOT.h"
 #include "TFile.h"
 #include "TNtuple.h"
 #include "TTree.h"
 #include "TH1D.h"
#endif

class GdScintEventMessenger;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void WriteRecoilToFile(G4int);
    void SetThreshold(G4int);
    void SetFileName(const char*);

  public:
         G4double      particleke;

  private:
	 void writeHits(const G4Event*);
  G4int                       GdScintCollID;
  G4int n_hits;
 // G4int GdScintCollection;
  G4int totEvents; //total number of events processed
  G4int fThreshold;

  GdScintEventMessenger* fGdScintEvMess;



  G4THitsMap<G4double>* GetHitsCollection(G4int hcID, const G4Event* evt) const;
};

#endif








