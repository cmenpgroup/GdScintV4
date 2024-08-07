#ifndef GdScintSD_h
#define GdScintSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "GdScintHit.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class DetectorConstruction;
class G4Step;
class G4HCofThisEvent;

class GdScintSD : public G4VSensitiveDetector {

   public:
//      GdScintSD(G4String);
     GdScintSD(G4String, DetectorConstruction*);
     ~GdScintSD();
  
     void Initialize(G4HCofThisEvent*);
     G4bool ProcessHits(G4Step*,G4TouchableHistory*);
     void EndOfEvent(G4HCofThisEvent*);
     void clear();
     void DrawScint();
     void PrintScint();
  
  private:
  
     GdScintHitsCollection* GdScintCollection;
     DetectorConstruction* GdDetector;
     G4int HitID;
    
};

#endif
