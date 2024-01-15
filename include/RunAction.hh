//
/// \file GdScint/src/RunAction.hh
/// \brief Definition of the GdScint::RunAction class
#ifndef GdScintRunAction_h
#define GdScintRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
  ~RunAction();  // is virtual in B4

  public:
    void BeginOfRunAction(const G4Run*);   // is virtual in B4
    void EndOfRunAction(const G4Run*);    // is virtual in B4
};
#endif
