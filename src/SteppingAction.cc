//
/// \file GdScint/src/SteppingAction.cc
/// \brief Implementation of the GdScint::SteppingAction class

#include "SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "globals.hh"
#include "G4StepPoint.hh"
#include "G4ReactionProduct.hh"


//For outputting to ascii:
#include <iomanip>
#include <fstream>

SteppingAction::SteppingAction()
{}

SteppingAction::~SteppingAction()
{}

G4VProcess* SteppingAction::GetCurrentProcess()
{
  return fpSteppingManager->GetfCurrentProcess();
}
 void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

}
