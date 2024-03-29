#ifndef GdScintEventMessenger_h
#define GdScintEventMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class EventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class GdScintEventMessenger: public G4UImessenger
{
  public:
    GdScintEventMessenger(EventAction* );
   ~GdScintEventMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  EventAction* GdEvAct;
   
    G4UIdirectory*             GdEADir;
    G4UIcmdWithAString*        SetFileCmd;
    G4UIcmdWithAnInteger*        SetThresholdCmd;
};
#endif
