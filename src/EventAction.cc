
#ifdef G4UI_USE_ROOT
  #include "G4UIRoot.hh"
#endif

#include "SteppingAction.hh"



#include "EventAction.hh"

#include "GdScintHit.hh"
#include "RunAction.hh" 
#include "GdScintEventMessenger.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Run.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "GdScintSD.hh"

#include <iomanip>
#include <fstream>


#include "G4THitsMap.hh"
    G4int verboseLevel= 1;  
    G4int dataFiles = 1;    
    G4int Ntuples = 0;      
EventAction::EventAction()
{
    GdScintCollID=-1;
    totEvents=0;
    fGdScintEvMess=new GdScintEventMessenger(this);    
}

EventAction::~EventAction()
{
    delete fGdScintEvMess;
}
//
G4THitsMap<G4double>*
EventAction::GetHitsCollection(G4int hcID, const G4Event* evt) const
{
            //
    auto GdScintHitsCollection = static_cast<G4THitsMap<G4double>*>(evt->GetHCofThisEvent()->GetHC(hcID));
    if(! GdScintHitsCollection) { G4cout<<"yikes"<<G4endl;
    }
    return GdScintHitsCollection;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{

 //      totEvents = 0;
	G4int evtNb;
	evtNb = evt->GetEventID();
	G4cout<<"!!!!!!!!!!! Begin of event: "<<evtNb<<G4endl;

}

void EventAction::EndOfEventAction(const G4Event* evt)
{ 
	G4cout << "Total Events: " << G4endl;
    G4cout << totEvents << G4endl;
    totEvents++;
    
    G4int evtNb = evt->GetEventID();
        
    GdScintHitsCollection*          ScintHC    = NULL;
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    
    G4int           n_hit = 0;
    G4double        ScintTotalEnergy = 0;
    G4ThreeVector   ScintPosition;
    G4String        ScintParticle;
    G4double        ScintEdep = 0;
    G4String        ScintParticleType;
    G4ThreeVector   ScintMomentum;
    G4ThreeVector   ScintMomentumDir;
    G4String        ScintVertexVolumeName;
    G4ThreeVector   ScintVertexPosition;
    G4ThreeVector   ScintVertexMomDir;
    G4String        ScintNuclearProcess;

if (HCE) {
	ScintHC = (GdScintHitsCollection*)(HCE->GetHC(GdScintCollID));
    if ( ScintHC )
    { G4cout<<"ScintHC :)"<<G4endl;
       // n_hit = ScintHC->entries(); 
	//G4cout<<"n_hit: "<<n_hit<<G4endl;
       // if( verboseLevel > 1 ) {
        //	G4cout<<"Total Number Of Step on this Hit "<<n_hit<<G4endl;
       // }

       // for (G4int i=0;i<n_hit;i++)
       for (G4int i=0;i++;)
        {
            ScintEdep              = (*ScintHC)[i]-> GetScintEdep(); 
            ScintTotalEnergy += ScintEdep;        
       
	} 

        G4cout<<"Edep in scintillator = "<<ScintTotalEnergy/MeV<<" MeV"<<G4endl;

		std::ofstream theFile6("scintillator.dat", std::ios::app);
		theFile6
		<<evtNb<<"\t"
		<<ScintTotalEnergy/MeV<<G4endl;
    }
    if ( ! ScintHC)
    { G4cout<<"no scint HC :("<<G4endl;
 return;	    
    }
    } 
if(! HCE) { G4cout<<"No HCE :("<<G4endl; }
    if ( G4VVisManager::GetConcreteInstance() )
    {
        G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
        G4int n_trajectories = 0;

        if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
        for (G4int i=0; i<n_trajectories; i++)
        { G4Trajectory* trj = (G4Trajectory*)
	                                ((*(evt->GetTrajectoryContainer()))[i]);
         trj->DrawTrajectory();
        }
    }

}
void EventAction::SetThreshold(G4int th){}

void EventAction::SetFileName(char const*){}
