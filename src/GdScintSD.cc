#include "GdScintSD.hh"
#include "DetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4HEVector.hh"
#include "G4String.hh"
#include "EventAction.hh"
#include "GdScintHit.hh"
GdScintSD::GdScintSD(G4String name, DetectorConstruction* GdSD)  	
: G4VSensitiveDetector(name),GdDetector(GdSD) {

       	G4String HCname="GdScintCollection";
  collectionName.insert(HCname);
} 
GdScintSD::~GdScintSD()
{ ; }

void GdScintSD::Initialize(G4HCofThisEvent* HCE)
{

G4cout<<"initialize sd"<<G4endl;	
 GdScintCollection = new GdScintHitsCollection(GetName(), collectionName[0]); 
  HitID = -1;
  verboseLevel=2;
/*    
    static G4int HCID = -1;
    if(HCID<0){ HCID = G4SDManager::GetSDMpointer()->GetCollectionID(GdScintCollection); }
//if(HCID<0) HCID = GetCollectionID(0);
//GdScintCollection = new GdScintHitsCollection("myGdScintSD", collectionName[0]);
    HCE->AddHitsCollection(HCID,GdScintCollection);
    G4cout<<"HCID: "<<HCID<<G4endl;
    G4int nHits = GdScintCollection->entries();
    if (verboseLevel>=1) {
        G4cout << "     Number of Scint hits: " << nHits << G4endl;
    }
*/
}

G4bool GdScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) //R0hist
{
       G4cout<<"HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"<<G4endl;	
////Get the Scint hit
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
    G4VPhysicalVolume* physVol = theTouchable->GetVolume();
    G4int ScintNumber = 0;
    ScintNumber= physVol->GetCopyNo();
	G4double id = G4Threading::G4GetThreadId();
	G4double edep = aStep->GetTotalEnergyDeposit();
        //tally->AccumEdep(id, edep);	
//	if (edep==0) return false;
	
		
    GdScintHit* aScintHit = new GdScintHit();

    aScintHit->SetScintPosition(aStep->GetPostStepPoint()->GetPosition());
    aScintHit->SetScintTime(aStep->GetPostStepPoint()->GetGlobalTime());
    aScintHit->SetScintsTime(aStep->GetPreStepPoint()->GetGlobalTime());
    aScintHit->SetScintParticle(aStep->GetTrack()->GetDefinition()->GetParticleName());
    aScintHit->SetScintEnergy(aStep->GetTrack()->GetKineticEnergy());
    aScintHit->SetScintEdep(aStep->GetTotalEnergyDeposit());
    aScintHit->SetScintParticleType(aStep->GetTrack()->GetDefinition()->GetParticleType());
    aScintHit->SetScintTrackLength(aStep->GetStepLength());
    aScintHit->SetScintMomentum(aStep->GetTrack()->GetMomentum());
    aScintHit->SetScintMomentumDir(aStep->GetTrack()->GetMomentumDirection());
    aScintHit->SetScintVertexVolName(aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName());
    aScintHit->SetScintTrackID(aStep->GetTrack()->GetTrackID());
    aScintHit->SetScintParentID(aStep->GetTrack()->GetParentID());
    aScintHit->SetScintKineticEnergy(aStep->GetTrack()->GetKineticEnergy());
    aScintHit->SetScintCurrentStepNumber(aStep->GetTrack()->GetCurrentStepNumber());
    aScintHit->SetScintVertexPosition(aStep->GetTrack()->GetVertexPosition());
    aScintHit->SetScintVelocity(aStep->GetTrack()->GetVelocity());
    aScintHit->SetScintVertexMomentumDirection(aStep->GetTrack()->GetVertexMomentumDirection());
    aScintHit->SetScintVertexKineticEnergy(aStep->GetTrack()->GetVertexKineticEnergy());
    aScintHit->SetScintMass(aStep->GetTrack()->GetDynamicParticle()->GetMass());
    
    if (aStep->GetTrack()->GetCreatorProcess()) {
        aScintHit->SetScintNuclearProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
    } else {
        aScintHit->NoProcess();
    }

    // hit is nuclear recoil
    if ( aStep->GetTrack()->GetCreatorProcess() ) {
        G4String part = aStep->GetTrack()->GetDefinition()->GetParticleType();
        G4String proc = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
        
        if ( proc ) {
            if ( part == "nucleus" ){
                if ( (proc == "pi+Elastic") || (proc == "pi-Elastic") || (proc == "k+Elastic") || (proc == "k-Elastic") || (proc == "k0LElastic") ||(proc == "k0SElastic") || (proc == "protonElastic") || (proc == "neutronElastic") ) {
                    aScintHit->SetRecoil();
                }
            }
        }        
    }

    // hit is inelastic nuclear recoil
    if( aStep->GetTrack()->GetCreatorProcess() ) {
        G4String part = aStep->GetTrack()->GetDefinition()->GetParticleType();
        G4String proc = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
        
        if( part != "gamma" ) {
            if( (proc == "pi+Inelastic") || (proc == "pi-Inelastic") || (proc == "k+Inelastic") || (proc == "k-Inelastic") || (proc == "k0LInelastic") ||(proc == "k0SInelastic") || (proc == "protonInelastic") || (proc == "neutronInelastic")){
                aScintHit->SetInelasticRecoil();
            }
        }
    }

 G4cout << "TESTING AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << G4endl; 
HitID = GdScintCollection->insert(aScintHit);
 
   return true;

//aScintHit->AddEdep(edep);
G4cout<<"edep: "<<edep<<G4endl;

}

void GdScintSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	
//    if (verboseLevel>=1)
//        GdScintCollection->PrintScintHits();
// PrintAllHits() instead?
  



	    static G4int HCID = -1;
    if(HCID<0){ HCID = G4SDManager::GetSDMpointer()->GetCollectionID(GdScintCollection); }
//if(HCID<0) HCID = GetCollectionID(0);
//GdScintCollection = new GdScintHitsCollection("myGdScintSD", collectionName[0]);
    HCE->AddHitsCollection(HCID,GdScintCollection);
    G4cout<<"HCID HCID HCID HCID: "<<HCID<<G4endl;
    G4int nHits = GdScintCollection->entries();
    if (verboseLevel>=1) {
        G4cout << "     Number of Scint hits: " << nHits << G4endl;
    }

}

void GdScintSD::clear()
{;}

void GdScintSD::DrawScint()
{;}

void GdScintSD::PrintScint()
{;} 
