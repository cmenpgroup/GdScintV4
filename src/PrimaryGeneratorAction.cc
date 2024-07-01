//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file GdScint/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the GdScint::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
//#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4double inch = 0.0254*m;
G4double pi = 4*std::atan(1.);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 2;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
 /* G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
  fParticleGun->SetParticleEnergy(6.*MeV);
*/
G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="neutron");G4double inch = 0.0254*m;
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(10.*MeV);
  G4double xm = 2*G4UniformRand() - 1;
 G4double ym = 2*G4UniformRand() - 1;
 G4double zm = 2*G4UniformRand()- 1;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xm,ym,zm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{


  G4double Cfz0 = 0;
  G4double CfSizeRMax= 0;
  G4double CfSizeDPhi = 0;
 if (!fCalifornium_Rod)
  {
    G4LogicalVolume* CfLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Californium");
    if ( CfLV ) fCalifornium_Rod = dynamic_cast<G4Tubs*>(CfLV->GetSolid());
  }

  if ( fCalifornium_Rod ) {
    Cfz0 = fCalifornium_Rod->GetZHalfLength();
    CfSizeRMax = fCalifornium_Rod->GetOuterRadius(); 
   // CfSPhi = fCalifornium_Rod->GetSPhi()*2.;
    CfSizeDPhi = fCalifornium_Rod-> GetDeltaPhiAngle();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Cf detector volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("GdScintPrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }


  
  G4double SizeRMax = CfSizeRMax * (G4UniformRand()-0.5) -7*inch;

  //G4double SizeDPhi = CfSizeDPhi * (G4UniformRand()-.5) +3*inch;

  G4double x0 = CfSizeRMax * cos(G4UniformRand() * 2*pi) +3*inch;
  G4double y0 = CfSizeRMax *  sin(G4UniformRand() *2*pi) -7*inch;
  G4double z0 =  Cfz0 * (G4UniformRand()-.5) -3.5*inch;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


