//
/// \file GdScint/src/DetectorConstruction.cc
/// \brief Implementation of the GdScint::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4SDManager.hh"
#include "GdScintSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{
  SDman = G4SDManager::GetSDMpointer();
  myGdScintSD = NULL; 
}

DetectorConstruction::~DetectorConstruction()
{ ; }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  //////////////////////////////////////////////////////////////
  //UNIT                                                     //
  /////////////////////////////////////////////////////////////

  G4double inch = 0.0254*m;
  G4double foot = 0.3048*m;
  G4double pi = 4*std::atan(1.);

  ////************************************************************************
  //// MATERIALS: Data taken from http://www.webelements.com/

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, fractionmass;
  G4double pressure, temperature;
  G4String name, symbol;
  G4int ncomponents, natoms;

  //// Define all elements used

  a = 1.00794*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1, a);

  a = 12.0107*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6, a);

//  a = 14.0*g/mole;
//  G4Element* elC14 = new G4Element(name="Carbon14", symbol="C14", z=6, a);

//  a = 14.0067*g/mole;
//  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7, a);

  a = 15.9994*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8, a);

//  a=18.998*g/mole;
//  G4Element* elF = new G4Element(name="Flourine"  ,symbol="F" , z= 9., a);

//  a = 22.9897*g/mole;
//  G4Element* elNa = new G4Element(name="Sodium", symbol="Na", z=11, a);

//  a = 24.305*g/mole;
//  G4Element* elMg = new G4Element(name="Magnesium", symbol="Mg", z=12, a);

  a=26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium"  ,symbol="Al" , z= 13., a);

  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14, a);

  a = 30.9738*g/mole;
  G4Element* elP = new G4Element(name="Phosphorus", symbol="P", z=15, a);

  a = 32.066*g/mole;
  G4Element* elS = new G4Element(name="Sulphur", symbol="S", z=16, a);

//  a = 35.453*g/mole;
//  G4Element* elCl = new G4Element(name="Chlorine", symbol="Cl", z=17, a);

  a = 39.95*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18, a);

//  a=39.964*g/mole;
//  G4Element* elK = new G4Element(name="Potassium"  ,symbol="K" , z= 19, a);

//  a=40.078*g/mole;
//  G4Element* elCa = new G4Element(name="Calcium",symbol="Ca" , z= 20., a);

  a = 51.9961*g/mole;
  G4Element* elCr = new G4Element(name="Chromium", symbol="Cr", z=24, a);

  a = 54.938*g/mole;
  G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z=25, a);

  a = 55.845*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26, a);

  a = 58.6934*g/mole;
  G4Element* elNi = new G4Element(name="Nickel", symbol="Ni", z=28, a);

//  a = 63.546*g/mole;
//  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z=29, a);

  a = 92.9064*g/mole;
  G4Element* elNb = new G4Element(name="Niobium", symbol="Nb", z=41, a);

//  a = 126.9*g/mole;
//  G4Element* elI = new G4Element(name="Iodine", symbol="I", z=53, a);

//  a=157.25*g/mole;
//  G4Element* elGd = new G4Element(name="Gadolinium"  ,symbol="Gd" , z= 64, a);

    a = 207.2*g/mole;
    G4Element* elPb = new G4Element(name="Lead", symbol="Pb", z=82, a);
    
    a = 251*g/mole;
    G4Element* elCf = new G4Element(name="Californium" ,symbol="Cf" , z = 98, a);


  // Temperature of DRIFT lab
  temperature = 293.15*kelvin;

  // Create Argon gas for the Oxy GEM detector
  density =  0.001784*g/cm3; // at STP
  pressure = 1.0*atmosphere;
  G4Material* Ar_gas = new G4Material(name="Argon_gas", density, ncomponents=1,
                          kStateGas, temperature, pressure);
  Ar_gas->AddElement(elAr, natoms=1);

  // Create CO2 gas at STP
  density = 0.001964*g/cm3; // at STP
  pressure = 1.0*atmosphere;
  G4Material* CO2 = new G4Material(name="CO2_gas", density, ncomponents=2,
                                             kStateGas, temperature);
  CO2->AddElement(elC, natoms=1);
  CO2->AddElement(elO, natoms=2);

  // 1 atm Ar-CO2 mixture
  density = 0.5*(0.001784 + 0.001964)*g/cm3; //
  pressure= 1.0*atmosphere; //
  G4Material*gem_gas = new G4Material(name="gas_GEM", density, ncomponents=2, kStateGas, temperature);
  gem_gas->AddMaterial(Ar_gas, fractionmass=50.0*perCent);
  gem_gas->AddMaterial(CO2, fractionmass=50.0*perCent);

  // Create stainless steel
  density = 8.027*g/cm3;
  G4Material* SS = new G4Material(name="Stainless_steel", density, ncomponents=9);
  SS->AddElement(elFe, fractionmass=67.445*perCent);
  SS->AddElement(elCr, fractionmass=18*perCent);
  SS->AddElement(elNi, fractionmass=11*perCent);
  SS->AddElement(elMn, fractionmass=2*perCent);
  SS->AddElement(elSi, fractionmass=1*perCent);
  SS->AddElement(elNb, fractionmass=0.4*perCent);
  SS->AddElement(elS, fractionmass=0.08*perCent);
  SS->AddElement(elP, fractionmass=0.045*perCent);
  SS->AddElement(elS, fractionmass= 0.03*perCent);

  //  Mylar
  G4Material* Mylar = new G4Material("Mylar", density = 1.39*g/cm3, ncomponents = 3);
  Mylar->AddElement(elC, natoms = 10);
  Mylar->AddElement(elH, natoms = 8);
  Mylar->AddElement(elO, natoms =4);

  //Create lead -Aaron Edit
  G4Material* Lead = new G4Material(name="Lead", density = 11.4*g/cm3, ncomponents=1);
  Lead->AddElement(elPb, natoms=1);

// create acrylic
  G4Material* Acrylic = new G4Material("Acrylic", density = .94*g/cm3, ncomponents =3);
  Acrylic->AddElement(elC, natoms = 5);
  Acrylic->AddElement(elH, natoms = 8);
  Acrylic->AddElement(elO, natoms = 2);

  G4Material* Aluminium = new G4Material("Aluminium", density =2.70*g/cm3 , ncomponents=1);
  Aluminium->AddElement(elAl, natoms=1);
//Create Californium
 
  G4Material* Californium = new G4Material("Californium", density =15.1*g/cm3 , ncomponents=1);
  Californium->AddElement(elCf, natoms=1);
  // Envelope parameters
  //
  G4double env_sizeXY = 100*cm, env_sizeZ = 100*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  // World volume
    G4double world_sizeZ = 1.2*m;
    G4double world_sizeXY = 1.2*m;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

    //----------------------------
    // Creation of the volumes
    //----------------------------

      /////////////////////////////////////////////////////////
      //World volume
      ////////////////////////////////////////////////////////

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



/*  //


  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",                    // its name
    0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
*/

  // OXY CERN Triple GEM
  //
  G4double oxy_gem_sizeZ = 10.*cm;   // full z length
  G4double oxy_gem_sizeY = 1.5*cm;		//full y length
  G4double oxy_gem_sizeX = 10.*cm;  // full x length
  G4double oxy_gem_posX = 3.*inch;     // x position with respect to world origin
  G4double oxy_gem_posY = -3.8*inch;     // y position with respect to world origin
  G4double oxy_gem_posZ = -5 - 0.5*oxy_gem_sizeZ;   // z position with respect to world origin
  G4Material* oxy_gem_mat = gem_gas; // GEM material

  G4Box* solidGEM = new G4Box("GEM",0.5*oxy_gem_sizeX,0.5*oxy_gem_sizeY, 0.5*oxy_gem_sizeZ);
  G4LogicalVolume* logicGEM = new G4LogicalVolume(solidGEM,oxy_gem_mat,"GEM");
  new G4PVPlacement(0,
                    G4ThreeVector(oxy_gem_posX,oxy_gem_posY,oxy_gem_posZ),
                    logicGEM,
                    "GEM",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);

  fScoringVolume = logicGEM;

//
//make Cf Rod


G4double Californium_Rod_sizeZ = 1*inch;      //  1/2 z length
  G4double Californium_Rod_sizeRmin = 0*inch;   // inner radius
  G4double Californium_Rod_Radius = .2*inch;   // outer radius
  G4double Californium_Rod_sizeSPhi = 0.0;        // starting angle phi in radians
  G4double Californium_Rod_sizeDPhi = 2.0*pi;     // ending angle phi in radians


  G4double Californium_Rod_posX = 3.*inch;         // x position with respect to world origin
  G4double Californium_Rod_posY = -7*inch;         // y position with respect to world origin
  G4double Californium_Rod_posZ = -3.5*inch;         // z position with respect to world origin
  G4Material* Californium_Rod_mat = Lead;          // stainless steel material


  G4VSolid* solidCalifornium_Rod;
	solidCalifornium_Rod =
	  new G4Tubs(			"Californium_Rod",
							Californium_Rod_sizeRmin,
							Californium_Rod_Radius,
							Californium_Rod_sizeZ, // used to be Californium_RodThickness
							0.*deg,
							360.*deg);
  G4LogicalVolume* logicCalifornium_Rod =
    new G4LogicalVolume(solidCalifornium_Rod,            //its solid
                        Californium_Rod_mat,             //its material
                        "Californium");       //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Californium_Rod_posX,Californium_Rod_posY,Californium_Rod_posZ), //position in the world volume
                    logicCalifornium_Rod,             //its logical volume
                    "Californium",                //its name
                    logicWorld,              //its world  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //
  // Vacuum Vessel - cylinder
  //
  G4double VesselOut_sizeZ = 0.5*3.0*foot;      //  1/2 z length
  G4double VesselOut_sizeRmin = 0.0*foot;   // inner radius
  G4double VesselOut_sizeRmax = 9.875*inch;   // outer radius
  G4double VesselOut_sizeSPhi = 0.0;        // starting angle phi in radians
  G4double VesselOut_sizeDPhi = 2.0*pi;     // ending angle phi in radians

  G4double VesselIn_sizeZ = VesselOut_sizeZ - 0.375*inch; // outer Vessel length minus wall thickness
  G4double VesselIn_sizeRmin = VesselOut_sizeRmin;   // inner radius
  G4double VesselIn_sizeRmax = VesselOut_sizeRmax - 0.375*inch;   // outer radius minus wall thickness
  G4double VesselIn_sizeSPhi = VesselOut_sizeSPhi;     // starting angle phi in radians
  G4double VesselIn_sizeDPhi = VesselOut_sizeDPhi;     // ending angle phi in radians

  G4double Vessel_posX = 1.*cm;         // x position with respect to world origin
  G4double Vessel_posY = 0.*cm;         // y position with respect to world origin
  G4double Vessel_posZ = 0.*cm;         // z position with respect to world origin
  G4Material* Vessel_mat = SS;          // stainless steel material

  G4Tubs* VesselOut = new G4Tubs("VesselOut",
                                VesselOut_sizeRmin,
                                VesselOut_sizeRmax,
                                VesselOut_sizeZ,
                                VesselOut_sizeSPhi,
                                VesselOut_sizeDPhi);

  G4Tubs* VesselIn = new G4Tubs("VesselIn",
                                VesselIn_sizeRmin,
                                VesselIn_sizeRmax,
                                VesselIn_sizeZ,
                                VesselIn_sizeSPhi,
                                VesselIn_sizeDPhi);

  G4SubtractionSolid* solidVessel = new G4SubtractionSolid("Vessel", VesselOut, VesselIn);

  G4LogicalVolume* logicVessel =
    new G4LogicalVolume(solidVessel,            //its solid
                        Vessel_mat,             //its material
                        "Vessel");       //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Vessel_posX,Vessel_posY,Vessel_posZ), //position in the world volume
                    logicVessel,             //its logical volume
                    "Vessel",                //its name
                    logicWorld,              //its world  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
//
 //plexiglass
  G4double Plexiglass_Sheet_sizeZ = 2.75*.5*foot;   // full z length
  G4double Plexiglass_Sheet_sizeY = .25*inch;		//full y length
  G4double Plexiglass_Sheet_sizeX = 15.65*.5*inch;  // full x length
  G4double Plexiglass_Sheet_posX = 0.*inch;     // x position with respect to world origin
  G4double Plexiglass_Sheet_posY = -4.5*inch;     // y position with respect to world origin
  G4double Plexiglass_Sheet_posZ = 0.*inch ;   // z position with respect to world origin
  G4Material* Plexiglass_Sheet_mat = Acrylic; // Plexiglass material

  G4Box* solidPlexiglass = new G4Box("Plexiglass",Plexiglass_Sheet_sizeX,Plexiglass_Sheet_sizeY,  Plexiglass_Sheet_sizeZ);
  G4LogicalVolume* logicPlexiglass = new G4LogicalVolume(solidPlexiglass,Plexiglass_Sheet_mat,"Plexiglass");
  new G4PVPlacement(0,
                    G4ThreeVector(Plexiglass_Sheet_posX,Plexiglass_Sheet_posY,Plexiglass_Sheet_posZ),
                    logicPlexiglass,
                    "Plexiglass",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);
  //always return the physical World

//
//Lead Brick
//

  G4double Lead_brick_sizeZ = 7.675*inch;   // full z length
  G4double Lead_brick_sizeY = 1.8*inch;		//full y length
  G4double Lead_brick_sizeX = 3.8*inch;  // full x length
  G4double Lead_brick_posX = -1*inch;     // x position with respect to world origin
  G4double Lead_brick_posY = -8*inch;     // y position with respect to world origin
  G4double Lead_brick_posZ = -4*inch;   // z position with respect to world origin
  G4Material* Lead_brick_mat = Lead; // GEM material

  G4Box* solidLead_brick = new G4Box("Lead_brick",0.5*Lead_brick_sizeX,0.5*Lead_brick_sizeY, 0.5*Lead_brick_sizeZ);
  G4LogicalVolume* logicLead_brick = new G4LogicalVolume(solidLead_brick,Lead_brick_mat,"Lead_brick");
  new G4PVPlacement(0,
                    G4ThreeVector(Lead_brick_posX,Lead_brick_posY,Lead_brick_posZ),
                    logicLead_brick,
                    "Lead_brick",
                    logicWorld,
                    false,
                    0,
                    checkOverlaps);
// lead pig
 G4double PigOut_sizeZ = 3.25*.5*inch;      //  1/2 z length
  G4double PigOut_sizeRmin = .84375*inch;   // inner radius
  G4double PigOut_sizeRmax = 1.96875*inch;   // outer radius
  G4double PigOut_sizeSPhi = 0.0;        // starting angle phi in radians
  G4double PigOut_sizeDPhi = 2.0*pi;     // ending angle phi in radians

  G4double PigIn_sizeZ = PigOut_sizeZ - 1.*inch; // outer Pig length minus wall thickness
  G4double PigIn_sizeRmin = 0.;   // inner radius
  G4double PigIn_sizeRmax = PigOut_sizeRmin ;   // outer radius minus wall thickness
  G4double PigIn_sizeSPhi = PigOut_sizeSPhi;     // starting angle phi in radians
  G4double PigIn_sizeDPhi = PigOut_sizeDPhi;     // ending angle phi in radians

  G4double Pig_posX = 3*inch;         // x position with respect to world origin
  G4double Pig_posY = -6.9*inch;         // y position with respect to world origin
  G4double Pig_posZ = -3.5*inch;         // z position with respect to world origin
  G4Material* Pig_mat = Lead;          // Lead material



  G4Tubs* PigOut = new G4Tubs("PigOut",
                                PigOut_sizeRmin,
                                PigOut_sizeRmax,
                                PigOut_sizeZ,
                                PigOut_sizeSPhi,
                                PigOut_sizeDPhi);

  G4Tubs* PigIn = new G4Tubs("PigIn",
                                PigIn_sizeRmin,
                                PigIn_sizeRmax,
                                PigIn_sizeZ,
                                PigIn_sizeSPhi,
                                PigIn_sizeDPhi);

  G4SubtractionSolid* solidPig = new G4SubtractionSolid("Pig", PigOut , PigIn);

  G4LogicalVolume* logicPig =
    new G4LogicalVolume(solidPig,            //its solid
                        Pig_mat,             //its material
                        "Pig");       //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Pig_posX,Pig_posY,Pig_posZ), //position in the world volume
                    logicPig,             //its logical volume
                    "Pig",                //its name
                    logicWorld,              //its world  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

// Pigend

 G4double Pigend_sizeZ = 1.875*.5*inch;      //  1/2 z length
 G4double Pigend_sizeRmin = 0*inch;   // inner radius
  G4double PigendRadius = 1.96875*inch;   // outer radius
  G4double Pigend_sizeSPhi = 0.0;        // starting angle phi in radians
  G4double Pigend_sizeDPhi = 2.0*pi;     // ending angle phi in radians


  G4double Pigend_posX = 3*inch;         // x position with respect to world origin
  G4double Pigend_posY = -6.9*inch;         // y position with respect to world origin
  G4double Pigend_posZ = -0.915*inch;         // z position with respect to world origin
  G4Material* Pigend_mat = Lead;          // Lead material

  G4VSolid* solidPigend;
	solidPigend =
	  new G4Tubs(			"Pigend",
							Pigend_sizeRmin,
							PigendRadius,
							Pigend_sizeZ, // used to be PigendThickness
							0.*deg,
							360.*deg);
  G4LogicalVolume* logicPigend =
    new G4LogicalVolume(solidPigend,            //its solid
                        Pigend_mat,             //its material
                        "Pigend");       //its name
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Pigend_posX,Pigend_posY,Pigend_posZ), //position in the world volume
                    logicPigend,             //its logical volume
                    "Pigend",                //its name
                    logicWorld,              //its world  volume b
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

//
//Pig cap

 G4double Pigcap_sizeZ = .46875*inch;      //  1/2 z length
 G4double Pigcap_sizeRmin = 0*inch;   // inner radius
  G4double PigcapRadius = 1.96875*inch;   // outer radius
  G4double Pigcap_sizeSPhi = 0.0;        // starting angle phi in radians
  G4double Pigcap_sizeDPhi = 2.0*pi;     // caping angle phi in radians


  G4double Pigcap_posX = 3*inch;         // x position with respect to world origin
  G4double Pigcap_posY = -6.9*inch;         // y position with respect to world origin
  G4double Pigcap_posZ = -5.59375*inch;         // z position with respect to world origin
  G4Material* Pigcap_mat = Lead;          // Lead material

  G4VSolid* solidPigcap;
	solidPigcap =
	  new G4Tubs(			"Pigcap",
							Pigcap_sizeRmin,
							PigcapRadius,
							Pigcap_sizeZ, // used to be PigcapThickness
							0.*deg,
							360.*deg);
  G4LogicalVolume* logicPigcap =
    new G4LogicalVolume(solidPigcap,            //its solid
                        Pigcap_mat,             //its material
                        "Pigcap");       //its name
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Pigcap_posX,Pigcap_posY,Pigcap_posZ), //position in the world volume
                    logicPigcap,             //its logical volume
                    "Pigcap",                //its name
                    logicWorld,              //its world  volume b
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    /////////////////////////
    //SENSITIVE DETECTORS  //
    /////////////////////////

  // // G4SDManager* SDman = G4SDManager::GetSDMpointer();
   // G4SDManager::GetSDMpointer()->AddNewDetector("GdScintSD");    
////    if(!myGdScintSD)
   // {
    ////    myGdScintSD = new GdScintSD( "GdScintSD", this  );
        ////SDman -> AddNewDetector( myGdScintSD );
   // } 
   //// logicGEM -> SetSensitiveDetector ( myGdScintSD );

  //
  //always return the physical World
  //
  return physWorld;
}
void DetectorConstruction::ConstructSDandField()
{
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	myGdScintSD = new GdScintSD("GdScintSD", this);
	SDman ->AddNewDetector(myGdScintSD);
	logicGEM -> SetSensitiveDetector(myGdScintSD);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
