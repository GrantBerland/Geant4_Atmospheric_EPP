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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "F03FieldSetup.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"

#include "G4AutoDelete.hh"
#include "G4SDManager.hh"

#include "DetectorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fAtmosphereFilename("MSISE00_atmosphere.csv"),
  fDetectorMessenger(),
  fTableSize(0),
  fLogicWorld(0)
{
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();


  // Material: Vacuum
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  G4Material* low_density_material = new G4Material("Low_dens",
              1.0 , 1.01*g/mole, 1.0E-10*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-12*pascal );

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World
  G4double world_sizeXY = 1000.0*km;
  G4double world_sizeZ  = 1000.0*km;
  // G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Tubs* solidWorld = new G4Tubs("World",     //its name
       0.,  			// inner radius
       world_sizeXY,  	        // outer radius
       0.5*world_sizeZ,   	// z half length
       0.,			// starting phi
       360.*deg);		// segment angle

  fLogicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_material,     //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fLogicWorld,           //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // msisAtmosTable contents:
  //
  // [0] - alt [km] 
  // [1] - O   * 
  // [2] - N2  *
  // [3] - O2  *
  // [4] - Total mass den [kg/m^3] 
  // [5] - Temp [K] 
  // [6] - He  *
  // [7] - Ar  *
  // [8] - H   *
  // [9] - N   *
  //
  // * all species' mass density in [kg/m^3]

  fTableSize = GetMSIStableSize(fAtmosphereFilename);

  // Cast to const for table instantiation
  unsigned const int tableSize = fTableSize;
  G4double msisAtmosTable[tableSize][10];
 
  // Populate array with MSIS atmosphere table
  GetMSIStable(msisAtmosTable,      // array to populate 
	       fAtmosphereFilename, // filename of table
	       tableSize);          // int of table size

  // Atmospheric material definitions
  G4Element* O = new G4Element("Oxygen",
		  	       "O",
		               8.,
			       16.0*g/mole);
  
  G4Element* N = new G4Element("Nitrogen",
		  	       "N",
		               7.,
			       14.0*g/mole);

  G4Element* Ar = new G4Element("Argon",
		  	       "Ar",
		               18.,
			       39.948*g/mole);
  
  G4Element* He = new G4Element("Helium",
		  	       "He",
		               2.,
			       4.0*g/mole);
  
  G4Element* H = new G4Element("Hydrogen",
		  	       "H",
		               1.,
			       1.0078*g/mole);
  


  std::cout << "Table size: " << tableSize << std::endl;

  // Layers are the size needed to fill the 1000 km column
  G4double layerThickness = ( int(1000 / tableSize) ) * km - 0. * m;
  G4double layerLocation;
  
  G4Tubs* atmosphereLayer = new G4Tubs("AtmosphereLayer",
		       	       		0.,
		  			world_sizeXY-1.*mm,
		       	       		0.5*layerThickness,
					0.,
					360.*deg);
  G4Material*      layerMaterial;
  G4Material*      diOxygen;
  G4Material*      diNitrogen;
  G4LogicalVolume* logicLayer;
  G4double         pressure;
  G4double         R_gas_constant_air = 287.;  // J/kg-K air
  G4double         totalAtmosMass;
  // approximate minimum number density [cm^-3] that Geant will tolerate
  G4double         zeroThreshold = 1e-21; 
  G4int            nComponents;
  for(int i=0; i<fTableSize; i++)
  {
     // Ideal gas law for atmospheric pressure
     // P [Pa] = R [J/kg-K air] * rho [kg/m^3] * T [K]
     pressure = R_gas_constant_air * 
	        (msisAtmosTable[i][2]) * 
	         msisAtmosTable[i][4];
     // TODO: replace R gas constant with layer mean mass! 

     nComponents    = 0;
     totalAtmosMass = 0;
     
     // atomic oxygen
     if(msisAtmosTable[i][1] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][1];}
     
     // di-Nitrogen
     if(msisAtmosTable[i][2] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][2];
        diNitrogen = new G4Material("Dinitrogen-Layer"+std::to_string(i),
		                    msisAtmosTable[i][2]*kg/m3,
				    1,
				    kStateGas,
				    msisAtmosTable[i][5]*kelvin,
				    pressure*pascal);
             diNitrogen->AddElement(N, 2);}
     else{ diNitrogen = nullptr; }
     
     // di-Oxygen
     if(msisAtmosTable[i][3] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][3];
             diOxygen = new G4Material("Dioxgyen-Layer"+std::to_string(i),
		                    msisAtmosTable[i][3]*kg/m3,
				    1,
				    kStateGas,
				    msisAtmosTable[i][4]*kelvin,
				    pressure*pascal);
             diOxygen->AddElement(O, 2);}
     
     // Helium
     if(msisAtmosTable[i][6] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][6];}
     
     // Argon
     if(msisAtmosTable[i][7] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][7];}
     
     // Hydrogen
     if(msisAtmosTable[i][8] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][8];}
     
     // Nitrogen
     if(msisAtmosTable[i][9] > zeroThreshold){ 
	     nComponents++;
      	     totalAtmosMass += msisAtmosTable[i][9];}
     

     // test material
     //layerMaterial = low_density_material;
     
     layerMaterial = new G4Material("AirLayer"+std::to_string(i), // name
		     totalAtmosMass*kg/m3,       // density
		     nComponents,                // number of components
		     kStateGas,                  // state
		     msisAtmosTable[i][4]*kelvin,// temperature
		     pressure*pascal); 	         // pressure
    
     // alt [km], O, N2, O2, Mass den [g/cm^3], Temp [K], He, Ar, H, N
    
     // O
     if(msisAtmosTable[i][1] > zeroThreshold) 
     {
	layerMaterial->AddElement(O, 
		       msisAtmosTable[i][1]/totalAtmosMass);
     }
     // N2
     if(msisAtmosTable[i][2] > zeroThreshold) 
     {
     	layerMaterial->AddMaterial(diNitrogen, // material, mass fraction
		     msisAtmosTable[i][2]/totalAtmosMass); 
     }
     // O2
     if(msisAtmosTable[i][3] > zeroThreshold) 
     {
     	layerMaterial->AddMaterial(diOxygen, // material, mass fraction
		     msisAtmosTable[i][3]/totalAtmosMass); 
     }
     // He 
     if(msisAtmosTable[i][6] > zeroThreshold) 
     {
     	layerMaterial->AddElement(He, // material, mass fraction
		     msisAtmosTable[i][6]/totalAtmosMass); 
     }
     // Ar 
     if(msisAtmosTable[i][7] > zeroThreshold) 
     {
     	layerMaterial->AddElement(Ar, // material, mass fraction
		     msisAtmosTable[i][7]/totalAtmosMass); 
     }
     // H 
     if(msisAtmosTable[i][8] > zeroThreshold) 
     {
     	layerMaterial->AddElement(H, // material, mass fraction
		     msisAtmosTable[i][8]/totalAtmosMass); 
     }
     // N
     if(msisAtmosTable[i][9] > zeroThreshold) 
     {
     	layerMaterial->AddElement(N, // material, mass fraction
		     msisAtmosTable[i][9]/totalAtmosMass); 
     }

     logicLayer = new G4LogicalVolume(atmosphereLayer,
		                   layerMaterial,
				   "AtmosphereLayer"+std::to_string(i));

     
     layerLocation = (i-fTableSize/2)*layerThickness + layerThickness/2.;
     
     new G4PVPlacement(0,
		       G4ThreeVector(0.,0.,layerLocation),
		       logicLayer,
		       "AtmosphereLayer"+std::to_string(i),
		       fLogicWorld,
		       false,
		       i,
		       checkOverlaps);

  }

  // always return the physical World
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{

  if(!fEmFieldSetup.Get())
  {    
    F03FieldSetup* emFieldSetup = new F03FieldSetup();    
    fEmFieldSetup.Put(emFieldSetup);
      
    G4AutoDelete::Register(emFieldSetup);
  }
    fLogicWorld->SetFieldManager(fEmFieldSetup.Get()->
		        GetGlobalFieldManager(), true); 
}

void DetectorConstruction::GetMSIStable(G4double tableEntry[][10],
					G4String filename,
					unsigned int tableSize)
{
  
  G4int counter_word = 0;
  std::ifstream filePtr;
  std::string line;

  // OPEN
  filePtr.open(filename, std::ios::in);

  // Check that file is open/accesible, throw error if not
  if(!filePtr.is_open())
  { 
    G4String errorMessage = "Could not open file ";
    errorMessage += filename;
    errorMessage += " !\n";
    throw std::runtime_error(errorMessage);
  }


  // Fill data array with data
  for(unsigned int lineIndex=0; lineIndex<tableSize; lineIndex++)
  {
    // Get line
    filePtr >> line;

    // Instantiate stringstream from line
    std::stringstream s_ptr(line); 
    
    counter_word = 0;
    while(s_ptr.good())
    {
      // Parse line into words delimited by commas
      G4String word;
      getline(s_ptr, word, ',');
      
      // Convert to double and assign table entry
      tableEntry[lineIndex][counter_word] = std::stod(word);

      counter_word++;
    }
  }

  // CLOSE
  filePtr.close();

}

G4int DetectorConstruction::GetMSIStableSize(G4String filename)
{

  // OPEN
  std::ifstream filePtr(filename, std::ios::in);
  G4String line;
  G4int counter = 0;

  while(getline(filePtr, line)) {counter++;}

  // CLOSE
  filePtr.close();

  return counter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
