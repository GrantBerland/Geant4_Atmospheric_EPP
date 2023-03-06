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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
// #include "DetectorAnalysis.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"

// Initialize autolock for multiple threads writing into a single file
namespace{ G4Mutex aMutex=G4MUTEX_INITIALIZER; } 

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fEnergyThreshold_keV(0.),
  fWindowAlt(500.),
  fDataCollectionType(0),
  fPhotonFilename(),
  fSteppingMessenger()
{

  fSteppingMessenger = new SteppingActionMessenger(this);

  switch(fDataCollectionType)
  {
    
    case(0):
      G4cout << "Energy deposition being recorded...";
      fRunAction->fEnergyHist_1->InitializeHistogram();
      G4cout << "Histogram initialized!" << G4endl;
      break;
    
    case(1):
      G4cout << "Particle trajectory being recorded..." << G4endl;
      break;
    
    case(2):
      G4cout << "Particle backscatter flux being recorded..." << G4endl;
      break;

    case(3):
      G4cout << "Photon data being recorded..." << G4endl;
      break;
   
    case(4):
      G4cout << "Radiation and Ionization being recorded..." << G4endl;
      fRunAction->fEnergyHist_1->InitializeHistogram();
      fRunAction->fEnergyHist2D_1->InitializeHistogram();
      fRunAction->fEnergyHist_2->InitializeHistogram();
      fRunAction->fEnergyHist2D_2->InitializeHistogram();
      break;

    default:
      throw std::invalid_argument("No data being recorded, exiting...");
      break;
  }

}

SteppingAction::~SteppingAction()
{
  delete fSteppingMessenger;
}


void SteppingAction::UserSteppingAction(const G4Step* step)
{

  G4Track* track = step->GetTrack();
 
  switch(fDataCollectionType)
  {
    
    case(0):  // Collects energy deposition per altitude
    { 
    	// Gets energy delta of particle over step length
    	const G4double energyDep = step->GetPreStepPoint()->GetKineticEnergy() - 
		step->GetPostStepPoint()->GetKineticEnergy();
    	
	if(energyDep > fEnergyThreshold_keV*keV)
    	{
	  // Gets altitude of particle
      	  G4ThreeVector position = track->GetPosition();
      	  G4double      zPos     = position.z();
      
          // Adds energy deposition to vector owned by RunAction, which is
          // written to a results file per simulation run
      	  G4int altitudeAddress = std::floor(500. + zPos/km);
      
	  // Thread lock this so only one thread can deposit energy into
	  // the histogram at a time. Unlocks when l goes out of scope.
	  if(altitudeAddress > 0 && altitudeAddress < 1000) 
	  {
	    LogEnergy(altitudeAddress, energyDep/keV);
	  }
        }
      break;
    }
    case(1):  // Collects particle trajectory (warning, lots of data!)
    {
      G4String particleName = 
	  track->GetDynamicParticle()->GetDefinition()->GetParticleName();
      
      if(particleName == "e-")
      {

        const G4double partEnergy = 
		step->GetPreStepPoint()->GetKineticEnergy();	
        G4ThreeVector position = track->GetPosition();
        G4double pos_array[4] = { position.x()/m, 
	      		          position.y()/m, 
			          position.z()/m,
				  partEnergy/keV};
        // Writes 3D position vector to results file
	// owned by RunAction
        fRunAction->fEnergyHist_1->WriteDirectlyToFile("part_traj.txt", 
			                             pos_array,
				sizeof(pos_array)/sizeof(*pos_array));
      
      }
      if(particleName == "gamma")
      {

        const G4double partEnergy = 
		step->GetPreStepPoint()->GetKineticEnergy();	
        G4ThreeVector position = track->GetPosition();
        G4double pos_array[4] = { position.x()/m, 
	      		          position.y()/m, 
			          position.z()/m,
				  partEnergy/keV};
        // Writes 3D position vector to results file
	// owned by RunAction
        fRunAction->fEnergyHist_1->WriteDirectlyToFile(fPhotonFilename, 
			                             pos_array,
				sizeof(pos_array)/sizeof(*pos_array));
      
      }
      break;
    }
    case(2): // Electron and photon backscatter tracking WIP
    {  
      G4String particleName = 
	  track->GetDynamicParticle()->GetDefinition()->GetParticleName();
      
      if(particleName == "e-" || particleName == "gamma")
      {
        
	G4ThreeVector momentumDirection = track->GetMomentumDirection();

	if(momentumDirection.z() > 0)
	{

	}

      }
      
      break;
    }

    case(3): // Photon tracking at 500 km altitude for precipitation inversion method
    {

      // Check if particle is a photon (most restrictive check)
      if(track->GetDynamicParticle()->GetDefinition()->GetParticleName() == "gamma")
      {
	
	// If particle is a photon, check that its altitude is at or greater than AEPEX altitude      
	// (less restrictive check)
	const G4ThreeVector position = track->GetPosition();
        
	if(position.z()/km > 0.) // AEPEX altitude of +500 km ASL
	{

	   // Lock scope to stop threads from overwriting data in same file
	   G4AutoLock lock(&aMutex); 

	   const G4ThreeVector momentumDirection = track->GetMomentumDirection();

	   const G4double partEnergy =  step->GetPreStepPoint()->GetKineticEnergy();

	   // Write data to file
	   std::ofstream dataFile;
	   dataFile.open(fPhotonFilename, std::ios_base::app);

	   dataFile << position.x()/m << ',' 
		    << position.y()/m << ','
		    << position.z()/m << ','
		    << momentumDirection.x() * partEnergy/keV << ','
		    << momentumDirection.y() * partEnergy/keV << ','
		    << momentumDirection.z() * partEnergy/keV << '\n'; 

	   dataFile.close();

	   // Kill photon after data collection to save processing time
	   track->SetTrackStatus(fStopAndKill);
	}
      }

      break;
    }


    case(4): // Radiation and ionization histograms 
    {

      G4int flag = 0;
      G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
      if( std::strcmp(particleName, "e-") == 0 ) // electron case
      {      
	flag = 1;
      }
      else if( std::strcmp(particleName, "gamma") == 0 ) // photon case
      {
	flag = 2;
      }
      else if( std::strcmp(particleName, "e+") == 0 ) // positron case
      {      
	flag = 1;
      }
      else
      {
	// If simulating high energy physics (~ GeV electrons), nuclear interactions are likely
	// and this line will be hit. Can change to warning or logging behavior in that case.
	throw std::runtime_error("New particle alert!! Particle: " + particleName);
      }




      if(flag > 0)
      {

        // Adds energy deposition to vector owned by RunAction, which is
        // written to a results file per simulation run
    	
	// Gets energy delta of particle over step length
    	G4double energyBefore = track->GetWeight() * step->GetPreStepPoint()->GetKineticEnergy(); 
    	G4double energyAfter  = track->GetWeight() * step->GetPostStepPoint()->GetKineticEnergy();
	
	//G4double energyDep  = energyBefore - energyAfter;

	// This line shouldn't be hit but it's a good check either way 
	if( energyBefore < energyAfter) throw std::runtime_error("Particle gained energy!");

	G4double energyDep = track->GetWeight() * step->GetTotalEnergyDeposit();
	
	// Gets altitude of particle
      	G4double zPos = track->GetPosition().z();
     
	// Rounds altitude to nearest kilometer
 	G4int altitudeAddress = std::round(500. + zPos/km);
	  
	
	// Check for valid altitude address
	if(altitudeAddress >= 0 && altitudeAddress < 500) 
	{
	  // Arguments:
	  // altitude, energy lost, current energy before loss calculation, particle type
	  LogEnergyToSpecificHistogram(altitudeAddress, energyDep, energyBefore, flag);
	}
        else if(altitudeAddress >= 500) 
            // escape outside simulation, top boundary 
	{
	  // Arguments:
	  // altitude, current energy, current energy, particle type
	  LogEnergyToSpecificHistogram(500, energyAfter, energyAfter, flag);

          // Log pitch angle of backscattered electrons with energy
          if(flag == 1) 
          {
	      const G4ThreeVector mom = track->GetMomentumDirection();
              
              LogPitchAngle(energyAfter, mom.x(), mom.y(), mom.z());

          }
	  track->SetTrackStatus(fStopAndKill);
	}
        else if(altitudeAddress <= 0)
            // escape outside simulation, bottom boundary 
        {  
          // Arguments:
	  // altitude, current energy, current energy, particle type
	  LogEnergyToSpecificHistogram(0, energyAfter, energyAfter, flag);
	  track->SetTrackStatus(fStopAndKill);
	}
        //track->GetNextVolume() == nullptr
        else
            // escape outside simulation, side boundary 
        {
          // Arguments:
	  // altitude, current energy, current energy, particle type
	  LogEnergyToSpecificHistogram(altitudeAddress, energyAfter, energyAfter, flag);
	  track->SetTrackStatus(fStopAndKill);
	}
      }

      break;
    }

    default: 
      throw std::runtime_error("Enter a valid data collection type!");
      break;
  
  }

}

void SteppingAction::LogEnergy(G4int histogramAddress, G4double energy)
{

  G4AutoLock lock(&aMutex);

  fRunAction->fEnergyHist_1->AddCountToBin(histogramAddress, energy/keV);

}

void SteppingAction::LogEnergyToSpecificHistogram(G4int histogramAddress, G4double entry1, G4double entry2, G4int whichHistogram)
{

  G4AutoLock lock(&aMutex);

  switch(whichHistogram)
  {
    case(1):
      
      // Electron and positron case
      fRunAction->fEnergyHist_1->AddCountToBin(histogramAddress, entry1/keV);
      fRunAction->fEnergyHist2D_1->AddCountTo2DHistogram(histogramAddress, entry2/keV);
      break;
    
    case(2):

      // Photon case
      fRunAction->fEnergyHist_2->AddCountToBin(histogramAddress, entry1/keV);
      fRunAction->fEnergyHist2D_2->AddCountTo2DHistogram(histogramAddress, entry2/keV);
      break;
    
    default:
      throw std::runtime_error("Enter a valid histogram selection!");
      break;
  }

}
void SteppingAction::LogPitchAngle(G4double energyAfter, G4double px, G4double py, G4double pz)
{

  G4AutoLock lock(&aMutex);

  std::ofstream backscatter_file;

  backscatter_file.open("bs_" + fRunAction->GetScoringFilename(), std::ios_base::app);

  // Encode energy into momentum direction to save some write time
  backscatter_file << px * energyAfter/keV << ','    
                   << py * energyAfter/keV << ',' 
                   << pz * energyAfter/keV << '\n';

  backscatter_file.close();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
