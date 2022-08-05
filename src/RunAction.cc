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
// $Id: RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4Electron.hh"

#include "myHistogram.hh"
#include "RunActionMessenger.hh"

#include <fstream>

RunAction::RunAction()
: G4UserRunAction(),
  fRunActionMessenger(),
  fHistogramFileName()
{

  fWarningEnergy   =   0.0 * keV;  // Arbitrary 
  fImportantEnergy =   0.0 * keV;  // Arbitrary 
  fNumberOfTrials  =   300;  // Arbitrary

  fRunActionMessenger     = new RunActionMessenger(this); 

  fEnergyHist_1               = new myHistogram(); // 1000 km in 1 km bins
  fEnergyHist2D_1             = new myHistogram(std::log10(0.250), std::log10(10000), 101); 
  							// [250 eV, 10 MeV] in 100 bins
  fEnergyHist_2               = new myHistogram(); // 1000 km in 1 km bins
  fEnergyHist2D_2             = new myHistogram(std::log10(0.250), std::log10(10000), 101); 
  							// [250 eV, 10 MeV] in 100 bins

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fEnergyHist_1;
  delete fEnergyHist2D_1;
  delete fEnergyHist_2;
  delete fEnergyHist2D_2;
  delete fRunActionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  ChangeLooperParameters( G4Electron::Definition() );
}

void RunAction::ChangeLooperParameters(const G4ParticleDefinition* particleDef )
{
  if( particleDef == nullptr )
     particleDef = G4Electron::Definition();
  auto transportPair= findTransportation( particleDef );
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if( transport != nullptr )
  { 
     // Change the values of the looping particle parameters of Transportation 
     if( fWarningEnergy   >= 0.0 )
        transport->SetThresholdWarningEnergy(  fWarningEnergy ); 
     if( fImportantEnergy >= 0.0 )
        transport->SetThresholdImportantEnergy(  fImportantEnergy ); 
     if( fNumberOfTrials  > 0 )
        transport->SetThresholdTrials( fNumberOfTrials );
  }
  else if( coupledTransport != nullptr )
  { 
     // Change the values for Coupled Transport
     if( fWarningEnergy   >= 0.0 )
        coupledTransport->SetThresholdWarningEnergy(  fWarningEnergy ); 
     if( fImportantEnergy >= 0.0 )
        coupledTransport->SetThresholdImportantEnergy(  fImportantEnergy ); 
     if( fNumberOfTrials  > 0 )
        coupledTransport->SetThresholdTrials( fNumberOfTrials );
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run*)
{

  G4cout << "Writing results to histogram...";
  fEnergyHist_1->WriteHistogramToFile("electron_dep_" + fHistogramFileName);
  fEnergyHist2D_1->Write2DHistogram("electron_ene_"   + fHistogramFileName);
  fEnergyHist_2->WriteHistogramToFile("photon_dep_"   + fHistogramFileName);   
  fEnergyHist2D_2->Write2DHistogram("photon_ene_"     + fHistogramFileName);   
  G4cout << "complete!" << G4endl;

}


std::pair<G4Transportation*, G4CoupledTransportation*> RunAction::findTransportation( 
				  const G4ParticleDefinition* particleDef,
                                  bool reportError )
{
  const auto *partPM=  particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport=
     dynamic_cast<G4CoupledTransportation*>(partTransport);

  if( reportError && !transport && !coupledTransport )
  {
     G4cerr << "Unable to find Transportation process for particle type "
            << particleDef->GetParticleName()
            << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
            << G4endl;
  }
  
  return
     std::make_pair( transport, coupledTransport );
         // <G4Transportation*, G4CoupledTransportation*>  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
