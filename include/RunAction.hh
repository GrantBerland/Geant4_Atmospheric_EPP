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
// $Id: RunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <vector>
#include "myHistogram.hh"

// Choose your fighter:
// #include "g4root.hh"
#include "g4xml.hh"
// #include "g4csv.hh"

class G4Run;
class SteppingAction;
class myHistogram;

class RunActionMessenger;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    inline void AddPhotonEnergy(G4double ene);
    inline void AddPhotonAltitude(G4double alt);
    inline void AddElectronEnergy(G4double ene);
    inline void AddElectronAltitude(G4double alt);

    void SetHistFileName(G4String name){fHistogramFileName=name;};

  public:
    myHistogram           *fEnergyHist;
  private:
    RunActionMessenger*    fRunActionMessenger;
    G4String 		   fHistogramFileName;
    std::vector<G4double> *fPhotonEnergyVector;
    std::vector<G4double> *fPhotonAltitudeVector;
    std::vector<G4double> *fElectronEnergyVector;
    std::vector<G4double> *fElectronAltitudeVector;
    
};

// TODO: make this 1 or 2 methods with particle type selection, etc.
inline void RunAction::AddPhotonEnergy(G4double ene)
{
  fPhotonEnergyVector->push_back(ene);
}

inline void RunAction::AddPhotonAltitude(G4double alt)
{
  fPhotonAltitudeVector->push_back(alt);
}

inline void RunAction::AddElectronAltitude(G4double alt)
{
  fElectronAltitudeVector->push_back(alt);
}

inline void RunAction::AddElectronEnergy(G4double ene)
{
  fElectronEnergyVector->push_back(ene);
}

#endif
