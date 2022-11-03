

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim)
  :G4UImessenger(),
   fPrimaryGenerator(prim) 
{
  fPrimDir = new G4UIdirectory("/energy/");
  fPrimDir->SetGuidance("Select particle distributions.");

  fcmd = new G4UIcmdWithAnInteger("/energy/setEnergyDistributionType",this);
  fcmd->SetParameterName("0 - Exponential, 1 - Monoenergetic",true);
  fcmd->SetDefaultValue(0);
  fcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fcmd2 = new G4UIcmdWithAnInteger("/energy/setPitchAngleDistributionType",this);
  fcmd2->SetParameterName("TODO: write options in",true);
  fcmd2->SetDefaultValue(0);
  fcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);

  fcmd3 = new G4UIcmdWithAnInteger("/energy/setSourceType",this);
  fcmd3->SetParameterName("0: electrons, 1: solar X-ray spectra, 2: CXB",true);
  fcmd3->SetDefaultValue(0);
  fcmd3->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fDcmd = new G4UIcmdWithADouble("/energy/setFoldingEnergy",this);
  fDcmd->SetParameterName("Folding or Mono Energy [keV]",true);
  fDcmd->SetDefaultValue(0);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fDcmd2 = new G4UIcmdWithADouble("/energy/setMaximumPitchAngle",this);
  fDcmd2->SetParameterName("Maximum pitch angle used in distributions [deg]",true);
  fDcmd2->SetDefaultValue(50.);
  fDcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDcmd3 = new G4UIcmdWithADouble("/energy/particleStartingAltitude",this);
  fDcmd3->SetParameterName("Particle initial altitude [km]",true);
  fDcmd3->SetDefaultValue(1000.);
  fDcmd3->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fcmd2;
  delete fcmd3;
  delete fDcmd;
  delete fDcmd2;
  delete fDcmd3;
}
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{

  if(command == fcmd){
    fPrimaryGenerator->SetEnergyDistribution(std::stoi(newValue));
  }    	  
  
  if(command == fcmd2){
    fPrimaryGenerator->SetPitchAngleDistribution(std::stoi(newValue));
  }    	  

  if(command == fDcmd){
    fPrimaryGenerator->SetEnergy(std::stod(newValue));
  }
  
  if(command == fDcmd2){
    fPrimaryGenerator->SetMaxPitchAngle(std::stod(newValue));
  }
  
  if(command == fDcmd3){
    fPrimaryGenerator->SetPartInitialAlt(std::stod(newValue));
  }

  if(command == fcmd3){
    fPrimaryGenerator->SetSourceInputType(std::stoi(newValue));
  }
}
