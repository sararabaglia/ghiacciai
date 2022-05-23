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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int t = 0;

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:G4UserSteppingAction(),detector(det),eventAct(evt)
{
  first = true;
  lvol_world = lvol_calorimeter = lvol_module = lvol_layer = lvol_fiber = lvol_mountain = 0;
  
  
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void SteppingAction::UserSteppingAction(const G4Step* step )
{ 
 //some initialisation
 // 
 if (first) {
   lvol_world  = detector->GetLvolWorld();
   lvol_calorimeter = detector->GetLvolCalorimeter();
   lvol_module = detector->GetLvolModule();   
   lvol_layer  = detector->GetLvolLayer();
   lvol_fiber  = detector->GetLvolFiber();
   lvol_mountain = detector->GetLvolMountain();
   first = false;   
 }
 
 //if no edep, return
 //
 G4double edep = step->GetTotalEnergyDeposit();
 ///if (edep == 0.) return;
 
 //G4ThreeVector pre_position = step -> GetPreStepPoint()->GetPosition();
 //G4ThreeVector post_position = step -> GetPostStepPoint()->GetPosition();
 G4Track* track = step->GetTrack();
 //G4ThreeVector positions = track->GetPosition();
 G4String particle_name = track-> GetParticleDefinition()->GetParticleName();
 //G4double kenergy = track -> GetKineticEnergy();
 
 //t++;
 //G4cout << positions << " " << t << G4endl;
 
 //G4double pre_x = positions.x();
 //G4double pre_y = positions.y();
 //G4double pre_z = positions.z();
 //locate point in geometry
 //
 G4int iModule = 0;
 G4int iLayer  = 0;
 G4int iFiber  = 0;
 G4int iDetector  = 0;
 //G4String * particle_name;
 //G4double * kEnergy;
 
 G4ThreeVector MomentumDirection(0.,0.,0.);
 G4ThreeVector HitPoint(0.,0.,0.);
 G4double kEnergy = 0;
   
 G4TouchableHandle touch1 = step->GetPreStepPoint()->GetTouchableHandle(); 
 G4LogicalVolume* lvol = touch1->GetVolume()->GetLogicalVolume();
 
 
 
 //G4String * particle_name;
 //G4double * kEnergy;
 
 /*if (nbtrk) {
        for (size_t lp=0; lp<(*secondary).size(); lp++) {
        	
        	particle_name[lp] = (*secondary)[lp]->GetDefinition()->GetParticleName();
        	kEnergy[lp] = (*secondary)[lp]->GetKineticEnergy();
        }
 
 	//eventAct -> StepInformation (particle_name, kEnergy, nbtrk);
 }
 */

      //G4cout << "primo punto " << t << G4endl;
      
      if (lvol == lvol_mountain) {
      if (particle_name == "mu-"){
      	if (t == 0){
      	MomentumDirection = step->GetPreStepPoint()->GetMomentumDirection();
      	HitPoint = step->GetPreStepPoint()->GetPosition();
      	kEnergy = step->GetPreStepPoint()->GetKineticEnergy();
      	t++;
      	//G4cout << "eccomi qua mont" << G4endl;
      	//G4cout << "mont punto " << t << G4endl;
      	}
      	}
      	}
      else {
      if (particle_name == "mu-"){
      	if (t == 1)
      	{
      		MomentumDirection = step->GetPreStepPoint()->GetMomentumDirection();
      		HitPoint = step->GetPreStepPoint()->GetPosition();
      		kEnergy = step->GetPreStepPoint()->GetKineticEnergy();
      		t=0;
      		//G4cout << "no mont punto " << t << G4endl;
      	}
      }
      }
      //G4cout << "quarto punto " << t << G4endl;
      //kEnergy = step->GetPreStepPoint()->GetKineticEnergy();
      //if (kEnergy == 0) t=0;
      
      if (lvol == lvol_world) { if (particle_name == "mu-"){eventAct->ExtraInformation	(MomentumDirection, HitPoint, kEnergy, t);}
      			       return;
      			       }
 else if (lvol == lvol_calorimeter) {iDetector = touch1->GetCopyNumber(0);}
 else if (lvol == lvol_module) { iModule = touch1->GetCopyNumber(0);
 				 iDetector = touch1->GetCopyNumber(1);}
 else if (lvol == lvol_layer)  { iLayer  = touch1->GetCopyNumber(0);
                                 iModule = touch1->GetCopyNumber(1);
                                 iDetector = touch1->GetCopyNumber(2);}
 else if (lvol == lvol_fiber)  { iFiber  = touch1->GetCopyNumber(0);
                                 iLayer  = touch1->GetCopyNumber(1);
                                 iModule = touch1->GetCopyNumber(2);
                                 iDetector = touch1->GetCopyNumber(3);}

 /*t++;
 const std::vector<const G4Track*>* secondary 
                                    = step->GetSecondaryInCurrentStep();
      size_t nbtrk = (*secondary).size();
      if (nbtrk) {
        G4cout << "\n    :----- List of secondaries "<< t << "----------------" << G4endl;
        G4cout.precision(4);
        for (size_t lp=0; lp<(*secondary).size(); lp++) {
          G4cout << "   "
                 << std::setw(13)                 
                 << (*secondary)[lp]->GetDefinition()->GetParticleName()
                 << ":  energy ="
                 << std::setw(6)
                 << G4BestUnit((*secondary)[lp]->GetKineticEnergy(),"Energy")
                 << "  time ="
                 << std::setw(6)
                 << G4BestUnit((*secondary)[lp]->GetGlobalTime(),"Time");
          G4cout << G4endl;
        }
              
        G4cout << "    :------------------------------------------\n" << G4endl;
      }*/
 // sum edep
 //
 //G4cout << "sesto punto " << t << G4endl;
 if (particle_name == "mu-") eventAct->ExtraInformation(MomentumDirection, HitPoint, kEnergy, t);
 eventAct->SumDeStep(iDetector, iModule, iLayer, iFiber, edep);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 G4Material* material = aStep->GetTrack()->GetMaterial();
 G4double birk1       = material->GetIonisation()->GetBirksConstant();
 G4double destep      = aStep->GetTotalEnergyDeposit();
 G4double stepl       = aStep->GetStepLength();  
 G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 //
 G4double response = destep;
 if (birk1*destep*stepl*charge != 0.)
   {
     response = destep/(1. + birk1*destep/stepl);
   }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

