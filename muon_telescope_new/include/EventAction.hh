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

#ifndef EventAction_h
#define EventAction_h 1
#include "G4Types.hh"
#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <map>

class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:  
    EventAction(DetectorConstruction*, PrimaryGeneratorAction*, G4String);
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SumDeStep(G4int, G4int, G4int, G4int, G4double);
	
    void WriteFibers(const G4Event*);
    
    void ExtraInformation (G4ThreeVector, G4ThreeVector, G4double, G4int);
			         	    

  private:  
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    G4String identity;
	
	G4int nbOfModules, nbOfLayers, kLayerMax, nbDetectors;     
    std::vector<G4double>   EtotLayer;
    std::vector<G4double>   EvisLayer;
	
	G4double EtotCalor;
	G4double EvisCalor;
	
	std::map<G4int, G4double> EvisFiber;
	//G4String* secondary_type;
	//G4double* secondary_kEnergy;
	//size_t n_secondary;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    