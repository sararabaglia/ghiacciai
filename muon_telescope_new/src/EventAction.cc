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

#include "EventAction.hh"
#include "SteppingAction.hh"

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det,PrimaryGeneratorAction* prim, G4String Id)
:detector(det), primary(prim), identity(Id)
{
  nbDetectors = detector->GetNbDetector();
  nbOfModules = detector->GetNbModules();	 	
  nbOfLayers  = detector->GetNbLayers();
  kLayerMax = nbOfModules*nbOfLayers + 1;
  
  EtotCalor = EvisCalor = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int state = 0;
G4int hitmont = 0;
G4int event_number = 0;

void EventAction::BeginOfEventAction(const G4Event*)
{
  EtotLayer.resize(kLayerMax);
  EvisLayer.resize(kLayerMax);			
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k] = EvisLayer[k] = 0.0;
  }
  EtotCalor = EvisCalor = 0.;
  EvisFiber.clear();
  event_number++;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void EventAction::SumDeStep(G4int iDetector, G4int iModule, G4int iLayer, G4int iFiber,
                            G4double deStep)
{
  if (iModule > 0) EtotCalor += deStep;
  		
  G4int kLayer = 0; G4int kFiber = 0; G4int kDetector = 0;
  if (iDetector > 0) kDetector = iDetector;
  if (iLayer > 0) {
	kLayer = (iModule-1)*nbOfLayers + iLayer;
	EtotLayer[kLayer] += deStep;
  }
  
  if (iFiber > 0) {
	EvisLayer[kLayer] += deStep;
	EvisCalor += deStep;
    kFiber = 1000*kLayer + iFiber;
	EvisFiber[kFiber] += deStep;	  	
  }
  
  /*std::ofstream File("secondary_particle.txt", std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
  
  G4int code = 0;
  if (particle_name == "mu-") code = 1;
  if (particle_name == "e+") code = 2;
  if (particle_name == "e-") code = 3;
  if (particle_name == "gamma") code = 4;
   
   File << event_number << " " << code << " " << pre_x << " " << pre_y << " " << pre_z<< " " << kenergy << G4endl;
   
   File.setf(mode,std::ios::floatfield);
   File.precision(prec);
   */
}

void EventAction::ExtraInformation(G4ThreeVector MomentumDirection, 				    G4ThreeVector HitPoint, G4double kEnergy, G4int t){
  std::ofstream File("mountain_intersection.txt", std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
   //G4cout << t << G4endl;
   if (t == 0){
   	if (state == 1) {
   	//G4cout << "ci entro qua" << G4endl;
   	File << " " << kEnergy << " " << MomentumDirection << " " << HitPoint << G4endl; state = 0;
   	//File << " " << kEnergy; state = 0;
   	}

   }
   else if (t == 1) {
   	if (state == 0) {File << event_number-1 << " " << hitmont+1 << " " << kEnergy << " " << MomentumDirection << " " << HitPoint; 
   	state++; hitmont++;
   	}
   }
   else G4cout << "c'è qualcosa che non va" << G4endl;
   
   
   File.setf(mode,std::ios::floatfield);
   File.precision(prec);

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*evt)
{
  //pass informations to Run
  //
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  	
  for (G4int k=0; k<kLayerMax; k++) {
     run->SumEvents_1(k,EtotLayer[k],EvisLayer[k]);   
  }
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->FillH1(1,EtotCalor);
  analysisManager->FillH1(2,EvisCalor);
  
  G4double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
  G4double Eleak = Ebeam - EtotCalor;
  run->SumEvents_2(EtotCalor,EvisCalor,Eleak);
  
  std::ofstream File("mountain_intersection.txt", std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
  
  //if (hitmont == 0) File << event_number-1 << " " << hitmont << G4endl;
  //else {
  //G4cout << "evento " << event_number-1 << " ++++++++++++++++++++++++ " << G4endl;
  	if (state == 1) File << " " << 0.0 << G4endl;
  	//}
   File.setf(mode,std::ios::floatfield);
   File.precision(prec);
  
  state = 0;
  hitmont = 0;
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     G4int kFiber = it->first;
	 G4int iFiber = kFiber%1000;
     G4double Evis = it->second;
	 analysisManager->FillH1(5,iFiber+0.5,Evis);
  }
    
  //write fired fibers on a file
  //
  WriteFibers(evt); 
}

/*void EventAction::StepInformation(G4String* p_name, G4double* kEnergy, size_t size)
{
  secondary_type = p_name;
  secondary_kEnergy = kEnergy;
  n_secondary = size;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        


void EventAction::WriteFibers(const G4Event* evt)
{
  // event is appended on a file
  //
  G4String name = G4AnalysisManager::Instance()->GetFileName();
  G4String fileName = name + "-" + identity + "_fibers.txt";
  
  std::ofstream File(fileName, std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
    
  //write event number  
  //
  //File << evt->GetEventID() << G4endl;
  
  
  //gun particle informations
  //
  G4ParticleGun* gun = primary->GetParticleGun();
  G4double ekin = gun->GetParticleEnergy();
  G4ThreeVector direction = gun->GetParticleMomentumDirection();
  G4ThreeVector position  = gun->GetParticlePosition();
  G4String particle = gun->GetParticleDefinition()->GetParticleName();
  //File << ekin << " " << direction << " " << position << G4endl;  
  //File << particle << G4endl;
  //if (EvisFiber.size()!=0){
  	File << evt->GetEventID() << " ";
  	File << ekin << " " << direction.x() << " " << direction.y() << " " << direction.z() << " " << position.x() << " " << position.y() <<  " " << position.z() << " ";
  
  
  //write fibers
  //
  //File << EvisFiber.size() << G4endl;
  	File << EvisFiber.size() << " ";
  //
  	std::map<G4int,G4double>::iterator it;         
  	for (it = EvisFiber.begin(); it != EvisFiber.end(); it++) {
     	  G4int kFiber = it->first;
     //G4double Evis = it->second;
     //File << " " << std::setw(7) << kFiber << " "<< std::setw(10) << Evis
            //<< G4endl;
     	  File << kFiber << " " ;
  	}
  
  /*if (n_secondary){
  G4int n = (G4int)n_secondary;
  for (G4int h = 0; h < n; h++)
  {
  	File << secondary_type[h] << " " << secondary_kEnergy[h] << G4endl;
  }}*/
           
  	File << G4endl;
  //}  
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

