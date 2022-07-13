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

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include <math.h>

#define _USE_MATH_DEFINES

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:Detector(det)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  //beam = 0.*mm;
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematic()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4String particleName;
      G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="mu-");
      particleGun->SetParticleDefinition(particle);
      G4double i = G4UniformRand();
      G4double j = G4UniformRand();
      G4double theta = 0.5*M_PI*i + (M_PI*0.25);
      G4double phi  = acos(2*cos(0.25*M_PI)*j -cos(0.25*M_PI));
      G4double normDirY   = cos(phi)* sin(theta);
      G4double normDirX   = sin(phi) * sin(theta);
      G4double normDirZ   = cos(theta);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(normDirX, normDirY,normDirZ)); // inclinati
      //particleGun->SetParticleMomentumDirection(G4ThreeVector(1, 0, 0)); // tutti dritti
      particleGun->SetParticleEnergy(1.*GeV);
      G4double positionX = 0.5*(Detector->GetWorldSizeX())+60;
      G4double positionY = ((Detector->GetCalorSizeYZ())*G4UniformRand())-0.5*(Detector-> GetCalorSizeYZ());
      G4double positionZ = ((Detector->GetCalorSizeYZ())*G4UniformRand())-0.5*(Detector-> GetCalorSizeYZ());
      particleGun->SetParticlePosition(G4ThreeVector(positionX,positionY,positionZ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //this function is called at the begining of event
  //
  //randomize beam, if requested.
  //
  /*if (beam > 0.) 
    {

      G4ThreeVector position = particleGun->GetParticlePosition();    
      G4double maxYZ = 0.49*(Detector->GetCalorSizeYZ());
      G4double x0 = position.x();
      G4double y0 = position.y() + (G4UniformRand()-0.5)*beam;
      G4double z0 = position.z() + (G4UniformRand()-0.5)*beam;
      if (std::abs(y0) > maxYZ) y0 = maxYZ;
      if (std::abs(z0) > maxYZ) z0 = maxYZ;      
      particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticlePosition(position);      
    }*/
  //else  {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4String particleName;
      G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="mu-");
      particleGun->SetParticleDefinition(particle);
      G4double i = G4UniformRand();
      G4double j = G4UniformRand();
      //G4double theta = acos(((1-0.5*sqrt(2))*j+(0.5*sqrt(2))));
      //G4double phi  = 2*M_PI*i;
      G4double theta = 20*deg;
      G4double phi  = 50*deg;
      G4double normDirX   = cos(phi)* sin(theta);
      G4double normDirY   = sin(phi) * sin(theta);
      G4double normDirZ   = cos(theta);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-normDirX, -normDirY, -normDirZ)); // inclinati
      //particleGun->SetParticleMomentumDirection(G4ThreeVector(0, -0.5, -0.5*sqrt(3))); // tutti dritti
      particleGun->SetParticleEnergy(10.*GeV);
      G4double positionZ = 0.5*(Detector->GetWorldSizeX())-40;
      G4double positionY = ((Detector->GetCalorSizeYZ())*G4UniformRand())-0.5*(Detector-> GetCalorSizeYZ());
      G4double positionX = ((Detector->GetCalorSizeYZ())*G4UniformRand())-0.5*(Detector-> GetCalorSizeYZ());
      particleGun->SetParticlePosition(G4ThreeVector(0.5*sqrt(2)*positionZ*cos(55*deg),0.5*sqrt(2)*positionZ*sin(55*deg),0.5*sqrt(2)*positionZ));
      //particleGun->SetParticlePosition(G4ThreeVector(0.,0., positionZ));
      particleGun->GeneratePrimaryVertex(anEvent);
      //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

