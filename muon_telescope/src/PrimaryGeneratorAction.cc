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
      G4double seed = time(NULL);
      G4Random::setTheSeed(seed);
      G4double i = G4UniformRand();
      G4double j = G4UniformRand();
      G4double theta = acos(sqrt(i));
      G4double phi  = 2*M_PI*j;
      G4double normDirX   = cos(phi) * sin(theta);
      G4double normDirY   = sin(phi)* sin(theta);
      G4double normDirZ   = -cos(theta);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(normDirX, normDirY,normDirZ)); // inclinati
      //particleGun->SetParticleMomentumDirection(G4ThreeVector(1, 0, 0)); // tutti dritti
      particleGun->SetParticleEnergy(10.*GeV);
      
      G4double l = G4UniformRand();
      G4double k = G4UniformRand();
      G4double theta2 = acos(k);
      G4double phi2  = 2*M_PI*l;
      
      G4double positionX = (Detector->GetWorldSizeR())*sin(theta2)*cos(phi2);
      G4double positionY = (Detector->GetWorldSizeR())*sin(theta2)*sin(phi2);
      G4double positionZ = (Detector->GetWorldSizeR())*cos(theta2);
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
      //G4double seed = time(NULL);
      //G4Random::setTheSeed(seed);
      G4double i = G4UniformRand();
      G4double j = G4UniformRand();
      G4double theta = acos(sqrt(i));
      G4double phi  = 2*M_PI*j;
      G4double normDirX   = cos(phi) * sin(theta);
      G4double normDirY   = sin(phi)* sin(theta);
      G4double normDirZ   = -cos(theta);
      //particleGun->SetParticleMomentumDirection(G4ThreeVector(-4.914e-01, -4.646e-01,-7.367e-01));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(normDirX, normDirY,normDirZ)); // inclinati
      //particleGun->SetParticleMomentumDirection(G4ThreeVector(-1, 0., 0.)); // tutti dritti
      
      G4double E = G4UniformRand();
      //G4double energy = (ttt+1)*(0.5*GeV)+400*GeV;
      G4double energy = (pow(E, -(1/2.7)))*400*GeV;
      particleGun->SetParticleEnergy(energy);
      
      G4double l = G4UniformRand();
      G4double k = G4UniformRand();
      G4double theta2 = acos(k);
      G4double phi2  = 2*M_PI*l;
      
      G4double positionX = (Detector->GetWorldSizeR())*sin(theta2)*cos(phi2);
      G4double positionY = (Detector->GetWorldSizeR())*sin(theta2)*sin(phi2);
      G4double positionZ = (Detector->GetWorldSizeR())*cos(theta2);
      //particleGun->SetParticlePosition(G4ThreeVector(8.474e+04,1.486e+05,1.274e+06));
      particleGun->SetParticlePosition(G4ThreeVector(positionX,positionY,positionZ));
      //particleGun->SetParticlePosition(G4ThreeVector((0.4*Detector->GetWorldSizeR()),0.,50+(0.5*(0.8*Detector->GetWorldSizeR())/1.2)));
      particleGun->GeneratePrimaryVertex(anEvent);
      //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

