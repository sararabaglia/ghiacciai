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

#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4EllipticalCone.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fiberMat(0),lvol_subfiber(0),lvol_fiber(0), absorberMat(0),lvol_layer(0),
 moduleMat(0),lvol_module(0), calorimeterMat(0),lvol_calorimeter(0), mountainMat(0), lvol_mountain(0), 
 iceMat(0), lvol_ice(0), lvol_ice_top(0), worldMat(0), pvol_world(0), defaultMat(0)
{
  // materials
  DefineMaterials();
  
  // default parameter values of calorimeter
  //
  subFiberDiameter    = (3)*mm;
  fiberDiameter       = (9)*mm; 	//1.08*mm
  nbOfFibers          = 110;		//490
  distanceInterFibers = (9)*mm;	//1.35*mm
  layerThickness      = (9)*mm;	//1.68*mm  
  milledLayer         = (180.00)*mm;    //1.40*mm ?
  nbOfLayers          = 2;		    //10
  nbOfModules         = 5;		    //9
  nbDetectors	      = 2600*m;		//12693 per circa 2km 8250 per 1.3km
     
  fiberLength         = (nbOfFibers)*distanceInterFibers;	//662.175*mm
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // define Elements
  //
  G4Element* H  = new G4Element("Hydrogen","H", 1,  1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  "C", 6, 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen","N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O", 8, 16.00*g/mole);
  G4Element* Si  = new G4Element("Silicon",  "Si", 14, 28.09*g/mole);
  G4Element* Al  = new G4Element("Aluminium",  "Al", 13, 27.00*g/mole);
  G4Element* Fe  = new G4Element("Iron",  "Fe", 26, 55.85*g/mole);
  G4Element* Ca  = new G4Element("Calcium",  "Ca", 20, 40.08*g/mole);
  G4Element* Na  = new G4Element("Sodium",  "Na", 11, 22.99*g/mole);
  G4Element* K  = new G4Element("Potassium",  "K", 19, 39.10*g/mole);
  G4Element* Mg  = new G4Element("Magnesium",  "Mg", 12, 24.31*g/mole);

  G4int natoms, ncomponents;
  G4double density, massfraction;				     

  // Lead
  //
  G4Material* Pb =   
  new G4Material("Lead", 82., 207.20*g/mole, density= 0.98*11.20*g/cm3);

  // Scintillator
  //
  G4Material* Sci = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
  
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // Air
  //
  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, massfraction=70*perCent);
  Air->AddElement(O, massfraction=30.*perCent);
  
  //Rock
  //
  G4Material* Rock = 
  new G4Material("Rock", density= 2.8*g/cm3, ncomponents=8);
  Rock->AddElement(O, massfraction=48.1*perCent);
  Rock->AddElement(Si, massfraction=27.7*perCent);
  Rock->AddElement(Al, massfraction=8.1*perCent);
  Rock->AddElement(Fe, massfraction=5.*perCent);
  Rock->AddElement(Ca, massfraction=3.6*perCent);
  Rock->AddElement(Na, massfraction=2.8*perCent);
  Rock->AddElement(K, massfraction=2.6*perCent);
  Rock->AddElement(Mg, massfraction=2.1*perCent);

  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Vacuum =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  
  //Water
  //
  G4Material* Water = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  Water -> AddElement(H, natoms=2);
  Water -> AddElement(O, natoms=1);
  
  //Ice
  G4Material* Ice = 
  new G4Material("Ice", density= 0.917*g/cm3, ncomponents=2);
  Ice -> AddElement(H, natoms=2);
  Ice -> AddElement(O, natoms=1);
  
  //attribute materials
  //
  defaultMat     = Air;  
  fiberMat       = Sci;
  absorberMat    = Water;
  moduleMat      = defaultMat;
  calorimeterMat = defaultMat;
  worldMat       = Vacuum;
  mountainMat	 = Rock;
  iceMat	 = Ice;

  // print table
  //      
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  //subfibers
  //
  G4Tubs*
  svol_subfiber = new G4Tubs("subfiber",		//name
                         0*mm, 0.5*subFiberDiameter,	//r1, r2
			 0.5*fiberLength,		//half-length 
			 0., twopi);			//theta1, theta2
			 
  lvol_subfiber = new G4LogicalVolume(svol_subfiber,	//solid
                                   fiberMat,		//material
                                   "subfiber");		//name
  
  // fibers
  //
  G4Tubs*
  svol_fiber = new G4Tubs("fiber",			//name
                         0*mm, 0.5*fiberDiameter,	//r1, r2
			 0.5*fiberLength,		//half-length 
			 0., twopi);			//theta1, theta2
			 
  lvol_fiber = new G4LogicalVolume(svol_fiber,		//solid
                                   absorberMat,		//material
                                   "fiber");		//name
				   
  
  //put subfibers within fiber
  //
  G4double Xcenter = 0.;
  G4double Ycenter = -subFiberDiameter;
  G4int n_subfiber = 0;
  new G4PVPlacement(0,		   			//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber);
  
  Xcenter = 0.;
  Ycenter = 0.;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+1);
  
  Xcenter = 0.;
  Ycenter = +subFiberDiameter;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+2);
  
  Xcenter = -2.58;
  Ycenter = -0.5*subFiberDiameter;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+3);
  
  Xcenter = -2.58;
  Ycenter = +0.5*subFiberDiameter;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+4);
  
  Xcenter = +2.58;
  Ycenter = -0.5*subFiberDiameter;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+5);
  
  Xcenter = +2.58;
  Ycenter = +0.5*subFiberDiameter;
  new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_subfiber,     		//logical volume	
                      "subfiber",	   			//name
                      lvol_fiber,        		//mother
                      false,             		//no boulean operat
                      n_subfiber+6);
  
  // layer
  //
  G4double sizeX = layerThickness;
  G4double sizeY = distanceInterFibers*nbOfFibers;
  G4double sizeZ = fiberLength;
  
  G4Box*      
  svol_layer = new G4Box("layer",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_layer = new G4LogicalVolume(svol_layer,		//solid
                                   absorberMat,		//material
                                   "layer");		//name

  // put fibers within layer
  //
  Xcenter = 0.;
  Ycenter = -0.5*(sizeY + distanceInterFibers);		//sommo distanceInterFibers perchè nel loop prima di piazzare la prima fibra, incremento di distanceInterFibers
  
  for (G4int k=0; k<nbOfFibers; k++) {
    Ycenter += distanceInterFibers;
    new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_fiber,     		   	//logical volume	
                      "fiber",	   			//name
                      lvol_layer,        		//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number

  }
				   
  // modules
  //
  moduleThickness = layerThickness*nbOfLayers + milledLayer;       
  sizeX = moduleThickness;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_module = new G4Box("module",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_module = new G4LogicalVolume(svol_module,	//solid
                                   defaultMat,		//material
                                   "module");		//name

  // put layers within module
  //
  Xcenter   = -0.5*moduleThickness - 0.5*layerThickness;
  //Xcenter = -0.5*(nbOfLayers+1)*layerThickness;
  //Ycenter =  0.25*distanceInterFibers;
  
  for (G4int k=0; k<nbOfLayers; k++) {
    Xcenter += layerThickness;
    //Ycenter  = - Ycenter;
    G4RotationMatrix rotm;                    //rotation matrix to place modules    
    if ((k+1)%2 == 0) rotm.rotateX(90*deg);
	G4Transform3D transform(rotm, G4ThreeVector(Xcenter,0.,0.));
    new G4PVPlacement(transform,    //position
                      lvol_layer,     		   	//logical volume	
                      "layer",	   			//name
                      lvol_module,        		//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number

  }			   				   

  // calorimeter
  //
  calorThickness = moduleThickness*nbOfModules;
  sizeX = calorThickness;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_calorimeter = new G4Box("calorimeter",		//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_calorimeter = new G4LogicalVolume(svol_calorimeter,	//solid
                                   calorimeterMat,		//material
                                   "calorimeter");		//name  

  // put modules inside calorimeter
  //  
  Xcenter = -0.5*(calorThickness + moduleThickness);
  

  for (G4int k=0; k<nbOfModules; k++) {
    Xcenter += moduleThickness;		    
    new G4PVPlacement(0,
    		      G4ThreeVector(Xcenter,0.,0.), //rotation+position
                      lvol_module,	     		//logical volume	
                      "module", 	   		    //name
                      lvol_calorimeter,        	//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number
  }
  
  
  //G4double R = 10*m;
  //mountain
  /*G4EllipticalCone*      
  svol_mountain = new G4EllipticalCone("mountain",		//name
                  2/R*0.5, 2/R*0.5, R*0.5,R*0.5);	//size


  lvol_mountain = new G4LogicalVolume(svol_mountain,	//solid
                                   calorimeterMat,		//material
                                   "mountain");		//name
  */
  
  //MONTAGNA + GHIACCIO
  
  G4Cons*      
  svol_mountain = new G4Cons("mountain", 0.,		//name
                  nbDetectors*0.2, 0., 0., nbDetectors*0.3, 0., twopi);

  lvol_mountain = new G4LogicalVolume(svol_mountain,	//solid
                                   mountainMat,		//material
                                   "mountain");
   
    
  /*G4Cons*      
  svol_ice = new G4Cons("ice", R*0.24,		//name
                  R*0.30, R*0.04, R*0.05, R*0.25, 0., twopi);

  lvol_ice = new G4LogicalVolume(svol_ice,	//solid
                                   iceMat,		//material
                                   "ice");
   
  
  G4Cons*      
  svol_ice_top = new G4Cons("ice_top", 0.,		//name
                  R*0.05, 0., 0., R*0.05, 0., twopi);

  lvol_ice_top = new G4LogicalVolume(svol_ice_top,	//solid
                                   iceMat,		//material
                                   "ice_top");
  
  new G4PVPlacement(0,
    		      G4ThreeVector(0.,0., R*0.05), //rotation+position
                      lvol_ice,	     		//logical volume	
                      "ice", 	   		    //name
                      lvol_mountain,        	//mother
                      false,             		//no boulean operat
                      0);
  
  
  new G4PVPlacement(0,
    		      G4ThreeVector(0.,0., 2*R*0.15+R*0.05), //rotation+position
                      lvol_ice_top,	     		//logical volume	
                      "ice_top", 	   		    //name
                      lvol_mountain,        	//mother
                      false,             		//no boulean operat
                      0);
  
  */
  //MONTAGNA sfera per mondo
  sizeX = nbDetectors;                            
  G4Sphere*      
  svol_ice = new G4Sphere("mountain",			//name
                  0, sizeX, 0*deg, 110*deg, 0*deg, 90*deg);	//size

  lvol_ice = new G4LogicalVolume(svol_ice,		//solid
                                   worldMat,		//material
                                   "mountain");
  
  // world
  //
  
  sizeX = nbDetectors;
 
  
  worldSizeR = sizeX;
  
  G4Sphere*      
  svol_world = new G4Sphere("world",			//name
                  0, sizeX, 0*deg, 110*deg, 0*deg, 90*deg);	//size

  lvol_world = new G4LogicalVolume(svol_world,		//solid
                                   worldMat,		//material
                                   "world");		//name 
				    
  pvol_world = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 lvol_world,		//logical volume
                                 "world",		//name
                                 0,			//mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //put calorimeter in world
  G4RotationMatrix rotm1;
  rotm1.rotateY(135*deg);
  rotm1.rotateZ(55*deg);
  
  G4Transform3D transform1(rotm1, G4ThreeVector(sqrt(2)*calorThickness*cos(55*deg),sqrt(2)*calorThickness*sin(55*deg),0.5*sqrt(2)*calorThickness));
  
  new G4PVPlacement(transform1,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    1);

  /*for (G4int i = 0; i<nbDetectors; i++)
  {
  	
  	if (i!=0) rotm1.rotateZ(alpha);
  	xx = apotemaCorretta*cos(i*alpha);
  	yy = apotemaCorretta*sin(i*alpha);
  	G4Transform3D transform1(rotm1, G4ThreeVector(xx,yy,(sqrt(2)*calorThickness/2)+10));
  	
  	new G4PVPlacement(transform1,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    i+1);
  	
  }*/
  
  //MONTAGNA
  new G4PVPlacement(0,
    		      G4ThreeVector(0.6*(nbDetectors)*cos(55*deg),0.6*(nbDetectors)*sin(55*deg), 0.3*nbDetectors), //rotation+position
                      lvol_mountain,	     		//logical volume	
                      "mountain", 	   		    //name
                      lvol_world,        	//mother
                      false,             		//no boulean operat
                      0);
  
  /*new G4PVPlacement(0,
    		      G4ThreeVector(0.,0., 0.), //rotation+position
                      lvol_ice,	     		//logical volume	
                      "mountain", 	   		    //name
                      lvol_world,        	//mother
                      false,             		//no boulean operat
                      0);
  */
  /*
  rotm1.rotateZ(0*deg);
  G4Transform3D transform2(rotm1, G4ThreeVector(sqrt(2)*calorThickness,0,sqrt(2)*calorThickness/2));
  new G4PVPlacement(transform2,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    0);				//copy number
		    				 
  G4RotationMatrix rotm2;
  rotm2.rotateY(45*deg);
  rotm2.rotateZ(180*deg);
  
  new G4PVPlacement(transform1,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    1);
  */
  /*G4RotationMatrix rotm3;
  rotm3.rotateY(-45*deg);
  rotm3.rotateZ(0*deg);
  G4Transform3D transform3(rotm3, G4ThreeVector(0,0,calorThickness/2+10));
  new G4PVPlacement(transform3,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    1);*/
  /*
  G4RotationMatrix rotm4;
  rotm4.rotateY(45*deg);
  rotm4.rotateZ(270*deg);
  G4Transform3D transform4(rotm4, G4ThreeVector(0,-sqrt(2)*calorThickness,sqrt(2)*calorThickness/2));
  new G4PVPlacement(transform4,		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    1);*/
  
  PrintCalorParameters();
  
  // Visualization attributes
  //
  lvol_fiber->SetVisAttributes (G4VisAttributes::GetInvisible());  
  lvol_layer->SetVisAttributes (G4VisAttributes::GetInvisible());
  lvol_world->SetVisAttributes (G4VisAttributes::GetInvisible());
    
  //always return the physical World
  //
  return pvol_world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UnitsTable.hh"

void DetectorConstruction::PrintCalorParameters()
{
  /*G4cout << "\n-------------------------------------------------------------"
     << "\n ---> The calorimeter is " << nbOfModules << " Modules"
     << "\n ---> A Module is " << nbOfLayers << " Layers + 1 milled Layer";
     
  //G4cout  
     //<< "\n ---> A Layer is " << G4BestUnit(layerThickness,"Length")  
     //<< " thickness of " << absorberMat->GetName();    
     
  G4cout 
     << "\n ---> A Layer includes " << nbOfFibers << " fibers";
  
  G4cout 
     << "\n ---> A Fiber includes " << 7 << " subfibers of " 
     << fiberMat->GetName();
     
  G4cout 
     << "\n      ---> diameter : " << G4BestUnit(fiberDiameter,"Length")
     << "\n      ---> length   : " << G4BestUnit(fiberLength,"Length")
     << "\n      ---> distance : " << G4BestUnit(distanceInterFibers,"Length");
     
  G4cout  
     << "\n ---> The milled Layer is " << G4BestUnit(milledLayer,"Length")  
     << " thickness of " << absorberMat->GetName();
     
  G4cout 
   << "\n\n ---> Module thickness " << G4BestUnit(moduleThickness,"Length");
  
  G4cout 
   << "\n\n ---> Total calor thickness " << G4BestUnit(calorThickness,"Length")
   <<   "\n      Tranverse size        " << G4BestUnit(fiberLength,"Length");

  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
