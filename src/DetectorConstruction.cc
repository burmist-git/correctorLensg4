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

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"

#include "globals.hh"

#include "math.h"

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //bool overlapsChecking = false;
  bool overlapsChecking = true;
  bool buildSmallBoxIncenter = false;
  //bool buildSingleSiPMArray = true;
  //bool buildSingleSiPMArray = false;
  
  //     
  // World
  //
  G4double world_sizeXY = 100*cm;
  G4double world_sizeZ  = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,world_mat,"World");                                   
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                     //no rotation
						   G4ThreeVector(),       //at (0,0,0)
						   logicWorld,            //its logical volume
						   "World",               //its name
						   0,                     //its mother  volume
						   false,                 //no boolean operation
						   0,                     //copy number
						   overlapsChecking);     //overlaps checking
  ////////////////////////////////////////////
  //load_fit_parameters_correction_lens_profile_function("../ellipse_Arc_Splines/lens_profile/lens_profile_n1.4_long_pol8_sym_par.dat");
  load_fit_parameters_correction_lens_profile_function("/home/burmist/home2/work/POEMMA/geant4/terzinag4_ana/hist_v02.00.00b/lens_n_reff_1.4.dat");
  //lens_profile_function(0);
  ////////////////////////////////////////////
  //
  G4double lens_ring_R_min;
  G4double lens_ring_R_max;
  G4double cons_ring_R_max;
  G4double cons_ring_R_min;
  G4double cons_ring_dR = 1.0;
  //
  G4double cons_ring_R_min_not_subt;
  G4double cons_ring_R_max_not_subt;
  G4double cons_ring_h_not_subt;
  G4double cons_ring_dR_not_subt = 0.0;
  //
  G4double cons_ring_h;
  bool cons_sub_not_add = true; 
  //
  int nn_ring = 100;
  G4double dl_ring = lens_profile_xmax/nn_ring;
  G4double lens_ring_openning = 360.0;//deg
  //
  G4double x0_ring;
  G4double x2_ring;
  G4double y0_ring;
  G4double y2_ring;
  G4double hcons_ring;
  G4double htubs_ring;
  //
  lens_thickness = -999.0;
  for(int i = 0;i<nn_ring;i++){
    G4double xval = lens_profile_xmax/(nn_ring-1)*i;
    G4double yval = lens_profile_function(xval);
    if(yval>lens_thickness)
      lens_thickness = yval;
  }
  
  //
  // Lens envelope 
  //
  G4double lens_envelope_x0 = 0.0;
  G4double lens_envelope_y0 = 0.0;
  G4double lens_envelope_z0 = 300.0;
  G4double lens_envelope_dl = 0.3;
  G4double lens_envelope_R_min = 0.0;
  G4double lens_envelope_R_max = (lens_profile_xmax - lens_profile_xmin)/2.0 + lens_envelope_dl;
  G4double lens_envelope_thickness = lens_thickness + lens_envelope_dl;
  //
  //G4cout<<"lens_envelope_R_max     "<<lens_envelope_R_max<<G4endl
  //	<<"lens_envelope_thickness "<<lens_envelope_thickness<<G4endl;
  //
  G4Tubs *lens_envelope_solid = new G4Tubs("lens_envelope_solid",       //name
					   lens_envelope_R_min,         //pRMin
					   lens_envelope_R_max,         //pRMax
					   lens_envelope_thickness/2.0, //pDz
					   0,                           //pSPhi
					   twopi);                      //pDPhi
  G4LogicalVolume* lens_envelope_logic = new G4LogicalVolume(lens_envelope_solid,world_mat,"lens_envelope_logic");                        
  new G4PVPlacement(0,                               //no rotation
		    G4ThreeVector(lens_envelope_x0,
				  lens_envelope_y0,
				  lens_envelope_z0), //at (0,0,0)
		    lens_envelope_logic,             //its logical volume
		    "lens_envelope_phys",            //its name
		    logicWorld,                      //its mother  volume
		    false,                           //no boolean operation
		    0,                               //copy number
		    overlapsChecking);               //overlaps checking
  //
  //
  //

  //G4cout<<"lens_thickness "<<lens_thickness<<G4endl;
  //
  for(int i = 0;i<nn_ring;i++){
    lens_ring_R_min = dl_ring*i;
    lens_ring_R_max = dl_ring*(i+1);
    x0_ring = lens_ring_R_min;
    x2_ring = lens_ring_R_max;
    y0_ring = lens_profile_function(x0_ring);
    y2_ring = lens_profile_function(x2_ring);
    //
    cons_sub_not_add = get_cons_ring_par(x0_ring, y0_ring, x2_ring, y2_ring, cons_ring_dR, cons_ring_R_min, cons_ring_R_max, cons_ring_h);
    if(cons_sub_not_add){
      //
      htubs_ring = y2_ring;
      hcons_ring = y2_ring - y0_ring;
      //G4cout<<"x0_ring = "<<x0_ring<<G4endl
      //    <<"y0_ring = "<<y0_ring<<G4endl
      //    <<"x2_ring = "<<x2_ring<<G4endl
      //    <<"y2_ring = "<<y2_ring<<G4endl;
      G4Tubs *lens_ring_solid = new G4Tubs("lens_ring_solid",               //name
					   lens_ring_R_min,                 //pRMin
					   lens_ring_R_max,                 //pRMax
					   htubs_ring/2.0,                  //pDz
					   0,                               //pSPhi
					   twopi/360.0*lens_ring_openning); //pDPhi    
      G4Cons *cons_solid = new G4Cons("cons_solid",       //pName
				      0.0,                //pRmin1,
				      cons_ring_R_max,    //pRmax1,
				      0.0,                //pRmin2,
				      cons_ring_R_min,    //pRmax2,
				      cons_ring_h/2.0,    //pDz,
				      0.0,                //pSPhi,
				      twopi/360.0*lens_ring_openning); //pDPhi
      //
      G4RotationMatrix* rotMatrix = new G4RotationMatrix();
      //G4ThreeVector transVector(0.0,0.0,(-y2/2.0 + y0/2.0 - hcons_ring/2.0));
      G4ThreeVector transVector(0.0,0.0, -htubs_ring/2.0 + hcons_ring/2.0);
      G4SubtractionSolid *ring_m_cons_solid = new G4SubtractionSolid("ring_m_cons_solid",lens_ring_solid,cons_solid,rotMatrix,transVector);
      G4LogicalVolume *ring_m_cons_logical = new G4LogicalVolume(ring_m_cons_solid,world_mat,"ring_m_cons_logical");
      G4ThreeVector lensRingPosition(0.0,0.0,-htubs_ring/2.0 + lens_envelope_thickness/2.0);
      new G4PVPlacement(0,                      //rotation
			lensRingPosition,       //position
			ring_m_cons_logical,    //its logical volume
			"ring_m_cons_physical", //its name
			lens_envelope_logic,    //its mother  volume
			false,                  //no boolean operation
			0,                      //copy number
			overlapsChecking);      //overlaps checking
      //
      //G4RotationMatrix Ra;
      //G4ThreeVector Ta;
      //G4Transform3D Tr;
      //G4LogicalVolume *lens_ring_logical = new G4LogicalVolume(lens_ring_solid,world_mat,"lens_ring_logical");
      //lens_ring_logical->GetName(); 
      //Ta.setX(0);
      //Ta.setY(0);
      //Ta.setZ(-htubs_ring/2.0);
      //Ra.rotateX(0.0);
      //Tr = G4Transform3D(Ra, Ta);
      //new G4PVPlacement(Tr,                   //Transformation
      //lens_ring_logical,    //its logical volume                                 
      //"lens_ring_physical", //its name
      //logicWorld,           //its mother  volume
      //false,                //no boolean operation
      //0);                   //copy number
      //Ra.rotateX(-0.0);
      ////
      //G4LogicalVolume *cons_logical = new G4LogicalVolume(cons_solid,world_mat,"cons_test_logical");
      //cons_logical->GetName(); 
      //Ta.setX(0);
      //Ta.setY(0);
      //Ta.setZ(-(htubs_ring - hcons_ring) - hcons_ring/2.0);
      //Ra.rotateX(0.0);
      //Tr = G4Transform3D(Ra, Ta);
      //new G4PVPlacement(Tr,                   //Transformation
      //cons_logical,    //its logical volume                                 
      //"cons_physical", //its name
      //logicWorld,           //its mother  volume
      //false,                //no boolean operation
      //0);                   //copy number
      //Ra.rotateX(-0.0);
    }
    else{
      get_cons_ring_par(x0_ring, y0_ring, x2_ring, y2_ring, cons_ring_dR_not_subt, cons_ring_R_min_not_subt, cons_ring_R_max_not_subt, cons_ring_h_not_subt);
      htubs_ring = y0_ring;
      hcons_ring = y0_ring - y2_ring;
      G4Tubs *lens_ring_solid = new G4Tubs("lens_ring_solid",  //name
					   lens_ring_R_min,    //pRMin
					   lens_ring_R_max,    //pRMax
					   htubs_ring/2.0,          //pDz
					   0,                  //pSPhi
					   twopi/360.0*lens_ring_openning); //pDPhi    
      G4LogicalVolume *lens_ring_logical = new G4LogicalVolume(lens_ring_solid,world_mat,"lens_ring_logical");
      lens_ring_logical->GetName(); 
      G4ThreeVector lensRingPosition(0.0,0.0,-htubs_ring/2.0 + lens_envelope_thickness/2.0);
      new G4PVPlacement(0,                    //rotation
			lensRingPosition,     //position
			lens_ring_logical,    //its logical volume
			"lens_ring_physical", //its name
			lens_envelope_logic,  //its mother  volume
			false,                //no boolean operation
			0,                    //copy number
			overlapsChecking);    //overlaps checking
      //
      G4Cons *cons_solid = new G4Cons("cons_solid",             //pName
				      cons_ring_R_min_not_subt, //pRmin1,
				      cons_ring_R_min_not_subt, //pRmax1,
				      cons_ring_R_min_not_subt, //pRmin2,
				      cons_ring_R_max_not_subt, //pRmax2,
				      cons_ring_h_not_subt/2.0, //pDz,
				      0.0,                      //pSPhi,
				      twopi/360.0*lens_ring_openning);       //pDPhi
      //
      G4LogicalVolume *cons_logical = new G4LogicalVolume(cons_solid,world_mat,"cons_logical");
      cons_logical->GetName(); 
      G4ThreeVector consRingPosition(0.0,0.0,-htubs_ring - cons_ring_h_not_subt/2.0 + lens_envelope_thickness/2.0);
      new G4PVPlacement(0,                   //rotation
			consRingPosition,    //position
			cons_logical,        //its logical volume
			"cons_physical",     //its name
			lens_envelope_logic, //its mother  volume
			false,               //no boolean operation
			0,                   //copy number
			overlapsChecking);   //overlaps checking
    }
  }
  //
  
  /*
  //
  G4Cons *cons_test_solid = new G4Cons("cons_test_solid",  //pName
				       0.0,                //pRmin1,
				       lens_ring_R_max,    //pRmax1,
				       0.0,                //pRmin2,
				       lens_ring_R_max-30, //pRmax2,
				       lens_thickness/2.0  //pDz,
				       0.0,                //pSPhi,
				       twopi);             //pDPhi
  G4LogicalVolume *cons_test_logical = new G4LogicalVolume(cons_test_solid,world_mat,"cons_test_solid");
  cons_test_logical->GetName(); 
  Ta.setX(0);
  Ta.setY(0);
  Ta.setZ(0);
  Ra.rotateX(0.0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *cons_test_physical = new G4PVPlacement(Tr,                   //Transformation
							    cons_test_logical,    //its logical volume                                 
							    "cons_test_physical", //its name
							    logicWorld,           //its mother  volume
							    false,                //no boolean operation
							    0);                   //copy number
  Ra.rotateX(-0.0);
  */
  
  //
  // Small box for orientation 
  //
  G4VSolid *boxsmall_solid = new G4Box("boxsmall_solid", 1.0*cm, 3.0*cm, 10.0*cm);
  G4LogicalVolume *boxsmall_logical = new G4LogicalVolume(boxsmall_solid,world_mat,"boxsmall_solid");
  new G4PVPlacement(0,                      //no rotation
		    G4ThreeVector(90*cm/2.0,90*cm/2.0,35*cm),       //at (0,0,0)
		    boxsmall_logical,       //its logical volume
		    "World",                //its name
		    logicWorld,             //its mother  volume
		    false,                  //no boolean operation
		    0,                      //copy number
		    overlapsChecking);      //overlaps checking

  //
  // Small box incenter 
  //
  G4VSolid *boxsmall_centre_solid = new G4Box("boxsmall_centre_solid", 1.0*mm, 1.0*mm, 1.0*mm);
  G4LogicalVolume *boxsmall_centre_logical = new G4LogicalVolume(boxsmall_centre_solid,world_mat,"boxsmall_centre_solid");
  if(buildSmallBoxIncenter)
    new G4PVPlacement(0,                       //no rotation
		      G4ThreeVector(),         //at (0,0,0)
		      boxsmall_centre_logical, //its logical volume
		      "World",                 //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      overlapsChecking);       //overlaps checking
  
  return physWorld;
}

void DetectorConstruction::load_fit_parameters_correction_lens_profile_function(std::string dataFileIn){
  std::ifstream fileIn(dataFileIn.c_str());
  G4double par;
  if (fileIn.is_open()){
    fileIn >> lens_profile_xmin;
    fileIn >> lens_profile_xmax;
    for(int i = 0 ;i<npar_pol8_symmetric;i++){
      fileIn >> par;
      par_pol8_symmetric[i] = par;
    }
    fileIn.close();
  }
  else G4cout<<"Unable to open file"<<G4endl;
}

G4double DetectorConstruction::lens_profile_function(G4double x){
  if(x < lens_profile_xmin || x > lens_profile_xmax )
    return -999.0;
  return par_pol8_symmetric[0] +
    par_pol8_symmetric[1]*pow(x,2) +
    par_pol8_symmetric[2]*pow(x,4) +
    par_pol8_symmetric[3]*pow(x,6) +
    par_pol8_symmetric[4]*pow(x,8);
}

bool DetectorConstruction::get_cons_ring_par(G4double x0, G4double y0, G4double x2, G4double y2, G4double dr, G4double &rMin, G4double &rMax, G4double &conh){
  if( (x2-x0)<0.0){
    G4cout<<"ERROR --> (x2-x0) < 0.0 "<<G4endl
	  <<"                  x2 = "<<x2<<G4endl
      	  <<"                  x0 = "<<x0<<G4endl;
    assert(0);
  }
  double k = (y2 - y0)/(x2-x0);
  double b = y0 - k*x0;
  rMin = x0 - dr;
  rMax = x2 + dr;
  double y0_dr = k*(x0 - dr) + b;
  double y2_dr = k*(x2 + dr) + b;
  conh = abs(y2_dr - y0_dr);
  if(k>0)
    return true;
  return false;
}
