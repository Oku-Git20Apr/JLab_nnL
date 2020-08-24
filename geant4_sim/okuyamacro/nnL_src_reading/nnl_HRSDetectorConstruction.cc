/*
  HRSDetectorConstruction.cc

  geometry for HRS Spectrometer

  fujita

  modifyed by KNS
*/

#include "HRSDetectorConstruction.hh"
#include "HRSGlobalSize.hh"

//=== Sensitive Detector ===//
#include "G4SDManager.hh"
#include "HRSTargetSD.hh"
#include "HRSVDSD.hh"
#include "EDC1SD.hh"
#include "EDC2SD.hh"
#include "HRSEH1SD.hh"
#include "HRSEH2SD.hh"
//=== Material Define ===//
#include "MaterialList.hh"

//=== Geometry Define ===//
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

//=== Beam Line Element ===//
#include "Target.hh"
#include "SplitterMagnet.hh"
#include "SeptumMagnet.hh"
#include "APEXSeptum.hh"
#include "Q1Magnet.hh"
#include "Q2Magnet.hh"
#include "Q3Magnet.hh"
#include "Q4Magnet.hh"
#include "DMagnet.hh"
#include "VirtualDetector.hh"
#include "VDetectorPB.hh"
#include "VDetectorQ1.hh"
#include "VDetectorQ2.hh"
#include "EDC1.hh"
#include "EDC2.hh"
#include "EHodo1.hh"
#include "EHodo2.hh"
#include "Collimator.hh"
#include "Collimator1.hh"
#include "CollimatorBox.hh"
//#include "TargetChamber.hh"
#include "SieveSlit.hh"


//=== ElectoMagnetic Field ===//
#include "HRSField.hh"
#include "MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"

//=== GDML parser include ===// by Aida 20171006
//
#include "G4GDMLParser.hh"
#include "HKSComponent.hh"

#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>

#include "HRSParamMan.hh"
//#include "HRSPrimaryGeneratorAction.hh"


//=== Program Start ===//
///////////////////////////////////////////////////
HRSDetectorConstruction::HRSDetectorConstruction()
: mList_(0), EMField_(0)
  ///////////////////////////////////////////////////
{
}

////////////////////////////////////////////////////
HRSDetectorConstruction::~HRSDetectorConstruction()
///////////////////////////////////////////////////
{
  delete mList_;
  delete EMField_;
}

/////////////////////////////////////////////////////////
MaterialList *HRSDetectorConstruction::DefineMaterials()
/////////////////////////////////////////////////////////
{
  MaterialList *ml = new MaterialList;

  return ml;
}

////////////////////////////////////////////////////////
G4VPhysicalVolume *HRSDetectorConstruction::Construct()
///////////////////////////////////////////////////////
{
  mList_ = DefineMaterials();

  G4VPhysicalVolume *world = ConstructPayload();
  //	G4GDMLParser parser;
  //	parser.Write("./gdml/World.gdml", world, false);

  return world;
}

///////////////////////////////////////////////////////
HRSField *HRSDetectorConstruction::MakeDetectorField()
///////////////////////////////////////////////////////
{
  return new HRSField();
}

///////////////////////////////////////////////////////////////
G4VPhysicalVolume *HRSDetectorConstruction::ConstructPayload()
//////////////////////////////////////////////////////////////
{ 
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();

  // ****** Parameters ****** //

  ////// APEX Septum /////
  const G4double E_Theta = 13.2*degree; //default
  G4double PIVOTtoQ1en = 169.*cm; // from pivot to Q1 entrance (fixed parameter)
  paramMan->SetE_Theta(E_Theta);
  G4RotationMatrix *rotHRS = new  G4RotationMatrix();
  rotHRS->rotateY(-E_Theta);
  G4SDManager *SDManager = G4SDManager::GetSDMpointer();
  if( !EMField_ ) EMField_ = MakeDetectorField();
  else            EMField_->cleanupElemList();
  const G4ThreeVector PIVOTpos(0.0*CLHEP::cm, 0.0*CLHEP::cm, 0.0*CLHEP::cm);
  G4AffineTransform *HRSTrans = new G4AffineTransform( rotHRS, PIVOTpos );
  
  G4cout<< "E_Theta [deg] = " << E_Theta/CLHEP::degree <<G4endl;
  G4cout<< "PIVOTpos [cm] = " << PIVOTpos.x()/CLHEP::cm << " "  << PIVOTpos.z()/CLHEP::cm <<G4endl;
  
  //=== World ===//
  G4Box *worldSolid = new G4Box( "World",WorldSizeX,
				 WorldSizeY,WorldSizeZ );
  G4LogicalVolume *worldLV =
    new G4LogicalVolume( worldSolid, mList_->Vacuum, "World LV");
  G4VPhysicalVolume *world = 
    new G4PVPlacement( 0, G4ThreeVector( 0.*CLHEP::cm, 0.*CLHEP::cm, 0.*CLHEP::cm ),
		       worldLV, "World", 0, false, 0);
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);

  G4Colour colourWorld(0., 1., 1.); // Cyan   
  //  G4VisAttributes *worldVisAtt = new G4VisAttributes(true, colourWorld);
  //  worldLV->SetVisAttributes(worldVisAtt); 

#if 1
  /////pivot /////suzuki
  G4Box *pivotSolid = new G4Box("piv_test", 100, 1, 1); //+-xyz
  G4LogicalVolume *pivotLV = 
    new G4LogicalVolume(pivotSolid, mList_->Vacuum, "pivotLV");
  G4VPhysicalVolume *pivotPV = 
    new G4PVPlacement( 0, G4ThreeVector( 0.*CLHEP::cm, 0.*CLHEP::cm, 0.*CLHEP::cm ), pivotLV, "pivotPV", worldLV, false, 0);

  G4VisAttributes *pivotVisAtt = new G4VisAttributes(true, colourWorld);
  pivotLV->SetVisAttributes(pivotVisAtt); 

  /////pivot2 /////suzuki
  G4VPhysicalVolume *pivot2PV = 
    new G4PVPlacement( 0, G4ThreeVector( 0.*CLHEP::cm, 0.*CLHEP::cm, 12.5*CLHEP::cm ), pivotLV, "pivot2PV", worldLV, false, 0);

  /////pivot3 /////suzuki
  G4VPhysicalVolume *pivot3PV = 
    new G4PVPlacement( 0, G4ThreeVector( 0.*CLHEP::cm, 0.*CLHEP::cm, -12.5*CLHEP::cm ), pivotLV, "pivot3PV", worldLV, false, 0);

#endif

#if 0

  ////// HALL A NIM VALUE (Q1 is old. Use Place() for Q1 placement)
  const G4double Q1length = 800.*CLHEP::mm; 
  const G4double Q2length = 1800.*CLHEP::mm; 
  const G4double Dlength = 6597.3*CLHEP::mm; // arc length
  const G4double Q3length = 1800.*CLHEP::mm; 

  const G4double Q1extoQ2en = 1250.*CLHEP::mm; // from Q1 exit to Q2 entrance
  const G4double Q2extoDen = 4420.*CLHEP::mm; // frpm Q2 exit to D entrance

#endif

#if 1

  // suzuki 20200107
  ////// HALL A NIM VALUE (Q1 is SOS Q1 and the center of SOS Q1 is same as old Q1. Use GDMLPlace() for Q1 Placement.)

  const G4double Q1length = 935.228*CLHEP::mm;
  PIVOTtoQ1en = (1690 - ( 935.228 - 800.0 )/2.0)*CLHEP::mm;
  const G4double Q2length = 1800.*CLHEP::mm; 
  const G4double Dlength = 6597.3*CLHEP::mm; // arc length
  const G4double Q3length = 1800.*CLHEP::mm; 

  const G4double Q1extoQ2en = (1250. - ( 935.228 - 800.0 )/2.0) *CLHEP::mm; // from Q1 exit to Q2 entrance
  const G4double Q2extoDen = 4420.*CLHEP::mm; // frpm Q2 exit to D entrance

#endif

  
#if 0

  //////  (suzuki 20200105)
  // Q1 entrance face is same to NIM Q1 but Q1-Q2 length is modified. Use GDMLPlace() for Q1 Placement.
  const G4double Q1length = 935.228*CLHEP::mm; 
  const G4double Q2length = 1800.*CLHEP::mm; 
  const G4double Dlength = 6597.3*CLHEP::mm; // arc length
  const G4double Q3length = 1800.*CLHEP::mm; 

  const G4double Q1extoQ2en = 1114.772*CLHEP::mm; // from Q1 exit to Q2 entrance
  const G4double Q2extoDen = 4420.*CLHEP::mm; // frpm Q2 exit to D entrance

#endif

  const G4double Q3entoVDC1 = 3570.*CLHEP::mm; //NIM VALUE
  //  const G4double Q3entoVDC1 = 4000.*CLHEP::mm; //test

  //==== SEP + HRS ========== // aida
  /* "VT" means the HRS spectrometer pivot */
  //  const G4double VTtoVD0 = 0.0*CLHEP::mm + 20.*CLHEP::cm; 
  const G4double VTtoVD3 = PIVOTtoQ1en; // from pivot to Q1 entrance
  const G4double VTtoQ1  = VTtoVD3 + Q1length/2.; 
  const G4double VTtoVD4 = VTtoVD3 + Q1length; // Q1 exit 
  const G4double VTtoVD5 = VTtoVD4 + Q1extoQ2en; // Q2 entrance 
  const G4double VTtoQ2  = VTtoVD5 + Q2length/2.; 
  const G4double VTtoVD6 = VTtoVD5 + Q2length; // Q2 exit
  const G4double VTtoVD7 = VTtoVD6 + Q2extoDen; // D entrance
  const G4double VTtoD   = VTtoVD7; // D entrance
  //  const G4double VTtoDen = VTtoVD5 + 0.05*CLHEP::cm + 330*CLHEP::cm ;
  const G4double R = 2.*840.*CLHEP::cm*sin(22.5*CLHEP::degree);
  const G4double VD8y = R*sin(22.5*CLHEP::degree); 
  const G4double VD8z = R*cos(22.5*CLHEP::degree)+VTtoVD7; 
  const G4double VD8toVD9 = 1500.*CLHEP::mm;
  const G4double VD8toQ3 = VD8toVD9 + Q3length/2.; 
  const G4double VD8toVD10 = VD8toVD9 + Q3length; 
  const G4double VDCGap = 335*CLHEP::mm;
  const G4double VD8toVD11 = VD8toVD10 + Q3entoVDC1 + VDCGap/2.;//org
  const G4double VD8toVD12 = VD8toVD10 + Q3entoVDC1; // aida 
  const G4double VD8toVD13 = VD8toVD10 + Q3entoVDC1 + VDCGap; // aida
  const G4double VD8toS2 = VD8toVD12 + 3179.*CLHEP::mm; // suzuki (db_L.s2.dat)

  paramMan->SetpVD8(VD8y,VD8z);

  //=== Target ===//
  const G4double TargetThickness = paramMan->GetTargetThickness();
  std::string TargetMaterial = paramMan->GetTargetMaterial();
  G4Material *Target = mList_->Vacuum;
  
  if( TargetMaterial == "Vanadium"){
    Target = mList_->V51;
  }
  else if( TargetMaterial == "Yttrium"){
    Target = mList_->Y89;
  }
  else if( TargetMaterial == "Lead"){
    Target = mList_->Pb208;
  }
  else if( TargetMaterial == "Carbon"){
    Target = mList_->C12;
  }
  else if( TargetMaterial == "Silicon"){
    Target = mList_->Si28;
  }
  else if( TargetMaterial == "Polyethylene"){
    Target = mList_->CH2;
  }
  else if( TargetMaterial == "Lithium"){
    Target = mList_->Li7;
  }
  else if( TargetMaterial == "Beryllium"){
    Target = mList_->Be9;
  }
  else if( TargetMaterial == "Boron"){
    Target = mList_->B10;
  }
  else if( TargetMaterial == "Calcium"){
    Target = mList_->Ca40;
  }
  else if( TargetMaterial == "Water"){
    Target = mList_->Water;
  }
  else if( TargetMaterial == "Cromium"){
    Target = mList_->Cr52;
  }
  else if( TargetMaterial == "Vacuum"){
    Target = mList_->Vacuum;
  }
  else {
    G4cout << "You choose strange target '" 
	   << TargetMaterial 
	   << "' !! I don't know such a target..." 
	   <<G4endl;
    G4cout << " --> No target will be used." <<G4endl;
    Target = mList_->Vacuum;
  }
  
  const G4double TargetDensity = Target->GetDensity();
  const G4double TargetSizeX = 10.*CLHEP::cm; //org 1cm
  const G4double TargetSizeY = 10.*CLHEP::cm;
  const G4double TargetSizeZ = (TargetThickness)/(TargetDensity/(CLHEP::g/cm3)*1000)*CLHEP::cm;
  G4RotationMatrix* rotTarget = new G4RotationMatrix();//Rotaion of target
  G4Box *targetSolid = new G4Box( "Target",
				  TargetSizeX/2,TargetSizeY/2,TargetSizeZ/2 );
  G4LogicalVolume *targetLV = 
    new G4LogicalVolume( targetSolid, Target, "Target LV");

  ///////////////////
  // Visualisation //
  ///////////////////
  G4Colour colourtarget(1., 0., 0.); // red
  G4VisAttributes* targetVisAtt = new G4VisAttributes(true, colourtarget);
  targetLV-> SetVisAttributes(targetVisAtt);

  if (TargetMaterial != "Vacuum"){
    //new G4PVPlacement( rotSplitter, 
    new G4PVPlacement( rotTarget, 
		       G4ThreeVector(0*cm, 0*cm, 0*cm),   
		       targetLV, "Target" , worldLV , false, 0);
  }
    
  //  int flag = paramMan->GetEventFlag();
  //  int EDCFlag = paramMan->GetEDCFlag();

#if 0 
  //=== Virtual Detector for Beamprofile ===//
  G4ThreeVector posVD0( 0.*cm, 0.*cm, VTtoVD0);
  //  posVD0.rotate(CentTheta, G4ThreeVector(0,1,0)); 
  const G4double VD0SizeX = 150.*CLHEP::cm;
  const G4double VD0SizeY = 150.*CLHEP::cm;
  const G4double VD0SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD0 = new  G4RotationMatrix();
  VirtualDetector *VD0 =
    new VirtualDetector( "VD0", posVD0, rotVD0, 
			 VD0SizeX, VD0SizeY, VD0SizeZ, 0 );
  VD0->SetMaterials( mList_->Vacuum );
  VD0->Place( world );
#endif
	
#if 1 //Q1
  //=== Virtual Detector for Q1 Entrance ===//
  //  const G4ThreeVector preposVD3( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD3 - 0.05*CLHEP::cm ); //default
  const G4ThreeVector preposVD3( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD3 - 0.05*CLHEP::cm ); //temp
  const G4ThreeVector posVD3 = HRSTrans->TransformPoint(preposVD3);
  const G4double VD3SizeX = 60.*CLHEP::cm;
  const G4double VD3SizeY = 60.*CLHEP::cm;
  const G4double VD3SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD3 = new  G4RotationMatrix();
  rotVD3->rotateY(-E_Theta);
  VirtualDetector *VD3 =
    new VirtualDetector( "VD3", posVD3, rotVD3, 
			 VD3SizeX, VD3SizeY, VD3SizeZ, 3 );
  VD3->SetMaterials( mList_->Vacuum );
  VD3->Place( world );
  
  //=== Q1  Collimator ======= 
  G4RotationMatrix *rotcoll1 = new  G4RotationMatrix();
  rotcoll1->rotateZ(-90*degree);
  Collimator *Coll1 = new Collimator( "Collimator" , posVD3, rotcoll1);
  Coll1->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  Coll1->SetRadius(15.*CLHEP::cm);
  //  //Coll1->Place( world );
  

  //=== Quadrupole 1 ===//
  std::string Q1map = paramMan->GetQ1FieldMap();
  const G4ThreeVector preposQ1( 0*CLHEP::cm, 0*CLHEP::cm, VTtoQ1 );
  const G4ThreeVector posQ1 = HRSTrans->TransformPoint(preposQ1);
  G4RotationMatrix *rotQ1 = new G4RotationMatrix();
  rotQ1->rotateY(-E_Theta);

  Q1Magnet *Q1 = 
    new Q1Magnet( "Q1 Magnet", posQ1, rotQ1, Q1map );
  Q1->SetMaterials( mList_->Fe, mList_->Vacuum );
  //Q1->Place( world ); //default tube
  Q1->GDMLPlace( world ); // SOS_Q1
  EMField_->AddElement( Q1 );
  paramMan->SetpQ1(posQ1);

    //  DMagnet *Dtemp = new DMagnet( "Dtemp", posQ1, rotQ1, Q1map, true);
    //  Dtemp->SetMaterials( mList_->Fe, mList_->Vacuum );
    //  Dtemp->GDMLYokePlace( world );
    //  EMField_->AddElement( Dtemp );


  //=== Virtual Detector for Q1 Exit ===//
  //  const G4ThreeVector preposVD4( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD4 + 0.05*CLHEP::cm ); //default
  const G4ThreeVector preposVD4( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD4 + 0.05*CLHEP::cm ); //temp
  const G4ThreeVector posVD4 = HRSTrans->TransformPoint(preposVD4);
  const G4double VD4SizeX = 60.*cm;
  const G4double VD4SizeY = 60.*cm;
  const G4double VD4SizeZ = 0.1*cm;
  G4RotationMatrix *rotVD4 = new  G4RotationMatrix();
  rotVD4->rotateY(-E_Theta);
  VirtualDetector *VD4 =
    new VirtualDetector( "VD4", posVD4, rotVD4, 
			 VD4SizeX, VD4SizeY, VD4SizeZ, 4 );
  VD4->SetMaterials( mList_->Vacuum );
  VD4->Place( world );
#endif
  
#if 1 //Q2
  //=== Collimator for Q2 Entrance ===//
  const G4ThreeVector preposCollQ2i( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD5 - 2.*CLHEP::cm );
  const G4ThreeVector posCollQ2i = HRSTrans->TransformPoint(preposCollQ2i);
  G4RotationMatrix *rotCollQ2i = new  G4RotationMatrix();
  rotCollQ2i->rotateY(-E_Theta);
  Collimator *CollQ2i = new Collimator( "Collimator" , posCollQ2i, rotCollQ2i);
  CollQ2i->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  CollQ2i->SetRadius(60./2.*CLHEP::cm);
  CollQ2i->SetColBoxSize(150.*CLHEP::cm, 150.*CLHEP::cm, 1.*CLHEP::cm);
  //  CollQ2i->Place( world );
  
  //=== Virtual Detector for Q2 Entrance ===//
  const G4ThreeVector preposVD5( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD5 - 0.05*CLHEP::cm );
  const G4ThreeVector posVD5 = HRSTrans->TransformPoint(preposVD5);
  const G4double VD5SizeX = 150.*CLHEP::cm;
  const G4double VD5SizeY = 150.*CLHEP::cm;
  const G4double VD5SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD5 = new  G4RotationMatrix();
  rotVD5->rotateY(-E_Theta);
  VirtualDetector *VD5 =
    new VirtualDetector( "VD5", posVD5, rotVD5, 
			 VD5SizeX, VD5SizeY, VD5SizeZ, 5 );
  VD5->SetMaterials( mList_->Vacuum );
  VD5->Place( world );

  //=== Quadrupole 2 ===//
  std::string Q2map = paramMan->GetQ2FieldMap();
  const G4ThreeVector preposQ2( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoQ2 );
  const G4ThreeVector posQ2 = HRSTrans->TransformPoint(preposQ2);
  G4RotationMatrix *rotQ2 = new G4RotationMatrix();
  rotQ2->rotateY(-E_Theta);
  Q2Magnet *Q2 = 
    new Q2Magnet( "Q2 Magnet", posQ2, rotQ2, Q2map );
  Q2->SetMaterials( mList_->Fe, mList_->Vacuum );
  Q2->Place( world );
  EMField_->AddElement( Q2 );
  paramMan->SetpQ2(posQ2);

  //=== Virtual Detector for Q2 Exit ===//
  const G4ThreeVector preposVD6( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD6 + 0.05*CLHEP::cm );
  const G4ThreeVector posVD6 = HRSTrans->TransformPoint(preposVD6);
  const G4double VD6SizeX = 150.*cm;
  const G4double VD6SizeY = 150.*cm;
  const G4double VD6SizeZ = 0.1*cm;
  G4RotationMatrix *rotVD6 = new  G4RotationMatrix();
  rotVD6->rotateY(-E_Theta);
  VirtualDetector *VD6 =
    new VirtualDetector( "VD6", posVD6, rotVD6, 
			 VD6SizeX, VD6SizeY, VD6SizeZ, 6 );
  VD6->SetMaterials( mList_->Vacuum );
  VD6->Place( world );

  //=== Q2 Collimator ======= by fujita
  G4RotationMatrix *rotcoll2 = new  G4RotationMatrix();
  //rotcoll2->rotateZ(-90*degree);
  Collimator *Coll2 = new Collimator( "Collimator" , posVD5, rotcoll2);
  Coll2->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  Coll2->SetRadius(30.*CLHEP::cm);
  //Coll2->Place( world );
#endif

#if 1 //Dipole
  //=== Collimator for Dipole Septum Entrance ===// aida
  const G4ThreeVector preposCollDi(	0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD7 - 30.*CLHEP::cm );
  G4RotationMatrix *rotCollDi = new  G4RotationMatrix();
  rotCollDi->rotateY(-E_Theta);
  const G4ThreeVector posCollDi = HRSTrans->TransformPoint(preposCollDi);
  CollimatorBox *CollDi = new CollimatorBox( "D Collimator forward" , posCollDi, rotCollDi);
  CollDi->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  CollDi->SetColBoxSize(300.*CLHEP::cm, 300.*CLHEP::cm, 5.*CLHEP::cm);
  CollDi->SetHoleSize(40.*CLHEP::cm, 150.*CLHEP::cm, 5.1*CLHEP::cm);
  //  CollDi->Place( world );
  
  //=== Virtual Detector for  Dipole Entrance ===//
  //  const G4ThreeVector preposVD7( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD7 - 0.05*CLHEP::cm -offset);
  const G4ThreeVector preposVD7( 0.*CLHEP::cm, 0.*CLHEP::cm, VTtoVD7 - 25*CLHEP::cm);
  const G4ThreeVector posVD7 = HRSTrans->TransformPoint(preposVD7);
  const G4double VD7SizeX = 200.*CLHEP::cm;
  const G4double VD7SizeY = 200.*CLHEP::cm;
  const G4double VD7SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD7 = new  G4RotationMatrix();
  rotVD7->rotateY(-E_Theta);
  //  rotVD7->rotateX(-30*CLHEP::degree);
  VirtualDetector *VD7 =
    new VirtualDetector( "VD7", posVD7, rotVD7, 
			 VD7SizeX, VD7SizeY, VD7SizeZ, 7 );
  VD7->SetMaterials( mList_->Vacuum );
  VD7->Place( world );
  
  //=== Dipole Magnet ==//
  const G4double r = 840.*CLHEP::cm; //central radius
  const G4double theta = 45*CLHEP::degree; //
  const G4ThreeVector preposD( 0.*CLHEP::cm, r , VTtoD ); //fujita
  const G4ThreeVector posD = HRSTrans->TransformPoint(preposD);
  const G4double posDz = posD.z() ;
  G4RotationMatrix *rotDmag = new  G4RotationMatrix(); 
  rotDmag->rotateZ(90.*CLHEP::degree);
  rotDmag->rotateX(90.*CLHEP::degree);
  rotDmag->rotateX(E_Theta);
  rotDmag->rotateZ(22.5*CLHEP::degree); // aida
  std::string Dmap = paramMan->GetDFieldMap();
  DMagnet *DYoke = new DMagnet( "DYoke", posD, rotDmag, Dmap, true);
  //  DMagnet *DCoil = new DMagnet( "DCoil", posD, rotDmag , "param/test.table", false);
  //  DMagnet *DYoke = new DMagnet( "DYoke", G4ThreeVector(0.,0.,0), 0 ,Dmap, true);
  DYoke->SetMaterials( mList_->Fe, mList_->Vacuum );
  //  DCoil->SetMaterials( mList_->Fe, mList_->Vacuum );
  //  DYoke->Place( world );
  DYoke->GDMLYokePlace( world );
  //  DYoke->GDMLCoilPlace( world );
  EMField_->AddElement( DYoke );
  paramMan->SetpVD7(r,posDz);

  //=== Virtual Detector for Dipole Exit ===//
  const G4ThreeVector preposVD8( 0.*cm, VD8y + 50.*CLHEP::cm*sin(45.*CLHEP::degree) , VD8z + 50.*CLHEP::cm*cos(45.*CLHEP::degree) ); //rotate
  const G4ThreeVector posVD8 = HRSTrans->TransformPoint(preposVD8);
  const G4double VD8SizeX = 200.*CLHEP::cm;
  const G4double VD8SizeY = 200.*CLHEP::cm;
  const G4double VD8SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD8 = new  G4RotationMatrix();
  rotVD8->rotateY(-E_Theta);
  rotVD8->rotateX(45*CLHEP::degree);
  VirtualDetector *VD8 =
    new VirtualDetector( "VD8", posVD8, rotVD8, 
			 VD8SizeX, VD8SizeY, VD8SizeZ, 8 );
  VD8->SetMaterials( mList_->Vacuum );
  VD8->Place( world );

  //  //=== Dipole Collimator ======= by fujita
  //  G4RotationMatrix *rotColl3 = new  G4RotationMatrix();
  //  rotColl3->rotateX(-30*CLHEP::degree);
  //  Collimator1 *Coll3 = new Collimator1( "Collimator1" , posVD7, rotColl3);
  //  Coll3->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  //  	Coll3->SetRadius(25.*CLHEP::cm);
  //  //Coll3->Place( world );
#endif

#if 1 //Q3
  //=== Collimator for Q3 Entrance ===//
  const G4ThreeVector preposCollQ3i( 0.*CLHEP::cm, VD8y+VD8toVD9*sin(theta)-5.*CLHEP::cm,VD8z+VD8toVD9*cos(theta)-5.*CLHEP::cm );
  const G4ThreeVector posCollQ3i = HRSTrans->TransformPoint(preposCollQ3i);
  G4RotationMatrix *rotCollQ3i = new  G4RotationMatrix();
  rotCollQ3i->rotateY(-E_Theta);
  rotCollQ3i->rotateX(45.*CLHEP::degree);
  Collimator *CollQ3i = new Collimator( "Collimator" , posCollQ3i, rotCollQ3i);
  CollQ3i->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  CollQ3i->SetRadius(60./2.*CLHEP::cm);
  CollQ3i->SetColBoxSize(150.*CLHEP::cm, 150.*CLHEP::cm, 1.*CLHEP::cm);
  //  CollQ3i->Place( world );
  
  //=== Virtual Detector for Q3 Entrance ===//
  const G4ThreeVector preposVD9( 0.*CLHEP::cm, VD8y+VD8toVD9*sin(theta)-0.03*CLHEP::cm,VD8z+VD8toVD9*cos(theta)-0.03*CLHEP::cm ); //rotate fujita
  const G4ThreeVector posVD9 = HRSTrans->TransformPoint(preposVD9);
  paramMan->SetpVD9(preposVD9.y(),preposVD9.z());
  const G4double VD9SizeX = 150.*CLHEP::cm;
  const G4double VD9SizeY = 150.*CLHEP::cm;
  const G4double VD9SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD9 = new  G4RotationMatrix();
  rotVD9->rotateY(-E_Theta);
  rotVD9->rotateX(45.*CLHEP::degree);
  VirtualDetector *VD9 =
    new VirtualDetector( "VD9", posVD9, rotVD9, 
			 VD9SizeX, VD9SizeY, VD9SizeZ, 9 );
  VD9->SetMaterials( mList_->Vacuum );
  VD9->Place( world );
  
  //=== Quadrupole 3 ===//
  std::string Q3map = paramMan->GetQ2FieldMap();
  const G4ThreeVector preposQ3( 0, VD8y+VD8toQ3*sin(theta),VD8z+VD8toQ3*cos(theta) ); //rotate
  const G4ThreeVector posQ3 = HRSTrans->TransformPoint(preposQ3);
  G4RotationMatrix *rotQ3 = new G4RotationMatrix();
  rotQ3->rotateY(-E_Theta);
  rotQ3->rotateX(45*CLHEP::degree); 
  Q3Magnet *Q3 = 
    new Q3Magnet( "Q3 Magnet", posQ3, rotQ3, Q3map );
  Q3->SetMaterials( mList_->Fe, mList_->Vacuum );
  Q3->Place( world );
  EMField_->AddElement( Q3 );
  
  //=== Virtual Detector for Q3 Exit ===//
  const G4ThreeVector preposVD10( 0.*CLHEP::cm, VD8y+VD8toVD10*sin(theta)+0.03*CLHEP::cm,
				  VD8z+VD8toVD10*cos(theta)+0.03*CLHEP::cm ); //rotate fujita
  const G4ThreeVector posVD10 = HRSTrans->TransformPoint(preposVD10);
  const G4double VD10SizeX = 150.*CLHEP::cm;
  const G4double VD10SizeY = 150.*CLHEP::cm;
  const G4double VD10SizeZ = 0.1*CLHEP::cm;
  paramMan->SetpVD10(preposVD10.y(),preposVD10.z());
  G4RotationMatrix *rotVD10 = new  G4RotationMatrix();
  rotVD10->rotateY(-E_Theta);
  rotVD10->rotateX(45*CLHEP::degree); //rotate
  VirtualDetector *VD10 =
    new VirtualDetector( "VD10", posVD10, rotVD10, 
			 VD10SizeX, VD10SizeY, VD10SizeZ, 10 );
  VD10->SetMaterials( mList_->Vacuum );
  VD10->Place( world );

  //=== Q3 Collimator ======= by fujita
  G4RotationMatrix *rotColl4 = new  G4RotationMatrix();
  rotColl4->rotateX(45*CLHEP::degree);
  Collimator *Coll4 = new Collimator( "Collimator" , posVD9, rotColl4);
  Coll4->SetMaterials( mList_->Heavymet, mList_->Vacuum );
  Coll4->SetRadius(30.*CLHEP::cm);
  //Coll4->Place( world );
#endif 

#if 1 //VD11 FP
  //=== Virtual Detector for RP ===//
  const G4double preposVD11y = VD8y + VD8toVD11*sin(theta);
  const G4double preposVD11z = VD8z + VD8toVD11*cos(theta);
  const G4ThreeVector preposVD11( 0 , preposVD11y, preposVD11z); //rotate
  const G4ThreeVector posVD11 = HRSTrans->TransformPoint(preposVD11);
  const G4double VD11SizeX = 500.*CLHEP::cm;
  const G4double VD11SizeY = 500.*CLHEP::cm;
  const G4double VD11SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD11 = new  G4RotationMatrix();
  rotVD11->rotateY(-E_Theta);
  rotVD11->rotateX(90.*CLHEP::degree); //org rotate
  VirtualDetector *VD11 =
    new VirtualDetector( "VD11", posVD11, rotVD11, 
			 VD11SizeX, VD11SizeY, VD11SizeZ, 11 );
  VD11->SetMaterials( mList_->Vacuum );
  VD11->Place( world );
#endif

#if 1 //VD12 (VDC1) // aida 20180105
  //=== Virtual Detector for VDC1 ===//
  const G4double preposVD12y = VD8y + VD8toVD12*sin(theta);
  const G4double preposVD12z = VD8z + VD8toVD12*cos(theta);
  const G4ThreeVector preposVD12( 0, preposVD12y, preposVD12z); //rotate
  const G4ThreeVector posVD12 = HRSTrans->TransformPoint(preposVD12);
  //  const G4double VD12SizeX = 28.8*CLHEP::cm; // ref.) J.Alcorn et al., NIM 522 (2004)
  //  const G4double VD12SizeY = 211.8*CLHEP::cm;
  const G4double VD12SizeX = 500*CLHEP::cm;
  const G4double VD12SizeY = 500*CLHEP::cm;
  const G4double VD12SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD12 = new  G4RotationMatrix();
  rotVD12->rotateY(-E_Theta);
  rotVD12->rotateX(90*CLHEP::degree); // org. rotate
  VirtualDetector *VD12 =
    new VirtualDetector( "VD12", posVD12, rotVD12, 
			 VD12SizeX, VD12SizeY, VD12SizeZ, 12 );
  VD12->SetMaterials( mList_->Vacuum );
  VD12->Place( world );
#endif

#if 1 //VD13 (VDC2) // aida 20180105
  //=== Virtual Detector for VDC2 ===//
  const G4double preposVD13y = VD8y + VD8toVD13*sin(theta);
  const G4double preposVD13z = VD8z + VD8toVD13*cos(theta);
  const G4ThreeVector preposVD13( 0, preposVD13y, preposVD13z); //rotate
  const G4ThreeVector posVD13 = HRSTrans->TransformPoint(preposVD13);
  const G4double VD13SizeX = 28.8*CLHEP::cm; // ref.) J.Alcorn et al., NIM 522 (2004)
  const G4double VD13SizeY = 211.8*CLHEP::cm;
  const G4double VD13SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD13 = new  G4RotationMatrix();
  rotVD13->rotateY(-E_Theta);
  rotVD13->rotateX(90*CLHEP::degree); // org. rotate
  VirtualDetector *VD13 =
    new VirtualDetector( "VD13", posVD13, rotVD13, 
			 VD13SizeX, VD13SizeY, VD13SizeZ, 13 );
  VD13->SetMaterials( mList_->Vacuum );
  VD13->Place( world );
#endif

#if 1
  //VD15 (S2) // suzuki 20191012
  // === Virtual Detector for S2 === //
  G4double preposVD15y = VD8y + VD8toS2*sin(theta);
  G4double preposVD15z = VD8z + VD8toS2*cos(theta);
  const G4double VD15shift = 12.1*CLHEP::cm; // db_L.s2.dat
  preposVD15y += -VD15shift*sin(theta);
  preposVD15z += VD15shift*cos(theta);
  const G4ThreeVector preposVD15(0, preposVD15y, preposVD15z); //rotate
  const G4ThreeVector posVD15 = HRSTrans->TransformPoint(preposVD15);
  const G4double VD15SizeX = 43.2*CLHEP::cm;
  const G4double VD15SizeY = 223.6*CLHEP::cm;
  const G4double VD15SizeZ = 0.1*CLHEP::cm;
  G4RotationMatrix *rotVD15 = new G4RotationMatrix();
  rotVD15->rotateY(-E_Theta);
  rotVD15->rotateX(45*CLHEP::degree); 
  VirtualDetector *VD15 = 
    new VirtualDetector( "VD15", posVD15, rotVD15,
			 VD15SizeX, VD15SizeY, VD15SizeZ, 15);
  VD15->SetMaterials( mList_->Vacuum );
  VD15->Place( world );

#endif


  /////////////////////////
  // Sensitive Detectors //
  /////////////////////////
  
  //=== Virtual Detector ===//
  HRSVDSD *VDSD = new HRSVDSD( "/HRS/VD" );
  SDManager->AddNewDetector( VDSD );
#if 1
  //  VD0->GetDetectorLV()->SetSensitiveDetector( VDSD );
  //  VD1->GetDetectorLV()->SetSensitiveDetector( VDSD );
  //  VD2->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD3->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD4->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD5->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD6->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD7->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD8->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD9->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD10->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD11->GetDetectorLV()->SetSensitiveDetector( VDSD ); // reference plane
  VD12->GetDetectorLV()->SetSensitiveDetector( VDSD ); // VDC1
  VD13->GetDetectorLV()->SetSensitiveDetector( VDSD ); // VDC2
  //  VD14->GetDetectorLV()->SetSensitiveDetector( VDSD );
  VD15->GetDetectorLV()->SetSensitiveDetector( VDSD ); // S2
  //  VDSD->SetElement( 0, VD0 );
  //  VDSD->SetElement( 1, VD1 );
  //  VDSD->SetElement( 2, VD2 );
  VDSD->SetElement( 3, VD3 );
  VDSD->SetElement( 4, VD4 );
  VDSD->SetElement( 5, VD5 );
  VDSD->SetElement( 6, VD6 );
  VDSD->SetElement( 7, VD7 );
  VDSD->SetElement( 8, VD8 );
  VDSD->SetElement( 9, VD9 );
  VDSD->SetElement( 10, VD10 );
  VDSD->SetElement( 11, VD11 );
  VDSD->SetElement( 12, VD12 );
  VDSD->SetElement( 13, VD13 );
  //  VDSD->SetElement( 14, VD14 );
  VDSD->SetElement( 15, VD15 );
#endif
  
  ///////////////////
  // Visualisation //
  ///////////////////
  
  //=== World ===//
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);

  //=== Splitter Yoke ===//
  
  //=== Q1 Yoke ===//
  
  //=== Q2 Yoke ===//
  
  //  //=== Target ===//
  //    G4Colour colourTarget(0., 1., 1.); // cyan
  //    G4VisAttributes *targetVisAtt = new G4VisAttributes(true, colourTarget);
  //    targetLV->SetVisAttributes(targetVisAtt);
  
  
  ///////////////////////
  // Special User Cuts //
  ///////////////////////
  //  G4UserLimits  *userLimitsYoke =
  //    new G4UserLimits( 5.*CLHEP::cm, 30.*CLHEP::cm, 100.*ns, 10.*MeV, 5.*CLHEP::cm);
  //  splitter->SetUserLimits( userLimitsYoke );
  //  Q1->SetUserLimits( userLimitsYoke );
  //  Q2->SetUserLimits( userLimitsYoke );
  //  Dipole->SetUserLimits( userLimitsYoke );
  //  CExtension->SetUserLimits( userLimitsYoke );
  //  targetC->SetUserLimits( userLimitsYoke );
  
  return world;
}

