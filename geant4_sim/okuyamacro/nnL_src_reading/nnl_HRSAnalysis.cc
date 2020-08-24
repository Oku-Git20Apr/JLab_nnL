/*
  HRSAnalysis.cc

  D.Kawama
*/

#include "HRSAnalysis.hh"
#include "HRSTargetHit.hh"
#include "HRSVDHit.hh"
#include "EDC1Hit.hh"
#include "EDC2Hit.hh"
#include "HRSEH1Hit.hh"
#include "HRSEH2Hit.hh"

#include "HRSTrajectory.hh"
#include "HRSSteppingAction.hh"
#include "RootHelper.hh"

#include "DMagnet.hh"

#include "HRSParamMan.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4String.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

#include "Randomize.hh"
#include <CLHEP/Random/RandGaussQ.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>

//#ifdef USE_CSTREAM
#include <cstdio>
//#endif
G4double PI = 4*atan(1);

G4double Mom;
G4double Theta;
G4double Phi;
G4double Xpt,Ypt;
G4double Xt,Yt,Zt;

G4double X[16],Y[16],Z[16];
G4double Xg[16], Yg[16], Zg[16]; /////suzuki
G4double Xp[16],Yp[16];
G4double Px[16],Py[16];
G4double Xv[16],Yv[16],Zv[16];
G4double P[16];
G4double L[16];
G4int trackID[16];
G4int charge[16];
G4int particleID[16];

G4double EDC1X[10],EDC1Y[10],EDC1Z[10];
G4double EDC1Xp[10],EDC1Yp[10];
G4double EDC1P[10] ;
G4double EDCFPX,EDCFPY,EDCFPZ;
G4double EDCFPXp,EDCFPYp;
G4double EDCFPP,EDCFPL ;
G4double EDC2X[6],EDC2Y[6],EDC2Z[6];
G4double EDC2Xp[6],EDC2Yp[6];
G4double EDC2P[6];

G4double EH1X[29],EH1Y[29],EH1Z[29];
G4double EH1Xp[29],EH1Yp[29];
G4double EH1P[29],EH1A[29],EH1T[29];

G4double EH2X[29],EH2Y[29],EH2Z[29];
G4double EH2Xp[29],EH2Yp[29];
G4double EH2P[29],EH2A[29],EH2T[29];

G4double XG,YG,ZG;

G4bool VDTrig=false;
G4bool VtxTrig=false;
G4bool SPLTrig=false;
G4bool PBTrig=false;
G4bool EDCTrig=false;
G4bool EDC1Trig=false;
G4bool EDC2Trig=false;
G4bool EHTrig=false;
G4bool EH1Trig=false;
G4bool EH2Trig=false;
G4bool EDetTrig=false;
G4double pBeam,pVD;

///////////////////////////////////////
HRSAnalysis::HRSAnalysis()
  : fTriggered(false), fActive_(true), fOutput(false)
    ///////////////////////////////////////
{
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  filename_ = paramMan->GetROOTFileName();
  defineHistograms();
  start = time(NULL);
  time(&start);
}

///////////////////////////
HRSAnalysis::~HRSAnalysis()
///////////////////////////
{
  SaveFile();
}


//////////////////////////////////////////////////
void HRSAnalysis::BeginOfRun( const G4Run *)
//////////////////////////////////////////////////
{
  trigNum=0;
  
}

////////////////////////////////////////////////
void HRSAnalysis::EndOfRun( const G4Run *aRun )
////////////////////////////////////////////////
{
  G4cout << G4endl;
  G4cout << " Events generated = " << aRun->GetNumberOfEvent() << G4endl;
  G4cout << " Events triggered = " << trigNum << G4endl;
  G4cout << G4endl;
}

/////////////////////////////////////////////////////////
void HRSAnalysis::BeginOfEvent( const G4Event *)
/////////////////////////////////////////////////////////
{
  //G4cout << "Begin Of Event" <<G4endl;
}

////////////////////////////////////////////////////////////////
void HRSAnalysis::PrimaryGeneration( const G4ThreeVector &pos,
				     const G4ThreeVector &mom,
				     const G4int ID,
				     const G4double z)
////////////////////////////////////////////////////////////////
{
  //G4cout << "Primary Generation" <<G4endl;
  gPos_ = pos; 
  gMom_ = mom;
  evID = ID;
  //	G4cout<< "evID = " << evID <<G4endl;
  zraster = z;
}
////////////////////////////////////////////////////////////////
void HRSAnalysis::PrimaryGeneration( const G4ThreeVector &pos,
				     const G4ThreeVector &mom,
				     const G4int ID,
				     const G4double z,
				     const G4double thetap, G4double phip)
////////////////////////////////////////////////////////////////
{
  //G4cout << "Primary Generation" <<G4endl;
  gPos_ = pos; 
  gMom_ = mom;
  evID = ID;
  //	G4cout<< "evID = " << evID <<G4endl;
  zraster = z;
  thetap_ = thetap;
  phip_ = phip;
}

////////////////////////////////////////////////////////
void HRSAnalysis::EndOfEvent( const G4Event *anEvent  )
////////////////////////////////////////////////////////
{
  char anatime[100];
  G4int eventID = anEvent->GetEventID();
  if( eventID%1000 == 0 ){
    end = time(NULL);
    time(&end);
    sprintf( anatime,"%.0f sec",difftime(end,start) );
    G4cout << "(SEP-HRS) Event#, Trig# : " << std::setw(8) << eventID << ", "
	   << std::setw(8) << trigNum 
	   << "  ( " << anatime << " )" <<G4endl;
  }
  fTriggered = false;
  NormalAnalysis( anEvent );
  //G4cout << "End Of Event" <<G4endl;
}

///////////////////////////////////////////////////////////
void HRSAnalysis::NormalAnalysis( const G4Event *anEvent )
///////////////////////////////////////////////////////////
{
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  TTree *tree1 = dynamic_cast<TTree *>(gFile->Get("tree1"));
  InitializeEvent();
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  int flag = paramMan->GetEventFlag();
  int EDCFlag = paramMan->GetEDCFlag();
  int FillFlag = paramMan->GetFillFlag();
  G4double E_Theta = paramMan->GetRotAngle()*CLHEP::degree;
  G4double FPAng = paramMan->GetFPAngle()*CLHEP::degree;
  pBeam =  paramMan->GetpBeam0();
  pVD   =  paramMan->GetpBeam1();
  
  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  //  G4int colIdTarget = SDMan->GetCollectionID( "TargetCollection" );
  G4int colIdVD[16] = {0};
  G4double pt  = -1000., pxt = -1000.;
  G4double pyt = -1000., pzt = -1000.;
  G4double theta = -1000., phi   = -1000.;
  G4double xpt = -1000., ypt = -1000.;
  G4double xt = -1000., yt = -1000., zt = -1000.;
  //	G4double z0 = -1000.0;
  G4double VDEdep = 0.;
  G4double DCEdep = 0.;
  G4double EH1Edep = 0.;
  G4double EH2Edep = 0.;
  G4int nHits = 0;
  G4double TargetEdep = 0.0;
  G4int VDhit[16] = {0};
  HRSVDHitsCollection *VDHC[16];

  //  //=== Target ===//
  //  HRSTargetHitsCollection *TargetHC; 
  //  TargetHC = dynamic_cast<HRSTargetHitsCollection *>(HCE->GetHC(colIdTarget));
  //  if( TargetHC ){
  //    nHits = TargetHC->entries();
  //    for( G4int i=0; i<nHits; i++ ){
  //      HRSTargetHit *aHit = (*TargetHC)[i];
  //      TargetEdep += aHit->GetEdep();
  //    }
  //  } // if( HRSTargetHC )
  
  if( flag==0 ){
    
    //=== Virtual Detector ===//
    for( int i=0 ; i<=15; ++i){ //fujita
      //		for( int i=1 ; i<=16; ++i){
      std::ostringstream id;
      id << i;
      G4String Name = G4String("VD")+id.str().c_str()+"-Collection";
      colIdVD[i] = SDMan-> GetCollectionID( Name );
      VDHC[i] = dynamic_cast<HRSVDHitsCollection *>(HCE->GetHC(colIdVD[i]));
      if( VDHC[i] ){
	G4int nHits = VDHC[i]->entries();
	G4ThreeVector lPos,gPos;
	G4ThreeVector MField;
	G4ThreeVector lMom;
	G4ThreeVector Vertex;
	//	G4cout << "i = " << i << G4endl; 
	//	G4cout << "nHits = " <<nHits << G4endl; 
	for( G4int j=0; j<nHits; j++ ){
	  HRSVDHit *aHit = (*VDHC[i])[j];
	  VDhit[i] = nHits;
	  Vertex = aHit->GetVertex();
	  Yv[i] = Vertex.y()/CLHEP::cm;
	  Xv[i] = (Vertex.x()*cos(E_Theta) + Vertex.z()*sin(E_Theta))/CLHEP::cm;
	  Zv[i] = (Vertex.z()*cos(E_Theta) - Vertex.x()*sin(E_Theta))/CLHEP::cm;
	  L[i] = aHit->GettLength()/CLHEP::cm;
	  VDEdep = aHit->GetEdep();
	  //For Global Coordinates
	  /*lPos = aHit->GetGPos();
	    X[i] = (lPos.x()*cos(E_Theta)
	    +lPos.z()*sin(E_Theta))/cm;
	    Z[i] = (lPos.z()*cos(E_Theta)
	    -lPos.x()*sin(E_Theta))/cm;*/
	  gPos = aHit->GetGPos();  //////suzuki from
	  Xg[i] = gPos.x()/CLHEP::cm; 
	  Yg[i] = gPos.y()/CLHEP::cm;
	  Zg[i] = gPos.z()/CLHEP::cm;
	  /////// suzuki to
	  //For Local Coordinates
	  lPos = aHit->GetPos();
	  X[i] = lPos.x()/CLHEP::cm;
	  Y[i] = lPos.y()/CLHEP::cm;
	  Z[i] = lPos.z()/CLHEP::cm;
	  //					gPos = aHit->GetGPos();
	  //					X[i] = gPos.x()/cm;
	  //					Y[i] = gPos.y()/cm;
	  //					Z[i] = gPos.z()/cm;
	  
	  //Don't comment out!!
	  
	  lMom = aHit->GetMom();
	  G4double pxVD = lMom.x();
	  G4double pyVD = lMom.y();
	  G4double pzVD = lMom.z();
	  P[i] = sqrt(pxVD*pxVD + pyVD*pyVD + pzVD*pzVD)/CLHEP::GeV;
	  /////	  G4cout << "P[" << i << "] = "<< P[i] << G4endl; /////suzuki
	  /*G4double pxVD = lMom.x()/GeV;
	    G4double pyVD = lMom.y()/GeV;
	    G4double pzVD = lMom.z()/GeV;
	    P[i] = sqrt(pxVD*pxVD+pyVD*pyVD+pzVD*pzVD);*/
	  Xp[i] = pxVD/pzVD;
	  Yp[i] = pyVD/pzVD;
	  Px[i] = pxVD/CLHEP::GeV; // added by fujita
	  Py[i] = pyVD/CLHEP::GeV;

	  trackID[i] = aHit->GetTrackID(); // trackID
	  charge[i] = aHit->GetCharge(); // charge
	  G4String pname = aHit->GetPname();
	  if( pname == "e-" || pname == "e+" || pname == "gamma" ){
		particleID[i] = 0;
	  }
	  else if( pname == "proton"){
		particleID[i] = 1;
	  }
	  else if( pname == "kaon+" || pname == "kaon-" || pname == "kaon0" ){
		particleID[i] = 2;
	  }
	  else if( pname == "pi+" || pname == "pi-" || pname == "pi0" ){
		particleID[i] = 3;
	  }
	  else{
		particleID[i] = -1;
	  }
	  
	  //					if ( i == 11 ) {
	  //						G4cout << "==========" << i << "========="<< G4endl;
	  //						G4cout << "x y z = " << X[i] << " " << Y[i] << " " << Z[i] << G4endl; //fujita
	  //						G4cout << "Xp Yp P = " << Xp[i] << " " << Yp[i] << " " << P[i] << G4endl; //fujita
	  //						G4cout << "P Px Py = " << P[i] << " " << Px[i] << " " << Py[i] << G4endl; //fujita
	  //					}
	  //					if ( i == 13 ) {
	  //	 				  G4cout << "x y z = " << X[i] << " " << Y[i] << " " << Z[i] << G4endl; //fujita
	  //	 				  G4cout << "Xp Yp P = " << Xp[i] << " " << Yp[i] << " " << P[i] << G4endl; //fujita
	  //	 				  getchar();
	  //	 				}
	  
	} // for i
      } // if( VDHC[i] )
    } // for j
    //   X[9]*=cos(FPAng);
    //   Xp[9]+=FPAng/CLHEP::degree*3.14/180;
    //
    
#if 0
    //===Drift Chamber===//
    
    if (EDCFlag!=0){
      //=== EDC1 ===//
      G4int colIdEDC1 = SDMan->GetCollectionID( "EDC1Collection" );
      EDC1HitsCollection *EDC1HC;
      EDC1HC = dynamic_cast<EDC1HitsCollection *>
	(HCE->GetHC(colIdEDC1));
      if(EDC1HC){
	G4int nHits = EDC1HC->entries();
	G4ThreeVector lPos;
	G4ThreeVector lMom;
	G4ThreeVector Vertex;
	G4int LID;
	EDC1Trig=true;
	for( G4int j=0; j<nHits; j++ ){
	  EDC1Hit *aHit = (*EDC1HC)[j];
	  DCEdep = aHit->GetEdep();
	  LID = aHit->GetLayerID();
	  if (LID==10){
	    EDCFPL = aHit->GettLength()/CLHEP::cm;
	    //For Local Coordinates
	    lPos = aHit->GetLPos();
	    EDCFPX = lPos.x()/CLHEP::cm;
	    EDCFPZ = lPos.z()/CLHEP::cm;
	    EDCFPY = lPos.y()/CLHEP::cm;
	    lMom = aHit->GetLMom();
	    G4double pxDC = lMom.x()/CLHEP::GeV;
	    G4double pyDC = lMom.y()/CLHEP::GeV;
	    G4double pzDC = lMom.z()/CLHEP::GeV;
	    EDCFPP = sqrt(pxDC*pxDC+pyDC*pyDC+pzDC*pzDC);
	    //EDCFPXp = atan2(pxDC, pzDC)/radian;
	    //EDCFPYp = atan2(pyDC, pzDC)/radian;
	    EDCFPXp = pxDC/pzDC;
	    EDCFPYp = pyDC/pzDC;
	  }
	  else {
	    //For Local Coordinates
	    lPos = aHit->GetLPos();
	    EDC1X[LID] = lPos.x()/CLHEP::cm;
	    EDC1Z[LID] = lPos.z()/CLHEP::cm;
	    EDC1Y[LID] = lPos.y()/CLHEP::cm;
	    lMom = aHit->GetLMom();
	    G4double pxDC = lMom.x()/CLHEP::GeV;
	    G4double pyDC = lMom.y()/CLHEP::GeV;
	    G4double pzDC = lMom.z()/CLHEP::GeV;
	    EDC1P[LID] = sqrt(pxDC*pxDC+pyDC*pyDC+pzDC*pzDC);
	    //EDC1Xp[LID] = atan2(pxDC, pzDC)/radian;
	    //EDC1Yp[LID] = atan2(pyDC, pzDC)/radian;
	    EDC1Xp[LID] = pxDC/pzDC;
	    EDC1Yp[LID] = pyDC/pzDC;
	  }
	}
      }
      //=== EDC2 ===//
      G4int colIdEDC2 = SDMan->GetCollectionID( "EDC2Collection" );
      EDC2HitsCollection *EDC2HC;
      EDC2HC = dynamic_cast<EDC2HitsCollection *>
	(HCE->GetHC(colIdEDC2));
      if(EDC2HC){
	G4int nHits = EDC2HC->entries();
	G4ThreeVector lPos;
	G4ThreeVector lMom;
	G4ThreeVector Vertex;
	G4int LID;
	EDC2Trig=true;
	for( G4int j=0; j<nHits; j++ ){
	  EDC2Hit *aHit = (*EDC2HC)[j];
	  DCEdep = aHit->GetEdep();
	  LID = aHit->GetLayerID();
	  //For Local Coordinates
	  lPos = aHit->GetLPos();
	  EDC2X[LID] = lPos.x()/CLHEP::cm;
	  EDC2Z[LID] = lPos.z()/CLHEP::cm;
	  EDC2Y[LID] = lPos.y()/CLHEP::cm;
	  lMom = aHit->GetLMom();
	  G4double pxDC = lMom.x()/CLHEP::GeV;
	  G4double pyDC = lMom.y()/CLHEP::GeV;
	  G4double pzDC = lMom.z()/CLHEP::GeV;
	  EDC2P[LID] = sqrt(pxDC*pxDC+pyDC*pyDC+pzDC*pzDC);
	  //EDC2Xp[LID] = atan2(pxDC, pzDC)/radian;
	  //EDC2Yp[LID] = atan2(pyDC, pzDC)/radian;
	  EDC2Xp[LID] = pxDC/pzDC;
	  EDC2Yp[LID] = pyDC/pzDC;
	}
      }
    }
#endif
#if 0
    //=== HRS HodoScope ===//
    //=== EH1 ===//
    G4int colIdEH1 = SDMan->GetCollectionID( "EH1Collection" );
    EH1HitsCollection *EH1HC;
    EH1HC = dynamic_cast<EH1HitsCollection *>
      (HCE->GetHC(colIdEH1));
    if(EH1HC){
      G4int nHits = EH1HC->entries();
      G4ThreeVector lPos;
      G4ThreeVector lMom;
      G4ThreeVector Vertex;
      G4int ELID;
      EH1Trig=true;
      if (nHits>0){
	for( G4int j=0; j<nHits; j++ ){
	  EH1Hit *aHit = (*EH1HC)[j];
	  ELID = aHit->GetLayerID();
	  EH1Edep = aHit->GetEdep();
	  EH1A[ELID] += EH1Edep;
	  G4double EHTime = aHit->GetTime()/CLHEP::ns;
	  EH1T[ELID] = EHTime;
	  //For Local Coordinates
	  lPos = aHit->GetLPos();
	  EH1X[ELID] = lPos.x()/CLHEP::cm;
	  EH1Z[ELID] = lPos.z()/CLHEP::cm;
	  EH1Y[ELID] = lPos.y()/CLHEP::cm;
	  lMom = aHit->GetLMom();
	  G4double pxDC = lMom.x()/CLHEP::GeV;
	  G4double pyDC = lMom.y()/CLHEP::GeV;
	  G4double pzDC = lMom.z()/CLHEP::GeV;
	  EH1P[ELID] = sqrt(pxDC*pxDC+pyDC*pyDC+pzDC*pzDC);
	  //EH1Xp[ELID] = atan2(pxDC, pzDC)/radian;
	  //EH1Yp[ELID] = atan2(pyDC, pzDC)/radian;
	  EH1Xp[ELID] = pxDC/pzDC;
	  EH1Yp[ELID] = pyDC/pzDC;
	}
      }
      else{
	for( G4int j=0; j<nHits; j++){
	  EH1A[j]=-1000.;
	}
      }
    }
    //=== EH2 ===//
    G4int colIdEH2 = SDMan->GetCollectionID( "EH2Collection" );
    EH2HitsCollection *EH2HC;
    EH2HC = dynamic_cast<EH2HitsCollection *>
      (HCE->GetHC(colIdEH2));
    if(EH2HC){
      G4int nHits = EH2HC->entries();
      G4ThreeVector lPos;
      G4ThreeVector lMom;
      G4ThreeVector Vertex;
      G4int ELID;
      EH2Trig=true;
      if (nHits>0){
	for( G4int j=0; j<nHits; j++ ){
	  EH2Hit *aHit = (*EH2HC)[j];
	  ELID = aHit->GetLayerID();
	  EH2Edep = aHit->GetEdep();
	  EH2A[ELID] += EH2Edep;
	  G4double EHTime = aHit->GetTime()/CLHEP::ns;
	  EH2T[ELID] = EHTime;
	  //For Local Coordinates
	  lPos = aHit->GetLPos();
	  EH2X[ELID] = lPos.x()/CLHEP::cm;
	  EH2Z[ELID] = lPos.z()/CLHEP::cm;
	  EH2Y[ELID] = lPos.y()/CLHEP::cm;
	  lMom = aHit->GetLMom();
	  G4double pxDC = lMom.x()/CLHEP::GeV;
	  G4double pyDC = lMom.y()/CLHEP::GeV;
	  G4double pzDC = lMom.z()/CLHEP::GeV;
	  EH2P[ELID] = sqrt(pxDC*pxDC+pyDC*pyDC+pzDC*pzDC);
	  //EH2Xp[ELID] = atan2(pxDC, pzDC)/radian;
	  //EH2Yp[ELID] = atan2(pyDC, pzDC)/radian;
	  EH2Xp[ELID] = pxDC/pzDC;
	  EH2Yp[ELID] = pyDC/pzDC;
	}
      }
      else{
	for( G4int j=0; j<nHits; j++){
	  EH2A[j]=-1000.;
	}
      }
    }
#endif 
  
  
  } // if( flag==0 )
  
  if ( EDC1X[0]>-1000.
       && EDC1X[1]>-1000.
       && EDC1X[2]>-1000.
       && EDC1X[3]>-1000.
       && EDC1X[4]>-1000.
       && EDC1X[5]>-1000.
       && EDC1X[6]>-1000.
       && EDC1X[7]>-1000.
       && EDC1X[8]>-1000.
       && EDC1X[9]>-1000.)
    {EDC1Trig=true;}
  else{EDC1Trig=false;}
  if ( EDC2X[0]>-1000.
       && EDC2X[1]>-1000.
       && EDC2X[2]>-1000.
       && EDC2X[3]>-1000.
       && EDC2X[4]>-1000.
       && EDC2X[5]>-1000.)
    {EDC2Trig=true;}
  else{EDC2Trig=false;}
  
  if( //VDhit[1]>0 && VDhit[2]>0 && // Sep
      //  if( // w/o septum
     VDhit[3]>0 && VDhit[4]>0 && // Q1
     VDhit[5]>0 && VDhit[6]>0 && // Q2
     VDhit[7]>0 && VDhit[8]>0 && // D
     VDhit[9]>0 && VDhit[10]>0 && // Q3
     VDhit[11]>0 // RP
      ){ 
    VDTrig = true; 
    trigNum++; // aida 20180105 (This number is shown in a terminal)
  } 
  //			fabs(X[7])<53./2.){ VDTrig = true;}
  //  if( VDTrig && Zv[7]<1.){ VtxTrig = true; }
  //  if( VDTrig ){ VtxTrig = true; } // VtxTrig is not be used
  //  if( VtxTrig && Y[13]<17.){ SPLTrig = true;trigNum++; }
  //  if( SPLTrig && EDC1Trig && EDC2Trig){ EDCTrig = true;}
  if( VDhit[12]>0 && VDhit[13]>0){ EDCTrig = true;} // VDC1,2
  //  else {EDCTrig=false;}
  //  if( VDhit[12]>0 ){ PBTrig = true; }
  //  if( EH1Trig && EH2Trig ){ EHTrig = true; }
  //  else{ EHTrig = false; }
  //  fTriggered = SPLTrig;
  if( VDhit[15] > 0 ){ EDetTrig = true;} //S2 trig
  fTriggered = VDTrig;
  
  xt = gPos_.x();
  yt = gPos_.y();
  zt = gPos_.z();
  
  /*
    xt = gPos_.x()*cos(E_Theta)+gPos_.z()*sin(E_Theta);
    yt = gPos_.y();
    zt = -gPos_.x()*sin(E_Theta)+gPos_.z()*cos(E_Theta);
  */
  
  //  pt  = gMom_.mag();
  //  pxt = gMom_.x()*cos(E_Theta)+gMom_.z()*sin(E_Theta);
  //  pyt = gMom_.y();
  //  pzt = -gMom_.x()*sin(E_Theta)+gMom_.z()*cos(E_Theta);
  pt  = gMom_.mag();
  pxt = gMom_.x();
  pyt = gMom_.y();
  pzt = gMom_.z();
  theta = acos( pzt/pt );
  phi  = atan2( pyt, pxt );
  xpt = pxt/pzt;
  ypt = pyt/pzt;
  
  
#if 0 // 
  G4ThreeVector Gpos(xt-35.34*CLHEP::mm,yt,zt);
  //xt = 0;    
  //yt = 0;    
  //zt = gPos_.mag();    //test (for genfromfile)
  zt = Gpos.mag() - 80*CLHEP::cm;    //test (for genUni)
  Mom=pt/CLHEP::GeV;
  Theta=theta/CLHEP::degree;
  Phi=phi/CLHEP::degree;
  G4double CentTheta = paramMan->GetCentTheta();
  //G4cout << "###### centtheta " << CentTheta << " " << xpt<< G4endl; 
  //Xpt=xpt - CentTheta;    //add by fujita (for HRS + SEP)
  //Xpt=xpt + 0.0959931;    //add by fujita (for HRS + SEP)
  //Xpt = xpt + 0.12 - CentTheta ;  //for seed file
  Xpt=xpt;  //org
  Ypt=ypt;
#endif
#if 0 // test by fujita
  //  const double PI = 4.*atan(1.);
  G4double sTheta = paramMan->GetCentTheta();//scattered angle
  sTheta = -sTheta*180/3.1415926*CLHEP::degree;
  double zlength = 80*CLHEP::cm;
  double z0 = (gPos_.z()+ zlength*cos(sTheta) )/sin(77.5*CLHEP::degree);
  double testx = zlength*sin(sTheta) - z0*CLHEP::cm*cos(77.5*CLHEP::degree);
  double testz = -zlength*cos(sTheta) + z0*CLHEP::cm*sin(77.5*CLHEP::degree);
  //G4cout << "XXXXXXXX " << testx << " " <<  gPos_.x() << G4endl;
  //getchar();
  xt = gPos_.x() -testx -35.34;
  yt = gPos_.y() ;
  zt = z0 ;
  G4ThreeVector Gpos(xt,yt,zt);
  G4double CentTheta = paramMan->GetCentTheta();
  Xpt= xpt +  7.0*PI/180. - CentTheta;  //for seed file
  //Xpt=xpt;  //org
  Ypt=ypt;

#endif
  Mom = pt/CLHEP::GeV;
  Theta = theta/CLHEP::radian;
  Phi = phi/CLHEP::radian;
  Xpt = xpt/CLHEP::radian;
  Ypt = ypt/CLHEP::radian;
  Xt = xt/CLHEP::cm;
  Yt = yt/CLHEP::cm;
  Zt = zt/CLHEP::cm;
  //	G4cout<< "========== target ========="<<G4endl;
  //	G4cout<< "Xp Yp P = " << Xpt << " " << Ypt << " " << Mom <<G4endl;
  //	G4cout<< "P Px Py = " << Mom << " " << pxt/CLHEP::GeV << " " << pyt/CLHEP::GeV <<G4endl;
  //	G4cout<<G4endl;
  //	G4cout << "nnnnnnnnnnnnnn " << Xt << " " << Yt << " " << Zt << G4endl;  
  if(FillFlag == 0){
    //    G4cout << "Before Fill : ParticleName = " << P[3] << G4endl;  //////suzuki
    tree->Fill();
  }
  else{
    if(fTriggered){
      tree->Fill();
    }
  }
  
}  

//////////////////////////////////////////
void HRSAnalysis::InitializeEvent( void )
//////////////////////////////////////////
{
  Mom = -1000.0;
  Theta = -1000.0;Phi = -1000.0;
  Xpt = -1000.0;Ypt = -1000.0;
  VDTrig=false;
  VtxTrig=false;
  SPLTrig=false;
  EDetTrig=false;
  for (int i=0;i<16;i++){
    X[i]=-1000,Y[i]=-1000,Z[i]=-1000;
    Xg[i]=-1000,Yg[i]=-1000,Zg[i]=-1000; /////suzuki
    Xp[i]=-1000,Yp[i]=-1000;
    Xv[i]=-1000,Yv[i]=-1000,Zv[i]=-1000;
    P[i]=-1000;
    L[i]=-1000;
	trackID[i] = -10;
	charge[i] = -10;
	particleID[i] = -10;
  }
		
  EDCTrig=false;
  EDCFPX=-1000,EDCFPY=-1000,EDCFPZ=-1000;
  EDCFPXp=-1000,EDCFPYp=-1000;
  EDCFPP=-1000, EDCFPL=-1000;
  for (int i=0;i<10;i++){
    EDC1X[i]=-1000,EDC1Y[i]=-1000,EDC1Z[i]=-1000;
    EDC1Xp[i]=-1000,EDC1Yp[i]=-1000;
    EDC1P[i]=-1000;
  }
  for (int i=0;i<6;i++){
    EDC2X[i]=-1000,EDC2Y[i]=-1000,EDC2Z[i]=-1000;
    EDC2Xp[i]=-1000,EDC2Yp[i]=-1000;
    EDC2P[i]=-1000;
  }
  for (int i=0;i<29;i++){
    EH1X[i]=-1000,EH1Y[i]=-1000,EH1Z[i]=-1000;
    EH1Xp[i]=-1000,EH1Yp[i]=-1000;
    EH1P[i]=-1000,EH1A[i]=0.,EH1T[i]=-1000;
    EH2X[i]=-1000,EH2Y[i]=-1000,EH2Z[i]=-1000;
    EH2Xp[i]=-1000,EH2Yp[i]=-1000;
    EH2P[i]=-1000,EH2A[i]=0.,EH2T[i]=-1000;
  } 
}

////////////////////////////////////////
void HRSAnalysis::SaveFile( void ) const
////////////////////////////////////////
{
  if( fActive_ )
    gFile->Write();
}

////////////////////////////////////////
void HRSAnalysis::Terminate( void ) const
////////////////////////////////////////
{
  if( fActive_ )
    gFile->Write();
  gFile->Close();
}

//////////////////////////////////////////////
void HRSAnalysis::outputData( const G4Event *)
//////////////////////////////////////////////
{
  if( fDummy ){
    fprintf( fDummy, "Kine: %5d %4.3f %4.3f %4.3f %7.3f %7.3f %7.3f \n",
	     evID, gPos_.x()/CLHEP::cm, gPos_.y()/CLHEP::cm, gPos_.z(),
	     gMom_.x()/CLHEP::MeV, gMom_.y()/CLHEP::MeV, gMom_.z()/CLHEP::MeV);
  }

}

///////////////////////////////////////////
void HRSAnalysis::defineHistograms( void )
///////////////////////////////////////////
{
  new TFile( filename_, "recreate","",0 );
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  int flag = paramMan->GetEventFlag();
  char sPEDC[100];
  char sXEDC[100];
  char sYEDC[100];
  char sZEDC[100];
  char sXpEDC[100];
  char sYpEDC[100];
  char sPEDCD[100];
  char sXEDCD[100];
  char sYEDCD[100];
  char sZEDCD[100];
  char sXpEDCD[100];
  char sYpEDCD[100];
  
  char sPEH[100];
  char sXEH[100];
  char sYEH[100];
  char sZEH[100];
  char sXpEH[100];
  char sYpEH[100];
  char sAEH[100];
  char sTEH[100];
  char sPEHD[100];
  char sXEHD[100];
  char sYEHD[100];
  char sZEHD[100];
  char sXpEHD[100];
  char sYpEHD[100];
  char sAEHD[100];
  char sTEHD[100];
  
  if( flag==0 ){ 
    HBTree("tree","tree of HRS");
    TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
    //    tree->Branch("evID",      &evID,   "evID/I");
    tree->Branch("EevID",      &evID,   "EevID/I"); // aida, 20171114
    tree->Branch("pBeam",     &pBeam, "pBeam/D");
    tree->Branch("pVD",       &pVD,    "pVD/D");
    tree->Branch("EMom",      &Mom,   "EMom/D");
    tree->Branch("ETheta",  &Theta, "ETheta/D");
    tree->Branch("EPhi",      &Phi,   "EPhi/D");
    tree->Branch("EThetaPre",  &thetap_, "EThetaPre/D"); // zenith angle with respect to the central ray, aida 20171127
    tree->Branch("EPhiPre",      &phip_, "EPhiPre/D"  ); // azimuthal angle with respect to the central ray
    tree->Branch("EXpt",      &Xpt,   "EXpt/D");
    tree->Branch("EYpt",      &Ypt,   "EYpt/D");
    tree->Branch("EXt",      &Xt,   "EXt/D");
    tree->Branch("EYt",      &Yt,   "EYt/D");
    tree->Branch("EZt",      &Zt,   "EZt/D");
    
    tree->Branch("EPBT",     &P[0],  "EPBT/D");//fujita , 1Oct2015 Before target
    tree->Branch("EPxBT",   &Px[0], "EPxBT/D");
    tree->Branch("EPyBT",   &Py[0], "EPyBT/D");
    tree->Branch("EXBT",     &X[0],  "EXBT/D");
    tree->Branch("EYBT",     &Y[0],  "EYBT/D");
    tree->Branch("EZBT",     &Z[0],  "EZBT/D");
    tree->Branch("EXpBT",   &Xp[0], "EXpBT/D");
    tree->Branch("EYpBT",   &Yp[0], "EYpBT/D");
    //  tree->Branch("EPAT",     &P[13],  "EPAT/D");//fujita , 1Oct2015 After target
    //  tree->Branch("EXAT",     &X[13],  "EXAT/D");
    //  tree->Branch("EYAT",     &Y[13],  "EYAT/D");
    //  tree->Branch("EZAT",     &Z[13],  "EZAT/D");
    //  tree->Branch("EXpAT",   &Xp[13], "EXpAT/D");
    //  tree->Branch("EYpAT",   &Yp[13], "EYpAT/D");
    //  tree->Branch("EPxAT",   &Px[13], "EPxAT/D");
    //  tree->Branch("EPyAT",   &Py[13], "EPyAT/D");
    
    tree->Branch("EPSepi",    &P[1],  "EPSepi/D");
    tree->Branch("EXSepi",    &X[1],  "EXSepi/D");
    tree->Branch("EYSepi",    &Y[1],  "EYSepi/D");
    tree->Branch("EZSepi",    &Z[1],  "EZSepi/D");
    tree->Branch("EXpSepi",  &Xp[1], "EXpSepi/D");
    tree->Branch("EYpSepi",  &Yp[1], "EYpSepi/D");
    tree->Branch("EPSepe",    &P[2],  "EPSepe/D");
    tree->Branch("EXSepe",    &X[2],  "EXSepe/D");
    tree->Branch("EYSepe",    &Y[2],  "EYSepe/D");
    tree->Branch("EZSepe",    &Z[2],  "EZSepe/D");
    tree->Branch("EXpSepe",  &Xp[2], "EXpSepe/D");
    tree->Branch("EYpSepe",  &Yp[2], "EYpSepe/D");
    tree->Branch("EPQ1i",    &P[3],  "EPQ1i/D");
    tree->Branch("EXQ1i",    &X[3],  "EXQ1i/D");
    tree->Branch("EYQ1i",    &Y[3],  "EYQ1i/D");
    tree->Branch("EZQ1i",    &Z[3],  "EZQ1i/D");
    tree->Branch("EXpQ1i",  &Xp[3], "EXpQ1i/D");
    tree->Branch("EYpQ1i",  &Yp[3], "EYpQ1i/D");
    tree->Branch("EPQ1e",    &P[4],  "EPQ1e/D");
    tree->Branch("EXQ1e",    &X[4],  "EXQ1e/D");
    tree->Branch("EYQ1e",    &Y[4],  "EYQ1e/D");
    tree->Branch("EZQ1e",    &Z[4],  "EZQ1e/D");
    tree->Branch("EXpQ1e",  &Xp[4], "EXpQ1e/D");
    tree->Branch("EYpQ1e",  &Yp[4], "EYpQ1e/D");
    tree->Branch("EPQ2i",    &P[5],  "EPQ2i/D");
    tree->Branch("EXQ2i",    &X[5],  "EXQ2i/D");
    tree->Branch("EYQ2i",    &Y[5],  "EYQ2i/D");
    tree->Branch("EZQ2i",    &Z[5],  "EZQ2i/D");
    tree->Branch("EXpQ2i",  &Xp[5], "EXpQ2i/D");
    tree->Branch("EYpQ2i",  &Yp[5], "EYpQ2i/D");
    tree->Branch("EPQ2e",    &P[6],  "EPQ2e/D");
    tree->Branch("EXQ2e",    &X[6],  "EXQ2e/D");
    tree->Branch("EYQ2e",    &Y[6],  "EYQ2e/D");
    tree->Branch("EZQ2e",    &Z[6],  "EZQ2e/D");
    tree->Branch("EXpQ2e",  &Xp[6], "EXpQ2e/D");
    tree->Branch("EYpQ2e",  &Yp[6], "EYpQ2e/D");
    tree->Branch("EPDi",     &P[7],   "EPDi/D");//toshi , 10Aug2012
    tree->Branch("EXDi",     &X[7],   "EXDi/D");
    tree->Branch("EYDi",     &Y[7],   "EYDi/D");
    tree->Branch("EZDi",     &Z[7],   "EZDi/D");
    tree->Branch("EXpDi",   &Xp[7],  "EXpDi/D");
    tree->Branch("EYpDi",   &Yp[7],  "EYpDi/D");
    tree->Branch("EPDe",     &P[8],   "EPDe/D");//toshi , 10Aug2012
    tree->Branch("EXDe",     &X[8],   "EXDe/D");
    tree->Branch("EYDe",     &Y[8],   "EYDe/D");
    tree->Branch("EZDe",     &Z[8],   "EZDe/D");
    tree->Branch("EXpDe",   &Xp[8],  "EXpDe/D");
    tree->Branch("EYpDe",   &Yp[8],  "EYpDe/D");
    tree->Branch("EPQ3i",    &P[9],  "EPQ3i/D");//fujita , 11Aug2014
    tree->Branch("EXQ3i",    &X[9],  "EXQ3i/D");
    tree->Branch("EYQ3i",    &Y[9],  "EYQ3i/D");
    tree->Branch("EZQ3i",    &Z[9],  "EZQ3i/D");
    tree->Branch("EXpQ3i",  &Xp[9], "EXpQ3i/D");
    tree->Branch("EYpQ3i",  &Yp[9], "EYpQ3i/D");
    tree->Branch("EPQ3e",    &P[10],  "EPQ3e/D");//fujita , 11Aug2014
    tree->Branch("EXQ3e",    &X[10],  "EXQ3e/D");
    tree->Branch("EYQ3e",    &Y[10],  "EYQ3e/D");
    tree->Branch("EZQ3e",    &Z[10],  "EZQ3e/D");
    tree->Branch("EXpQ3e",  &Xp[10], "EXpQ3e/D");
    tree->Branch("EYpQ3e",  &Yp[10], "EYpQ3e/D");
    tree->Branch("EPFP",     &P[11],  "EPFP/D"); // RP, fujita , 11Aug2014
    tree->Branch("EPxFP",   &Px[11], "EPxFP/D");
    tree->Branch("EPyFP",   &Py[11], "EPyFP/D");
    tree->Branch("EXFP",     &X[11],  "EXFP/D");
    tree->Branch("EYFP",     &Y[11],  "EYFP/D");
    tree->Branch("EZFP",     &Z[11],  "EZFP/D");
    tree->Branch("EXgFP",    &Xg[11], "EXgFP/D"); /////suzuki
    tree->Branch("EYgFP",    &Yg[11], "EYgFP/D");
    tree->Branch("EZgFP",    &Zg[11], "EZgFP/D");
    tree->Branch("EXpFP",   &Xp[11], "EXpFP/D");
    tree->Branch("EYpFP",   &Yp[11], "EYpFP/D");
	tree->Branch("EtrackIDFP", &trackID[11], "EtrackID/I");
	tree->Branch("EchargeFP", &charge[11], "Echarge/I");
	tree->Branch("EparticleIDFP", &particleID[11], "EparticleID/I");
    tree->Branch("EPFP1",     &P[12],  "EPFP1/D"); // VDC1, aida 20180105
    tree->Branch("EPxFP1",   &Px[12], "EPxFP1/D");
    tree->Branch("EPyFP1",   &Py[12], "EPyFP1/D");
    tree->Branch("EXFP1",     &X[12],  "EXFP1/D");
    tree->Branch("EYFP1",     &Y[12],  "EYFP1/D");
    tree->Branch("EZFP1",     &Z[12],  "EZFP1/D");
    tree->Branch("EXgFP1",    &Xg[12], "EXgFP1/D");  //////suzuki
    tree->Branch("EYgFP1",    &Yg[12], "EYgFP1/D");
    tree->Branch("EZgFP1",    &Zg[12], "EZgFP1/D");
    tree->Branch("EXpFP1",   &Xp[12], "EXpFP1/D");
    tree->Branch("EYpFP1",   &Yp[12], "EYpFP1/D");
	tree->Branch("EtrackIDFP1", &trackID[12], "EtrackID/I");
	tree->Branch("EchargeFP1", &charge[12], "Echarge/I");
	tree->Branch("EparticleIDFP1", &particleID[12], "EparticleID/I");
    tree->Branch("EPFP2",     &P[13],  "EPFP2/D"); // VDC2, aida 20180105
    tree->Branch("EPxFP2",   &Px[13], "EPxFP2/D");
    tree->Branch("EPyFP2",   &Py[13], "EPyFP2/D");
    tree->Branch("EXFP2",     &X[13],  "EXFP2/D");
    tree->Branch("EYFP2",     &Y[13],  "EYFP2/D");
    tree->Branch("EZFP2",     &Z[13],  "EZFP2/D");
    tree->Branch("EXgFP2",    &Xg[13], "EXgFP2/D");   /////suzuki
    tree->Branch("EYgFP2",    &Yg[13], "EYgFP2/D");
    tree->Branch("EZgFP2",    &Zg[13], "EZgFP2/D");
    tree->Branch("EXpFP2",   &Xp[13], "EXpFP2/D");
    tree->Branch("EYpFP2",   &Yp[13], "EYpFP2/D");
	tree->Branch("EtrackIDFP2", &trackID[13], "EtrackID/I");
	tree->Branch("EchargeFP2", &charge[13], "Echarge/I");
	tree->Branch("EparticleIDFP2", &particleID[13], "EparticleID/I");
    tree->Branch("EPS2",     &P[15],  "EPS2/D"); // S2, suzuki 20191012
    tree->Branch("EPxS2",   &Px[15], "EPxS2/D");
    tree->Branch("EPyS2",   &Py[15], "EPyS2/D");
    tree->Branch("EXS2",     &X[15],  "EXS2/D");
    tree->Branch("EYS2",     &Y[15],  "EYS2/D");
    tree->Branch("EZS2",     &Z[15],  "EZS2/D");
    //  tree->Branch("EPFP1",    &P[7],  "EPFP1/D");
    //  tree->Branch("ELFP1",    &L[7],  "ELFP1/D");
    //  tree->Branch("EXFP1",    &X[7],  "EXFP1/D");
    //  tree->Branch("EYFP1",    &Y[7],  "EYFP1/D");
    //  tree->Branch("EZFP1",    &Z[7],  "EZFP1/D");
    //  tree->Branch("EXpFP1",  &Xp[7], "EXpFP1/D");
    //  tree->Branch("EYpFP1",  &Yp[7], "EYpFP1/D");
    //  tree->Branch("EXvFP1",  &Xv[7], "EXvFP1/D");
    //  tree->Branch("EYvFP1",  &Yv[7], "EYvFP1/D");
    //  tree->Branch("EZvFP1",  &Zv[7], "EZvFP1/D");
    //  tree->Branch("EPFP2",    &P[8],  "EPFP2/D");
    //  tree->Branch("ELFP2",    &L[8],  "ELFP2/D");
    //  tree->Branch("EXFP2",    &X[8],  "EXFP2/D");
    //  tree->Branch("EYFP2",    &Y[8],  "EYFP2/D");
    //  tree->Branch("EZFP2",    &Z[8],  "EZFP2/D");
    //  tree->Branch("EXpFP2",  &Xp[8], "EXpFP2/D");
    //  tree->Branch("EYpFP2",  &Yp[8], "EYpFP2/D");
    //  tree->Branch("EXvFP2",  &Xv[8], "EXvFP2/D");
    //  tree->Branch("EYvFP2",  &Yv[8], "EYvFP2/D");
    //  tree->Branch("EZvFP2",  &Zv[8], "EZvFP2/D");
    //  tree->Branch("EPFPS3",    &P[9],  "EPFP3/D");
    //  tree->Branch("ELFP3",    &L[9],  "ELFP3/D");
    //  tree->Branch("EXFP3",    &X[9],  "EXFP3/D");
    //  tree->Branch("EYFP3",    &Y[9],  "EYFP3/D");
    //  tree->Branch("EZFP3",    &Z[9],  "EZFP3/D");
    //  tree->Branch("EXpFP3",  &Xp[9], "EXpFP3/D");
    //  tree->Branch("EYpFP3",  &Yp[9], "EYpFP3/D");
    //  tree->Branch("EXvFP3",  &Xv[9], "EXvFP3/D");
    //  tree->Branch("EYvFP3",  &Yv[9], "EYvFP3/D");
    //  tree->Branch("EZvFP3",  &Zv[9], "EZvFP3/D");
    //  tree->Branch("EPFP4",    &P[10],  "EPFP4/D");
    //  tree->Branch("ELFP4",    &L[10],  "ELFP4/D");
    //  tree->Branch("EXFP4",   &X[10],  "EXFP4/D");
    //  tree->Branch("EYFP4",   &Y[10],  "EYFP4/D");
    //  tree->Branch("EZFP4",   &Z[10],  "EZFP4/D");
    //  tree->Branch("EXpFP4", &Xp[10], "EXpFP4/D");
    //  tree->Branch("EYpFP4", &Yp[10], "EYpFP4/D");
    //  tree->Branch("EXvFP4", &Xv[10], "EXvFP4/D");
    //  tree->Branch("EYvFP4", &Yv[10], "EYvFP4/D");
    //  tree->Branch("EZvFP4", &Zv[10], "EZvFP4/D");
    //  tree->Branch("EPFP5",   &P[11],  "EPFP5/D");
    //  tree->Branch("ELFP5",   &L[11],  "ELFP5/D");
    //  tree->Branch("EXFP5",   &X[11],  "EXFP5/D");
    //  tree->Branch("EYFP5",   &Y[11],  "EYFP5/D");
    //  tree->Branch("EZFP5",   &Z[11],  "EZFP5/D");
    //  tree->Branch("EXpFP5", &Xp[11], "EXpFP5/D");
    //  tree->Branch("EYpFP5", &Yp[11], "EYpFP5/D");
    //  tree->Branch("EXvFP5", &Xv[11], "EXvFP5/D");
    //  tree->Branch("EYvFP5", &Yv[11], "EYvFP5/D");
    //  tree->Branch("EZvFP5", &Zv[11], "EZvFP5/D");
    //  tree->Branch("EXSS",   &X[12], "EXSS/D");
    //  tree->Branch("EYSS",   &Y[12], "EYSS/D");
    //  tree->Branch("EZSS",   &Z[12], "EZSS/D");
    //  tree->Branch("EXpSS",   &Xp[12], "EXpSS/D");
    //  tree->Branch("EYpSS",   &Yp[12], "EYpSS/D");
    //  tree->Branch("EXvSS", &Xv[12], "EXvSS/D");
    //  tree->Branch("EYvSS", &Yv[12], "EYvSS/D");
    //  tree->Branch("EZvSS", &Zv[12],   "EZvSS/D");
    //  tree->Branch("EXVC",   &X[13],  "EXVC/D");
    //  tree->Branch("EYVC",   &Y[13],  "EYVC/D");
    //  tree->Branch("EZVC",   &Z[13],  "EZVC/D");
    //  tree->Branch("EXpVC",   &Xp[13],  "EXpVC/D");
    //  tree->Branch("EYpVC",   &Yp[13],  "EYpVC/D");
    //  tree->Branch("EXvVC", &Xv[13], "EXvVC/D");
    //  tree->Branch("EYvVC", &Yv[13], "EYvVC/D");
    //  tree->Branch("EZvVC", &Zv[13], "EZvVC/D");
    //  tree->Branch("EPSPLP",   &P[14],  "EPSPLP/D");
    //  tree->Branch("EXSPLP",  &X[14], "EXSPLP/D");
    //  tree->Branch("EYSPLP",  &Y[14], "EYSPLP/D");
    //  tree->Branch("EXpSPLP",  &Xp[14], "EXpSPLP/D");
    //  tree->Branch("EYpSPLP",  &Yp[14], "EYpSPLP/D");
    //  tree->Branch("EXvSPLP",&Xv[14],"EXvSPLP/D");
    //  tree->Branch("EYvSPLP",&Yv[14],"EYvSPLP/D");
    //  tree->Branch("EZvSPLP",&Zv[14],"EZvSPLP/D");
    //  tree->Branch("EPSPLP",   &P[14],  "EPSPLP/D");
    /************* comment out by suzuki from  
     //===Drift Chamber===//
    tree->Branch("PEDCFP",   &EDCFPP,  "PEDCFP/D");
    tree->Branch("LEDCFP",   &EDCFPL,  "LEDCFP/D");
    tree->Branch("XEDCFP",   &EDCFPX,  "XEDCFP/D");
    tree->Branch("YEDCFP",   &EDCFPY,  "YEDCFP/D");
    tree->Branch("ZEDCFP",   &EDCFPZ,  "ZEDCFP/D");
    tree->Branch("XpEDCFP",   &EDCFPXp,  "XpEDCFP/D");
    tree->Branch("YpEDCFP",   &EDCFPYp,  "YpEDCFP/D");
    for (int i=0;i<10;i++){
      sprintf(sPEDC,"PEDC%d",i+1); 
      sprintf(sPEDCD,"PEDC%d/D",i+1); 
      sprintf(sXEDC,"XEDC%d",i+1); 
      sprintf(sXEDCD,"XEDC%d/D",i+1); 
      sprintf(sYEDC,"YEDC%d",i+1); 
      sprintf(sYEDCD,"YEDC%d/D",i+1); 
      sprintf(sZEDC,"ZEDC%d",i+1); 
      sprintf(sZEDCD,"ZEDC%d/D",i+1); 
      sprintf(sXpEDC,"XpEDC%d",i+1); 
      sprintf(sXpEDCD,"XpEDC%d/D",i+1); 
      sprintf(sYpEDC,"YpEDC%d",i+1); 
      sprintf(sYpEDCD,"YpEDC%d/D",i+1); 
      tree->Branch(sPEDC,    &EDC1P[i],   sPEDCD);
      tree->Branch(sXEDC,    &EDC1X[i],   sXEDCD);
      tree->Branch(sYEDC,    &EDC1Y[i],   sYEDCD);
      tree->Branch(sZEDC,    &EDC1Z[i],   sZEDCD);
      tree->Branch(sXpEDC,   &EDC1Xp[i],   sXpEDCD);
      tree->Branch(sYpEDC,   &EDC1Yp[i],   sYpEDCD);
    }
    for (int i=0;i<6;i++){
      sprintf(sPEDC,"PEDC%d",i+11); 
      sprintf(sPEDCD,"PEDC%d/D",i+11); 
      sprintf(sXEDC,"XEDC%d",i+11); 
      sprintf(sXEDCD,"XEDC%d/D",i+11); 
      sprintf(sYEDC,"YEDC%d",i+11); 
      sprintf(sYEDCD,"YEDC%d/D",i+11); 
      sprintf(sZEDC,"ZEDC%d",i+11); 
      sprintf(sZEDCD,"ZEDC%d/D",i+11); 
      sprintf(sXpEDC,"XpEDC%d",i+11); 
      sprintf(sXpEDCD,"XpEDC%d/D",i+11); 
      sprintf(sYpEDC,"YpEDC%d",i+11); 
      sprintf(sYpEDCD,"YpEDC%d/D",i+11); 
      tree->Branch(sPEDC,    &EDC2P[i],   sPEDCD);
      tree->Branch(sXEDC,    &EDC2X[i],   sXEDCD);
      tree->Branch(sYEDC,    &EDC2Y[i],   sYEDCD);
      tree->Branch(sZEDC,    &EDC2Z[i],   sZEDCD);
      tree->Branch(sXpEDC,   &EDC2Xp[i],   sXpEDCD);
      tree->Branch(sYpEDC,   &EDC2Yp[i],   sYpEDCD);
    }
    //===HRS HodoScope===//
    for (int i=0;i<29;i++){
      
      sprintf(sPEH,"PEH%d",100+i+1); 
      sprintf(sPEHD,"PEH%d/D",100+i+1); 
      sprintf(sXEH,"XEH%d",100+i+1); 
      sprintf(sXEHD,"XEH%d/D",100+i+1); 
      sprintf(sYEH,"YEH%d",100+i+1); 
      sprintf(sYEHD,"YEH%d/D",100+i+1); 
      sprintf(sZEH,"ZEH%d",100+i+1); 
      sprintf(sZEHD,"ZEH%d/D",100+i+1); 
      sprintf(sXpEH,"XpEH%d",100+i+1); 
      sprintf(sXpEHD,"XpEH%d/D",100+i+1); 
      sprintf(sYpEH,"YpEH%d",100+i+1); 
      sprintf(sYpEHD,"YpEH%d/D",100+i+1); 
      sprintf(sAEH,"AEH%d",100+i+1); 
      sprintf(sAEHD,"AEH%d/D",100+i+1); 
      sprintf(sTEH,"TEH%d",100+i+1); 
      sprintf(sTEHD,"TEH%d/D",100+i+1); 
      tree->Branch(sPEH,    &EH1P[i],   sPEHD);
      tree->Branch(sXEH,    &EH1X[i],   sXEHD);
      tree->Branch(sYEH,    &EH1Y[i],   sYEHD);
      tree->Branch(sZEH,    &EH1Z[i],   sZEHD);
      tree->Branch(sXpEH,   &EH1Xp[i],   sXpEHD);
      tree->Branch(sYpEH,   &EH1Yp[i],   sYpEHD);
      tree->Branch(sAEH,    &EH1A[i],   sAEHD);
      tree->Branch(sTEH,    &EH1T[i],   sTEHD);
      
      sprintf(sPEH,"PEH%d",200+i+1); 
      sprintf(sPEHD,"PEH%d/D",200+i+1); 
      sprintf(sXEH,"XEH%d",200+i+1); 
      sprintf(sXEHD,"XEH%d/D",200+i+1); 
      sprintf(sYEH,"YEH%d",200+i+1); 
      sprintf(sYEHD,"YEH%d/D",200+i+1); 
      sprintf(sZEH,"ZEH%d",200+i+1); 
      sprintf(sZEHD,"ZEH%d/D",200+i+1); 
      sprintf(sXpEH,"XpEH%d",200+i+1); 
      sprintf(sXpEHD,"XpEH%d/D",200+i+1); 
      sprintf(sYpEH,"YpEH%d",200+i+1); 
      sprintf(sYpEHD,"YpEH%d/D",200+i+1); 
      sprintf(sAEH,"AEH%d",200+i+1); 
      sprintf(sAEHD,"AEH%d/D",200+i+1); 
      sprintf(sTEH,"TEH%d",200+i+1); 
      sprintf(sTEHD,"TEH%d/D",200+i+1); 
      tree->Branch(sPEH,    &EH2P[i],   sPEHD);
      tree->Branch(sXEH,    &EH2X[i],   sXEHD);
      tree->Branch(sYEH,    &EH2Y[i],   sYEHD);
      tree->Branch(sZEH,    &EH2Z[i],   sZEHD);
      tree->Branch(sXpEH,   &EH2Xp[i],   sXpEHD);
      tree->Branch(sYpEH,   &EH2Yp[i],   sYpEHD);
      tree->Branch(sAEH,    &EH2A[i],   sAEHD);
      tree->Branch(sTEH,    &EH2T[i],   sTEHD);
    }
    */  //****************    comment out by suzuki to
    //===Trigger==// 
    tree->Branch("EVDTrig",&VDTrig,"EVDTrig/B");
    tree->Branch("EVtxTrig",&VtxTrig,"EVtxTrig/B");
    tree->Branch("ESPLTrig",&SPLTrig,"ESPLTrig/B");
    tree->Branch("EDCTrig",&EDCTrig,"EDCTrig/B");
    tree->Branch("PBTrig",&PBTrig,"PBTrig/B");
    tree->Branch("EHTrig",&EHTrig,"EHTrig/B");
    tree->Branch("EDetTrig",&EDetTrig,"EDetTrig/B");
    
  }//if(Event Flag)
}

/////////////////////////////////////////////////////
void HRSAnalysis::SetDataFile(const char *datafile )
/////////////////////////////////////////////////////
{
  DataFile_.open( datafile );
}

/////////////////////////////////////////////////////////////////
G4double HRSAnalysis::RandGauss( G4double center, G4double dev )
/////////////////////////////////////////////////////////////////
{
  G4double rand1 = G4UniformRand();
  G4double rand2 = G4UniformRand();

  G4double a = sqrt(-2.0*log(rand1)) * cos(2.0*M_PI*rand2);

  return dev*a + center;
}
