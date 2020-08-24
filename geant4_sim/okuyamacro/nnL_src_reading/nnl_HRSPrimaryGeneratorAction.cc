/*
  HRSPrimaryGeneratorAction.cc
  Modified by T. Gogami (Feb 5, 2017)
*/

#include "HRSPrimaryGeneratorAction.hh"
#include "HRSAnalysis.hh"
#include "HRSParamMan.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>

int gevID = 0;
G4ParticleDefinition* particle;
//Particle ID

/////////////////////////////////////////////////////////////
HRSPrimaryGeneratorAction::
HRSPrimaryGeneratorAction( HRSDetectorConstruction* , 
			   HRSAnalysis* analysisManager )
  :G4VUserPrimaryGeneratorAction(), anaMan(analysisManager)
/////////////////////////////////////////////////////////////
{
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  int pID = paramMan->GetParticleID();
  G4ParticleTable* particleTable 
                  = G4ParticleTable::GetParticleTable();
  if ( pID==1 ){
    particle = particleTable->FindParticle("e-");
  }
  else if (pID==2){
    particle = particleTable->FindParticle("kaon+");
  }
  else if (pID==3){
    particle = particleTable->FindParticle("pi+");
  }
  else{
    particle = particleTable->FindParticle("geantino");
  }
  dataFile=fopen( paramMan->GetInputDataName(), "r");
  //CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine);
  //CLHEP::HepRandom::setTheSeed(time(NULL),4);

}

/////////////////////////////////////////////////////////
HRSPrimaryGeneratorAction::~HRSPrimaryGeneratorAction()
////////////////////////////////////////////////////////
{
  if ( dataFile ) fclose( dataFile );
}

//Generator event action
//////////////////////////////////////////////////////////////////////
void HRSPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent )
////////////////////////////////////////////////////////////////////
{
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  int generateID = paramMan->GetGenerationID();
  if ( generateID == 0 ) {
    return GenParUni( anEvent );
  }
  else if ( generateID == 1 ) {
    return GenFromFile( anEvent );
  }
  else if ( generateID == 2 || generateID == 3 ){
	return GenFromDist( anEvent );
  }
}

///////////////////////////////////////////////////////////////
void HRSPrimaryGeneratorAction::GenerateBeam( G4Event* )
///////////////////////////////////////////////////////////////
{
/*  G4ThreeVector prodPoint(0.*cm, 0.*cm, -1.*cm);
  G4ThreeVector BeamMom(0.*GeV, 0.*GeV, 2.5*GeV);

  gunParticle-> SetParticleEnergy(0.0); // reset
  gunParticle-> SetParticleMomentum(BeamMom);
  gunParticle-> SetParticlePosition(prodPoint);
  gunParticle-> GeneratePrimaryVertex(anEvent);

  if( anaMan ){
    anaMan-> PrimaryGeneration(prodPoint, BeamMom);
  }*/
}

/////////////////////////////////////////////////////////////////////
void HRSPrimaryGeneratorAction::GenParUni( G4Event* anEvent )
/////////////////////////////////////////////////////////////////////
{
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  G4int evID = anEvent->GetEventID();
  gevID = evID;
  G4int StepMomFlag = paramMan->GetStepMomFlag();
  G4int StepPhiFlag = paramMan->GetStepPhiFlag();
  G4int StepHorFlag = paramMan->GetStepHorFlag();
  G4int StepVerFlag = paramMan->GetStepVerFlag();
  G4int ChangeTargetPosFlag = paramMan->GetChangeTargetPosFlag();
  
  G4double prepx,prepy,prepz; 
  G4double px,py,pz; 
  G4double E_Theta = paramMan->GetRotAngle()*CLHEP::degree;
  G4double E_Theta_gen = 13.2*CLHEP::degree; //"I want to use E_Theta but E_Theta is not well defined in this program so I set new E_theta" //default
  G4double CentMom = paramMan->GetCentMom();
  G4double AcptMom = paramMan->GetAcptMom();
  G4double CentTheta = paramMan->GetCentTheta();
  G4double AcptTheta = paramMan->GetAcptTheta();
  G4double CentPhi = paramMan->GetCentPhi();
  G4double AcptPhi;
	if((AcptPhi = paramMan->GetAcptPhi())==3.1415){
		AcptPhi = PI;
	}
  G4double AcptXp = paramMan->GetAcptXp();
  G4double AcptYp = paramMan->GetAcptYp();
  G4double xR = paramMan->GetRasterX();
  G4double yR = paramMan->GetRasterY();
  G4double zR = paramMan->GetRasterZ();
	G4double Mom = 0.;
	if(StepMomFlag==0){
		Mom = (CentMom + (G4UniformRand()-0.5)*2.0*AcptMom);
	}
	else if(StepMomFlag==1){
		const G4double minmom = 94.; // [%]
		const G4double maxmom = 106.; // [%]
		const G4double step = 1.;
		const int n = (int)((maxmom - minmom)/step) + 1;
		Mom = CentMom*(minmom + step*(evID%n))/100.; 
//		G4cout<< "Mom = " << Mom <<G4endl;
	}
	else if(StepMomFlag==2){
		const G4double minmom = 94.; // [%]
		const G4double maxmom = 106.; // [%]
		const G4double step = 1.;
		const int n = (int)((maxmom - minmom)/step) + 1; 
//		G4cout << "n = " << n << G4endl;
		Mom = CentMom*(minmom + step*(rand()%n))/100.; 
//		G4cout<< "Mom = " << Mom <<G4endl;
	}
  G4double xp = 0.;
  G4double yp = 0.;
  G4int RanSign = 0;
  G4double Theta = 0.;
  G4double Phi = 0.;
//  G4double min = cos(CentTheta+AcptTheta);
//  G4double max = cos(CentTheta-AcptTheta);
//  G4double da = (max-min)/2. ;
//  G4double a = min + da ;
  G4double min = cos(AcptTheta);
  G4double max = 1.;
  G4double da = (max-min)/2. ;
  G4double a = min + da ;
  G4ThreeVector prepvec;


  ////////// Define ANGLE
	if(StepHorFlag==0 && StepVerFlag==0){    
		while(1){
			Theta = acos( a + 2.*da*(G4UniformRand()-0.5) ) ;
			if(CentPhi!=0. || AcptPhi!=PI){
				G4cerr << "HRSPrimaryGeneratorAction: Phi Central should be 0. and Phi Width shoudle be 3.1415" << G4endl;
				exit( -1 );
			}
			if(StepPhiFlag==0){
				Phi = CentPhi + (2.*AcptPhi*(G4UniformRand() - 0.5));
				Theta = acos( a + 2.*da*(G4UniformRand()-0.5) ) ;
			}
			else{
				G4double stepsize = 45.*PI/180.; // [rad]
				Theta = AcptTheta;
				Phi = 0. + stepsize*evID;
			}
			prepvec.setRThetaPhi(Mom, Theta, Phi);
			G4RotationMatrix *rot = new  G4RotationMatrix();
			//			rot->rotateY(-CentTheta); //comment out by suzuki
			rot->rotateY(-E_Theta_gen); /////suzuki 
			G4AffineTransform *affine = new G4AffineTransform( rot, G4ThreeVector(0,0,0) );
			///// G4AffineTransform *affine = new G4AffineTransform( rot, G4ThreeVector(0,0,67.001) ); /////suzuki
 			G4ThreeVector pvec = affine->TransformPoint(prepvec);
			px = pvec.getX();
			py = pvec.getY();
			pz = pvec.getZ();
			double theta_ = acos(pz/pvec.getR());
//			G4cout << evID << " " << theta_ << G4endl;
			delete rot, affine;
//			if( StepPhiFlag==0 && (0.08<theta_ && theta_<0.17) ) break;
			break;
		} // while(1)
	}
	else if(StepHorFlag!=0 && StepVerFlag==0){
		const G4double minxp = -0.06; // [rad]
		const G4double maxxp = 0.06; // [rad]
		const G4double stepxp = 0.020;
		const int n = (int)((maxxp - minxp)/stepxp) + 1;
		if(StepHorFlag==1){
			xp = minxp + stepxp*(evID%n); 
//			G4cout<< "xp = " << xp <<G4endl;
		}
		else if(StepHorFlag==2){
			xp = minxp + stepxp*(rand()%n); 
		}
		prepz = Mom/sqrt(1. + xp*xp + yp*yp);
		prepx = prepz*xp;
		prepy = prepz*yp;
		prepvec.set(prepx, prepy, prepz);
		G4RotationMatrix *rot = new  G4RotationMatrix();
		//		rot->rotateY(-CentTheta);  /////comment out by suzuki
		rot->rotateY(E_Theta_gen); ///// suzuki  cneter = 13.2degree
		G4AffineTransform *affine = new G4AffineTransform( rot, G4ThreeVector(0,0,0) );
		G4ThreeVector pvec = affine->TransformPoint(prepvec);
		px = pvec.getX();
		py = pvec.getY();
		pz = pvec.getZ();
	}
	else if(StepHorFlag==0 && StepVerFlag!=0){
		const G4double minyp = -0.05; // [rad]
		const G4double maxyp = 0.05; // [rad]
		const G4double stepyp = 0.005;
		const int n = (int)((maxyp - minyp)/stepyp) + 1;
//		G4cout << "n = " << n << G4endl;
		if(StepVerFlag==1){
			yp = minyp + stepyp*(evID%n); 
//			G4cout << "yp = " << yp << G4endl;
//			getchar();
		}
		else if(StepVerFlag==2){
			yp = minyp + stepyp*(rand()%n); 
		}
		prepz = Mom/sqrt(1. + xp*xp + yp*yp);
		prepx = prepz*xp;
		prepy = prepz*yp;
		prepvec.set(prepx, prepy, prepz);
		G4RotationMatrix *rot = new  G4RotationMatrix();
		rot->rotateY(-CentTheta);
		G4AffineTransform *affine = new G4AffineTransform( rot, G4ThreeVector(0,0,0) );
		G4ThreeVector pvec = affine->TransformPoint(prepvec);
		px = pvec.getX();
		py = pvec.getY();
		pz = pvec.getZ();
	}
	else{
		G4cerr << "HRSPrimaryGeneratorAction: kouji tyuu" << G4endl;
		exit( -1 );
	}
	

	/////// Define POSITION 

	//////////////////////////////////////////////////////////////////////////////////////
	/////Setting for Acceptance study Changing distance pivot to Q1entrance //////////////
	/////                       Change PIVOTtoGEN                           //////////////        
	G4double VTtoPIVOT_z = 0;  /////suzuki for nnl
      
	G4double PIVOTtoGEN = 0.0; ///// Distance from PIVOT to Generation point  a value 0.0 means PIVOT  //suzuki
	G4ThreeVector prepos_gen(0, 0, PIVOTtoGEN);  //suzuki
	G4RotationMatrix *genRot = new G4RotationMatrix; //suzuki
	genRot->rotateY(-E_Theta_gen);
	G4AffineTransform *genAffine = new G4AffineTransform(genRot, G4ThreeVector(0, 0, VTtoPIVOT_z)); //suzuki
	
	G4ThreeVector pos_gen =genAffine->TransformPoint( prepos_gen ); //suzuki

	G4double x = pos_gen.x() + (G4UniformRand()-0.5)*xR*2.; // +-xR
	G4double y = pos_gen.y() + (G4UniformRand()-0.5)*yR*2.; // +-yR
	G4double z = pos_gen.z() + (G4UniformRand()-0.5)*zR*2.; // +-zR ///// suzuki for nnl
	//	std::cout << "debug  = (" << x << ", " << ", " << y << ", " << z << ")" << std::endl; //debug
	//	G4double z = VTtoPIVOT_z + (G4UniformRand()-0.5)*zR*2.; // +-zR ///// suzuki for nnl

	
	//  G4double z = 1.0*zR*2.;

/////////////////////////////////////////
#if 0
	x = 0.;
	y = 0.;
	z = 0.;
	px = 0.;
	py = 0.;
	pz = 3.0296;
//	G4cout<< "evID = " << evID <<G4endl;
//	if(evID%3==0){
//		px = 5.;
//		py = 0.;
//		pz = 0.;
//	}
//	else if(evID%3==1){
//		px = 0.;
//		py = 5.;
//		pz = 0.;
//	}
//	else if(evID%3==2){
//		px = 0.;
//		py = 0.;
//		pz = 5.;
//	}
#endif
/////////////////////////////////////////

	G4ThreeVector gPos, gMom;
	if(ChangeTargetPosFlag==0){
		gMom.set(px*CLHEP::GeV, py*CLHEP::GeV, pz*CLHEP::GeV);
		gPos.set(x*CLHEP::cm, y*CLHEP::cm, z*CLHEP::cm);
	}
/*
	else{
		if(ChangeTargetPosFlag==1){
			gPos = paramMan->GetpQ1();
		}
		else if(ChangeTargetPosFlag==2){
			gPos = paramMan->GetpQ2();
		}
		gMom.set(0.*CLHEP::GeV, 0.*CLHEP::GeV, Mom*CLHEP::GeV);
		G4ThreeVector XAxis(1., 0., 0.);
		gMom.rotateY(paramMan->GetE_Theta() - CentTheta + Theta);
		XAxis.rotateY(paramMan->GetE_Theta() - CentTheta + Theta);
		gMom.rotate(-Phi, XAxis); // phi means elevation angle
	}
*/
//	G4cout<< "gPos = " << gPos.x()/CLHEP::cm << ", " << gPos.y()/CLHEP::cm << ", " << gPos.z()/CLHEP::cm <<G4endl;
//  G4cout<< "gMom " << gMom.x()/CLHEP::GeV << " " << gMom.y()/CLHEP::GeV << " " << gMom.z()/CLHEP::GeV <<G4endl;
//	G4cout<< "E_Theta = " << paramMan->GetE_Theta() <<G4endl;
//  SetMom(gMom,gPos,anEvent,0,z);
//  SetMom(gMom,gPos,anEvent,evID,z);
  SetMom(gMom, gPos, anEvent, evID, z, Theta, Phi); // aida 20171209

  //G4ThreeVector gMomtmp(0.*GeV, 0.*GeV, 1.*GeV);
  //G4ThreeVector gPostmp(0.*mm, 0.*mm, 0.*mm);
  //SetMom(gMomtmp,gPostmp,anEvent,0,z);
}

////////////////////////////////////////////////////////////////////
void HRSPrimaryGeneratorAction::GenFromFile( G4Event *anEvent )
////////////////////////////////////////////////////////////////////
{
//  const double PI = 4.*atan(1.);
  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  G4double E_Theta = paramMan->GetRotAngle()*CLHEP::degree;
  G4double CentMom = paramMan->GetCentMom();
  G4double xR = paramMan->GetRasterX();
  G4double yR = paramMan->GetRasterY();
  G4double zR = paramMan->GetRasterZ();
  char str[MAXCHAR];
  fgets( str, MAXCHAR, dataFile );
  
  G4double px, py, pz;
  G4double Mom=0.;
  G4double prepx; 
  G4double prepy; 
  G4double prepz; 
  G4double x0, y0, z0;
  G4double t, dp,pe, prob;
  G4double xp, yp;
  G4double xp0, yp0;
  G4double xpe, ype;
  G4double xpk,ypk,pk;
  G4double pbeam0,pbeam1;
  G4int evID;

  //suzuki seed file  
  if( sscanf(str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &evID, &pe, &xpe, &ype, &pk, &xpk, &ypk, &x0, &y0, &z0) == 10 ){
    Mom = pe;
  }  
  else if (sscanf( str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	      &evID,  &x0, &xp0, &y0, &yp0, &pe, &z0, &xpk, &ypk, &pk, &pbeam1, &pbeam0)==12){ // [GeV], [cm], [mrad]
    Mom = pe;
    paramMan->SetpBeam(pbeam0,pbeam1);
  }
 else if (sscanf( str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	      &evID,  &x0, &xp0, &y0, &yp0, &pe, &z0, &xpk, &ypk, &pk)==10){
//    G4cout << x0 << " " << y0 << " " << z0 << " "  << G4endl;
    Mom = pe;
  }
  else if(sscanf( str, "%lf %lf %lf %lf %lf %lf %lf %lf",
		  &x0, &xp0, &y0, &yp0, &z0, &t, &dp, &prob)==8){
    Mom = CentMom*(1.+(dp/100.));
  }
  else{std::cout << "*****I can't read seed file ...*****"<< std::endl;}

  G4double CentTheta = paramMan->GetCentTheta();

//  xp = xp0/1000. -7.0*PI/180. + CentTheta ;//for seed file 
//	xp = xp0/1000. - CentTheta ; // Uniform
//  xp = xp0/1000.; // [rad]
//  yp = yp0/1000.; // [rad]
//  prepz = (Mom/sqrt(1+xp*xp+yp*yp));
//  prepz = Mom/sqrt(1. + xp*xp + yp*yp);
//  prepx = prepz*xp;
//  prepy = prepz*yp;
//  prepx = prepz*xp;
//  prepy = prepz*yp;
//  px = prepx*cos(E_Theta)-prepz*sin(E_Theta);
//  py = prepy;
//  pz = prepx*sin(E_Theta)+prepz*cos(E_Theta);
//	px = prepx;
//	py = prepy;
//	pz = prepz;
  pz = Mom/sqrt(1. + xpe*xpe + ype*ype);
  px = pz*xpe;
  py = pz*ype;
  // ~~~~~~~~ Original ~~~~~~~~~~~~~~~~
  G4double x = xR*(x0*cos(E_Theta)-z0*sin(E_Theta));
  G4double y = yR*y0;
  G4double z = zR*(x0*sin(E_Theta)+z0*cos(E_Theta));

  G4ParticleMomentum gMom( px*CLHEP::GeV, py*CLHEP::GeV, pz*CLHEP::GeV);
//  G4ThreeVector gPos( x*CLHEP::cm, y*CLHEP::cm, z*CLHEP::cm);
  G4ThreeVector gPos( x0*CLHEP::cm, y0*CLHEP::cm, z0*CLHEP::cm);  //org

//  G4ThreeVector gMom( 0*CLHEP::GeV, 0*CLHEP::GeV, 3.0*CLHEP::GeV);
//  G4ThreeVector gPos( 0*CLHEP::cm, 0*CLHEP::cm, 0*CLHEP::cm);  
//  G4cout << gMom.x()/CLHEP::GeV << gMom.y()/CLHEP::GeV << gMom.z()/CLHEP::GeV <<G4endl;
//  G4cout << "######## "<< gPos.x()/CLHEP::cm << gPos.y()/CLHEP::cm << gPos.z()/CLHEP::cm <<G4endl;

  SetMom(gMom,gPos,anEvent,evID,z0);
}




/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void HRSPrimaryGeneratorAction::GenFromDist( G4Event* anEvent )
/////////////////////////////////////////////////////////////////////
{

  HRSParamMan *paramMan = HRSParamMan::GetParamMan();
  G4int evID = anEvent->GetEventID();
  gevID = evID;
  G4double prepx,prepy,prepz; 
  G4double px,py,pz; 
  G4double E_Theta = paramMan->GetRotAngle()*CLHEP::degree;
  G4double E_Theta_gen = 13.2*CLHEP::degree; //"I want to use E_Theta but E_Theta is not well defined in this program so I set new E_theta" //default
  G4double CentMom = paramMan->GetCentMom();
  G4double AcptMom = paramMan->GetAcptMom();
  G4double CentTheta = paramMan->GetCentTheta();
  G4double AcptTheta = paramMan->GetAcptTheta();
  G4double CentPhi = paramMan->GetCentPhi();
  G4double AcptPhi;
	if((AcptPhi = paramMan->GetAcptPhi())==3.1415){
		AcptPhi = PI;
	}
  G4double AcptXp = paramMan->GetAcptXp();
  G4double AcptYp = paramMan->GetAcptYp();
  G4double xR = paramMan->GetRasterX();
  G4double yR = paramMan->GetRasterY();
  G4double zR = paramMan->GetRasterZ();
  G4double Mom = 0.;

  G4double xp = 0.;
  G4double yp = 0.;
  G4int RanSign = 0;
  G4double Theta = 0.;
  G4double Phi = 0.;
  G4double min = cos(AcptTheta);
  G4double max = 1.;
  G4double da = (max-min)/2. ;
  G4double a = min + da ;
  G4ThreeVector prepvec;



  /// DISTRIBUTION
  double mom_min = CentMom - AcptMom;
  double mom_max = CentMom + AcptMom;
  double th_min = (CentTheta - AcptTheta) * (180./TMath::Pi());
  double th_max = (CentTheta + AcptTheta) * (180./TMath::Pi());
  //  TF2* MyDist = new TF2( "MyDist", "my2Dfunc(x,par[0])", mom_min, mom_max, th_min, th_max, 8); // to be checkedw

  TF2* MyDist = new TF2( "MyDist", my2Dfunc, mom_min, mom_max, th_min, th_max, 8); // to be checkedw

  int generateID = paramMan->GetGenerationID();
  if( generateID == 2 ){ // H kine electron (small range fit)
    double setParams[8] = { 3.91*1e5, -5.64*1e2, 2.05*1e3, -1.85*1e1, 6.74*1e0, 3.67*1e-1, -1.62*1e1, 1.90*1e2};
    MyDist->SetParameters(setParams);
  }
  else if( generateID == 3 ){ // pi ( wide range fit)
    double setParams[8] = { 5.067*1e6, -5.895*1e3, 1.529*1e3, 1.217*1e1, 6.032*1e0, -4.809*1e-2, 1.060*1e0, -1.530*1e0};
    MyDist->SetParameters(setParams);
  }
  double funcx, funcy;
  double rnd_max = MyDist->GetMaximumXY(funcx, funcy);
  //  G4cout << "MyFuncMax: " << rnd_max << G4endl;
  
  while(1){
	// Mom
	Mom = (CentMom + (G4UniformRand()-0.5)*2.0*AcptMom);

	// Angle
	Theta = acos( a + 2.*da*(G4UniformRand()-0.5) ) ;
	if(CentPhi!=0. || AcptPhi!=PI){
	  G4cerr << "HRSPrimaryGeneratorAction: Phi Central should be 0. and Phi Width shoudle be 3.1415" << G4endl;
	  exit( -1 );
	}
	Phi = CentPhi + (2.*AcptPhi*(G4UniformRand() - 0.5));
	Theta = acos( a + 2.*da*(G4UniformRand()-0.5) ) ;

	prepvec.setRThetaPhi(Mom, Theta, Phi);
	G4RotationMatrix *rot = new  G4RotationMatrix();
	//			rot->rotateY(-CentTheta); //comment out by suzuki
	rot->rotateY(-E_Theta_gen); /////suzuki 
	G4AffineTransform *affine = new G4AffineTransform( rot, G4ThreeVector(0,0,0) );
	G4ThreeVector pvec = affine->TransformPoint(prepvec);
	px = pvec.getX();
	py = pvec.getY();
	pz = pvec.getZ();
	double theta_ = acos(pz/pvec.getR()) * (180.0/TMath::Pi()); //  azimuth angle

	// Check
	double rnd_tmp = G4UniformRand()*rnd_max;
	//	G4cout << "func:" << MyDist->Eval( Mom, theta_) << " rnd:" << rnd_tmp << G4endl;
	if( MyDist->Eval( Mom, theta_ ) > rnd_tmp ){
	  //	  G4cout << "OK" << G4endl;
	  delete rot, affine;
	  break;
	}
	delete rot, affine;
  }



	

	/////// Define POSITION 

	//////////////////////////////////////////////////////////////////////////////////////
	/////Setting for Acceptance study Changing distance pivot to Q1entrance //////////////
	/////                       Change PIVOTtoGEN                           //////////////        
	G4double VTtoPIVOT_z = 0;  /////suzuki for nnl
      
	G4double PIVOTtoGEN = 0.0; ///// Distance from PIVOT to Generation point  a value 0.0 means PIVOT  //suzuki
	G4ThreeVector prepos_gen(0, 0, PIVOTtoGEN);  //suzuki
	G4RotationMatrix *genRot = new G4RotationMatrix; //suzuki
	genRot->rotateY(-E_Theta_gen);
	G4AffineTransform *genAffine = new G4AffineTransform(genRot, G4ThreeVector(0, 0, VTtoPIVOT_z)); //suzuki
	
	G4ThreeVector pos_gen =genAffine->TransformPoint( prepos_gen ); //suzuki

	G4double x = pos_gen.x() + (G4UniformRand()-0.5)*xR*2.; // +-xR
	G4double y = pos_gen.y() + (G4UniformRand()-0.5)*yR*2.; // +-yR
	G4double z = pos_gen.z() + (G4UniformRand()-0.5)*zR*2.; // +-zR ///// suzuki for nnl


	G4ThreeVector gPos, gMom;
	gMom.set(px*CLHEP::GeV, py*CLHEP::GeV, pz*CLHEP::GeV);
	gPos.set(x*CLHEP::cm, y*CLHEP::cm, z*CLHEP::cm);
		
	SetMom(gMom, gPos, anEvent, evID, z, Theta, Phi); // aida 20171209


}


//////////////////////////////////////////
void HRSPrimaryGeneratorAction::SetMom
( G4ParticleMomentum gMom, G4ThreeVector gPos, G4Event* anEvent,G4int evID,G4double zraster)
//////////////////////////////////////////
{
  //G4cout << gMom.x()/CLHEP::GeV << gMom.y()/CLHEP::GeV << gMom.z()/CLHEP::GeV <<G4endl;
  //G4cout << "#############" << zraster << G4endl;
  G4int n_particle = 1;
  gunParticle = new G4ParticleGun(n_particle);
  gunParticle-> SetParticleEnergy(0.0*CLHEP::GeV); //reset
  gunParticle-> SetParticleDefinition(particle);
  gunParticle-> SetParticleMomentum(gMom);
  gunParticle-> SetParticlePosition(gPos);
  gunParticle-> GeneratePrimaryVertex(anEvent);

  if( anaMan ){
    anaMan->PrimaryGeneration( gPos, gMom, evID ,zraster );
    //G4cout << gPos.x() << G4endl;
  }
  delete gunParticle;

}
//////////////////////////////////////////
void HRSPrimaryGeneratorAction::SetMom
( G4ParticleMomentum gMom, G4ThreeVector gPos, 
	G4Event* anEvent, G4int evID, G4double zraster,
	G4double theta, G4double phi)
//////////////////////////////////////////
{
  //G4cout << gMom.x()/CLHEP::GeV << gMom.y()/CLHEP::GeV << gMom.z()/CLHEP::GeV <<G4endl;
  //G4cout << "#############" << zraster << G4endl;
  G4int n_particle = 1;
  gunParticle = new G4ParticleGun(n_particle);
  gunParticle-> SetParticleEnergy(0.0*CLHEP::GeV); //reset
  gunParticle-> SetParticleDefinition(particle);
  gunParticle-> SetParticleMomentum(gMom);
  gunParticle-> SetParticlePosition(gPos);
  gunParticle-> GeneratePrimaryVertex(anEvent);

  
  if( anaMan ){
    anaMan->PrimaryGeneration( gPos, gMom, evID ,zraster, theta, phi);
    //    G4cout << "gPos.x = " <<gPos.x()/CLHEP::cm << G4endl; /////suzuki
    //    G4cout << "gPos.y = " <<gPos.y()/CLHEP::cm << G4endl;
    //    G4cout << "gPos.z = " <<gPos.z()/CLHEP::cm << G4endl;

  }
  delete gunParticle;

}

/////////////////////////////////////////////////////////////////////
double gauss2D( double *x, double *par )
/////////////////////////////////////////////////////////////////////
{
  double z1 = double( (x[0]-par[1])/par[2] );
  double z2 = double( (x[1]-par[3])/par[4] );
  return par[0]*exp(-0.5*(z1*z1+z2*z2));
}

/////////////////////////////////////////////////////////////////////
double my2Dfunc( double *x, double *par )
/////////////////////////////////////////////////////////////////////
{
  return gauss2D(x, &par[0]) + par[5]*x[1]*x[1] + par[6]*x[1] + par[7];
}

