#include "main.h"

using namespace std;

///////////////////////////////
//return momentum resolution //
//  made by suzuki           //
///////////////////////////////

#define Z_RASTER 0

#define MOM 1
#define XP 1
#define YP 1
#define ANG 1
#define Z 1

const int N_hist = 50;

ofstream file_Momfwhm("../fwhm/Momfwhm.dat");
ofstream file_Xpfwhm("../fwhm/Xpfwhm.dat");
ofstream file_Ypfwhm("../fwhm/Ypfwhm.dat");

const int MaxChar = 144;
const double SigToFWHM = 2.354820;

char name_Mmom0[200];
char name_Mxpt0[200];
char name_Mypt0[200];
char name_Mlen0[200];
char name_Mzt0[200];
char name_ifile[200];

void resolution( TObject** h_res, TFile* saveFile, _Info I )
{
  TRandom3* rnd = new TRandom3();
  
  strcpy(name_Mmom0 , I.name_read_Mmom);
  strcpy(name_Mxpt0 , I.name_read_Mxpt);
  strcpy(name_Mypt0 , I.name_read_Mypt);
  strcpy(name_Mlen0 , I.name_read_Mlen);
  strcpy(name_Mzt0 , I.name_read_Mzt);

  double **ParMom = new double*[5];
  double **ParXp = new double*[5];
  double **ParYp = new double*[5];
  double **ParZ = new double*[5];
  for(int i = 0; i < 5; i++){
    ParMom[i] = new double[I.nParam];
    ParXp[i] = new double[I.nParam];
    ParYp[i] = new double[I.nParam];
    ParZ[i] = new double[I.nParam];
  }


  TH1D *histMom_ori = new TH1D("histMom_ori", "histMom_ori", 500, -0.002, 0.002);
  TH1D *histXp_ori = new TH1D("histXp_ori", "histXp_ori", 500, -0.005, 0.005);
  TH1D *histYp_ori = new TH1D("histYp_ori", "histYp_ori", 500, -0.005, 0.005);
  TH1D *histAng_ori = new TH1D("histAng_ori", "histAng_ori", 500, -0.005, 0.005);
  TH1D *histZ_ori = new TH1D("histZ_ori", "histZ_ori", 500, -15.0, 15.0);

  char tmp_c[144] = "";
  char str[144] = "";
  
  ifstream fMom( name_Mmom0 );
  ifstream fXp( name_Mxpt0);
  ifstream fYp( name_Mypt0);
  ifstream fZ( name_Mzt0 );
  for(int i = 0; i < I.nParam; i++){
    fMom >> ParMom[0][i] >> ParMom[1][i] >> ParMom[2][i] >> ParMom[3][i] >> ParMom[4][i];
    fXp >> ParXp[0][i] >> ParXp[1][i] >> ParXp[2][i] >> ParXp[3][i] >> ParXp[4][i];
    fYp >> ParYp[0][i] >> ParYp[1][i] >> ParYp[2][i] >> ParYp[3][i] >> ParYp[4][i];
    fZ >> ParZ[0][i] >> ParZ[1][i] >> ParZ[2][i] >> ParZ[3][i] >> ParZ[4][i];
  }

  TChain *tree = new TChain("tree");
  
#if 0  //Read rootfile from paramfile
  tree->Add(name_ifile);
#endif
#if 1 //Read rootfile here
  tree->Add( I.readFile );
#endif

  cout << " Reading : " << I.readFile << endl;
    
  const int nEntry = tree->GetEntries();

  cout << "nEntry = " << nEntry << endl;
    
  _HRS_data D = {0};
  
  tree->SetBranchAddress("EMom", &D.Mom);
  tree->SetBranchAddress("EXpt", &D.Xpt);
  tree->SetBranchAddress("EYpt", &D.Ypt);
  tree->SetBranchAddress("EXFP", &D.XFP);
  tree->SetBranchAddress("EYFP", &D.YFP);
  tree->SetBranchAddress("EXpFP", &D.XpFP);
  tree->SetBranchAddress("EYpFP", &D.YpFP);
  tree->SetBranchAddress("EXSepi", &D.XSepi);
  tree->SetBranchAddress("EYSepi", &D.YSepi);
  tree->SetBranchAddress("EXSepe", &D.XSepe);
  tree->SetBranchAddress("EYSepe", &D.YSepe);
  tree->SetBranchAddress("EXQ1i", &D.XQ1i);
  tree->SetBranchAddress("EYQ1i", &D.YQ1i);
  tree->SetBranchAddress("EXQ2i", &D.XQ2i);
  tree->SetBranchAddress("EYQ2i", &D.YQ2i);
  tree->SetBranchAddress("EXQ3i", &D.XQ3i);
  tree->SetBranchAddress("EYQ3i", &D.YQ3i);
  tree->SetBranchAddress("ESPLTrig", &D.SPLTrig);
  tree->SetBranchAddress("EVDTrig", &D.EVDTrig);
  tree->SetBranchAddress("EDCTrig", &D.EDCTrig);
  //  tree->SetBranchAddress("EDetTrig", &D.EDetTrig);
  tree->SetBranchAddress("EZt", &D.zt0);
  if( I.pid != -1 ){
    tree->SetBranchAddress("EtrackIDFP1", &D.trackIDFP1);
    tree->SetBranchAddress("EchargeFP1", &D.chargeFP1);
    tree->SetBranchAddress("EparticleIDFP1", &D.particleID1);
  }

  
  
  for(int i = 1; i < nEntry; i++){ //default
    tree->GetEntry(i);
    
    double Mom_matrix = 0;
    double Mom_true = D.Mom;
    double Xp_matrix = 0;
    double Xp_true = D.Xpt;
    double Yp_matrix = 0;
    double Yp_true = D.Ypt;
    double Z_matrix = 0;
    double Z_true = D.zt0;

    double theta_true = getTheta( D.Mom, D.Xpt, D.Ypt );
    double theta_matrix = 0;
    
    if(i%100000 == 0){
      cout << "[ " << i << " ] Filled" << endl;
    }
    
    // ------ constraint to fill ----------------- //
    if( judgement( D, I ) ){
      
      if(I.resoflag){
	D.XFP = rnd->Gaus(D.XFP, I.xres); //cm //default 0.010
	D.YFP = rnd->Gaus(D.YFP, I.yres); //cm //default 0.010
	D.XpFP = rnd->Gaus(D.XpFP, I.xpres); //rad //default 5.0e-4
	D.YpFP = rnd->Gaus(D.YpFP, I.ypres); //rad //default 5.0e-4
      }
      
      double XFP_sca = (D.XFP - XFPm)/XFPr;
      double XpFP_sca = (D.XpFP - XpFPm)/XpFPr;
      double YFP_sca = (D.YFP - YFPm)/YFPr;
      double YpFP_sca = (D.YpFP - YpFPm)/YpFPr;
      
      for(int j = 0; j < I.nParam; j++){
	Mom_matrix  = Mom_matrix 
	  + ParMom[0][j] 
	  * pow(XFP_sca, ParMom[1][j])
	  * pow(XpFP_sca, ParMom[2][j])
	  * pow(YFP_sca, ParMom[3][j])
	  * pow(YpFP_sca, ParMom[4][j]);  
	
	Xp_matrix  = Xp_matrix 
	  + ParXp[0][j] 
	  * pow(XFP_sca, ParXp[1][j])
	  * pow(XpFP_sca, ParXp[2][j])
	  * pow(YFP_sca, ParXp[3][j])
	  * pow(YpFP_sca, ParXp[4][j]);  
	
	Yp_matrix  = Yp_matrix 
	  + ParYp[0][j] 
	  * pow(XFP_sca, ParYp[1][j])
	  * pow(XpFP_sca, ParYp[2][j])
	  * pow(YFP_sca, ParYp[3][j])
	  * pow(YpFP_sca, ParYp[4][j]);  

	Z_matrix  = Z_matrix 
	  + ParZ[0][j] 
	  * pow(XFP_sca, ParZ[1][j])
	  * pow(XpFP_sca, ParZ[2][j])
	  * pow(YFP_sca, ParZ[3][j])
	  * pow(YpFP_sca, ParZ[4][j]);  

      }
      
      Mom_matrix = Mom_matrix * Momr + Momm;
      Xp_matrix = Xp_matrix * Xptr + Xptm;
      Yp_matrix = Yp_matrix * Yptr + Yptm;
      Z_matrix = Z_matrix * ztr + ztm;

      theta_matrix = getTheta( Mom_matrix, Xp_matrix, Yp_matrix );
      
#if MOM
      histMom_ori->Fill( (Mom_true - Mom_matrix)/Mom_true );

#endif
#if XP	    
      histXp_ori->Fill( (Xp_true - Xp_matrix) );
#endif
#if YP	    
      histYp_ori->Fill( (Yp_true - Yp_matrix) );
#endif
#if ANG
      histAng_ori->Fill( theta_true - theta_matrix );
#endif
#if Z
      histZ_ori->Fill( Z_true - Z_matrix );
#endif
      
    }
  }
  
  string title = "Q1:" + to_string(I.Q1) + " Q2:" + to_string(I.Q2) + " Q3:" + to_string(I.Q3);

  
  ////////////////////////////////
  //// save histgram  ////////////
  ////////////////////////////////
  saveFile->cd();
  
#if MOM
  histMom_ori->GetXaxis()->SetMaxDigits(2);
  histMom_ori->GetXaxis()->SetNdivisions(515);
  histMom_ori->SetTitle( title.c_str() );
  h_res[0] = histMom_ori;
  h_res[0]->Write( ("histResMom" + to_string(I.i)).c_str() );
  delete histMom_ori;
#endif
  
#if XP
  histXp_ori->GetXaxis()->SetMaxDigits(2);
  histXp_ori->GetXaxis()->SetNdivisions(515);
  histXp_ori->SetTitle( title.c_str() );
  h_res[1] = histXp_ori;
  h_res[1]->Write( ("histResXp" + to_string(I.i)).c_str() );
  delete histXp_ori;
#endif
  
#if YP
  histYp_ori->GetXaxis()->SetMaxDigits(2);
  histYp_ori->GetXaxis()->SetNdivisions(515);
  histYp_ori->SetTitle( title.c_str() );    
  h_res[2] = histYp_ori;
  h_res[2]->Write( ("histResYp" + to_string(I.i)).c_str() );
  delete histYp_ori;
#endif
  
#if ANG
  histAng_ori->GetXaxis()->SetMaxDigits(2);
  histAng_ori->GetXaxis()->SetNdivisions(515);
  histAng_ori->SetTitle( title.c_str() );    
  h_res[3] = histAng_ori;
  h_res[3]->Write( ("histResAng" + to_string(I.i)).c_str() );
  delete histAng_ori;
#endif

#if Z
  histZ_ori->GetXaxis()->SetMaxDigits(2);
  histZ_ori->GetXaxis()->SetNdivisions(515);
  histZ_ori->SetTitle( title.c_str() );    
  h_res[4] = histZ_ori;
  h_res[4]->Write( ("histResZ" + to_string(I.i)).c_str() );
  delete histZ_ori;
#endif
  
  delete tree;

  for(int i = 0; i < 5; i++){
    delete[] ParMom[i];
    delete[] ParXp[i];
    delete[] ParYp[i];
    delete[] ParZ[i];
  }


  cout << "********************************" << endl;
  cout << "**** Resolution Calculated *****" << endl;
  cout << "********************************" << endl;


}
