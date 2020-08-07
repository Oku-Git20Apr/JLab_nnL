#include "main.h"

#define RP_LINE_DRAW 0
#define CANV2 0

#define pi 3.1415

double alpha = 1.0/137.0;

int acceptance( TObject** h_acc, TFile* saveFile, _Info I )
int acceptance()
{

  string d_name = "test";
  TChain *tree = new TChain("tree");

  tree->Add( I.readFile );

  const Int_t nEntry = tree->GetEntries();
  const Double_t SolidAngle = 61.1475; // mrad   5degree 23.9mrad   8degree  61.1475mrad  4.58deg  20.09mrad

  int N_hist = 10;

  _HRS_data D = {0};
  
  tree->SetBranchAddress("EMom", &D.Mom);  
  tree->SetBranchAddress("ETheta", &D.Theta);  
  tree->SetBranchAddress("EPhi", &D.Phi);  
  tree->SetBranchAddress("EXpt", &D.Xpt);
  tree->SetBranchAddress("EYpt", &D.Ypt);  
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
  tree->SetBranchAddress("EXFP", &D.XFP);
  tree->SetBranchAddress("EYFP", &D.YFP);
  tree->SetBranchAddress("EXpFP", &D.XpFP);
  tree->SetBranchAddress("EYpFP", &D.YpFP);
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

  
  const  int Nbin_mom = 50;

  const  int Nbin_theta = 30;

  double mom_min;
  double mom_max;

  cout << "Data Flag = " << I.dataflag  << endl;
  if( I.dataflag == 0 ){ //K40
    mom_min = 2.70;
    mom_max = 3.30;
  }
  else if( I.dataflag == 10){ // nnl Tkine
    mom_min = 1.95;
    mom_max = 2.45;
  }
  else if( I.dataflag == 11){ // nnl Hkine
    mom_min = 1.82;
    mom_max = 2.32;
  }
  else if( I.dataflag == 12){ // nnl kaon
    mom_min = 1.55;
    mom_max = 2.05;
  }


  TH1D *All_mom = new TH1D("All_mom", "All_mom", Nbin_mom, mom_min, mom_max);
  TH1D *All_theta = new TH1D("All_theta", "All_theta", Nbin_theta, 0.4, 0.1);
  TH2D *All = new TH2D("All", "All", 50, mom_min, mom_max, 30, 0.4, 0.1);
  TH1D *Acc_mom = new TH1D("Acc_mom", "Acc_mom", Nbin_mom, mom_min, mom_max);
  TH1D *Acc_theta = new TH1D("Acc_theta", "Acc_theta", Nbin_theta, 0.4, 0.1);
  TH2D *Acc = new TH2D("Acc", "Acc", 50, mom_min, mom_max, 30, 0.4, 0.1);

  /************************************
   ******** Loop Start  ***************
   ***********************************/
  for(int i = 0; i < nEntry; i++){
    tree->GetEntry(i);
    double theta_eK = atan(D.Xpt);

    if(
       judgement( D, I )
	)
      {
	Acc_mom->Fill(D.Mom);
	Acc_theta->Fill(theta_eK);
	Acc->Fill(D.Mom, theta_eK);
      }
    
    All_mom->Fill(D.Mom);
    All_theta->Fill(theta_eK);
    All->Fill(D.Mom, theta_eK);

    if( i % 100000 == 0){
      cout << "[ " << i << " ] Filled" << endl;
    }
  }
  
  double err_mom[Nbin_mom] = {0};
  for(int i = 0; i < Nbin_mom; i++){
    double err_Acc_mom = Acc_mom->GetBinContent(i);
    double err_All_mom = All_mom->GetBinContent(i);
    if(err_All_mom != 0 && err_Acc_mom != 0){
      err_mom[i] = sqrt( pow(1.0/sqrt(err_Acc_mom), 2) + pow(1.0/sqrt(err_All_mom), 2) );
    } else{
      err_mom[i] = 0;
    }
  }

  Acc_mom->Divide(Acc_mom, All_mom, SolidAngle);
  for(int i = 0; i < Nbin_mom; i++){
    err_mom[i] = err_mom[i] * Acc_mom->GetBinContent(i);
    Acc_mom->SetBinError(i, err_mom[i]);
  }
  Acc_mom->SetMarkerStyle(20);
  Acc_mom->SetMarkerSize(0.4);

  double err_theta[Nbin_theta] = {0};
  for(int i = 0; i < Nbin_theta; i++){
    double err_Acc_theta = Acc_theta->GetBinContent(i);
    double err_All_theta = All_theta->GetBinContent(i);
    if(err_All_theta != 0 && err_Acc_theta != 0){
      err_theta[i] = sqrt( pow(1.0/sqrt(err_Acc_theta), 2) + pow(1.0/sqrt(err_All_theta), 2) );
    } else{
      err_theta[i] = 0;
    }
  }

  Acc_theta->Divide(Acc_theta, All_theta, SolidAngle);
  for(int i = 0; i < Nbin_theta; i++){
    err_theta[i] = err_theta[i] * Acc_theta->GetBinContent(i);
    Acc_theta->SetBinError(i, err_theta[i]);
  }

  Acc_theta->SetMarkerStyle(20);
  Acc_theta->SetMarkerSize(0.4);
  
  Acc->Divide(Acc, All, SolidAngle);

  string title = "Q1:" + to_string(I.Q1) + " Q2:" + to_string(I.Q2) + " Q3:" + to_string(I.Q3);
  Acc_mom->SetTitle( title.c_str() );
  Acc_theta->SetTitle( title.c_str() );
  Acc->SetTitle( title.c_str() );

  h_acc[0] = Acc_mom;
  h_acc[1] = Acc_theta;
  h_acc[2] = Acc;

  saveFile->cd();
  h_acc[0]->Write( ("histAccMom" + to_string(I.i)).c_str() );
  h_acc[1]->Write( ("histAccTheta" + to_string(I.i)).c_str() );
  h_acc[2]->Write( ("histAcc" + to_string(I.i)).c_str() );
  
  delete All_mom;
  delete All_theta;
  delete All;
  
  delete Acc_mom;
  delete Acc_theta;
  delete Acc;

  cout << "*******************************************" << endl;
  cout << "********** Acceptance calculated **********" << endl;
  cout << "*******************************************" << endl;
  
  return 0;

}
