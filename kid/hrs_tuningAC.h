
#ifndef hrs_tuningAC_h
#define hrs_tuningAC_h 1

double s2f1_off(int i,char* ARM,char* MODE,int KINE);
const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "Setting.h"


int nth=3;//change, nth=0 originally
char* mode="H";
int kine=1;
double tdc_time=56.23;
bool ac2_min=true;
bool single_flag=false;

class tuningAC{

 public:
  tuningAC();
  ~tuningAC();
  void SetRunList(string ifname);
  void SetRun(string ifname);
  void SetRoot(string ofname);
  void SetBranch();
  void SetParam();  
  void MakeHist();
  void Fill();
  void Fitting();
  void Tuning(); 
  void Draw();
  void Print(string ofname);
  void Write();
  void Write_coin();  
  void Comment();
  TFile* fnew;
  public:
  Setting* set;



  
  //=== SetRunList ====//
  TChain* T;
  int ENum;

  //=== SetRoot =======//
  TTree* tnew;
  double mm_ac1[100][100];
  double mm_ac2[100][100];
  double fom_ac1[100];
  double fom_ac2[100];
  double mm_c;  
  double ct_c;
  //==== SetBranch ====//


   double RF1[100],LF1[100];
  double Rs0r_ac[100],Rs0l_ac[100],Ls0r_ac[100],Ls0l_ac[100];
  double Rs2r_ac[100],Rs2l_ac[100],Ls2r_ac[100],Ls2l_ac[100];
  double Rs0r_tc[100],Rs0l_tc[100],Ls0r_tc[100],Ls0l_tc[100];
  double Rs2r_tc[100],Rs2l_tc[100],Ls2r_tc[100],Ls2l_tc[100];
  double Ra1t[100],Ra1a[100],Ra1a_p[100],Ra1a_c[100],Ra1sum;
  double Ra2t[100],Ra2a[100],Ra2a_p[100],Ra2a_c[100],Ra2sum;
  double La1t[100],La1a[100],La1a_p[100],La1a_c[100],La1sum;
  double La2t[100],La2a[100],La2a_p[100],La2a_c[100],La2sum;
  double Rp[100],Rpx[100],Rpy[100],Lp[100],Lpx[100],Lpy[100];
  double Rth[100],Rph[100],Rx[100],Rvz[100],Lth[100],Lph[100],Lx[100],Lvz[100];
  double Rbeta[100],Lbeta[100];
  double rs2pathl[100],rs0pathl[100],rtrpathl[100];
  double ls2pathl[100],ls0pathl[100],ltrpathl[100];
  double trigger[100];
  double hallap;
  double Rs2tpads[100],Ls2tpads[100];
  double Rs2trpad[100],Ls2trpad[100];
  double Ru1_time[100];
  int NRu1_time;
  //---- Gogami root ---------//
  double ctime[100];
  double DRT5;
  //---- Toyama ana_Lambda ----//
  double Rz,Lz,Rpz,Lpz; 
  double tcoin_t;
  //------------------------//
 double mm; 
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rtof[16];
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off; 
 double Rs2_tcorr,Ls2_tcorr;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double tof_r,tof_l; 
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr;
 double ct_acc; 
 double acc;
 double ct;
 int ev;

 //===== SetParam ======//


 double ac1_adc[100],ac2_adc[100];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;
 double th1_max,th2_max,th2_min;
 double ac1_kcut,ac2_kcut_min,ac2_kcut_max;
 double th_ac2_t,th_ac2_b; 

 
 //===== Make Hist =======//
 TH1F* hmm;
 TH1F* hmm_acc;
 TH1F* hmm_p;
 TH1F* hRu1_time_c;
 TH1F* hRu1_time_s; 
 TH2F* hcoin_ac1[100];
 TH2F* hcoin_ac2[100];
 TH2F* hcoin_ac1_acc[100];
 TH2F* hcoin_ac2_acc[100];
 TH2F* hvdc1_ac1[100];
 TH2F* hvdc2_ac1[100];
 TH2F* hvdc1_ac2[100];
 TH2F* hvdc2_ac2[100];
 TH2F* hs0_ac1[100];
 TH2F* hs2_ac1[100];
 TH2F* hs0_ac2[100];
 TH2F* hs2_ac2[100];
 TH1F* hcoin_t1[100];
 TH1F* hcoin_t2[100];
 TH1F* hcoin_ac1_max[100];
 TH1F* hcoin_ac2_max[100];
 TH1F* hcoin_t3[100][100];
 TH1F* hcoin_t;
 TH1F* hcoin_tc;
 TH1F* hcoin_acc_ac1[100];
 TH1F* hcoin_acc_ac2[100];
 TH2F* ha1_a2;
 TH1F* hcoin_k;
 TH1F* hcoin_pi;
 TH1F* hcoin_p;
 TH2F* hcoin_ac1_all;
 TH2F* hcoin_ac2_all;
 TH2F* hmm_ac1[100];
 TH2F* hmm_ac2[100];
 TH2F* hmm_ac1_acc[100];
 TH2F* hmm_ac2_acc[100];
 TH2F* hct_a1a_c[24];
 TH2F* hct_a2a_c[26]; 
 //----- Fill -----// 
 TH1D* hmm_ac1_p[100][100];
 TH1D* hmm_ac2_p[100][100];
 TF1* facc[100][100][2];
 TF1* fpi[100][100][2];
 TF1* fk[100][100][2];
 TF1* fcoin[100][100][2];
 TF1* fp[100][100][2];
 TF1* fbg[100][100][2];
 TF1* fbg_s[100][100][2];
 TF1* fLam[100][100][2];
 TF1* fSig[100][100][2];
 TF1* fLam_p;
 TF1* fSig_p;   
 TH2F* hfom_ac[100][100];
 TH2F* hAC;
 TH2F* hAC2;
 
 int iter_ac1;
 int iter_ac2;
 int iter_max;
 TH1D* hcoin_ac1_p[100][100];
 TH1D* hcoin_ac2_p[100][100]; 
 TH1D* hcoin_ac1_all_p[100][100];
 TH1D* hcoin_ac2_all_p[100][100]; 
 TH1D* hcoin_ac1_acc_p[100][100];
 TH1D* hcoin_ac2_acc_p[100][100]; 
 TH1D* hmm_ac1_all_p[100][100];
 TH1D* hmm_ac2_all_p[100][100];
 TH1D* hmm_ac1_acc_p[100][100];
 TH1D* hmm_ac2_acc_p[100][100];

 TGraphErrors* gsum_pi_ac1[100][100];
 TGraphErrors* gsum_p_ac1[100][100];
 TGraphErrors* gsum_k_ac1[100][100];
 TGraphErrors* grate_k_ac1[100][100];
 TGraphErrors* grate_p_ac1[100][100];
 TGraphErrors* grate_pi_ac1[100][100];
 TGraphErrors* gsum_pi_ac2[100][100];
 TGraphErrors* gsum_p_ac2[100][100];
 TGraphErrors* gsum_k_ac2[100][100];
 TGraphErrors* grate_k_ac2[100][100];
 TGraphErrors* grate_pi_ac2[100][100];
 TGraphErrors* grate_p_ac2[100][100];
 TGraphErrors* gSN_k_ac1[100][100];
 TGraphErrors* gSN_k_ac2[100][100];
 TGraphErrors* gfom_ac1[100];
 TGraphErrors* gfom_ac2[100];
 TGraphErrors* gfom;
 TGraphErrors* gmm_SN_ac1[100];
 TGraphErrors* gmm_SN_ac2[100];
 TGraphErrors* gmm_S_ac1[100];
 TGraphErrors* gmm_S_ac2[100];
 TGraphErrors* gmm_ac2[100]; 
 TGraphErrors* gmm_ac1[100]; 
 TGraphErrors* gL_ac1[100];
 TGraphErrors* gL_ac2[100];  
 TGraphErrors* gS_ac1[100];
 TGraphErrors* gS_ac2[100];  
 TGraphErrors* gL_eff_ac1[100];
 TGraphErrors* gL_eff_ac2[100];  
 TGraphErrors* gS_eff_ac1[100];
 TGraphErrors* gS_eff_ac2[100];  
 TGraphErrors* gS_SN_ac1[100];
 TGraphErrors* gS_SN_ac2[100];  
 TGraphErrors* gL_SN_ac1[100];
 TGraphErrors* gL_SN_ac2[100];  
 TGraphErrors* gL_N_ac1[100];
 TGraphErrors* gL_N_ac2[100];
 TGraphErrors* gL_FOM_ac1[100];    
 TGraphErrors* gL_FOM_ac2[100];  


 TF1* facc_t1def[100][100];
 TF1* fpi_t1def[100][100];
 TF1* fk_t1def[100][100];
 TF1* fcoin_t1def[100][100];
 TF1* fp_t1def[100][100];
 TF1* facc_t2def[100][100];
 TF1* fpi_t2def[100][100];
 TF1* fk_t2def[100][100];
 TF1* fcoin_t2def[100][100];
 TF1* fp_t2def[100][100];
 TF1* facc_t3def[100][100];
 TF1* fpi_t3def[100][100];
 TF1* fk_t3def[100][100];
 TF1* fcoin_t3def[100][100];
 TF1* fp_t3def[100][100];
 TF1* fcoin_t1[100][100];
 TF1* fcoin_t2[100][100];  
 TF1* fcoin_t3[100][100]; 
 TF1* facc_kc;
 TF1* fk_kc;
 TF1* fpi_pic;
 TF1* fp_pc;

 //----- Tuning hist ----//
 TH1F* hcoin_fom;
 TH1F* hcoin_acc;
 TH1F* hmm_fom;
 TH1F* hmm_fom_acc;
 TH1F* hmm_fom_p;
 TF1* fbg_L;
 TF1* fbg_S; 
 TF1*fL_p;    
 TF1* fS_p;  
 TF1* fL_all;
 TF1* fS_all;  
 TF1* fS_fom;
 TF1* fL_fom;
 TF1* fL_fom_bg;
 TF1* fS_fom_bg;
 TF1* fL_fom_p;
 TF1* fS_fom_p;
 TH1F* hcoin_fom_p;
 TF1* fk_fom; 
 //----- paremters ----//
 double bin_vdc,min_vdc,max_vdc;
 double min_s0,max_s0;
 int bin_s0;
 double min_s2,max_s2;
 int bin_s2;
 double bin_coin;
 double bin_coin_c;
 int bin_beta; 
 int bin_adc;
 int bin_ac1;
 int bin_ac2;
 double min_mm,max_mm,bin_mm;

 //===== Fill ========//
 
 //--- Coin Offset -----//
 double pathl_off,s2_offset,coin_offset;
 //----- Cut Parameters ----------//
 double coin_cutmin=-248;
 double coin_cutmax=-244; 
 double rpathl_cutmin=28.7;
 double rpathl_cutmax=29.4;
 double lpathl_cutmin=28.6;
 double lpathl_cutmax=29.2;
 double rbeta_cutmin=0.0;
 double rbeta_cutmax=1.0;
 double lbeta_cutmin=0.9;
 double lbeta_cutmax=1.0;
 double Rvz_cutmin=-0.1;
 double Rvz_cutmax= 0.1;
 double Lvz_cutmin=-0.1;
 double Lvz_cutmax= 0.1;
 double Rx_cutmin= -0.4;
 double Rx_cutmax= 0.4;
 //-------------------------------//
 bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig,cut_track,cut_s0;


 //===== Fitting =========//


 //--- Parameters -----//
  double bg_min,bg_max;
  double bgs_min,bgs_max;
  double Lfom[100][100][2],Sfom[100][100][2];
  double L0_err[100][100][2],L1_err[100][100][2],L2_err[100][100][2];
  double S0_err[100][100][2],S1_err[100][100][2],S2_err[100][100][2];
  double nL_err[100][100][2],nS_err[100][100][2];
  double bgL_ac1[100][100], bgL_ac2[100][100],bgS_ac1[100][100], bgS_ac2[100][100];
  double totL_ac1[100][100], totL_ac2[100][100],totS_ac1[100][100], totS_ac2[100][100];
  double nL[100][100][2],sigL[100][100][2],meanL[100][100][2];
 double nS[100][100][2],sigS[100][100][2],meanS[100][100][2];
 double kmin[100][100][2],kmax[100][100][2];
 double inte_ktot[100][100][2], inte_ksig[100][100][2];
 double p0_acc[100][100][2], p1_acc[100][100][2];
 double n_p[100][100][2],sig_p[100][100][2],mean_p[100][100][2];
 double n_pi[100][100][2],sig_pi[100][100][2],mean_pi[100][100][2];
 double n_k[100][100][2],sig_k[100][100][2],mean_k[100][100][2];
 int bin_ac1_adc[100][100],bin_min_ac1,bin_max_ac1,bin_ac2_adc[100][100],bin_max_ac2,bin_min_ac2;
 double sum_k[100][100][2],sum_p[100][100][2],sum_pi[100][100][2]; 
 double sum_k_err[100][100][2],sum_p_err[100][100][2],sum_pi_err[100][100][2]; 
 double inte_acc[100][100][2];
 double th_ac1[100],th_ac2[100];
 int bin_th_ac1[100][100],bin_th_ac2[100][100]; 
 double nk[100][100][100][100][2],npi[100][100][100][100][2],np[100][100][100][100][2];
 double max_nk[100][100][2],max_npi[100][100][2],max_np[100][100][2];
 double n_p_err[100][100][2],n_pi_err[100][100][2],n_k_err[100][100][2];
 double FOM_ac1[100][100],FOM_ac2[100][100];
 double max_fom_ac1,max_fom_ac2;
 int fom_th1,fom_th2;
 double nLam_ac1,nLam_ac2,SNLam_ac1,SNLam_ac2;
 int fom_max_th2,fom_max_th1;
 double FOM_max_ac1[100],FOM_max_ac2[100],FOM_th1[100],FOM_th2[100]; 
 
 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_sig_k,def_mean_k,def_acc;
 double def_num_k,def_num_p,def_num_pi,def_acc_k,def_acc_pi,def_acc_p;

 double def_t1_k[100][100],def_t1_pi[100][100],def_t1_p[100][100],def_t1_acc[100][100];
 double def_t1_k_err[100][100],def_t1_pi_err[100][100],def_t1_p_err[100][100],def_t1_acc_err[100][100];
 double t1sig_k[100][100],t1sig_p[100][100],t1sig_pi[100][100],t1mean_p[100][100],t1mean_k[100][100],t1mean_pi[100][100];
 double t1sum_k[100],t1sum_pi[100],t1sum_p[100];
 double t1sum_k_err[100],t1sum_pi_err[100],t1sum_p_err[100];
 double def_t2_k[100][100],def_t2_pi[100][100],def_t2_p[100][100],def_t2_acc[100][100];
 double t2sig_k[100][100],t2sig_p[100][100],t2sig_pi[100][100],t2mean_p[100][100],t2mean_k[100][100],t2mean_pi[100][100];
 double def_t2_k_err[100][100],def_t2_pi_err[100][100],def_t2_p_err[100][100],def_t2_acc_err[100][100];
 double t2sum_k[100],t2sum_pi[100],t2sum_p[100];
 double t2sum_k_err[100],t2sum_pi_err[100],t2sum_p_err[100];
 double def_t3_k[100][100],def_t3_pi[100][100],def_t3_p[100][100],def_t3_acc[100][100];
 double t3sig_k[100][100],t3sig_p[100][100],t3sig_pi[100][100],t3mean_p[100][100],t3mean_k[100][100],t3mean_pi[100][100];
 double t3sum_k[100][100],t3sum_pi[100][100],t3sum_p[100][100];
 double emp[100];

 double rate_k[100][100][2],rate_p[100][100][2],rate_pi[100][100][2];
 double rate_k_err[100][100][2],rate_p_err[100][100][2],rate_pi_err[100][100][2];
 double sum_acc[100][100][2];
 double max_SN_ac1[100],max_SN_ac2[100];
 int SN_ac1[100],SN_ac2[100];
 double bg_0[100][100][2],bg_1[100][100][2],bg_2[100][100][2];
 double bg_s0[100][100][2],bg_s1[100][100][2],bg_s2[100][100][2];
 double L0[100][100][2],L1[100][100][2],L2[100][100][2];
 double S0[100][100][2],S1[100][100][2],S2[100][100][2];
 double sum_k_max=1250.;
 double fom_max=0.0;
 


 //====== Tuning ============//
 bool ac2_up,ac2_down,ac2_flag; 
 double Lbg_fom[3],Sbg_fom[3];
 double Lam_p[3],Sig_p[3];
 double NL_err,NS_err;
 double Lam_p_err[3],Sig_p_err[3];

  double pbg[3];  
  double pbg_S[3];
  double pL[3],pL_err[3];
  double pS[3],pS_err[3];
  double sum_L,sum_S;

 double all_L;
 double all_S;
 double bg_L;
 double bg_S;
 double sumk_fom;
 double meank_fom;
 double sigk_fom;
 double sk_fom;
 double sk_fom_ct;
 double nk_fom;
 double snk_fom;
 double fom;

 //====== Draw ===========//
 TCanvas* c9;
 TCanvas* c10;
 TCanvas* c11;
 TCanvas* c12;
 TCanvas* c13; 
 TCanvas* c14;
 TCanvas* c15;
 TCanvas* c16;
 TCanvas* c17;
 TCanvas* c18;
 TCanvas* c19;
 TCanvas* c20;
 TCanvas* c21;
 TCanvas* c22;
 TCanvas* c23;
 TCanvas* c24;
 TCanvas* c25; 
 TCanvas* c30;
 TCanvas* c31;
 TCanvas* c32;
 TCanvas* c33;
 TCanvas* c34;


   
 
};
///////////////////////////////////////////////////////////////
tuningAC::tuningAC(){
  set= new Setting();
  set->Initialize();

  cout<<"====================================="<<endl;
  cout<<"===== ANALYSIS Infomation ==========="<<endl;
  cout<<"====================================="<<endl;
  cout<<"MODE : "<<mode<<endl;
  cout<<"F1TDC type : "<<kine<<endl;
  cout<<"AC2 Min mode: "<<ac2_min<<endl;

}

tuningAC::~tuningAC(){}

//////////////////////////////////////////////////////////////

void tuningAC::SetRun(string ifname){


  if(mode=="G"||mode=="T" ){T=new TChain("tree"); }
  else {T=new TChain("T");}
  T->Add(ifname.c_str());
  ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 
}

///////////////////////////////////////////////////////////////
void tuningAC::SetRunList(string ifname){
  cout<<"Set Run List "<<endl;
  
  if(mode=="G"||mode=="T" ){T=new TChain("tree"); }
  else {T=new TChain("T");}

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //   cout<<buf<<endl;
  }
  ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 
}
///////////////////////////////////////////////////////////////////////////

void tuningAC::SetRoot(string ofname){

  fnew = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew =new TTree("T",ofname.c_str());
  tnew = T->CloneTree(0);
    
  tnew->Branch("mm_c",&mm_c,"mm_c/D");
  tnew->Branch("ct_c",&ct_c,"ct_c/D");  
  /*

  tnew->Branch("mm_ac1",mm_ac1,"mm_ac1[100][100]/D");
  tnew->Branch("mm_ac2",mm_ac2,"mm_ac2[100][100]/D");
  tnew->Branch("fom_ac1",fom_ac1,"fom_ac1[100]/D");
  tnew->Branch("fom_ac2",fom_ac2,"fom_ac2[100]/D");  
  tnew->Branch("th_ac1",th_ac1,"th_ac1[100]/D");
  tnew->Branch("th_ac2",th_ac2,"th_ac2[100]/D");  
  */

  
}

////////////////////////////////////////////////////////////////////////////
void tuningAC::SetBranch(){

 T->SetBranchStatus("*",0);  
 
//------ Right Arm -------------//

 T->SetBranchStatus("RTDC.F1FirstHit",1);
 T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
 T->SetBranchStatus("R.s2.t_pads",1);
 T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
 T->SetBranchStatus("R.s2.trpad",1);
 T->SetBranchAddress("R.s2.trpad",Rs2trpad);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c); 
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
// T->SetBranchStatus("R.a1.asum_p",1);
// T->SetBranchAddress("R.a1.asum_p",&Ra1sum);
// T->SetBranchStatus("R.a2.asum_p",1);
// T->SetBranchAddress("R.a2.asum_p",&Ra2sum);
//change 
 T->SetBranchStatus("R.vdc.u1.time",1);
 T->SetBranchAddress("R.vdc.u1.time",Ru1_time);
 T->SetBranchStatus("Ndata.R.vdc.u1.time",1);
 T->SetBranchAddress("Ndata.R.vdc.u1.time",&NRu1_time);
 // path length//
 T->SetBranchStatus("R.s2.trpath",1); 
 T->SetBranchAddress("R.s2.trpath",rs2pathl); 
 T->SetBranchStatus("R.tr.pathl",1);  
 T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // Target positon information //
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm a---------------//
 T->SetBranchStatus("LTDC.F1FirstHit",1);
 T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
 T->SetBranchStatus("L.s2.t_pads",1);
 T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
 T->SetBranchStatus("L.s2.trpad",1);
 T->SetBranchAddress("L.s2.trpad",Ls2trpad);
  // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);

 if(mode=="T"){
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&tcoin_t);
 T->SetBranchStatus("ac1_sum",1);
 T->SetBranchAddress("ac1_sum",&Ra1sum);
 T->SetBranchStatus("ac2_sum",1);
 T->SetBranchAddress("ac2_sum",&Ra2sum);
 T->SetBranchStatus("mm",1);
 T->SetBranchAddress("mm",&mm);
 T->SetBranchStatus("Lz",1);    
 T->SetBranchAddress("Lz",&Lz);
 T->SetBranchStatus("Rz",1);    
 T->SetBranchAddress("Rz",&Rz);
 T->SetBranchStatus("Rp",1);
 T->SetBranchAddress("Rp",&Rpz);
 T->SetBranchStatus("Lp",1);
 T->SetBranchAddress("Lp",&Lpz);
 T->SetBranchStatus("ct_acc",1);
 T->SetBranchAddress("ct_acc",&ct_acc);
 }



}


///////////////////////////////////////////////////////
void tuningAC::SetParam(){

  //Paremater Setting //
 // AC set ACP//
   if(mode=="H"){
   min_coin=-10;
   max_coin=20.0;
   min_coin_c=-10;
   max_coin_c=20.0;

 //  min_coin_c=-100;
 //  max_coin_c=1000.0;
 min_ac1=0.0;
 max_ac1=5000.;
 min_ac2=0.0;
 max_ac2=20000.;
 min_adc=-500.0;
 max_adc=20000.;
//=== AC Threshold variable ===//
 th1_max=500.;
 ac1_adc[0]=th1_max;
 ac1_adc[1]=300.;
 ac1_adc[2]=100.;


 ac2_adc[0]=min_ac2;
 ac2_adc[1]=1000.;
 ac2_adc[2]=3000.;
 th2_max=6000.;

 //---Kaon Cut ----//
 ac1_kcut=100.;
 ac2_kcut_min=1000.;
 ac2_kcut_max=5000.;



 //ACT //
   }else if(mode=="T"){
   min_coin=-20;
   max_coin=10.0;
   min_coin_c=-20;
   max_coin_c=10.0;

 //==== NPE =======//
 min_ac1=0.0;
 max_ac1=30.0;
 //min_ac2=0.0;
 min_ac2=0.0;
 // max_ac2=100.;
 max_ac2=100.;
 th1_max=3.0;
 //th1_max=max_ac1;
 ac1_adc[0]=3;
 //ac1_adc[1]=th1_max;
 //ac1_adc[2]=th1_max;
 ac1_adc[1]=5;
 // ac1_adc[2]=max_ac1;
 ac1_adc[2]=7;



 
 if(ac2_min){

   //==== ACparam ===//
   
   //===========================//
   //== AC2 Upper threshold ====//
   //===========================//
 ac1_adc[0]=0.48;
 ac1_adc[1]=1.0;
 ac1_adc[2]=1.5;

 ac2_adc[0]=1.0;
 ac2_adc[1]=2.0;
 ac2_adc[2]=3.0;

 ac2_kcut_max=12.0;
 //---- AC1 -----//
 min_ac1=0.0; // AC1 minmum
 th1_max=10.0; // AC1 maximum 

 min_ac1=0.0;
 max_ac1=10.0;
 min_ac2=0.0; 
 max_ac2=50.0;  
 //--- AC2 -----//

 // min_ac2=0.0;

  th2_max=10.;
  th_ac2_t=10.0; // FIll Cut of hmm_ac1 (Upper cut)
  //  th_ac2_t=max_ac2;
  //  th_ac2_b=min_ac2;
 /*
 for(int i=0;i<nth;i++){
   ac1_adc[i]=th1_max-(th1_max-min_ac1)/nth*i;
   ac2_adc[i]=min_ac2+(th2_max-min_ac2)/nth*i;
 }
 */
 
 

 
 }else{

   th2_max= 20.;
   max_ac2=20.;
   //th2_max=max_ac2;
 min_ac2=0.0;
 ac2_adc[2]=max_ac2;
 ac2_adc[1]=20.0;
 ac2_adc[0]=15.0;
 th2_min=0.0;
 min_ac2=0.0;
 th_ac2_b=4.2;
 ac1_adc[2]=1.2;

 }

 
 

 //---Kaon Cut ----//
 ac1_kcut=1.0;
 ac2_kcut_min=3.0;
 ac2_kcut_max=7.0;

 /*
 ac2_adc[0]=max_ac2;
 ac2_adc[1]=6000.;
 ac2_adc[2]=5000.;
 */
 //----------------//
   }else if(mode=="G"){
 min_coin=-20;
 max_coin=20.0;
 min_coin_c=-20;
 max_coin_c=20.0;
 min_ac1=0.0;
 max_ac1=30.;
 min_ac2=0.0;
 max_ac2=50.;
 min_adc=-5.0;
 max_adc=20.;
//==== AC Threshold variable (ACTH)===//
 ac1_adc[0]=max_ac1;
 ac1_adc[1]=1.3;
 ac1_adc[2]=1.0;
 ac2_adc[0]=min_ac2;
 ac2_adc[1]=5.;
 ac2_adc[2]=7;
}


}
////////////////////////////////////////////////////////////

void tuningAC::MakeHist(){


 min_vdc=-0.2e-6;
 max_vdc= 1.2e-6;
 bin_vdc=(max_vdc-min_vdc)/tdc_time*1.0e6;
 bin_vdc=(int)bin_vdc;
 min_s0=-10;
 max_s0=10000;
 bin_s0=int(max_s0-min_s0);


 min_s2=-10;
 max_s2=5000;
 bin_s2=max_s2-min_s2;
        bin_coin=(int)(max_coin-min_coin)/tdc_time;
        bin_coin_c=(int)(max_coin_c-min_coin_c)/tdc_time;;
        bin_beta=6000;
	bin_adc=(int)max_adc-min_adc;
	bin_ac1=(int)(max_ac1-min_ac1)*3; 
	bin_ac2=(int)(max_ac2-min_ac2)*3; 
  min_mm=0.5;
  max_mm=1.5;
  bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
  bin_mm=(int)bin_mm;


  // iter_ac1=40; //NIter
 iter_ac1=30; //Iter
 iter_max=iter_ac1;

  
  if(mode=="T"){
  bin_ac1=1000.; 
  bin_ac2=(max_ac2-min_ac2)*100; 
}


 hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",(int)bin_coin,min_coin,max_coin);
 hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);


 ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(ha1_a2,"","","");

 hcoin_k=new TH1F("hcoin_k","Coincidence time w/ Correction Kaon Cut  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 hcoin_pi=new TH1F("hcoin_pi","Coincidence time w/ Correction Pion  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 hcoin_p=new TH1F("hcoin_p","Coincidence time w/ Correction Proton  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 hcoin_ac1_all=new TH2F("hcoin_ac1_all","Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1_all,"","","");
 hcoin_ac2_all=new TH2F("hcoin_ac2_all","Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2_all,"","","");


 for(int i=0;i<nth;i++){
   hmm_ac1[i]=new TH2F(Form("hmm_ac1[%d]",i),"Mass AC1 Cut",bin_mm,min_mm,max_mm,bin_ac1,min_ac1,max_ac1);
   hmm_ac2[i]=new TH2F(Form("hmm_ac2[%d]",i),"Mass AC2 Cut",bin_mm,min_mm,max_mm,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hmm_ac1[i],"Mass AC1 Cut","Mass [GeV/c^2]","AC1 NPE");
 set->SetTH2(hmm_ac2[i],"Mass AC2 Cut","Mass [GeV/c^2]","AC2 NPE");
 hmm_ac1_acc[i]=new TH2F(Form("hmm_ac1_acc[%d]",i),"Mass AC1 ACC Cut",bin_mm,min_mm,max_mm,bin_ac1,min_ac1,max_ac1);
   hmm_ac2_acc[i]=new TH2F(Form("hmm_ac2_acc[%d]",i),"Mass AC2 ACC Cut",bin_mm,min_mm,max_mm,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hmm_ac1_acc[i],"Mass AC1 ACC Cut","Mass [GeV/c^2]","Counts /2 MeV");
 set->SetTH2(hmm_ac2_acc[i],"Mass AC2 ACC Cut","Mass [GeV/c^2]","Counts /2 MeV");


 hcoin_ac1[i]=new TH2F(Form("hcoin_ac1[%d]",i),"Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1[i],"","","");
 hcoin_ac2[i]=new TH2F(Form("hcoin_ac2[%d]",i),"Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2[i],"","","");
 hcoin_ac1_acc[i]=new TH2F(Form("hcoin_ac1_acc[%d]",i),"Coinc time (ACC) vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1_acc[i],"","","");
 hcoin_ac2_acc[i]=new TH2F(Form("hcoin_ac2_acc[%d]",i),"Coinc time (ACC) vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2_acc[i],"","","");
 hvdc1_ac1[i]=new TH2F(Form("hvdc1_ac1[%d]",i),"VDC1 Efficincy vs AC1 threshold",bin_vdc,min_vdc,max_vdc,bin_ac1,min_ac1,max_ac1);
 hvdc1_ac2[i]=new TH2F(Form("hvdc1_ac2[%d]",i),"VDC1 Efficincy vs AC2 threshold",bin_vdc,min_vdc,max_vdc,bin_ac2,min_ac2,max_ac2);
 hvdc2_ac1[i]=new TH2F(Form("hvdc2_ac1[%d]",i),"VDC2 Efficincy vs AC1 threshold",bin_vdc,min_vdc,max_vdc,bin_ac1,min_ac1,max_ac1);
 hvdc2_ac2[i]=new TH2F(Form("hvdc2_ac2[%d]",i),"VDC2 Efficincy vs AC2 threshold",bin_vdc,min_vdc,max_vdc,bin_ac2,min_ac2,max_ac2);
 hs0_ac1[i]=new TH2F(Form("hs0_ac1[%d]",i),"S0 Efficiency vs AC1 threshold",bin_s0,min_s0,max_s0,bin_ac1,min_ac1,max_ac1);
 hs0_ac2[i]=new TH2F(Form("hs0_ac2[%d]",i),"S0 Efficiency vs AC2 threshold",bin_s0,min_s0,max_s0,bin_ac2,min_ac2,max_ac2);
 hs2_ac1[i]=new TH2F(Form("hs2_ac1[%d]",i),"S2 Efficiency vs AC1 threshold",bin_s2,min_s2,max_s2,bin_ac1,min_ac1,max_ac1);
 hs2_ac2[i]=new TH2F(Form("hs2_ac2[%d]",i),"S2 Efficiency vs AC2 threshold",bin_s2,min_s2,max_s2,bin_ac2,min_ac2,max_ac2);
 hcoin_t1[i]=new TH1F(Form("hcoin_t1[%d]",i), Form("Coincidence 0<AC1<%lf cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
  set->SetTH1(hcoin_t1[i],"Coincidence time ","","");
 hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence %lf<AC2<%lf  cut",ac1_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
   set->SetTH1(hcoin_t2[i],"Coincidence time ","","");
 hcoin_ac1_max[i]=new TH1F(Form("hcoin_ac1_max[%d]",i), Form("Coincidence time %lf<AC2<%lf cut",ac2_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_ac1_max[i],"","","");
 hcoin_ac2_max[i]=new TH1F(Form("hcoin_ac2_max[%d]",i), Form("Coincidence time AC1<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_ac2_max[i],"","","");
 hcoin_acc_ac1[i]=new TH1F(Form("hcoin_acc_ac1[%d]",i),"Coincidence time ACC BG ",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_acc_ac1[i],"","","");
 hcoin_acc_ac2[i]=new TH1F(Form("hcoin_acc_ac2[%d]",i),"Coincidence time ACC BG ",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_acc_ac2[i],"","","");
 for(int j=0;j<nth;j++){
  hcoin_t3[i][j]=new TH1F(Form("hcoin_t3[%d][%d]",i,j),Form("Coincidence time S2R-S2L[ns] ac1_adc< %lf, ac2_adc<%lf; Coin-Time [ns];Counts "
							    ,ac1_adc[i],ac2_adc[j]),bin_coin_c,min_coin_c,max_coin_c);
  set->SetTH1(hcoin_t3[i][j],"","","");
  hfom_ac[i][j]=new TH2F(Form("hfom_ac[%d][%d]",i,j),"",iter_ac1,min_ac1,th1_max,iter_ac1,th2_min,th2_max);
  set->SetTH2(hfom_ac[i][j],"","","");


 }
 }





 
 hAC=new TH2F("hAC","AC threshold",
	      iter_ac1+1,min_ac1-(th1_max-min_ac1)/(2*iter_ac1),th1_max+(th1_max-min_ac1)/(2*iter_ac1)
	      ,nth,ac2_adc[0]-(ac2_adc[nth-1]-ac2_adc[0])/(2.0*nth),ac2_adc[nth-1]+(ac2_adc[nth-1]-ac2_adc[0])/(2.0*nth));
 set->SetTH2(hAC,"AC threshold","AC1","AC2");

 hAC2=new TH2F("hAC2","AC threshold",
	      iter_ac1+1,min_ac1-(th1_max-min_ac1)/(2*iter_ac1),th1_max+(th1_max-min_ac1)/(2*iter_ac1)
	      ,nth,ac2_adc[0]-(ac2_adc[nth-1]-ac2_adc[0])/(2.0*nth),ac2_adc[nth-1]+(ac2_adc[nth-1]-ac2_adc[0])/(2.0*nth));

 set->SetTH2(hAC2,"AC threshold","AC1","AC2"); 



 // TH2F* hct_a1a_c[24];
 // TH2F* hct_a2a_c[26]; 

 for(int i=0;i<24;i++){

   hct_a1a_c[i]=new TH2F(Form("hct_a1a_c_%i",i),"",bin_coin_c,min_coin_c,max_coin_c,2000,-100,10000);
   set->SetTH2(hct_a1a_c[i],Form("coin-time vs AC1 seg %d ADC hist",i),"coin time [ns]","ADC [ch]");
 }







 
 gfom=new TGraphErrors();
 set->SetGr(gfom,"","","");
 for(int k=0;k<nth;k++){
   
   gL_FOM_ac1[k]=new TGraphErrors();
   gL_FOM_ac1[k]->SetName(Form("gL_FOM_ac1_%d",k));
   gL_FOM_ac2[k]=new TGraphErrors();
   gL_FOM_ac2[k]->SetName(Form("gL_FOM_ac2_%d",k));   
   gL_N_ac1[k]=new TGraphErrors();
   gL_N_ac1[k]->SetName(Form("gL_N_ac1_%d",k));
   gL_N_ac2[k]=new TGraphErrors();
   gL_N_ac2[k]->SetName(Form("gL_N_ac2_%d",k));   
   gS_SN_ac1[k]=new TGraphErrors();
   gS_SN_ac1[k]->SetName(Form("gA_SN_ac1_%d",k));
   gS_SN_ac2[k]=new TGraphErrors();
   gS_SN_ac2[k]->SetName(Form("gA_SN_ac2_%d",k));   
   gL_SN_ac1[k]=new TGraphErrors();
   gL_SN_ac1[k]->SetName(Form("gL_SN_ac1_%d",k));   
   gL_SN_ac2[k]=new TGraphErrors();
   gL_SN_ac2[k]->SetName(Form("gL_SN_ac2_%d",k));      
   gS_eff_ac1[k]=new TGraphErrors();
   gS_eff_ac1[k]->SetName(Form("gS_eff_ac1_%d",k));   
   gS_eff_ac2[k]=new TGraphErrors();
   gS_eff_ac2[k]->SetName(Form("gS_eff_ac2_%d",k));      
   gL_eff_ac1[k]=new TGraphErrors();
   gL_eff_ac1[k]->SetName(Form("gL_eff_ac1_%d",k));   
   gL_eff_ac2[k]=new TGraphErrors();
   gL_eff_ac2[k]->SetName(Form("gL_eff_ac2_%d",k));      
   gL_ac1[k]=new TGraphErrors();
   gL_ac1[k]->SetName(Form("gL_ac1_%d",k));
   gL_ac2[k]=new TGraphErrors();
   gL_ac2[k]->SetName(Form("gL_ac2_%d",k));      
   gS_ac1[k]=new TGraphErrors();
   gS_ac1[k]->SetName(Form("gS_ac1_%d",k));   
   gS_ac2[k]=new TGraphErrors();
   gS_ac2[k]->SetName(Form("gS_ac2_%d",k));      
   gmm_ac1[k]=new TGraphErrors();

   gmm_ac2[k]=new TGraphErrors();
   gmm_SN_ac1[k]=new TGraphErrors();
   gmm_SN_ac2[k]=new TGraphErrors();
   gmm_S_ac1[k]=new TGraphErrors();
   gmm_S_ac2[k]=new TGraphErrors();
   gfom_ac1[k]=new TGraphErrors();
   gfom_ac2[k]=new TGraphErrors();
   gfom_ac1[k]->SetTitle(Form("FOM Estimation Graph %lf<AC2<%lf;AC1 ADC ;FOM ",ac2_adc[k],th2_max));
   gfom_ac2[k]->SetTitle(Form("FOM Estimation Graph AC1<%lf;AC2 ADC;FOM",ac1_adc[k]));
   gfom_ac1[k]->SetMarkerStyle(21);
   gfom_ac1[k]->SetMarkerColor(kRed);
   gfom_ac1[k]->SetMarkerSize(0.5);
   gfom_ac2[k]->SetMarkerStyle(21);
   gfom_ac2[k]->SetMarkerColor(kRed);
   gfom_ac2[k]->SetMarkerSize(0.5);
   set->SetGr(gfom_ac1[k],"","","");
   set->SetGr(gfom_ac2[k],"","",""); 
   set->SetGr(gmm_SN_ac1[k],"MM AC1 Prot","AC1 th","S/N");
   set->SetGr(gmm_SN_ac2[k],"MM AC2 Prot","AC2 th ","S/N");
   set->SetGr(gmm_S_ac1[k],"MM AC1 Prot","AC1 th","Lambda Survival Rate");
   set->SetGr(gmm_S_ac2[k],"MM AC2 Prot","AC2 th ","Lambda Survival Rate");
   set->SetGr(gL_ac1[k],"Lambda AC1 Prot","AC1 th","Lambda counts");
   set->SetGr(gL_ac2[k],"Lambda AC2 Prot","AC2 th ","Lambda counts");
   set->SetGr(gL_eff_ac1[k],"Lambda AC1 Kaon Survival rate Prot","AC1 th","Kaon Survival Rate");
   set->SetGr(gL_eff_ac2[k],"Lambda AC2 Kaon Survival Rate Prot","AC2 th ","Kaon Survival Rate");
   set->SetGr(gS_ac1[k],"Sigma AC1 Prot","AC1 th","Sigma counts");
   set->SetGr(gS_ac2[k],"Sigma AC2 Prot","AC2 th ","Sigma counts");
   set->SetGr(gS_eff_ac1[k],"Sigma AC1 Survival Rrate Prot","AC1 th","Kaon Survival Rate");
   set->SetGr(gS_eff_ac2[k],"Sigma AC2 Survival Rrate Prot","AC2 th ","Kaon Survival Rate");
   set->SetGr(gS_SN_ac1[k],"Sigma AC1 S/N Prot","AC1 th","S/N Rate");
   set->SetGr(gS_SN_ac2[k],"Sigma AC2 S/N Prot","AC2 th ","S/N Rate");
   set->SetGr(gL_SN_ac1[k],"Lambda AC1 S/N Prot","AC1 th","S/N Rate");
   set->SetGr(gL_SN_ac2[k],"Lambda AC2 S/N Prot","AC2 th ","S/N Rate");
   set->SetGr(gL_N_ac1[k],"Lambda AC1 BG  Prot","AC1 th","S/N Rate");
   set->SetGr(gL_N_ac2[k],"Lambda AC2 BG  Prot","AC2 th ","S/N Rate");
   set->SetGr(gL_FOM_ac1[k],"Lambda AC1 FOM Prot","AC1 th","FOM");
   set->SetGr(gL_FOM_ac2[k],"Lambda AC2 FOM Prot","AC2 th ","FOM");
 for(int j=0;j<nth;j++){

 gSN_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(gSN_k_ac1[j][k],"","","");
 gSN_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(gSN_k_ac2[j][k],"","","");
 gsum_pi_ac1[j][k]=new TGraphErrors();
 gsum_p_ac1[j][k]=new TGraphErrors();
 gsum_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(gsum_k_ac1[j][k],"","","");
 grate_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_k_ac1[j][k],"","","");
 grate_p_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_p_ac1[j][k],"","","");
 grate_pi_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_pi_ac1[j][k],"","","");
 gsum_pi_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_pi_ac2[j][k],"","","");
 gsum_p_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_pi_ac2[j][k],"","","");
 gsum_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_k_ac2[j][k],"","","");
 grate_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_k_ac2[j][k],"","","");
 grate_pi_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_pi_ac2[j][k],"","","");
 grate_p_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_p_ac2[j][k],"","","");
 //--- TGraphErrors Setting ----//
 gsum_pi_ac1[j][k]->SetTitle(Form("Sum of Pion (AC2<%lf);AC1 ADC-th [ch];Pion Events [Counts]",ac2_adc[k]));
 gsum_p_ac1[j][k]->SetTitle(Form("SUM of Proton (AC2<%lf);AC1 ADC-th [ch];Proton Events [Counts]",ac2_adc[k]));
 gsum_k_ac1[j][k]->SetTitle(Form("SUM of Kaon (AC2<%lf);AC1 ADC-th [ch];Kaon Events [Counts]",ac2_adc[k]));
 grate_k_ac1[j][k]->SetTitle(Form("Kaoin Survival rate ((AC2<%lf));AC1 ADC-th [ch]; Survival rate ",ac2_adc[k]));
 grate_pi_ac1[j][k]->SetTitle(Form("Pion Survuval rate (AC<%lf);AC1 ADC-th [ch]; Surival rate ",ac2_adc[k]));
 grate_p_ac1[j][k]->SetTitle("Proton Survival rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Survival rate "); 
 gSN_k_ac1[j][k]->SetTitle(Form("Kaon S/N rate AC1 [%d][%d]",j,k));
 gsum_pi_ac1[j][k]->SetMarkerStyle(21);
 gsum_pi_ac1[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac1[j][k]->SetMarkerSize(0.5);
 gsum_p_ac1[j][k]->SetMarkerStyle(21);
 gsum_p_ac1[j][k]->SetMarkerColor(kRed);
 gsum_p_ac1[j][k]->SetMarkerSize(0.5);
 gsum_k_ac1[j][k]->SetMarkerStyle(21);
 gsum_k_ac1[j][k]->SetMarkerColor(kRed);
 gsum_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_k_ac1[j][k]->SetMarkerStyle(21);
 grate_k_ac1[j][k]->SetMarkerColor(kBlue);
 grate_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_p_ac1[j][k]->SetMarkerStyle(21);
 grate_p_ac1[j][k]->SetMarkerColor(kBlue);
 grate_p_ac1[j][k]->SetMarkerSize(0.5);
 grate_pi_ac1[j][k]->SetMarkerStyle(21);
 grate_pi_ac1[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac1[j][k]->SetMarkerSize(0.5);
 gSN_k_ac1[j][k]->SetMarkerStyle(21);
 gSN_k_ac1[j][k]->SetMarkerColor(kRed);
 gSN_k_ac1[j][k]->SetMarkerSize(0.5);
 gsum_pi_ac2[j][k]->SetTitle("SUM of Pion vs AC2 Threshold;AC2 ADC-th [ch];Pion Events [Counts]");
 gsum_p_ac2[j][k]->SetTitle("SUM of Proton vs AC2 Threshold;AC2 ADC-th [ch];Proton Events [Counts]");
 gsum_k_ac2[j][k]->SetTitle("SUM of Kaon vs AC2 Threshold;AC2 ADC-th [ch];Kaon Events [Counts]");
 grate_k_ac2[j][k]->SetTitle("Kaoin Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_pi_ac2[j][k]->SetTitle("Pion Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_p_ac2[j][k]->SetTitle("Proton Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 gSN_k_ac2[j][k]->SetTitle(Form("Kaon S/N rate AC2 [%d][%d]",j,k));
 gsum_pi_ac2[j][k]->SetMarkerStyle(21);
 gsum_pi_ac2[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac2[j][k]->SetMarkerSize(0.5);
 gsum_p_ac2[j][k]->SetMarkerStyle(21);
 gsum_p_ac2[j][k]->SetMarkerColor(kRed);
 gsum_p_ac2[j][k]->SetMarkerSize(0.5);
 gsum_k_ac2[j][k]->SetMarkerStyle(21);
 gsum_k_ac2[j][k]->SetMarkerColor(kRed);
 gsum_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_k_ac2[j][k]->SetMarkerStyle(21);
 grate_k_ac2[j][k]->SetMarkerColor(kBlue);
 grate_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_p_ac2[j][k]->SetMarkerStyle(21);
 grate_p_ac2[j][k]->SetMarkerColor(kBlue);
 grate_p_ac2[j][k]->SetMarkerSize(0.5);
 grate_pi_ac2[j][k]->SetMarkerStyle(21);
 grate_pi_ac2[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac2[j][k]->SetMarkerSize(0.5);
 gSN_k_ac2[j][k]->SetMarkerStyle(21);
 gSN_k_ac2[j][k]->SetMarkerColor(kRed);
 gSN_k_ac2[j][k]->SetMarkerSize(0.5);
 }
 }



 


 //===== Fitting =======//

  bg_min=1.0;
  bg_max=1.15;
  bgs_min=1.19;
  bgs_max=1.23;

 
 facc_kc=new TF1("facc_kc","[0]",min_coin_c,max_coin_c);
 facc_kc->SetNpx(2000);
 fk_kc=new TF1("fk_kc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc->SetNpx(2000);
 fpi_pic=new TF1("fpi_pic","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic->SetNpx(2000);
 fp_pc=new TF1("fp_pc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc->SetNpx(2000);
 
  fLam_p=new TF1("fLam_p","gausn(0)",min_mm,max_mm);
      fLam_p->SetNpx(2000);
      fLam_p->SetLineColor(4);
      fLam_p->SetFillColor(4);
      fLam_p->SetFillStyle(3001);  
  fSig_p=new TF1("fSig_p","gausn(0)",min_mm,max_mm);
      fSig_p->SetNpx(2000);
      fSig_p->SetLineColor(4);
      fSig_p->SetFillColor(4);
      fSig_p->SetFillStyle(3001);  


      //======= Tuning =======//
 
 hcoin_fom=new TH1F("hcoin_fom","hcoin_fom",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_fom,"","","");
 hcoin_acc=new TH1F("hcoin_acc","hcoin_acc",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_acc,"","","");

 hmm=new TH1F("hmm","hcoin",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm,"Mass w/o AC tuning","Mass [GeV]","Counts/2 MeV"); 
      hmm->GetXaxis()->SetRangeUser(1.0,1.3);
      hmm->GetYaxis()->SetRangeUser(0.0,450.0);
      
 hmm_acc=new TH1F("hmm_acc","hcoin_fom",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm_acc,"ACC Mass w/o AC tuning ","Mass [GeV]","Counts/2 MeV"); 
 hmm_p=new TH1F("hmm_p","hmm_p",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm_p,"Mass peak w/o AC tuning","Mass [GeV]","Counts/2 MeV"); 


      
 hmm_fom=new TH1F("hmm_fom","hcoin_fom",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm_fom,"Mass w/ AC tuning","Mass [GeV]","Counts/2 MeV"); 
      hmm_fom->GetXaxis()->SetRangeUser(1.0,1.3);
      hmm_fom->GetYaxis()->SetRangeUser(0.0,450.0);

 hmm_fom_acc=new TH1F("hmm_fom_acc","hcoin_fom",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm_fom_acc,"ACC Mass w/ AC tuning ","Mass [GeV]","Counts/2 MeV"); 
 hmm_fom_p=new TH1F("hmm_fom_p","hmm_fom_p",bin_mm,min_mm,max_mm);
 set->SetTH1(hmm_fom_p,"Mass peak w/ AC tuning","Mass [GeV]","Counts/2 MeV"); 


      fL_p=new TF1("fL_p","gausn(0)",0.5,1.5);  
      fL_p->SetNpx(2000);
      fL_p->SetLineColor(4);
      fL_p->SetFillColor(4);
      fL_p->SetFillStyle(3001);
      fS_p=new TF1("fS_p","gausn(0)",0.5,1.5);  
      fS_p->SetNpx(2000);
      fS_p->SetLineColor(4);
      fS_p->SetFillColor(4);
      fS_p->SetFillStyle(3001);


 fL_fom=new TF1("fL_fom","gausn(0)+pol1(3)",1.05,1.15);
	fL_fom->SetNpx(2000);	

 fS_fom=new TF1("fS_fom","gausn(0)+pol1(3)",1.15,1.25);
	fS_fom->SetNpx(2000);

	
 fL_fom_bg=new TF1("fL_fom_bg","pol1(0)",min_mm,max_mm);
 fS_fom_bg=new TF1("fS_fom_bg","pol1(0)",min_mm,max_mm);

 fL_fom_p=new TF1("fL_fom_p","gausn(0)",min_mm,max_mm);
      fL_fom_p->SetNpx(2000);
      fL_fom_p->SetLineColor(4);
      fL_fom_p->SetFillColor(4);
      fL_fom_p->SetFillStyle(3001);
 fS_fom_p=new TF1("fS_fom_p","gausn(0)",min_mm,max_mm);
      fS_fom_p->SetNpx(2000);
      fS_fom_p->SetLineColor(4);
      fS_fom_p->SetFillColor(4);
      fS_fom_p->SetFillStyle(3001);

 hcoin_fom_p=new TH1F("hcoin_fom_p","hcoin_fom_p",bin_coin_c,min_coin_c,max_coin_c);
 fk_fom=new TF1("fk_fom","gausn(0)",min_coin_c,max_coin_c);
      fk_fom->SetLineColor(4);
      fk_fom->SetFillColor(4);
      fk_fom->SetFillStyle(3001);


      hRu1_time_c=new TH1F("hRu1_time_c","RHRS VDC U1 TDC coin trigger w/ PID & Path length & z cut ",400,-0.1e-6,0.5e-6);
      set->SetTH1(hRu1_time_c,"RHRS VDC U1 TDC coin trigger w/ PID & Path length & z cut","TDC [sec]","Counts");
      hRu1_time_s=new TH1F("hRu1_time_s","RHRS VDC U1 TDC T4 trigger w/ PID & Path length & z cut ",400,-0.1e-6,0.5e-6);
      set->SetTH1(hRu1_time_s,"RHRS VDC U1 TDC T4 trigger w/ PID & Path length & z cut","TDC [sec]","Counts");      
      
}




//////////////////////////////////////////////////////////////////////////

void tuningAC::Fill(){

 //--- Coin Offset -----//
 pathl_off=0.0;
 pathl_off=0.0;
 if(mode=="H"&&kine==2){
 coin_offset=498.;
 s2_offset=-489.0;
 pathl_off=-485.5; 

 }else if(mode=="H" && kine==1){
 pathl_off=-498.+30.-3.0+0.5;
 s2_offset=-500.0+25.;
 coin_offset=-41.35+498.;
}


 int ev=0;
 double ct;
 for(int k=0;k<ENum;k++){
   T->GetEntry(k);

   if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<ENum<<endl; 
     ev=ev+1;
   }


 pe_=Lp[0];//*sqrt(1+pow(Lth[0],2)+pow(Lph[0],2));
 pk=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 ppi=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 pe=4.313; // GeV   //hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Epi=sqrt(pow(ppi,2)+pow(mpi,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rbeta=pk/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 

 
 if(mode=="H"){
 Rs2_off=s2f1_off(Rs2pads,"R",mode,kine);
 Ls2_off=s2f1_off(Ls2pads,"L",mode,kine);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
  coin_t=tof_r-tof_l-coin_offset-s2_offset; //coin time
  coin_tc=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction 
  }

//====== Cut condition ========================// 
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_track=false;
   cut_s0=false;
   coin_trig=false;

   if(mode=="H" || mode=="G"){
     Rz=Rvz[0];
     Lz=Lvz[0];
    }
   // if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
  
 
   if(Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax)cut_vz=true;

   if(mode=="G"){
     coin_t=ctime[0];
     coin_tc=ctime[0];
     if(coin_t>-100)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true;
  }else if(mode=="T"){
     coin_tc=tcoin_t;
     coin_t=tcoin_t;
  if(-80.<coin_tc && coin_tc<80.)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true; 
  }else if(mode=="H"){
 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
  if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true;  
   }

   

 //=======================================//
  


   //==========================================//
   //========= Fill Hist =====================//
   //========================================//


   if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){
 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2
 hcoin_ac1_all->Fill(coin_tc,Ra1sum);
 hcoin_ac2_all->Fill(coin_tc,Ra2sum);

 
 //========== W/o AC cut Missing mass ==============//

 if(coin_tc<1.0 && -1.0<coin_tc) hmm->Fill(mm);
 if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
	    hmm_acc->Fill(mm);
	    ct=-2222.0;
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	}
 }
   
 //==============================================//

 
 //--------- with AC1 Cut ---------------// 
  for(int j=0;j<nth;j++){
      cut_ac1=false;
    //cut_ac1=true;
   if(Ra1sum<ac1_adc[j])cut_ac1=true;
   //if(Ra1sum<0.04)cut_ac1=true;

   if(ac2_min && cut_ac1){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    hcoin_ac2_max[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
   if(coin_tc<1.0 && -1.0<coin_tc)hmm_ac2[j]->Fill(mm,Ra2sum);
   //------- ACC BG -------------------------------------------------------------//
     

   if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac2[j]->Fill(ct);
                 hcoin_acc_ac2[j]->Fill(ct-20);
		 hcoin_ac2_acc[j]->Fill(ct,Ra2sum);
                 hmm_ac2_acc[j]->Fill(mm,Ra2sum);
		 
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	}
	
   }

   
    //----------------------------------------------------------------------------//	  

   }else if(ac2_min==0 && cut_ac1){
	    // && Ra2sum>th2_min){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    hcoin_ac2_max[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    if(coin_tc<1.0 && -1.0<coin_tc)hmm_ac2[j]->Fill(mm,Ra2sum);
 //------- ACC BG -------------------------------------------------------------//
     
  if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac2[j]->Fill(ct);
                 hcoin_acc_ac2[j]->Fill(ct-20);
		 hcoin_ac2_acc[j]->Fill(ct,Ra2sum);
		 hmm_ac2_acc[j]->Fill(mm,Ra2sum);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	 }
     }
    //----------------------------------------------------------------------------//	
   }  
  }

 
 
   
   
   
  //-------with AC2 Cut --------------------//
   for(int i=0;i<nth;i++){
     cut_ac2=false;
     //    if(ac2_min && Ra2sum>ac2_adc[i] && Ra2sum<th_ac2_t)cut_ac2=true;
     //    if(ac2_min==0 && Ra2sum>th_ac2_b && Ra2sum<ac2_adc[i])cut_ac2=true;

     
     //      if(ac2_min && Ra2sum>ac2_adc[i] && Ra2sum<th_ac2_t)cut_ac2=true;
     //      if(ac2_min==0 && Ra2sum>th_ac2_b && Ra2sum<ac2_adc[i])cut_ac2=true;

     //-------- Coin time study -----------//
     if(ac2_min &&Ra2sum>ac2_adc[i] &&Ra2sum< ac2_kcut_max)cut_ac2=true;
     //   if(ac2_min==0 && Ra2sum<ac2_adc[i] &&Ra2sum>th2_min)cut_ac2=true;
     // if(10.0>Ra2sum && Ra2sum>3.0)cut_ac2=true;
    if(cut_ac2 && Ra1sum<th1_max){
    hcoin_t1[i]->Fill(coin_tc); // AC2 cut
    hcoin_ac1[i]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  B
    hcoin_ac1_max[i]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  
    if(coin_tc<1.0 && -1.0<coin_tc)hmm_ac1[i]->Fill(mm,Ra1sum);


 //------- ACC BG -------------------------------------------------------------//
 
   if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
        ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac1[i]->Fill(ct);
                 hcoin_acc_ac1[i]->Fill(ct-20);
		 hcoin_ac1_acc[i]->Fill(ct,Ra1sum);
		 hmm_ac1_acc[i]->Fill(mm,Ra1sum);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.0;}
	 }
    }
 
    //----------------------------------------------------------------------------//

    }
   }
   }



 //-------with AC1 && AC2 Cut --------------------//
  cut_ac1=false;
  cut_ac2=false; 
   for(int k=0;k<nth;k++){
   for(int j=0;j<nth;j++){ 
  if(Ra1sum<ac1_adc[j])cut_ac1=true;
  if(Ra2sum<ac2_adc[k])cut_ac2=true;
  if(cut_ac1 && cut_ac2)hcoin_t3[j][k]->Fill(coin_tc); //AC1 & AC2 cut
  } 
   }
   //--------- Kaon Cut Hist ---------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum<ac1_kcut 
    && ac2_kcut_min<Ra2sum && Ra2sum<ac2_kcut_max)hcoin_k->Fill(coin_tc);
   //--------- Pion Cut Hist ---------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut 
    && Ra2sum>ac2_kcut_max){
   hcoin_pi->Fill(coin_tc);
   if(NRu1_time==5)hRu1_time_c->Fill(Ru1_time[2]);
   
 }
  if(cut_Rs2 && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut 
    && Ra2sum>ac2_kcut_max && NRu1_time==5)hRu1_time_s->Fill(Ru1_time[2]);

  //-------- Proton Cut Hist --------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum<ac1_kcut 
    && Ra2sum<ac2_kcut_min)hcoin_p->Fill(coin_tc);

 }

 for(int i=0;i<nth;i++){
 hcoin_acc_ac1[i]->Scale(20./100.);
 hcoin_acc_ac2[i]->Scale(20./100.);
 hmm_ac1_acc[i]->Scale(20./100.);
 hmm_ac2_acc[i]->Scale(20./100.);
 hcoin_ac1_acc[i]->Scale(20./100.);
 hcoin_ac2_acc[i]->Scale(20./100.);
 }

 


}


////////////////////////////////////////////////////////////////////////////////////

void tuningAC::Fitting(){


 if(mode=="H" && kine==2){
 def_sig_p=0.852; def_mean_p=0.0;
 def_sig_pi=0.443; def_mean_pi=11;
 def_sig_k=0.644; def_mean_k=8.0;
 def_acc=27.7;
 }else if(mode=="H" && kine==1){
 def_sig_p=0.852; def_mean_p=0.0;
 def_sig_pi=0.443; def_mean_pi=11;
 def_sig_k=0.644; def_mean_k=8.;
 def_acc=27.7;
 }else if(mode=="G"){
 def_sig_p=0.852; def_mean_p=11.06;
 def_sig_pi=0.4; def_mean_pi=0.0;
 def_sig_k=0.644; def_mean_k=3.16;
 def_acc=22.7;
 }else if(mode=="T"){
 def_sig_p=0.841; def_mean_p=-8.21;
 def_sig_pi=0.275; def_mean_pi=3.049;
 def_sig_k=0.454; def_mean_k=-0.095;
 def_acc=22.7;
}
    
 sum_k_max=1250.;
 fom_max=0.0;

 for(int i=0;i<nth;i++){
 max_SN_ac1[i]=0.0;
 max_SN_ac2[i]=0.0;
 SN_ac1[i]=0;
 SN_ac2[i]=0;
 }


 for(int i=0;i<iter_max;i++)emp[i]=0.0;

cout<<"defined parameters in SR analysis"<<endl;
 
         //----------------------------------//


 fpi_pic->SetParameter(1,def_mean_pi);
 fpi_pic->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi_pic->SetParameter(2,def_sig_pi);
 fpi_pic->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi_pic->SetParameter(3,def_acc);

 fp_pc->SetParameter(1,def_mean_p);
 fp_pc->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp_pc->SetParameter(2,def_sig_p);
 fp_pc->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
 fp_pc->SetParameter(3,def_acc);

 hcoin_k->Fit("facc_kc","Rq","",min_coin_c,min_coin_c+3.0);
 def_acc_k=facc_kc->GetParameter(0);

 fk_kc->SetParameter(1,def_mean_k);
 fk_kc->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk_kc->SetParameter(2,def_sig_k);
 fk_kc->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk_kc->FixParameter(3,def_acc_k);

 hcoin_k->Fit("fk_kc","Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_pi->Fit("fpi_pic","Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_p->Fit("fp_pc","Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);

 def_num_k=fk_kc->GetParameter(0);
 def_mean_k=fk_kc->GetParameter(1);
 def_sig_k=fk_kc->GetParameter(2);

 def_num_p=fp_pc->GetParameter(0);
 def_mean_p=fp_pc->GetParameter(1);
 def_sig_p=fp_pc->GetParameter(2);

 def_num_pi=fpi_pic->GetParameter(0);
 def_mean_pi=fpi_pic->GetParameter(1);
 def_sig_pi=fpi_pic->GetParameter(2);


 //============= SRA ======================//

     for(int th1=0;th1<nth;th1++){
     for(int th2=0;th2<nth;th2++){
     if( th1==th2){
     for(int i=0;i<iter_ac1;i++){

    th_ac1[i]=th1_max-(th1_max-min_ac1)/iter_ac1*i;
    th_ac2[i]=min_ac2+(th2_max-min_ac2)/iter_ac1*i;


    
    //    cout<<"th1: "<<th1+1<<"/"<<nth<<":th2: "<<th2+1<<"/"<<nth<<": "
    //    <<"i "<<i<<"/"<<iter_ac1<<endl;
    cout<<"====== th: "<<th1+1<<"/"<<nth<<"   i: "<<i+1<<" / "<<iter_ac1<<" ========"<<endl;
    cout<<"th_ac1:"<<th_ac1[i]<<"/"<<th1_max<<"    th_ac2: "<<th_ac2[i]<<"/"<<th2_max<<endl;


 bin_th_ac1[i][th2]=hcoin_ac1[th2]->GetYaxis()->FindBin(th_ac1[i]);//hcoin_ac1[th2] is ac2 fixed threshold
 bin_min_ac1=hcoin_ac1[th2]->GetYaxis()->FindBin(min_ac1);
 hcoin_ac1_all_p[i][th2]=hcoin_ac1[th2]->ProjectionX(Form("hcoin_all_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_all_p[i][th2]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_all_p[i][th2],"","","");
 hcoin_ac1_acc_p[i][th2]=hcoin_ac1_acc[th2]->ProjectionX(Form("hcoin_ac1_acc_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_acc_p[i][th2]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_acc_p[i][th2],"","","");
 hcoin_ac1_acc_p[i][th2]->Scale(20./100.);
 hcoin_ac1_p[i][th2]=hcoin_ac1[th2]->ProjectionX(Form("hcoin_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_p[i][th2]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_p[i][th2],"","","");
//change
  hcoin_ac1_p[i][th2]->Add(hcoin_ac1_all_p[i][th2], hcoin_ac1_acc_p[i][th2],1.0,-1.0);
  // hcoin_ac1_p[i][th2]->Add(hcoin_ac1_acc_p[i][th2],-1.0);






 hmm_ac1_acc_p[i][th2]=hmm_ac1_acc[th2]->ProjectionX(Form("hcoin_ac1_acc_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hmm_ac1_acc_p[i][th2]->SetTitle(Form("Missing Mass (ACC) AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_acc_p[i][th2],"","","");
 hmm_ac1_acc_p[i][th2]->Scale(20./100.);

 hmm_ac1_all_p[i][th2]=hmm_ac1[th2]->ProjectionX(Form("hmm_ac1_all_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hmm_ac1_all_p[i][th2]->SetTitle(Form("Missing Mass (All) AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hmm_ac1_all_p[i][th2],Form("hmm_ac1_all_p[%d][%d]",i,th2),"","");

 hmm_ac1_p[i][th2]=hmm_ac1[th2]->ProjectionX(Form("hmm_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 set->SetTH1(hmm_ac1_p[i][th2],Form("hmm_ac1_p[%d][%d]",i,th2),"","");
 hmm_ac1_p[i][th2]->Add(hmm_ac1_acc_p[i][th2],-1.0);



 //===================================//
 //===== AC2 Upper threshold =========//
 //===================================//
  

  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));

  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  //hcoin_ac2_acc_p[i][th1]->Scale(20./100.);

  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");

  //============ MM AC2 =======================//
  //hmm_ac2_acc_p[fom_th2][fom_max_th1]->Draw("same");
  hmm_ac2_acc_p[i][th1]=hmm_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
						     //bin_min_ac2,bin_th_ac2[i][th1]);
  // hmm_ac2_acc_p[i][th1]->SetTitle(Form("Missing Mass (ACC) AC2<%lf & %lf<AC2<%lf",th_ac2[i],ac2_adc[th1],th1_max));
 set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
 //hmm_ac2_acc_p[i][th1]->Scale(20./100.);
 hmm_ac2_all_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
 set->SetTH1(hmm_ac2_all_p[i][th1],Form("hmm_ac2_all_p[%d][%d]",i,th1),"","");
 hmm_ac2_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
 set->SetTH1(hmm_ac2_p[i][th1],Form("hmm_ac2_p[%d][%d]",i,th1),"","");
 hmm_ac2_p[i][th1]->Add(hmm_ac2_acc_p[i][th1],-1.0);
 





 if(ac2_min){

   /*
  th_ac2[i]=min_ac2+(th2_max-min_ac2)/iter_ac1*i;

  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));

  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  hcoin_ac2_acc_p[i][th1]->Scale(20./100.);

  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  //hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
  hcoin_ac2_p[i][th1]->Add(hcoin_ac2_acc_p[i][th1],-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");

   


   
  //============ MM AC2 =======================//
  //hmm_ac2_acc_p[fom_th2][fom_max_th1]->Draw("same");
  hmm_ac2_acc_p[i][th1]=hmm_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
						     //bin_min_ac2,bin_th_ac2[i][th1]);
  // hmm_ac2_acc_p[i][th1]->SetTitle(Form("Missing Mass (ACC) AC2<%lf & %lf<AC2<%lf",th_ac2[i],ac2_adc[th1],th1_max));
 set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
 //hmm_ac2_acc_p[i][th1]->Scale(20./100.);
 hmm_ac2_all_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
					
 set->SetTH1(hmm_ac2_all_p[i][th1],Form("hmm_ac2_all_p[%d][%d]",i,th1),"","");

 hmm_ac2_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
 set->SetTH1(hmm_ac2_p[i][th1],Form("hmm_ac2_p[%d][%d]",i,th1),"","");
 hmm_ac2_p[i][th1]->Add(hmm_ac2_acc_p[i][th1],-1.0);

   */
   


 }else{

  th_ac2[i]= th2_max-(th2_max-min_ac2)/iter_ac1*i;
  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]);  
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],bin_min_ac2,th_ac2[i]));
  hcoin_ac2_acc_p[i][th1]->Scale(20./100.);
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],bin_min_ac2,bin_th_ac2[i][th1])); 
  //hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
hcoin_ac2_p[i][th1]->Add(hcoin_ac2_acc_p[i][th1],-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");


  hmm_ac2_acc_p[i][th1]=hmm_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]);

 set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
 //hmm_ac2_acc_p[i][th1]->Scale(20./100.);
 hmm_ac2_all_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_all_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); 
					
 set->SetTH1(hmm_ac2_all_p[i][th1],Form("hmm_ac2_all_p[%d][%d]",i,th1),"","");

 hmm_ac2_p[i][th1]=hmm_ac2[th1]->ProjectionX(Form("hmm_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]);
 set->SetTH1(hmm_ac2_p[i][th1],Form("hmm_ac2_p[%d][%d]",i,th1),"","");
 hmm_ac2_p[i][th1]->Add(hmm_ac2_acc_p[i][th1],-1.0);




 }

    

 // hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); //hcoin_ac2[th1] is ac1 fixed threshold

 //--- Initial Parameters -----------//

  facc[i][th2][0]=new TF1(Form("facc[%d][%d][0]",i,th2),"[0]",min_coin_c,max_coin_c);
  facc[i][th2][0]->SetNpx(2000);
  hcoin_ac1_p[i][th2]->Fit(Form("facc[%d][%d][0]",i,th2),"Rq","",min_coin_c,min_coin_c+3.);
  p0_acc[i][th2][0]=facc[i][th2][0]->GetParameter(0);
  facc[i][th1][1]=new TF1(Form("facc[%d][%d][1]",i,th1),"[0]",min_coin_c,max_coin_c);
  facc[i][th1][1]->SetNpx(2000);  
  hcoin_ac2_p[i][th1]->Fit(Form("facc[%d][%d][1]",i,th1),"Rq","",min_coin_c,min_coin_c+3.);
  p0_acc[i][th1][1]=facc[i][th1][1]->GetParameter(0);

  //Lambda BG  AC1
  //fbg[i][th2][0]=new TF1(Form("fbg[%d][%d][0]",i,th2),"pol2(0)",min_mm,max_mm);
  fbg[i][th2][0]=new TF1(Form("fbg[%d][%d][0]",i,th2),"pol1(0)",min_mm,max_mm);
  hmm_ac1_all_p[i][th2]->Fit(Form("fbg[%d][%d][0]",i,th2),"Rq","",bg_min,bg_max);
  bg_0[i][th2][0]=fbg[i][th2][0]->GetParameter(0);
  bg_1[i][th2][0]=fbg[i][th2][0]->GetParameter(1);
//  bg_2[i][th2][0]=fbg[i][th2][0]->GetParameter(2);
  
  //Sigma BG AC1   
  fbg_s[i][th2][0]=new TF1(Form("fbg_s[%d][%d][0]",i,th2),"pol1(0)",min_mm,max_mm);
  hmm_ac1_all_p[i][th2]->Fit(Form("fbg_s[%d][%d][0]",i,th2),"Rq","",bgs_min,bgs_max);
  bg_s0[i][th2][0]=fbg_s[i][th2][0]->GetParameter(0);
  bg_s1[i][th2][0]=fbg_s[i][th2][0]->GetParameter(1);
  //  bg_s2[i][th2][0]=fbg_s[i][th2][0]->GetParameter(2);
 
  //Lambda BG AC2
  //fbg[i][th1][1]=new TF1(Form("fbg[%d][%d][1]",i,th1),"pol2(0)",min_mm,max_mm);
  fbg[i][th1][1]=new TF1(Form("fbg[%d][%d][1]",i,th1),"pol1(0)",min_mm,max_mm);
  hmm_ac2_all_p[i][th1]->Fit(Form("fbg[%d][%d][1]",i,th1),"Rq","",bg_min,bg_max);
  bg_0[i][th1][1]=fbg[i][th1][1]->GetParameter(0);
  bg_1[i][th2][1]=fbg[i][th1][1]->GetParameter(1);
  //Sigma BG AC2
  fbg_s[i][th1][1]=new TF1(Form("fbg_s[%d][%d][1]",i,th1),"pol1(0)",min_mm,max_mm);
  hmm_ac2_all_p[i][th1]->Fit(Form("fbg_s[%d][%d][1]",i,th1),"Rq","",1.17,1.23);
			     //,bgs_min,bgs_max);
  bg_s0[i][th1][1]=fbg_s[i][th1][1]->GetParameter(0);
  bg_s1[i][th2][1]=fbg_s[i][th1][1]->GetParameter(1);
  //  bg_s2[i][th2][1]=fbg_s[i][th1][1]->GetParameter(2);
  //Lambda Peak AC1
  
  fLam[i][th2][0]=new TF1(Form("fLam[%d][%d][0]",i,th2),"gausn(0)+pol1(3)",min_mm,max_mm);
  fLam[i][th2][0]->SetParameter(1,1.116);
  //  fLam[i][th2][0]->SetParLimits(1,1.115,1.117);
  fLam[i][th2][0]->SetParameter(2,3.1e-03);
  //  fLam[i][th2][0]->SetParLimits(2,3e-03,3.5e-03);
  fLam[i][th2][0]->SetParameter(3,bg_0[i][th2][0]);
  fLam[i][th2][0]->SetParameter(4,bg_1[i][th2][0]);  
  //  fLam[i][th2][0]->FixParameter(3,bg_0[i][th2][0]);
  //  fLam[i][th2][0]->FixParameter(4,bg_1[i][th2][0]);
  hmm_ac1_all_p[i][th2]->Fit(Form("fLam[%d][%d][0]",i,th2),"Rq","",1.10,1.13);
  L0[i][th2][0]=fLam[i][th2][0]->GetParameter(0);
  L1[i][th2][0]=fLam[i][th2][0]->GetParameter(1);
  L2[i][th2][0]=fLam[i][th2][0]->GetParameter(2);


  nL[i][th2][0]=fLam[i][th2][0]->GetParameter(0);
  nL[i][th2][0]=  nL[i][th2][0]/0.002;
  nL_err[i][th2][0]=sqrt(nL[i][th2][0]);


    //fLam[i][th2][0]->GetParError(0);
  //  nL_err[i][th2][0]=  nL_err[i][th2][0]/0.002;
  meanL[i][th2][0]=fLam[i][th2][0]->GetParameter(1);
  sigL[i][th2][0]=fLam[i][th2][0]->GetParameter(2);
  
  // Sigma Peak AC1
  fSig[i][th2][0]=new TF1(Form("fSig[%d][%d][0]",i,th2),"gausn(0)+pol1(3)",min_mm,max_mm);

  fSig[i][th2][0]->SetParameter(0,1);  
  //  fSig[i][th2][0]->SetParLimits(0,0.0,100);
  fSig[i][th2][0]->SetParameter(1,1.20);
  //  fSig[i][th2][0]->SetParLimits(1,1.18,1.22);
  fSig[i][th2][0]->SetParameter(2,6e-3);
  //  fSig[i][th2][0]->SetParLimits(2,5e-03,7.0e-03);
  fSig[i][th2][0]->SetParameter(3,bg_s0[i][th2][0]);
  fSig[i][th2][0]->SetParameter(4,bg_s1[i][th2][0]);  
  //  fSig[i][th2][0]->FixParameter(3,bg_s0[i][th2][0]);
  //  fSig[i][th2][0]->FixParameter(4,bg_s1[i][th2][0]);
  hmm_ac1_all_p[i][th2]->Fit(Form("fSig[%d][%d][0]",i,th2),"Rq","",bgs_min,bgs_max);
  
  S0[i][th2][0]=fSig[i][th2][0]->GetParameter(0);
  S1[i][th2][0]=fSig[i][th2][0]->GetParameter(1);
  S2[i][th2][0]=fSig[i][th2][0]->GetParameter(2);
  nS[i][th2][0]=fSig[i][th2][0]->GetParameter(0);
  nS[i][th2][0]=nS[i][th2][0]/0.002;
  meanS[i][th2][0]=fSig[i][th2][0]->GetParameter(1);
  sigS[i][th2][0]=fSig[i][th2][0]->GetParameter(2);


 



  
  //Lambda Peak AC2
  fLam[i][th1][1]=new TF1(Form("fLam[%d][%d][1]",i,th1),"gausn(0)+pol1(3)",min_mm,max_mm);
  fLam[i][th1][1]->SetParameter(1,1.116); 
  //  fLam[i][th1][1]->SetParLimits(1,1.115,1.117); 
  fLam[i][th1][1]->SetParameter(2,3.1e-3); 
  //  fLam[i][th1][1]->SetParLimits(2,3.0e-03,3.5e-03);
  fLam[i][th1][1]->SetParameter(3,bg_0[i][th1][1]);
  fLam[i][th1][1]->SetParameter(4,bg_1[i][th1][1]);  
  //  fLam[i][th1][1]->FixParameter(3,bg_0[i][th1][1]);
  //  fLam[i][th1][1]->FixParameter(4,bg_1[i][th1][1]);
  //  fLam[i][th1][1]->FixParameter(5,bg_2[i][th1][1]);
  hmm_ac2_all_p[i][th1]->Fit(Form("fLam[%d][%d][1]",i,th1),"Rq","",1.10,1.13);

  L0[i][th1][1]=fLam[i][th1][1]->GetParameter(0);
  L1[i][th1][1]=fLam[i][th1][1]->GetParameter(1);
  L2[i][th1][1]=fLam[i][th1][1]->GetParameter(2);

  nL[i][th1][1]=fLam[i][th1][1]->GetParameter(0);
  nL[i][th1][1]= nL[i][th1][1]/0.002;
  nL_err[i][th1][1]= sqrt(nL[i][th1][1]);
    //fLam[i][th1][1]->GetParError(0);
  //  nL_err[i][th1][1]=nL_err[i][th1][1]/0.002;
  meanL[i][th1][1]=fLam[i][th1][1]->GetParameter(1);
  sigL[i][th1][1]=fLam[i][th1][1]->GetParameter(2);

  // Sigma Peak AC2
  fSig[i][th1][1]=new TF1(Form("fSig[%d][%d][0]",i,th1),"gausn(0)+pol1(3)",min_mm,max_mm);
  fSig[i][th1][1]->SetParameter(0,1.0);
  //  fSig[i][th1][1]->SetParLimits(0,0.0,100);
  fSig[i][th1][1]->SetParameter(1,1.20);  
  //  fSig[i][th1][1]->SetParLimits(1,1.19,1.21);
  fSig[i][th1][1]->SetParameter(2,6.0e-3);
  //  fSig[i][th1][1]->SetParLimits(2,5e-03,1e-02);
  fSig[i][th1][1]->SetParameter(3,bg_s0[i][th1][1]);
  fSig[i][th1][1]->SetParameter(4,bg_s1[i][th1][1]);  
  //  fSig[i][th1][1]->FixParameter(3,bg_s0[i][th1][1]);
  //  fSig[i][th1][1]->FixParameter(4,bg_s1[i][th1][1]);
  //  fSig[i][th1][1]->FixParameter(5,bg_s2[i][th1][1]);
  hmm_ac2_all_p[i][th1]->Fit(Form("fSig[%d][%d][0]",i,th1),"Rq","",bgs_min,bgs_max);
  S0[i][th1][1]=fSig[i][th1][1]->GetParameter(0);
  S1[i][th1][1]=fSig[i][th1][1]->GetParameter(1);
  S2[i][th1][1]=fSig[i][th1][1]->GetParameter(2);




  nS[i][th1][1]=fSig[i][th1][1]->GetParameter(0);
  nS[i][th1][1]=nS[i][th1][1]/0.002;
  meanS[i][th1][1]=fSig[i][th1][1]->GetParameter(1);
  sigS[i][th1][1]=fSig[i][th1][1]->GetParameter(2);




  gL_ac1[th2]->SetPoint(i,th_ac1[i],nL[i][th2][0]);
  gL_ac1[th2]->SetPointError(i,0.0,nL_err[i][th2][0]);
  gL_ac2[th1]->SetPoint(i,th_ac2[i],nL[i][th1][1]);
  gL_ac2[th1]->SetPointError(i,0.0,nL_err[i][th1][1]);
  gL_eff_ac1[th2]->SetPoint(i,th_ac1[i],nL[i][th2][0]/nL[0][th2][0]);
  gL_eff_ac2[th1]->SetPoint(i,th_ac2[i],nL[i][th1][1]/nL[0][th1][1]);

  gL_eff_ac1[th2]->SetPointError(i,0.0,sqrt(1./nL[i][th2][0]*(nL[i][th2][0]/nL[0][th2][0])*(1.-nL[i][th2][0]/nL[0][th2][0])));
  gL_eff_ac2[th1]->SetPointError(i,0.0,sqrt(1./nL[i][th1][1]*(nL[i][th1][1]/nL[0][th1][1])*(1.-nL[i][th1][1]/nL[0][th1][1])));





  gS_ac1[th2]->SetPoint(i,th_ac1[i],nS[i][th2][0]);
  gS_ac1[th2]->SetPointError(i,0,nS_err[i][th2][0]);
  gS_ac2[th1]->SetPoint(i,th_ac2[i],nS[i][th1][1]);
  gS_ac2[th1]->SetPointError(i,0,nS_err[i][th1][1]);
  gS_eff_ac1[th2]->SetPoint(i,th_ac1[i],nS[i][th2][0]/nS[0][th2][0]);
  gS_eff_ac2[th1]->SetPoint(i,th_ac2[i],nS[i][th1][1]/nS[0][th1][1]);







   totL_ac1[i][th2]=hmm_ac1_all_p[i][th2]->Integral(hmm_ac1_all_p[i][th2]->FindBin(-3*sigL[i][th2][0]+meanL[i][th2][0])
  		      ,hmm_ac1_all_p[i][th2]->FindBin(+3*sigL[i][th2][0]+meanL[i][th2][0]));

   bgL_ac1[i][th2]=totL_ac1[i][th2]-nL[i][th2][0];


   //   Lfom[i][th2][0]=sqrt(nL[i][th2][0]*nL[i][th2][0]/bgL_ac1[i][th2]);
   Lfom[i][th2][0]=sqrt(nL[i][th2][0]*nL[i][th2][0]/totL_ac1[i][th2]);
   gL_SN_ac1[th2]->SetPoint(i,th_ac1[i],nL[i][th2][0]/bgL_ac1[i][th2]);
   gL_N_ac1[th2]->SetPoint(i,th_ac1[i],bgL_ac1[i][th2]);

    

   //   gL_FOM_ac1[th2]->SetPoint(i,th_ac1[i],sqrt(nL[i][th2][0]*nL[i][th2][0]/bgL_ac1[i][th2]));
   gL_FOM_ac1[th2]->SetPoint(i,th_ac1[i],Lfom[i][th2][0]);

   totL_ac2[i][th1]=hmm_ac2_all_p[i][th1]->Integral(hmm_ac2_all_p[i][th1]->FindBin(-3*sigL[i][th1][1]+meanL[i][th1][1])
					     ,hmm_ac2_all_p[i][th1]->FindBin(+3*sigL[i][th1][1]+meanL[i][th1][1]));
   bgL_ac2[i][th1]= totL_ac2[i][th1]-nL[i][th1][1];

   gL_SN_ac2[th1]->SetPoint(i,th_ac2[i],nL[i][th1][1]/bgL_ac2[i][th1]);
   gL_N_ac2[th1]->SetPoint(i,th_ac2[i],bgL_ac2[i][th1]);
   //   Lfom[i][th1][1]=sqrt(nL[i][th1][1]*nL[i][th1][1]/bgL_ac2[i][th1]);
   Lfom[i][th1][1]=sqrt(nL[i][th1][1]*nL[i][th1][1]/totL_ac2[i][th1]);
   gL_FOM_ac2[th1]->SetPoint(i,th_ac2[i],Lfom[i][th1][1]);

 //------- AC1 --------------// 

 fp[i][th2][0]=new TF1(Form("fp[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][th2][0]->SetNpx(2000);
 fpi[i][th2][0] =new TF1(Form("fpi[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fpi[i][th2][0]->SetNpx(2000);
 fk[i][th2][0]=new TF1(Form("fk[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][th2][0]->SetNpx(2000);


 fp[i][th2][0]->FixParameter(1,def_mean_p); 
 fp[i][th2][0]->FixParameter(2,def_sig_p); 
 
 // fpi[i][th2][0]->FixParameter(1,def_mean_pi); 
 //fpi[i][th2][0]->FixParameter(2,def_sig_pi); 
 //fk[i][th2][0]->FixParameter(1,def_mean_k); 
 //fk[i][th2][0]->FixParameter(2,def_sig_k);

 fpi[i][th2][0]->FixParameter(1,def_mean_pi); 
 fk[i][th2][0]->FixParameter(1,def_mean_k);
 fk[i][th2][0]->SetParameter(2,def_sig_k);  
 fpi[i][th2][0]->SetParLimits(2,0.95*def_sig_pi,1.05*def_sig_pi); 
 fk[i][th2][0]->SetParLimits(2,0.75*def_sig_k,1.25*def_sig_k);



//----- AC1 Fitting -----------//

 hcoin_ac1_p[i][th2]->Fit(Form("fp[%d][%d][0]",i,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th2][0]=fp[i][th2][0]->GetParameter(0);
 mean_p[i][th2][0]=fp[i][th2][0]->GetParameter(1);
 sig_p[i][th2][0]=fp[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fpi[%d][%d][0]",i,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(0);
 mean_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(1);
 sig_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fk[%d][%d][0]",i,th2),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th2][0]=fk[i][th2][0]->GetParameter(0);
 mean_k[i][th2][0]=fk[i][th2][0]->GetParameter(1);
 sig_k[i][th2][0]=fk[i][th2][0]->GetParameter(2);

 //------- Get Error Paramters ---//
 n_pi_err[i][th2][0]=fpi[i][th2][0]->GetParError(0); 
 n_k_err[i][th2][0]=fk[i][th2][0]->GetParError(0); 
 n_p_err[i][th2][0]=fp[i][th2][0]->GetParError(0); 

 //----- AC1 Coint Fitting ---------//

 fcoin[i][th2][0] =new TF1(Form("fcoin[%d][%d][0]",i,th2),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][th2][0]->SetNpx(2000);
 fcoin[i][th2][0]->SetTitle(Form("Coin-Time w AC cut  (AC1<%d ch && AC2>%d ch);Coin time [ns];Counts [1/56 ns]",th_ac1[i],ac2_adc));
 fcoin[i][th2][0]->SetParameters(n_pi[i][th2][0],mean_pi[i][th2][0],sig_pi[i][th2][0],n_k[i][th2][0],mean_k[i][th2][0],sig_k[i][th2][0],n_p[i][th2][0],mean_p[i][th2][0],sig_p[i][th2][0]);





 //------- AC2 -------------//

 fp[i][th1][1]=new TF1(Form("fp[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][th1][1]->SetNpx(2000);
 fpi[i][th1][1] =new TF1(Form("fpi[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fpi[i][th1][1]->SetNpx(2000);
 fk[i][th1][1]=new TF1(Form("fk[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][th1][1]->SetNpx(2000);
 
 fpi[i][th1][1]->FixParameter(1,def_mean_pi); 
 fpi[i][th1][1]->FixParameter(2,def_sig_pi);
 
 fk[i][th1][1]->FixParameter(1,def_mean_k); 
 // fk[i][th1][1]->FixParameter(2,def_sig_k); 
 fk[i][th1][1]->SetParLimits(2,0.75*def_sig_k,1.25*def_sig_k); 
 fp[i][th1][1]->FixParameter(1,def_mean_p); 
 fp[i][th1][1]->FixParameter(2,def_sig_p); 


  //----- AC2 Fitting -----------// 

 hcoin_ac2_p[i][th1]->Fit(Form("fp[%d][%d][1]",i,th1),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th1][1]=fp[i][th1][1]->GetParameter(0);
 mean_p[i][th1][1]=fp[i][th1][1]->GetParameter(1);
 sig_p[i][th1][1]=fp[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fpi[%d][%d][1]",i,th1),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(0);
 mean_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(1);
 sig_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fk[%d][%d][1]",i,th1),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th1][1]=fk[i][th1][1]->GetParameter(0);
 mean_k[i][th1][1]=fk[i][th1][1]->GetParameter(1);
 sig_k[i][th1][1]=fk[i][th1][1]->GetParameter(2);

 // n_k[i][th1][1]=hcoin_ac2_p[i][th1]->Integral(hcoin_ac2_p[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_ac2_p[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 //------- Get Error Paramters ---//
 n_pi_err[i][th2][1]=fpi[i][th1][1]->GetParError(0); 
 n_k_err[i][th2][1]=fk[i][th1][1]->GetParError(0); 
 n_p_err[i][th2][1]=fp[i][th1][1]->GetParError(0); 

 //----- AC2 Coint Fitting ---------//

 fcoin[i][th1][1] =new TF1(Form("fcoin[%d][%d][1]",i,th1),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][th1][1]->SetNpx(2000);
 fcoin[i][th1][1]->SetTitle(Form("Coin-Time w AC cut (AC1 Cut<%lf ch && AC2 Cut>%lf ch);Coin time [ns];Counts [1/56 ns]",ac1_adc,th_ac2[i]));
 fcoin[i][th1][1]->SetParameters(n_pi[i][th1][1],mean_pi[i][th1][1],sig_pi[i][th1][1],n_k[i][th1][1],mean_k[i][th1][1],sig_k[i][th1][1],n_p[i][th1][1],mean_p[i][th1][1],sig_p[i][th1][1]);

  
 sum_k_max=1250.;
  //---- AC1 ----//
 sum_k[i][th2][0]=n_k[i][th2][0]/tdc_time;
 //sum_k[0][th2][0]=1200; 
 
 sum_pi[i][th2][0]=n_pi[i][th2][0]/tdc_time;
 sum_p[i][th2][0]=n_p[i][th2][0]/tdc_time;
 //sum_acc[i][th2][0]=p0_acc[i][th2][0]*(6*sig_k[i][th2][0]/tdc_time);
 sum_acc[i][th2][0]=hcoin_ac1_acc_p[i][th2]->Integral(hcoin_ac1_acc_p[i][th2]->FindBin(-3*sig_k[i][th2][0]+mean_k[i][th2][0])
					              ,hcoin_ac1_acc_p[i][th2]->FindBin(+3*sig_k[i][th2][0]+mean_k[i][th2][0]));
 if( sum_acc[i][th2][0]<0.0)sum_acc[i][th2][0]=0.0;
 sum_k_err[i][th2][0]=n_k_err[i][th2][0]/tdc_time; 
 sum_pi_err[i][th2][0]=n_pi_err[i][th2][0]/tdc_time;
 sum_p_err[i][th2][0]=n_p_err[i][th2][0]/tdc_time;
 //---- AC2 ----//
 sum_k[i][th1][1]=n_k[i][th1][1]/tdc_time; 
 //sum_k[0][th1][1]=1200.;
 sum_pi[i][th1][1]=n_pi[i][th1][1]/tdc_time;
 sum_p[i][th1][1]=n_p[i][th1][1]/tdc_time;
 //sum_acc[i][th1][1]=p0_acc[i][th1][1]*(6*sig_k[i][th1][1]/tdc_time);
 sum_acc[i][th1][1]=hcoin_ac2_acc_p[i][th1]->Integral(hcoin_ac2_acc_p[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
			     ,hcoin_ac2_acc_p[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 if(sum_acc[i][th1][1]<0.0)sum_acc[i][th1][1]=0.0;
 sum_k_err[i][th1][1]=n_k_err[i][th1][1]/tdc_time; 
 sum_pi_err[i][th1][1]=n_pi_err[i][th1][1]/tdc_time;
 sum_p_err[i][th1][1]=n_p_err[i][th1][1]/tdc_time;
 
  
 //---- AC1 ----//
 if(sum_pi[i][th2][0]<1.05*sum_pi[0][th2][0]){
 gsum_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_pi[i][th2][0]);
 gsum_pi_ac1[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][0]);}
 if(sum_p[i][th2][0]<1.05*sum_p[0][th2][0]){
 gsum_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_p[i][th2][0]);
 gsum_p_ac1[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th2][0]);}
 if(sum_k[i][th2][0]<1.05*sum_k[0][th2][0]){
 gsum_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]);
 gsum_k_ac1[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th2][0]);
 gSN_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]/sum_acc[i][th2][0]);
 gSN_k_ac1[th1][th2]->SetPointError(i,emp[i],emp[i]);}
 //---- AC2 ----//
 if(sum_pi[i][th1][1]<1.05*sum_pi[0][th1][1]){
 gsum_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_pi[i][th1][1]);
 gsum_pi_ac2[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][1]);}
 if(sum_p[i][th1][1]<1.05*sum_p[0][th1][1]){
 gsum_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_p[i][th1][1]);
 gsum_p_ac2[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th1][1]);}
 if(sum_k[i][th1][1]<1.05*sum_k[0][th1][1]){
 gsum_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]);
 gsum_k_ac2[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th1][1]);
 gSN_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]/sum_acc[i][th1][1]);
 gSN_k_ac2[th1][th2]->SetPointError(i,emp[i],emp[i]);}
 


 //----------- Rate TGraph ---------------------//
  if(rate_k[i][th2][0]<1.1){
   rate_k[i][th2][0]=sum_k[i][th2][0]/sum_k[0][th2][0];
   //rate_k[i][th2][0]=sum_k[i][th2][0]/sum_k_max;
   rate_k_err[i][th2][0]=sqrt(1./sum_k[i][th2][0]*rate_k[i][th2][0]*(1.-rate_k[i][th2][0]));}
  else{rate_k[i][th2][0]=0.0;
       rate_k_err[i][th2][0]=0.0;}


 if(rate_k[i][th1][1]<1.1){
   rate_k[i][th1][1]=sum_k[i][th1][1]/sum_k[0][th1][1];
   //rate_k[i][th1][1]=sum_k[i][th1][1]/sum_k_max;
 rate_k_err[i][th1][1]=sqrt(1./sum_k[i][th1][1]*rate_k[i][th1][1]*(1.-rate_k[i][th1][1])); 
 }else{
 rate_k[i][th1][1]=0.0;
 rate_k_err[i][th1][1]=0.0;}

 rate_pi[i][th2][0]=sum_pi[i][th2][0]/sum_pi[0][th2][0];
 rate_p[i][th2][0]=sum_p[i][th2][0]/sum_p[0][th2][0];

 rate_pi[i][th1][1]=sum_pi[i][th1][1]/sum_pi[0][th1][1];
 rate_p[i][th1][1]=sum_p[i][th1][1]/sum_p[0][th1][1];

 rate_p_err[i][th2][0]=sqrt(1./sum_p[i][th2][0]*rate_p[i][th2][0]*(1.-rate_p[i][th2][0]));
 rate_pi_err[i][th2][0]=sqrt(1./sum_pi[i][th2][0]*rate_pi[i][th2][0]*(1.-rate_pi[i][th2][0]));

 rate_p_err[i][th1][1]=sqrt(1./sum_p[i][th1][1]*rate_p[i][th1][1]*(1.-rate_p[i][th1][1]));
 rate_pi_err[i][th1][1]=sqrt(1./sum_pi[i][th1][1]*rate_pi[i][th1][1]*(1.-rate_pi[i][th1][1]));
 //---- AC1 ----//
 if(0.0<rate_k[i][th2][0] && rate_k[i][th2][0]<1.05){
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_k[i][th2][0]);
 grate_k_ac1[th1][th2]->SetPointError(i,0,rate_k_err[i][th2][0]);
 }else{
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],0.0);
 grate_k_ac1[th1][th2]->SetPointError(i,0,0.0);
}

 if(rate_p[i][th2][0]<1.05){
 grate_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_p[i][th2][0]);
 grate_p_ac1[th1][th2]->SetPointError(i,0,rate_p_err[i][th2][0]);}
 if(rate_pi[i][th2][0]<1.05){
 grate_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_pi[i][th2][0]);
 grate_pi_ac1[th1][th2]->SetPointError(i,0,rate_pi_err[i][th2][0]);}



 //---- AC2 ----//

 if(0.0<rate_k[i][th1][1]&& rate_k[i][th1][1]<1.05){
  grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_k[i][th1][1]);
  grate_k_ac2[th1][th2]->SetPointError(i,0,rate_k_err[i][th1][1]);
 }else{
  grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],0);
  grate_k_ac2[th1][th2]->SetPointError(i,0,0);}

 if(rate_p[i][th1][1]<1.05){
  grate_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_p[i][th1][1]);
  grate_p_ac2[th1][th2]->SetPointError(i,0,rate_p_err[i][th1][1]);}
 if(rate_pi[i][th1][1]<1.05){  
 grate_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_pi[i][th1][1]);
 grate_pi_ac2[th1][th2]->SetPointError(i,0,rate_pi_err[i][th1][1]);}





//-----------------------------//

 if(n_k[i][th2][0]>max_nk[th1][th2][0])max_nk[th1][th2][0]=n_k[i][th2][0];   if(n_k[i][th1][1]>max_nk[th1][th2][1])max_nk[th1][th2][1]=n_k[i][th1][1];   
 if(n_p[i][th2][0]>max_np[th1][th2][0])max_np[th1][th2][0]=n_p[i][th2][0];   if(n_p[i][th1][1]>max_np[th1][th2][1])max_np[th1][th2][1]=n_p[i][th1][1];
 if(n_pi[i][th2][0]>max_npi[th1][th2][0])max_npi[th1][th2][0]=n_pi[i][th2][0];  if(n_pi[i][th1][1]>max_npi[th1][th2][1])max_npi[th1][th2][1]=n_pi[i][th1][1];
          

 if(max_SN_ac1[th2]<sum_k[i][th2][0]/sum_acc[i][th2][0] && sum_k[i][th2][0]/sum_acc[i][th2][0]<6.0){
   max_SN_ac1[th2]=sum_k[i][th2][0]/sum_acc[i][th2][0];
   SN_ac1[th2]=i;}

 if(max_SN_ac2[th1]<sum_k[i][th1][1]/sum_acc[i][th1][1] && sum_k[i][th1][1]/sum_acc[i][th1][1]<6.0 ){
   max_SN_ac2[th1]=sum_k[i][th1][1]/sum_acc[i][th1][1];
   SN_ac2[th1]=i;}
   

 // if(FOM_ac1[i][th2]<1.0)
  FOM_ac1[i][th2]=sqrt(pow(sum_k[i][th2][0],2)/sum_acc[i][th2][0]);
// if(FOM_ac2[i][th1]<1.0)
  FOM_ac2[i][th1]=sqrt(pow(sum_k[i][th1][1],2)/sum_acc[i][th1][1]);
  //    th_ac1[i]=th1_max-(th1_max-min_ac1)/iter_ac1*i;
 gfom_ac1[th2]->SetPoint(i,th_ac1[i],FOM_ac1[i][th2]);
 gfom_ac2[th1]->SetPoint(i,th_ac2[i],FOM_ac2[i][th1]);
 gfom->SetPoint(i,th_ac1[i],sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2)));
 hfom_ac[th1][th2]->Fill(th_ac1[i],th_ac2[i],sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2)));

  
 // hAC->Fill(th_ac1[i],ac2_adc[th2],FOM_ac1[i][th2]);
 hAC->Fill(th_ac1[i],ac2_adc[th2],Lfom[i][th2][0]);//Lambdda FOM
 hAC2->Fill(ac1_adc[1],th_ac2[i],Lfom[i][th1][0]);//Lambdda FOM
 if(max_fom_ac1<Lfom[i][th2][0]){
   max_fom_ac1=Lfom[i][th2][0];
   fom_th1=i;
  nLam_ac1=nL[i][th2][0];
  SNLam_ac1=nL[i][th2][0]/bgL_ac1[i][th2];
  fom_max_th2=th2;
 }
 if(max_fom_ac2<Lfom[i][th1][1]){
   max_fom_ac2=Lfom[i][th1][1];
   fom_th2=i;
  nLam_ac2=nL[i][th1][1];
  SNLam_ac2=nL[i][th1][1]/bgL_ac2[i][th1];
  fom_max_th1=th1;
}


 if(FOM_max_ac1[th2]<Lfom[i][th2][0]) FOM_max_ac1[th2]=Lfom[i][th2][0];
 if(FOM_max_ac2[th1]<Lfom[i][th1][1]) FOM_max_ac2[th1]=Lfom[i][th2][1];




if(fom_max<sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2))){
   fom_max=sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2));}

     }

     }
     }
     }

//change
//	string ofname_tmp = "../test.pdf";
//	c11=new TCanvas("c11","c11");
//	c11->cd();
//	//hmm_fom->Draw();
//	hcoin_ac1_all_p[0][0]->Draw();
//
//  cout<<"Print is starting"<<endl;
//  cout<<"pdf name: "<<ofname_tmp<<endl;
// c11->Print(Form("%s[",ofname_tmp.c_str()));
// c11->Print(Form("%s",ofname_tmp.c_str()));
//
//
////	hmm_fom->SetName("test");
////	hmm_fom->Write();
//	hcoin_ac1_all_p[0][0]->SetName("test");
//	hcoin_ac1_all_p[0][0]->Write();
}//end Fitting

//////////////////////////////////////////////////////////////////

void tuningAC::Tuning(){ 

  cout<<"==========================="<<endl;
  cout<<"========= Tuning =========="<<endl;
  cout<<"==========================="<<endl;
	//=========== FOMANA ========================//

	ac2_up=false;
        ac2_down=false;
	bool coin_flag=false;
   if(ac2_min)ac2_up=true;
   if(ac2_min==0)ac2_down=true;
	ev=0;

	bool mass_flag;
	bool z_cut;
	for(int k=0;k<ENum;k++){
	  T->GetEntry(k);
cout << " kkkkkkkkkk" << k << endl;


	  if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<ENum<<endl; 
     ev=ev+1;
   }

   mm_c=-100.0; //Initialization
   ct_c=-100;
   ac2_flag = false;
   coin_flag = false;
   z_cut=false;
   
   if(mode=="G"){
     coin_tc=ctime[0];
     if(coin_tc>-100)coin_trig=true;

   }else if(mode=="T"){
     coin_tc=tcoin_t;
  if(-80.<coin_tc && coin_tc<80.)coin_trig=true;
  }else if(mode=="H"){
     Ls2pads=(int)Ls2tpads[0];
     Rs2pads=(int)Rs2tpads[0];     
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
  if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true;  
   }

   
   if(mode=="H" || mode=="G"){
//     coin_tc=ctime[0];
     Rz=Rvz[0];
     Lz=Lvz[0];
    }
//     if(coin_tc>-100)coin_trig=true;

   
//cout << " Rz:Lz" << Rz << ":" << Lz << endl;
   if(Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax)z_cut=true;
   
//cout << " fom_th1:fom_th2" << fom_th1 << ":" << fom_th2 << endl;
//change  
//  th_ac2_t=10000.; // FIll Cut of hmm_ac1 (Upper cut)
   if(ac2_up && z_cut && Ra2sum<th_ac2_t)ac2_flag=true;

cout << " Ra2sum===" << Ra2sum << endl;
   if(ac2_down && z_cut && Ra2sum>th_ac2_b)ac2_flag=true;
//   if(ac2_up && z_cut
//      && Ra1sum<th_ac1[fom_th1] &&th_ac2[fom_th2]<Ra2sum && Ra2sum<th_ac2_t)ac2_flag=true;
//
//   if(ac2_down && z_cut
//       && Ra1sum<th_ac1[fom_th1] &&th_ac2[fom_th2]>Ra2sum && Ra2sum>th_ac2_b)ac2_flag=true;

   
   if(z_cut && coin_trig)ct_c=coin_tc;
   if(z_cut || ac2_flag){
cout << "ac2_flag" << ac2_flag <<endl;
cout << "z_cut" << z_cut <<endl;
cout << "coin_trig" << coin_trig <<endl;
   }
   if(ac2_flag && z_cut && coin_trig){
	    hcoin_fom->Fill(coin_tc);
   	    if(-1.0<coin_tc && coin_tc<1.0){
	      hmm_fom->Fill(mm);
	      mm_c=mm;
	    }
     //------- ACC BG -------------------------------------------------------------//     
     //------- ACC BG -------------------------------------------------------------//     
  if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
cout << "ACC BG" << k <<endl;
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc->Fill(ct);
		 hmm_fom_acc->Fill(mm);		
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	}
  }
    //----------------------------------------------------------------------------//	
   }

//cout << " kkk===" << k << endl;
   
   tnew->Fill();

   
	}



cout << "After fill" << endl;

	

	//========== w/o AC tuning hist =============//
	
	hmm_acc->Scale(20./100.);
	hmm_p->Add(hmm,hmm_acc,1.0,-1.0);

cout << "aaaaaaaaaaaaa" << endl;
cout << "w/o AC" << endl;
	
  //======== Lambda B.G. ============//
  fbg_L=new TF1("fbg_L","pol1(0)",1.02,1.3);  
  //  fbg_L->FixParameter(0,-1137.13); 
  //  fbg_L->FixParameter(1,1071.8);
  fbg_L->SetParameter(0,-1137.13); 
  fbg_L->SetParameter(1,1071.8);

  hmm_p->Fit("fbg_L","","",1.09,1.13);
  pbg[0]=fbg_L->GetParameter(0);
  pbg[1]=fbg_L->GetParameter(1);


cout << "aaaaaaaaaaaaa" << endl;
cout << "Lambda BG" << endl;
  //======== Sigma B.G. ==============//
  fbg_S=new TF1("fbg_S","pol1(0)",1.17,1.26);  
  //  fbg_S->FixParameter(0,448.634); 
  //  fbg_S->FixParameter(1,-352.65);
  fbg_S->SetParameter(0,448.634); 
  fbg_S->SetParameter(1,-352.65);
  hmm_p->Fit("fbg_S","","",1.17,1.26);
  pbg_S[0]=fbg_S->GetParameter(0);
  pbg_S[1]=fbg_S->GetParameter(1);

  
cout << "aaaaaaaaaaaaa" << endl;
cout << "Sigma BG" << endl;


  //======== Lambda Peak ============//
  double min_L=1.105;
  double max_L=1.13;

  fL_p->SetParameter(0,3.07090e+00); 
  fL_p->SetParameter(1,1.11546e+00);
  fL_p->SetParameter(2,3.6814e-03);
  hmm_p->Fit("fL_p","","",min_L,max_L);
  pL[0]=fL_p->GetParameter(0);
  pL[1]=fL_p->GetParameter(1); 
  pL[2]=fL_p->GetParameter(2);  

  pL_err[0]=fL_p->GetParError(0);
  pL_err[1]=fL_p->GetParError(1); 
  pL_err[2]=fL_p->GetParError(2);  

cout << "aaaaaaaaaaaaa" << endl;
cout << "Lambda peak" << endl;
  //======== Sigma Peak ============//
  double min_S=1.19;
  double max_S=1.21;


  fS_p->SetParameter(0,1.06770); 
  fS_p->SetParameter(1,1.19941);
  fS_p->SetParameter(2,4.9135e-3);
  hmm_p->Fit("fS_p","","",min_S,max_S);

  pS[0]=fS_p->GetParameter(0);
  pS[1]=fS_p->GetParameter(1); 
  pS[2]=fS_p->GetParameter(2);  

  pS_err[0]=fS_p->GetParError(0);
  pS_err[1]=fS_p->GetParError(1); 
  pS_err[2]=fS_p->GetParError(2);  


cout << "aaaaaaaaaaaaa" << endl;
cout << "Lambda peak" << endl;
  //======== Lambda Peak + B.G. =============//


  fL_all=new TF1("fL_all","pol1(0)+gausn(2)",min_L,max_L);
  fL_all->SetParameter(0,fbg_L->GetParameter(0)); 
  fL_all->SetParameter(1,fbg_L->GetParameter(1));
  fL_all->SetParameter(2,pL[0]); 
  fL_all->SetParameter(3,pL[1]);
  fL_all->SetParameter(4,pL[2]);

  hmm_p->Fit("fL_all","","",min_L,max_L);
  pbg[0]=fL_all->GetParameter(0);
  pbg[1]=fL_all->GetParameter(1);
  pL[0]=fL_all->GetParameter(2);
  pL[1]=fL_all->GetParameter(3); 
  pL[2]=fL_all->GetParameter(4);  
  pL_err[0]=fL_all->GetParError(2);
  pL_err[1]=fL_all->GetParError(3); 
  pL_err[2]=fL_all->GetParError(4);


 
  fL_p->SetParameters(pL[0],pL[1],pL[2]);
  fbg_L->SetParameters(pbg[0],pbg[1]);


  
  //======== Sigma Peak + B.G. =============//
  fS_all=new TF1("fS_all","pol1(0)+gausn(2)",min_S,max_S);
  fS_all->SetParameter(0,fbg_S->GetParameter(0)); 
  fS_all->SetParameter(1,fbg_S->GetParameter(1));
  fS_all->SetParameter(2,pS[0]); 
  fS_all->SetParameter(3,pS[1]);
  fS_all->SetParameter(4,pS[2]);

  hmm_p->Fit("fS_all","","",min_S,max_S);
  pbg_S[0]=fS_all->GetParameter(0);
  pbg_S[1]=fS_all->GetParameter(1);
  pS[0]=fS_all->GetParameter(2);
  pS[1]=fS_all->GetParameter(3); 
  pS[2]=fS_all->GetParameter(4);  
  pS_err[0]=fS_all->GetParError(2);
  pS_err[1]=fS_all->GetParError(3); 
  pS_err[2]=fS_all->GetParError(4);


 
  fS_p->SetParameters(pS[0],pS[1],pS[2]);
  fbg_S->SetParameters(pbg_S[0],pbg_S[1]);



  double bin_min_mmL,bin_max_mmL;
  bin_min_mmL=hmm->GetXaxis()->FindBin(pL[1]-3*pL[2]);
  bin_max_mmL=hmm->GetXaxis()->FindBin(pL[1]+3*pL[2]);
  sum_L=hmm->Integral(bin_min_mmL,bin_max_mmL);

  double bin_min_mmS,bin_max_mmS;
  bin_min_mmS=hmm->GetXaxis()->FindBin(pS[1]-3*pS[2]);
  bin_max_mmS=hmm->GetXaxis()->FindBin(pS[1]+3*pS[2]);
  sum_S=hmm->Integral(bin_min_mmS,bin_max_mmS);




  

	//======== w/ AC tuning hist ===============//

  
	hmm_fom_acc->Scale(2.0/100.);
	hmm_fom_p->Add(hmm_fom,hmm_fom_acc,1.0,-1.0);

	hmm_fom->Fit("fL_fom_bg","Rq","",1.08,1.15);
	hmm_fom->Fit("fS_fom_bg","Rq","",1.17,1.23);

	
	for(int i=0;i<3;i++){
	  Lbg_fom[i]=fL_fom_bg->GetParameter(i);
	  Sbg_fom[i]=fS_fom_bg->GetParameter(i);
	  fL_fom->FixParameter(i+3,Lbg_fom[i]);
	  fS_fom->FixParameter(i+3,Sbg_fom[i]);}

        fL_fom->SetParameter(0,1.0);
        fL_fom->SetParameter(1,1.115);
	fL_fom->SetParameter(2,1.5e-3);

        fS_fom->SetParLimits(0,0.0,10);
	fS_fom->SetParLimits(1,1.180,1.200);
	fS_fom->SetParameter(2,6.0e-3);

	hmm_fom->Fit("fL_fom","Rq","",1.1,1.15);
	hmm_fom->Fit("fS_fom","Rq","",1.17,1.23);

	NL_err=fL_fom->GetParError(0);
	NS_err=fS_fom->GetParError(0);

	for(int i=0;i<3;i++){
	 Lam_p[i]=fL_fom->GetParameter(i);
	 Sig_p[i]=fS_fom->GetParameter(i);
	 Lam_p_err[i]=fL_fom->GetParError(i);
	 Sig_p_err[i]=fS_fom->GetParError(i);
}



    all_L=hmm_fom->Integral(hmm_fom->FindBin(Lam_p[1]-3*Lam_p[2]),hmm_fom->FindBin(Lam_p[1]+3*Lam_p[2]));
    all_S=hmm_fom->Integral(hmm_fom->FindBin(Sig_p[1]-3*Sig_p[2]),hmm_fom->FindBin(Sig_p[1]+3*Sig_p[2]));


    bg_L=all_L-(Lam_p[0]/0.002);
    bg_S=all_S-(Sig_p[0]/0.002);



    //===== FOM Peak Fit =========//
    
	fL_fom_p->SetNpx(2000);	
	fL_fom_p->SetParameter(0,Lam_p[0]);
	fL_fom_p->SetParameter(1,Lam_p[1]);
	fL_fom_p->SetParameter(2,Lam_p[2]);

	fS_fom_p->SetNpx(2000);	
	fS_fom_p->SetParameter(0,Sig_p[0]);
	fS_fom_p->SetParameter(1,Sig_p[1]);
	fS_fom_p->SetParameter(2,Sig_p[2]);



	
	
        hcoin_acc->Scale(20./100.);
	hcoin_fom_p->Add(hcoin_fom,hcoin_acc,1.0,-1.0);
	set->SetTH1(hcoin_fom_p,"","","");


        fk_fom->SetNpx(2000);
        fk_fom->FixParameter(1,def_mean_k);
        fk_fom->FixParameter(2,def_sig_k);
	hcoin_fom_p->Fit("fk_fom","","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
	sumk_fom=fk_fom->GetParameter(0);
	meank_fom=fk_fom->GetParameter(1);
        sigk_fom=fk_fom->GetParameter(2);
	sumk_fom=sumk_fom/tdc_time;
        sk_fom=hcoin_fom_p->Integral(hcoin_fom_p->FindBin(-3*sigk_fom+meank_fom)
					    ,hcoin_fom_p->FindBin(3*sigk_fom+meank_fom));

        sk_fom_ct=hcoin_fom_p->Integral(hcoin_fom_p->FindBin(-1.0)
					    ,hcoin_fom_p->FindBin(1.0));


        nk_fom=hcoin_acc->Integral(hcoin_acc->FindBin(-3*sigk_fom+meank_fom)
					    ,hcoin_acc->FindBin(3*sigk_fom+meank_fom));

	snk_fom=sk_fom/nk_fom; 
	fom=sqrt(snk_fom*sumk_fom);



	cout<<"FOM ana is done"<<endl;

}

///////////////////////////////////////////////////////////////////////////////////
void tuningAC::Draw(){

//======== Draw TCanvas ==============//
  cout<<"start draw"<<endl;


 


  
      c11=new TCanvas("c11","c11");
      c11->cd();
      hcoin_fom->SetTitle("Coin time AC1 && AC2 Cut; coin time [ns];Counts/56 ns");
      hcoin_fom->SetTitleSize(0.06,"x");
      hcoin_fom->SetTitleSize(0.06,"y");
      hcoin_fom->SetLabelSize(0.04,"x");
      hcoin_fom->SetLabelSize(0.04,"y");
      hcoin_fom->GetXaxis()->SetTitleOffset(0.7);
      hcoin_fom->GetYaxis()->SetTitleOffset(0.7);
      hcoin_fom->GetXaxis()->CenterTitle();
      hcoin_fom->GetYaxis()->CenterTitle();
      hcoin_fom->SetLineColor(1);
      hcoin_fom->Draw();
      fk_fom->Draw("same");


      


 
      c14=new TCanvas("c14","c14");
      c14->cd(1);
      hmm_fom->Draw();
      fL_fom->SetLineColor(2);
      fS_fom->SetLineColor(4);
      fL_fom->Draw("same");
      fS_fom->Draw("same");

  cout<<"aaaaaaaaaaaaa"<<endl;

      int l=0;

   c15=new TCanvas("c15","c15");
   c15->Divide(4,4);
   for(int i=0;i<16;i++){
  cout<<"aaaaaaaaaaaaa"<< iter_ac1 << endl;
    if(iter_ac1<16) l=i;
    else if(iter_ac1<32)    l=2*i;
    else if(iter_ac1<48)    l=3*i;
    else l=4*i;
    if(iter_ac1<16) l=i;
         if(l>=iter_ac1) l=iter_ac1-1;
    c15->cd(i+1);

    hcoin_ac1_all_p[l][fom_max_th2]->Draw();
    fk[l][fom_max_th2][0]->SetLineColor(3);        
    fk[l][fom_max_th2][0]->SetFillStyle(3005);        
    fk[l][fom_max_th2][0]->SetFillColor(3);        
    fk[l][fom_max_th2][0]->Draw("same");    
   };

  cout<<"c15 is finished"<<endl;

   c16=new TCanvas("c16","c16");
   c16->Divide(4,4);
   for(int i=0;i<16;i++){
    if(iter_ac1<16) l=i;
    else if(iter_ac1<32)    l=2*i;
    else if(iter_ac1<48)    l=3*i;
    else l=4*i;
         if(l>=iter_ac1) l=iter_ac1-1;
    c16->cd(i+1);
    hcoin_ac2_all_p[l][fom_max_th1]->Draw();
    fk[l][fom_max_th1][1]->SetLineColor(3);        
    fk[l][fom_max_th1][1]->SetFillStyle(3005);        
    fk[l][fom_max_th1][1]->SetFillColor(3);            
    fk[l][fom_max_th1][1]->Draw("same");        
   };
   
  cout<<"c16 is finished"<<endl;

   gsum_k_ac1[0][0]->GetYaxis()->SetRangeUser(500.0,3000);
  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;
      c17=new TCanvas("c17","c17");
  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;
      c17->cd();
  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;
      for(int i=0;i<3;i++){
  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;
 gsum_k_ac1[i][i]->SetFillColor(i+1);
 gsum_k_ac1[i][i]->SetMarkerColor(i+1);
 gsum_k_ac1[i][i]->SetFillStyle(3005);

  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;
 if(i==0)gsum_k_ac1[i][i]->Draw("AP");
 else  gsum_k_ac1[i][i]->Draw("P");
  cout<<"aaaaaaaaaaaaaaaaaaaaaa"<<endl;

	   }

  cout<<"c17 is finished"<<endl;
      c18=new TCanvas("c18","c18");
      c18->cd();
      for(int i=0;i<3;i++){
 gsum_k_ac2[i][i]->SetFillColor(i+1);
 gsum_k_ac2[i][i]->SetMarkerColor(i+1);
 gsum_k_ac2[i][i]->SetFillStyle(3005);

 if(i==0)gsum_k_ac2[i][i]->Draw("AP");
 else  gsum_k_ac2[i][i]->Draw("P");
      }

      
  cout<<"c18 is finished"<<endl;
      
   c21=new TCanvas("c21","c21");
   c21->Divide(4,4);
   for(int i=0;i<16;i++){
    if(iter_ac1<16) l=i;
    else if(iter_ac1<32)    l=2*i;
    else if(iter_ac1<48)    l=3*i;
    else l=4*i;
         if(l>=iter_ac1) l=iter_ac1-1;
    c21->cd(i+1);
    hmm_ac1_all_p[l][fom_max_th2]->GetXaxis()->SetRangeUser(1.0,1.3);


    hmm_ac1_all_p[l][fom_max_th2]->Draw();
    //    fbg[l][fom_max_th2][0]->Draw("same");
    fLam[l][fom_max_th2][0]->SetLineColor(2);
    fLam[l][fom_max_th2][0]->SetRange(L1[l][fom_max_th2][0]-3*L2[l][fom_max_th2][0]
					  ,L1[l][fom_max_th2][0]+3*L2[l][fom_max_th2][0]);
    fLam[l][fom_max_th2][0]->Draw("same");
    fSig[l][fom_max_th2][0]->SetLineColor(4);
    fSig[l][fom_max_th2][0]->SetRange(S1[l][fom_max_th2][0]-3*S2[l][fom_max_th2][0]
					  ,S1[l][fom_max_th2][0]+3*S2[l][fom_max_th2][0]);    
    fSig[l][fom_max_th2][0]->Draw("same");

   };
      
  cout<<"c21 is finished"<<endl;


   c22=new TCanvas("c22","c22");
   c22->Divide(4,4);
   for(int i=0;i<16;i++){
     if(iter_ac1<16)l=i;
    else if(iter_ac1<32)    l=2*i;
    else if(iter_ac1<48)    l=3*i;
    else l=4*i;
     if(l>=iter_ac1) l=iter_ac1-1;
    c22->cd(i+1);
    hmm_ac2_all_p[l][fom_max_th1]->GetXaxis()->SetRangeUser(1.0,1.3);
    hmm_ac2_all_p[l][fom_max_th1]->Draw();
    //    fbg[l][fom_max_th1][1]->Draw("same");
    fLam[l][fom_max_th1][1]->SetLineColor(2);
    fLam[l][fom_max_th1][1]->SetRange(L1[l][fom_max_th1][1]-3*L2[l][fom_max_th1][1]
   				    ,L1[l][fom_max_th1][1]+3*L2[l][fom_max_th1][1]);    
    fLam[l][fom_max_th1][1]->Draw("same");
    fSig[l][fom_max_th1][1]->SetLineColor(4);
    fSig[l][fom_max_th1][1]->SetRange(S1[l][fom_max_th1][1]-3*S2[l][fom_max_th1][1]
   				    ,S1[l][fom_max_th1][1]+3*S2[l][fom_max_th1][1]);    
    fSig[l][fom_max_th1][1]->Draw("same");
   };

  cout<<"c22 is finished"<<endl;

   c23=new TCanvas("c23","c23");
   c23->Divide(2,1);
   c23->cd(1);
    hmm_ac1_all_p[iter_ac1-1][fom_max_th2]->Draw();
    fbg[iter_ac1-1][fom_max_th2][0]->Draw("same");
    fLam[iter_ac1-1][fom_max_th2][0]->SetRange(L1[iter_ac1-1][fom_max_th2][0]-3*L2[iter_ac1-1][fom_max_th2][0]
					      ,L1[iter_ac1-1][fom_max_th2][0]+3*L2[iter_ac1-1][fom_max_th2][0]);

    fLam[iter_ac1-1][fom_max_th2][0]->Draw("same");
    fSig[iter_ac1-1][fom_max_th2][0]->Draw("same");   
    c23->cd(2);
    hmm_ac2_all_p[iter_ac1-1][fom_max_th1]->Draw();
    fbg[iter_ac1-1][fom_max_th1][1]->Draw("same");
    fLam[iter_ac1-1][fom_max_th1][1]->SetRange(L1[iter_ac1-1][fom_max_th1][1]-3*L2[iter_ac1-1][fom_max_th1][1]
					      ,L1[iter_ac1-1][fom_max_th1][1]+3*L2[iter_ac1-1][fom_max_th1][1]);    
    fLam[iter_ac1-1][fom_max_th1][1]->Draw("same");
    fSig[iter_ac1-1][fom_max_th1][1]->Draw("same");
   
  cout<<"c23 is finished"<<endl;
    


      c24=new TCanvas("c24","c24");
      c24->cd(1);
 
      fL_fom->SetLineColor(2);
      fL_fom->Draw();
      hmm_fom->Draw("same"); 
  
  cout<<"c24 is finished"<<endl;


      c25=new TCanvas("c25","c25");
      c25->cd(1); 
      fS_fom->SetLineColor(2);
      fS_fom->Draw();
      hmm_fom->Draw("same"); 


  cout<<"c25 is finished"<<endl;
 
      c30=new TCanvas("c30","AC1 threshold");     
      c30->Divide(2,1);


	for(int i=0;i<nth;i++){
	  //	  c30->cd(i+1);
      
      gL_FOM_ac1[i]->SetFillColor(i+1);
      gL_FOM_ac1[i]->SetMarkerColor(i+1);
      gL_FOM_ac1[i]->SetFillStyle(3005);
      gL_FOM_ac2[i]->SetFillColor(i+1);
      gL_FOM_ac2[i]->SetMarkerColor(i+1);
      gL_FOM_ac2[i]->SetFillStyle(3005);
       if(i==0){
    c30->cd(1);
    gL_FOM_ac1[i]->Draw("AP");
    c30->cd(2);
    gL_FOM_ac2[i]->Draw("AP");
       }else{
	 c30->cd(1);
	 gL_FOM_ac1[i]->Draw("P");
	 c30->cd(2);
	 gL_FOM_ac2[i]->Draw("P");
       }
	}


  cout<<"c30 is finished"<<endl;
	
      c32=new TCanvas("c32","AC hist");
      c32->cd();
      hAC->Draw("colz");

  cout<<"c32 is finished"<<endl;

      c34=new TCanvas("c34","AC hist");
      c34->cd();
      hmm_ac1_all_p[fom_th1][fom_max_th2]->SetLineColor(1); hmm_ac2_all_p[fom_th2][fom_max_th1]->SetLineColor(4);
      hmm_fom->SetLineColor(2);      
      hmm_fom->Draw();
      hmm_ac1_all_p[fom_th1][fom_max_th2]->Draw("same"); 
      hmm_ac2_all_p[fom_th2][fom_max_th1]->Draw("same"); 


  cout<<"c34 is finished"<<endl;

       
      cout<<"drawing is done !!"<<endl;       



  
}


///////////////////////////////////////////////////////////////////


void tuningAC::Print(string ofname){
  cout<<"Print is starting"<<endl;
  cout<<"pdf name: "<<ofname<<endl;
 c11->Print(Form("%s[",ofname.c_str()));
 c11->Print(Form("%s",ofname.c_str()));
 c14->Print(Form("%s",ofname.c_str()));
 c15->Print(Form("%s",ofname.c_str()));
 c16->Print(Form("%s",ofname.c_str()));
 c17->Print(Form("%s",ofname.c_str()));
 c18->Print(Form("%s",ofname.c_str())); 
      cout<<"c18 is done"<<endl;       
 c21->Print(Form("%s",ofname.c_str()));  
 c22->Print(Form("%s",ofname.c_str())); 
 c23->Print(Form("%s",ofname.c_str()));
 c24->Print(Form("%s",ofname.c_str()));
 c25->Print(Form("%s",ofname.c_str())); 
 c30->Print(Form("%s",ofname.c_str()));
 c30->Print(Form("%s]",ofname.c_str()));
// c32->Print(Form("%s",ofname.c_str()));
// c32->Print(Form("%s]",ofname.c_str()));
 
    
 cout<<"Print is done "<<endl;
   


}


///////////////////////////////////////////////////////////////

void tuningAC::Write_coin(){


}

///////////////////////////////////////////////////////////////

void tuningAC::Write(){


 gL_ac1[fom_max_th2]->SetName(Form("gL_ac1_%d",fom_max_th2));
 gL_ac1[fom_max_th2]->Write();
 gS_ac1[fom_max_th2]->SetName(Form("gS_ac1_%d",fom_max_th2)); 
 gS_ac1[fom_max_th2]->Write();
 gL_FOM_ac1[fom_max_th2]->SetName(Form("gL_FOM_ac1_%d",fom_max_th2));
 gL_FOM_ac1[fom_max_th2]->Write();
 
 gL_ac2[fom_max_th1]->SetName(Form("gL_ac2_%d",fom_max_th1));  
 gL_ac2[fom_max_th1]->Write();
 gS_ac2[fom_max_th1]->SetName(Form("gS_ac2_%d",fom_max_th1)); 
 gS_ac2[fom_max_th1]->Write();  
 gL_FOM_ac2[fom_max_th1]->SetName(Form("gL_FOM_ac2_%d",fom_max_th1));
 gL_FOM_ac2[fom_max_th1]->Write();  

 for(int i=0;i<3;i++){
 gSN_k_ac1[i][i]->SetFillColor(i+1);
 gSN_k_ac1[i][i]->SetMarkerColor(i+1);
 gSN_k_ac1[i][i]->SetFillStyle(3005);
 gSN_k_ac2[i][i]->SetFillColor(i+1);
 gSN_k_ac2[i][i]->SetMarkerColor(i+1);
 gSN_k_ac2[i][i]->SetFillStyle(3005);
 // TGraphErrors* gsum_k_ac1[100][100];
 gsum_k_ac1[i][i]->SetFillColor(i+1);
 gsum_k_ac1[i][i]->SetMarkerColor(i+1);
 gsum_k_ac1[i][i]->SetFillStyle(3005);
 gsum_k_ac2[i][i]->SetFillColor(i+1);
 gsum_k_ac2[i][i]->SetMarkerColor(i+1);
 gsum_k_ac2[i][i]->SetFillStyle(3005); 
 grate_k_ac1[i][i]->SetFillColor(i+1);
 grate_k_ac1[i][i]->SetMarkerColor(i+1);
 grate_k_ac1[i][i]->SetFillStyle(3005);
 grate_k_ac2[i][i]->SetFillColor(i+1);
 grate_k_ac2[i][i]->SetMarkerColor(i+1);
 grate_k_ac2[i][i]->SetFillStyle(3005);
  //----- Pion -----------//
 grate_pi_ac1[i][i]->SetFillColor(i+1);
 grate_pi_ac1[i][i]->SetMarkerColor(i+1);
 grate_pi_ac1[i][i]->SetFillStyle(3005);
 grate_pi_ac2[i][i]->SetFillColor(i+1);
 grate_pi_ac2[i][i]->SetMarkerColor(i+1);
 grate_pi_ac2[i][i]->SetFillStyle(3005);  

  //----- Proton -----------//
 grate_p_ac1[i][i]->SetFillColor(i+1);
 grate_p_ac1[i][i]->SetMarkerColor(i+1);
 grate_p_ac1[i][i]->SetFillStyle(3005);
 grate_p_ac2[i][i]->SetFillColor(i+1);
 grate_p_ac2[i][i]->SetMarkerColor(i+1);
 grate_p_ac2[i][i]->SetFillStyle(3005);  
 
 
 gSN_k_ac1[i][i]->SetName(Form("gSN_k_ac1_%d",i));
 gSN_k_ac1[i][i]->Write(); 
 gSN_k_ac2[i][i]->SetName(Form("gSN_k_ac2_%d",i)); 
 gSN_k_ac2[i][i]->Write(); 
 grate_k_ac1[i][i]->SetName(Form("grate_k_ac1_%d",i));
 grate_k_ac1[i][i]->Write(); 
 grate_k_ac2[i][i]->SetName(Form("grate_k_ac2_%d",i)); 
 grate_k_ac2[i][i]->Write();
 gsum_k_ac1[i][i]->SetName(Form("gsum_k_ac1_%d",i));
 gsum_k_ac1[i][i]->Write(); 
 gsum_k_ac2[i][i]->SetName(Form("gsum_k_ac2_%d",i)); 
 gsum_k_ac2[i][i]->Write();

 //---- proton ----//
 grate_p_ac1[i][i]->SetName(Form("grate_p_ac1_%d",i));
 grate_p_ac1[i][i]->Write(); 
 grate_p_ac2[i][i]->SetName(Form("grate_p_ac2_%d",i)); 
 grate_p_ac2[i][i]->Write(); 

 //---- pion ----//
 grate_pi_ac1[i][i]->SetName(Form("grate_pi_ac1_%d",i));
 grate_pi_ac1[i][i]->Write(); 
 grate_pi_ac2[i][i]->SetName(Form("grate_pi_ac2_%d",i)); 
 grate_pi_ac2[i][i]->Write(); 
 
 
 hcoin_t1[i]->Write();
 hcoin_t2[i]->Write();
 
 }


 gSN_k_ac1[fom_max_th2][fom_max_th2]->SetName(Form("gSN_k_ac1_%d",fom_max_th2));
 gSN_k_ac1[fom_max_th2][fom_max_th2]->Write(); 
 gSN_k_ac2[fom_max_th1][fom_max_th1]->SetName(Form("gSN_k_ac2_%d",fom_max_th1)); 
 gSN_k_ac2[fom_max_th1][fom_max_th1]->Write(); 
 grate_k_ac1[fom_max_th2][fom_max_th2]->SetName(Form("grate_k_ac1_%d",fom_max_th2));
 grate_k_ac1[fom_max_th2][fom_max_th2]->Write(); 
 grate_k_ac2[fom_max_th1][fom_max_th1]->SetName(Form("grate_k_ac2_%d",fom_max_th1)); 
 grate_k_ac2[fom_max_th1][fom_max_th1]->Write(); 
 
 facc_kc->Write();
 fk_kc->Write();
 fpi_pic->Write();
 fp_pc->Write();
 set->SetTH1(hmm_ac1_all_p[fom_th1][fom_max_th2],"hmm_ac1_all_p","Mass [GeV]","Counts/2 MeV");
 hmm_ac1_all_p[fom_th1][fom_max_th2]->Write();
 set->SetTH1(hmm_ac2_all_p[fom_th2][fom_max_th1],"hmm_ac2_all_p","Mass [GeV]","Counts/2 MeV"); 
 hmm_ac2_all_p[fom_th2][fom_max_th1]->Write();
 hmm->Write();
 hmm_acc->Write(); 
 hmm_p->Write();
 fL_all->Write();
 fL_p->Write();
 fS_all->Write();
 fS_p->Write(); 
 
 hmm_fom->Write();
 hmm_fom_acc->Write();
 hmm_fom_p->Write();
 fL_fom->Write();
 fL_fom_p->Write();
 fS_fom->Write();
 fS_fom_p->Write(); 

 hmm_ac1_all_p[fom_th1][fom_max_th2]->Write();
 hmm_ac2_all_p[fom_th2][fom_max_th1]->Write();

 fLam[fom_th1][fom_max_th2][0]->Write();
 fSig[fom_th1][fom_max_th2][0]->Write();
 fLam_p->SetParameters(L0[fom_th1][fom_max_th2][0],L1[fom_th1][fom_max_th2][0],L2[fom_th1][fom_max_th2][0]);
 fSig_p->SetParameters(S0[fom_th1][fom_max_th2][0],S1[fom_th1][fom_max_th2][0],S2[fom_th1][fom_max_th2][0]);
 fLam_p->Write();
 fSig_p->Write(); 
 fk_fom->Write();
 
 tnew->Write();
 hAC->Write();
 hAC2->Write(); 
 hcoin_tc->Write();
 hcoin_fom->Write();
 hRu1_time_s->Write(); 
 hRu1_time_c->Write();

 
 fnew->Close();
}

///////////////////////////////////////////////////
///////////////////////////////


void tuningAC::Comment(){
   //======= COMMENT OUT ==================//
	cout<<"===================================="<<endl;
	cout<<"============= Comment =============="<<endl;
	cout<<"===================================="<<endl;	
	
	cout<<"==========AC1 analysis============="<<endl;
	cout<<Form("Lfom Max AC1[%d]:",fom_th1)<<max_fom_ac1<<endl;
	cout<<Form("N Lambda  Max AC1[%d]:",fom_th1)<<nLam_ac1<<endl;
	cout<<Form(" Lambda area Total Num  Max AC1[%d]:",fom_th1)<<totL_ac1[fom_th1][fom_max_th2]<<endl;
	cout<<Form(" Lambda BG  Max AC1[%d]:",fom_th1)<<bgL_ac1[fom_th1][fom_max_th2]<<endl;
	cout<<Form("SN Lambda  Max AC1[%d]:",fom_th1)<<SNLam_ac1<<endl;

	cout<<"==========AC2 analysis============="<<endl;
	cout<<Form("Lfom Max AC2[%d]:",fom_th2)<<max_fom_ac2<<endl;
        cout<<Form("Lfom Max AC2[%d]:",fom_th2)<<max_fom_ac2<<endl;
	cout<<Form("N Lambda  Max AC2[%d]:",fom_th2)<<nLam_ac2<<endl;
        cout<<Form("SN Lambda  Max AC2[%d]:",fom_th2)<<SNLam_ac2<<endl;
	
	cout<<"==========Tuning AC  analysis============="<<endl;
	cout<<"AC1 < "<<th_ac1[fom_th1]<<endl;
	if(ac2_min)cout<<th_ac2[fom_th2]<<" < AC2 < "<<th_ac2_t<<endl;
	if(ac2_min==0)cout<<th_ac2_b<<" < AC2 < "<<th_ac2[fom_th2]<<endl;
	cout<<"Lambda N:"<<Lam_p[0]/0.002<<endl;
	cout<<"Lambda N Error:"<<NL_err/0.002<<endl;
	cout<<"Lambda mean :"<<Lam_p[1]<<endl;
	cout<<"Lambda mean Error:"<<Lam_p_err[1]<<endl;
	cout<<"Lambda sigma:"<<Lam_p[2]<<endl;
	cout<<"Lambda sigma Error:"<<Lam_p_err[2]<<endl;
 	cout<<"BG in Lam: "<<bg_L<<endl; 
	cout<<"Lam S/N:"<<Lam_p[0]/0.002/bg_L<<endl;
	cout<<"Lam FOM: "<<sqrt(Lam_p[0]/0.002/(bg_L+Lam_p[0]/0.002)*Lam_p[0]/0.002)<<endl;

	cout<<"Sigma N:"<<Sig_p[0]/0.002<<endl;
	cout<<"Sigma N Error:"<<NL_err/0.002<<endl;
	cout<<"Sigma mean :"<<Sig_p[1]<<endl;
	cout<<"Sigma mean Error:"<<Sig_p_err[1]<<endl;
	cout<<"Sigma sigma:"<<Sig_p[2]<<endl;
	cout<<"Sigma sigma Error:"<<Sig_p_err[2]<<endl;
 	cout<<"BG in Sig: "<<bg_L<<endl; 
	cout<<"Sig S/N:"<<Sig_p[0]/0.002/bg_L<<endl;
	cout<<"Sig FOM: "<<sqrt(Sig_p[0]/0.002/(bg_L+Sig_p[0]/0.002)*Sig_p[0]/0.002)<<endl;
	//	sumk_fom=fk_fom->GetParameter(0);
	cout<<"Kaon Events: "<<sumk_fom<<endl;		
	cout<<"Kaon Efficiency: "<<(Lam_p[0]+Sig_p[0])/0.002/1250<<endl;

	cout<<"ac2_min: "<<ac2_min<<endl;


}



//==============================================//
//========== Defined Function ==================//
//=============================================//

double s2f1_off(int i,char* ARM,char* MODE, int KINE){


  double RS2_offset[16],LS2_offset[16];
  if(MODE=="H" && KINE==2){
 
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(MODE=="H" && KINE==1){
    
    //double  RS2_off_H1[16]={-16911.4,-16864.9,-16900,-16897.6,-16874.8,-16869.3,-16901.1,-16876.8,-16895.6,-16860.3,-16892.6,-16885,-16847.3,-16843.3,-16838.4,-16882.6};
    //double  LS2_off_H1[16]={-25336.9,-25385.7,-25367,-25392.2,-25391,-25386.3,-25422,-25428.9,-25415.2,-25425,-25438,-25381,-25394.4,-25417.5,-25432.8,-26082.1};

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,-16895.6,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}



#endif
