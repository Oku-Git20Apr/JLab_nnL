//const double c=299792458e-9;// [m/ns]
//const double mk=493.7e-3;// Kaon mass [GeV/c^2]
//const double me=0.511e-3;// electron mass [GeV/c^2] 
//const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
//const double mn=939.6e-3; // neutron mass [GeV/c^2]
//const double mpi=139.6e-3;// pion mass [GeV/c^2]
//#include <iostream>
//#include <fstream>
//#include <math.h>
//#include <string>
//#include <time.h>
//#include <stdio.h>
//#include <unistd.h>
//#include <sstream>
//using namespace std;
//#include "TApplication.h"
//#include "TH1F.h"
//#include "TH2F.h"
//#include "TF1.h"
//#include "TFile.h"
//#include "TLeaf.h"
//#include "TTree.h"
//#include "TCut.h"
//#include "TChain.h"
//#include "TCanvas.h"
//#include "TVector3.h"
//#include "TGraph.h"
//#include "TLine.h"
//#include "TLatex.h"
//#include "TText.h"
//#include "TStyle.h"
//#include "TROOT.h"
//#include "TGraphErrors.h"
//#include "TProfile.h"
//#include "TSystem.h"
//#include "TColor.h"
//#include "TPaveText.h"
//#include "TRandom.h"
#include "Setting.h"
#include "Param.h"
#include "ParamMan.h"
#include "Tree.h"
#include "define.h"
#include "TEfficiency.h"

struct TreeBranch{
  
  int z_cut,pid_cut,ct_cut;
  int nev,nrun;
  double missing_mass_MgL;
  double missing_mass_MgL_acc;
  double missing_mass, coin_time;
  double missing_mass_acc;
  double missing_mass_L;
  double missing_mass_nnL;
  double missing_mass_H3L;
  double missing_mass_cut;
  double missing_mass_Al;
  double missing_mass_Lb;
  double missing_mass_nnLb;
  double missing_mass_b;
  double missing_mass_Al_bg;
  double mm_tuned;
  double momR, momL;
  double momRz, momLz;
  double momRz_c, momLz_c;
  double zR, zL;
  double AC1_sum, AC2_sum;
  double AC1_npe_sum,AC2_npe_sum;
  double AC1_npe[24],AC2_npe[26];
  double yp_cor;
  double ctimecorR,ctimecorL;
  double ct_acc,ct_b,ct_c; 
  double ct_g,ct_gb;
  double Rs0ra_p,Rs0la_p,Rs0a_p;
  double Rs2ra_p[16],Rs2la_p[16],Rs2a_p[16];
  double Ls2ra_p[16],Ls2la_p[16],Ls2a_p[16];
  double Rvdc_u1,Rvdc_u2,Rvdc_v1,Rvdc_v2;
  double trig;
  double RXFP,RYFP,RXpFP,RYpFP;
  double RXt,RYt,RXpt,RYpt;
  double LXFP,LYFP,LXpFP,LYpFP;
  double LXt,LYt,LXpt,LYpt;
  double Lp[100],Rp[100],Bp;
  double Lp_c[100],Rp_c[100],Bp_c;  
  double dpe,dpe_[100],dpk[100];
  int Rs2_pad[100],Ls2_pad[100];
  double RS2T_F1[16],RS2B_F1[16],RS2T_ref,RS2B_ref,RS2T_F1_c[16],RS2B_F1_c[16],RS2T_F1_b[16],RS2B_F1_b[16];
  double LS2T_F1[16],LS2B_F1[16],LS2T_ref,LS2B_ref,LS2T_F1_c[16],LS2B_F1_c[16],LS2T_F1_b[16],LS2B_F1_b[16];
  double Rtof[100],Ltof[100];
  int ntrack_r,ntrack_l;
  double Rpathl,Lpathl,Rpathl_c,Lpathl_c;
  //int runnum;
};
static TreeBranch tr;
////////////////////////////////////
////SET PARAMETERS//////////////////
////////////////////////////////////
int nth=100;//change, nth=0 originally, max 99
char const* mode="H";
int kine=1;
//double tdc_time=56.23;
//double tdc_time=0.05623;//[ns/ch]
bool ac2_min = true;


class tuning : public Tree
{

 public:
  tuning();//constructer
  ~tuning();//destructer
  void SetRunList(string ifname);
	//pick up rootfiles from runlist in "ifname"
	
  void SetRun(string ifname);
	//pick up rootfiles from runlist in "ifname"
		//if single root file

  void SetRoot(string ifname);
	//if you want to save as a rootfile
		//you should use this with Write();

  void SetBranch();
	//This is my old version. 	
	//Now that I am using linked "Tree.cc", so this function has not been used anymore.

  void GetACParam();  
	//Reading AC parameters from offset_ac.dat,
	//in order to convert "ch" to "npe".

  void ReadParam(string pname);  
	//Reading parameters from "pname"
	
  void SetParam();  
	//Setting parameters such as temporary constants for tuning.

    void matrix(string mtparam);
    void MTParam_R();
    void MTParam_L();
    void MTParam_G();
    void MTP_mom();
	//Matrix elements are picked up form "mtparam"
	//These four MT functions are called in this matrix(mtparam) function.

  void MakeHist();
	//Define histograms and functions. (eg. TH1D, TF1, etc.)
	//This is also used to fill into a new tree named "tree_out".
	
  void Filling();
	//Loop analysis
	//Using GetEntry();
	
	
  void Fitting();
	//After making Cointime and MM histograms,
	//we fit those to evaluate something
	
  void ACtune(); 
	//Tuning of AC cut conditions.
	//I choose the best cut using Lambda/Sigma peak.
	
  void Draw();
  void Print(string ofname);
  void Write();
  void Write_coin();  
  void Comment();
    double  AC_npe(int nac, int seg, double adc);
    void Calib(int rt, int lt);
    double BG_Al(int events);
    void PathCalib(int rhit, int lhit);
    double Eloss(double yp,double z,char const* arm);
    void CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
    double CoinCalc_gogami(int RS2_seg, int LS2_seg, int rhit, int lhit);
  //    double CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  //double CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit);
 private:
  TFile* fnew;
  Setting* set;
  ParamMan *param;

 public:
    TH2D *h_rbay_rbax, *h_rbby_rbbx;
    TH2D *h_rby_rbx;
    TH3D *h3_fom;
    TH2D *hProjectionxy, *hProjectionyz, *hProjectionzx;
  TTree* tree_out;

////////////////////////////////////
////Define beforehand///////////////
////////////////////////////////////

  //=== SetRunList ====//
//  TChain* T;
//  int ENum;
//
////  //=== SetRoot =======//
  TTree* tnew;
  double mm_ac1[100][100];
  double mm_ac2[100][100];
  double fom_ac1[100];
  double fom_ac2[100];
  double mm_c;//  
  double ct_c;//these are also used in ACtune()
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
 //double ct;
 int ev;
//
// //===== SetParam ======//
//
//
 double ac1_adc[100],ac2l_adc[100],ac2u_adc[100],zver[100],zver_diff[100];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;
 double th1_max,th2_max,th2_min;
 double ac1_kcut,ac2_kcut_min,ac2_kcut_max;
 double th_ac2_t,th_ac2_b; 
//
// 
// //===== Make Hist =======//
 TH1F* hmm;
 TH1F* hmm_acc;
 TH1F* hmm_p;
 TH1F* npe_sum_a1;
 TH1F* npe_sum_a2;
 TH1F* h_pisr1;
 TH1F* h_ksr1;
 TH1F* h_psr1;
 TH1F* h_Lsr1;
 TH1F* h_Ssr1;
 TH1F* h_pisr2l;
 TH1F* h_ksr2l;
 TH1F* h_psr2l;
 TH1F* h_Lsr2l;
 TH1F* h_Ssr2l;
 TH1F* h_pisr2u;
 TH1F* h_ksr2u;
 TH1F* h_psr2u;
 TH1F* h_Lsr2u;
 TH1F* h_Ssr2u;
 TH1F* h_pitot1;
 TH1F* h_ktot1;
 TH1F* h_ptot1;
 TH1F* h_Ltot1;
 TH1F* h_Stot1;
 TH1F* h_pitot2l;
 TH1F* h_ktot2l;
 TH1F* h_ptot2l;
 TH1F* h_Ltot2l;
 TH1F* h_Stot2l;
 TH1F* h_pitot2u;
 TH1F* h_ktot2u;
 TH1F* h_ptot2u;
 TH1F* h_Ltot2u;
 TH1F* h_Stot2u;
 TH2F* h_pisr2d;
 TH2F* h_ksr2d;
 TH2F* h_psr2d;
 TH2F* h_Lsr2d;
 TH2F* h_Ssr2d;
 TH1F* h_ksr11;
 TH1F* h_ksr111;
 TH1F* h_ksr1111;
 TH1F* h_ktot11;
 TH1F* h_ktot111;
 TH1F* h_ktot1111;
// TH1F* hRu1_time_c;
// TH1F* hRu1_time_s; 
// TH2F* hcoin_ac1[100];
// TH2F* hcoin_ac2[100];
// TH2F* hcoin_ac1_acc[100];
// TH2F* hcoin_ac2_acc[100];
// TH2F* hvdc1_ac1[100];
// TH2F* hvdc2_ac1[100];
// TH2F* hvdc1_ac2[100];
// TH2F* hvdc2_ac2[100];
// TH2F* hs0_ac1[100];
// TH2F* hs2_ac1[100];
// TH2F* hs0_ac2[100];
// TH2F* hs2_ac2[100];
// TH1F* hcoin_t1[100];
 TH1F* hcoin_t2[100];
//------KID---------//
 TH1F* hcoin_k_ac1[100];
 TH1F* hcoin_k_ac2[100];
 TH1F* hcoin_k_fom[100][100][100];
 TH1F* hcoin_bg_fom[100][100][100];
 TH1F* hcoin_wo_bg_fom[100][100][100];
 TH1F* hcoin_k_fom_noAC1;
 TH1F* hcoin_bg_fom_noAC1;
 TH1F* hcoin_wo_bg_fom_noAC1;//SR(No AC1 Cut)
 TH1F* hcoin_k_fom_noAC2;
 TH1F* hcoin_bg_fom_noAC2;
 TH1F* hcoin_wo_bg_fom_noAC2;//SR(No AC2 Cut)
 TH1F* hcoin_k_fom_noAC;
 TH1F* hcoin_bg_fom_noAC;
 TH1F* hcoin_wo_bg_fom_noAC;//SR(No AC Cut, neither 1 nor 2)
 TH1F* hcoin_k_fom_noZ[100];
 TH1F* hcoin_bg_fom_noZ[100];
 TH1F* hcoin_wo_bg_fom_noZ[100];//SR(No AC Cut, neither 1 nor 2)
 TH1F* hmm_L_fom_noAC1;
 TH1F* hmm_L_fom_noAC2;
 TH1F* hmm_L_fom_noAC;
 TH1F* hmm_L_fom_noZ[100];
 TH1F* hmm_bg_fom_noAC1;
 TH1F* hmm_bg_fom_noAC2;
 TH1F* hmm_pibg_fom_noAC;
 TH1F* hmm_pibg_fom_noZ[100];
 TH1F* hmm_wo_bg_fom_noAC1;
 TH1F* hmm_wo_bg_fom_noAC2;
 TH1F* hmm_wo_bg_fom_noAC;
 TH1F* hmm_wo_bg_fom_noZ[100];
 TH1F* hmm_pi_fom_noAC;//MM if pion
 TH1F* hmm_bg_fom_noAC;
 TH1F* hmm_pi_fom_noZ[100];//MM if pion
 TH1F* hmm_bg_fom_noZ[100];
 TH1F* hmm_pi_wobg_fom_noAC;
 TH1F* hmm_pi_fom_no[100];//MM if pion
 TH1F* hmm_bg_fom_no[100];
 TH1F* hmm_pi_wobg_fom_noZ[100];
 TH2D* h_m2_mm;
 TH2D* h_m2_ac;

 TF1* fcoin_noAC1;
 TF1* fcoin_noAC2;
 TF1* fcoin_noAC;
 TF1* fcoin_noZ[100];
 TF1* fpi_noAC1;
 TF1* fk_noAC1;
 TF1* fp_noAC1;
 TF1* fpi_noAC2;
 TF1* fk_noAC2;
 TF1* fp_noAC2;
 TF1* fpi_noAC;
 TF1* fk_noAC;
 TF1* fp_noAC;
 TF1* fpi_noZ[100];
 TF1* fk_noZ[100];
 TF1* fp_noZ[100];
 TF1* fmmbg_noAC1;
 TF1* fmmbg_noAC2;
 TF1* fmmbg_noAC;
 TF1* fmmbg_noZ[100];
 TF1* fL_noAC1;
 TF1* fS_noAC1;
 TF1* fmm_noAC1;
 TF1* fL_noAC2;
 TF1* fS_noAC2;
 TF1* fmm_noAC2;
 TF1* fL_noAC;
 TF1* fS_noAC;
 TF1* fmm_noAC;
 TF1* fL_noZ[100];
 TF1* fS_noZ[100];
 TF1* fmm_noZ[100];
 TH1F* hmm_L_fom[100][100][100];
 TH1F* hmm_bg_fom[100][100][100];
 TH1F* hmm_wo_bg_fom[100][100][100];
 TH1F* hmm_pi_fom[100][100][100];
 TH1F* hmm_pibg_fom[100][100][100];
 TH1F* hmm_pi_wobg_fom[100][100][100];
// TH1F* hcoin_ac1_max[100];
// TH1F* hcoin_ac2_max[100];
// TH1F* hcoin_t3[100][100];
// TH1F* hcoin_t;
   TH1F* hcoin_tc;
// TH1F* hcoin_acc_ac1[100];
// TH1F* hcoin_acc_ac2[100];
// TH2F* ha1_a2;
 TH1F* hcoin_k;
 TH1F* hcoin_pi;
 TH1F* hcoin_p;
// TH2F* hcoin_ac1_all;
// TH2F* hcoin_ac2_all;
// TH2F* hmm_ac1[100];
// TH2F* hmm_ac2[100];
// TH2F* hmm_ac1_acc[100];
// TH2F* hmm_ac2_acc[100];
// TH2F* hct_a1a_c[24];
// TH2F* hct_a2a_c[26]; 
// //----- Fill -----// 
// TH1D* hmm_ac1_p[100][100];
// TH1D* hmm_ac2_p[100][100];
 TF1* facc[100][100][2];
 TF1* fpi[100][100][100];
 TF1* fk[100][100][100];
 TF1* fp[100][100][100];
 TF1* fcoin[100][100][100];
 TF1* fL[100][100][100];
 TF1* fS[100][100][100];
 TF1* fmm[100][100][100];
 TF1* fmmbg[100][100][100];
// TF1* fbg[100][100][2];
// TF1* fbg_s[100][100][2];
// TF1* fLam[100][100][2];
// TF1* fSig[100][100][2];
// TF1* fLam_p;
// TF1* fSig_p;   
// TH2F* hfom_ac[100][100];
// TH2F* hAC;
// TH2F* hAC2;
// 
 int iter_ac1;
 int iter_ac2;
// int iter_max;
// TH1D* hcoin_ac1_p[100][100];
// TH1D* hcoin_ac2_p[100][100]; 
// TH1D* hcoin_ac1_all_p[100][100];
// TH1D* hcoin_ac2_all_p[100][100]; 
// TH1D* hcoin_ac1_acc_p[100][100];
// TH1D* hcoin_ac2_acc_p[100][100]; 
// TH1D* hmm_ac1_all_p[100][100];
// TH1D* hmm_ac2_all_p[100][100];
// TH1D* hmm_ac1_acc_p[100][100];
// TH1D* hmm_ac2_acc_p[100][100];
//
// TGraphErrors* gsum_pi_ac1[100][100];
// TGraphErrors* gsum_p_ac1[100][100];
// TGraphErrors* gsum_k_ac1[100][100];
// TGraphErrors* grate_k_ac1[100][100];
// TGraphErrors* grate_p_ac1[100][100];
// TGraphErrors* grate_pi_ac1[100][100];
// TGraphErrors* gsum_pi_ac2[100][100];
// TGraphErrors* gsum_p_ac2[100][100];
// TGraphErrors* gsum_k_ac2[100][100];
// TGraphErrors* grate_k_ac2[100][100];
// TGraphErrors* grate_pi_ac2[100][100];
// TGraphErrors* grate_p_ac2[100][100];
// TGraphErrors* gSN_k_ac1[100][100];
// TGraphErrors* gSN_k_ac2[100][100];
// TGraphErrors* gfom_ac1[100];
// TGraphErrors* gfom_ac2[100];
// TGraphErrors* gfom;
// TGraphErrors* gmm_SN_ac1[100];
// TGraphErrors* gmm_SN_ac2[100];
// TGraphErrors* gmm_S_ac1[100];
// TGraphErrors* gmm_S_ac2[100];
// TGraphErrors* gmm_ac2[100]; 
// TGraphErrors* gmm_ac1[100]; 
// TGraphErrors* gL_ac1[100];
// TGraphErrors* gL_ac2[100];  
// TGraphErrors* gS_ac1[100];
// TGraphErrors* gS_ac2[100];  
// TGraphErrors* gL_eff_ac1[100];
// TGraphErrors* gL_eff_ac2[100];  
// TGraphErrors* gS_eff_ac1[100];
// TGraphErrors* gS_eff_ac2[100];  
// TGraphErrors* gS_SN_ac1[100];
// TGraphErrors* gS_SN_ac2[100];  
// TGraphErrors* gL_SN_ac1[100];
// TGraphErrors* gL_SN_ac2[100];  
// TGraphErrors* gL_N_ac1[100];
// TGraphErrors* gL_N_ac2[100];
// TGraphErrors* gL_FOM_ac1[100];    
// TGraphErrors* gL_FOM_ac2[100];  
 TGraphErrors* gcoin_pi_sr[100];  
 TGraphErrors* gcoin_k_sr[100];  
 TGraphErrors* gcoin_p_sr[100];  
//
//
// TF1* facc_t1def[100][100];
// TF1* fpi_t1def[100][100];
// TF1* fk_t1def[100][100];
// TF1* fcoin_t1def[100][100];
// TF1* fp_t1def[100][100];
// TF1* facc_t2def[100][100];
// TF1* fpi_t2def[100][100];
// TF1* fk_t2def[100][100];
// TF1* fcoin_t2def[100][100];
// TF1* fp_t2def[100][100];
// TF1* facc_t3def[100][100];
// TF1* fpi_t3def[100][100];
// TF1* fk_t3def[100][100];
// TF1* fcoin_t3def[100][100];
// TF1* fp_t3def[100][100];
// TF1* fcoin_t1[100][100];
// TF1* fcoin_t2[100][100];  
// TF1* fcoin_t3[100][100]; 
 TF1* facc_kc;
 TF1* fk_kc;
 TF1* fpi_pic;
 TF1* fp_pc;
 TF1* fac[100];
 TF1* fkk[100];
//
// //----- Tuning hist ----//
 TH1F* hcoin_fom;
 TH1F* hcoin_acc;
 TH1F* hmm_fom;//MM for FOM
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
//
// //===== Fill ========//
////from ana_Lambda.h line 116//////// 
//// LHRS ////
    TH1D *h_L_trig;
    TH1D* h_Rs2;
    TH1D* h_Ls2;
    TH1D *h_L_tr_n, *h_L_tr_ch2;
    TH1D *h_L_p, *h_L_pathl, *h_L_px, *h_L_py, *h_L_pz;
    TH1D *h_L_tgy, *h_L_tgth, *h_L_tgph;
    TH1D *h_L_vx, *h_L_vy, *h_L_vz;
    TH2D *h_L_y_x, *h_L_th_x, *h_L_ph_y;
    TH2D *h_L_tgph_tgth;

    TH1D *h_L_beta, *h_L_m2;
    TH2D *h_L_beta_p , *h_L_beta_m2;
    TH2D *h_L_dedx_p, *h_L_dedx_m2;
    TH1D *h_L_s0_dedx;
    TH2D *h_L_s0_dedx_x, *h_L_s0_beta_x;
    TH1D *h_L_s2_pad;
    TH1D *h_L_s2_dedx;
    TH2D *h_L_s2_dedx_x, *h_L_s2_beta_x;
    TH2D *h_L_s2_dedx_pad, *h_L_s2_beta_pad;

    TH1D *h_L_tgt;
    TH2D *h_L_s2pad_tgt;
    TH2D *h_L_p_tgt, *h_L_pathl_tgt, *h_L_tgy_tgt, *h_L_tgth_tgt, *h_L_tgph_tgt;
    TH2D *h_L_x_tgt, *h_L_y_tgt;

//// RHRS ////
    TH1D *h_R_trig;

    
    TH1D *h_R_tr_n, *h_R_tr_ch2;
    TH1D *h_R_p, *h_R_pathl, *h_R_px, *h_R_py, *h_R_pz;
    TH1D *h_R_tgy, *h_R_tgth, *h_R_tgph;
    TH1D *h_R_vx, *h_R_vy, *h_R_vz;
    TH2D *h_R_y_x, *h_R_th_x, *h_R_ph_y;
    TH2D *h_R_tgph_tgth;

    TH1D *h_R_beta, *h_R_m2;
    TH2D *h_R_beta_p , *h_R_beta_m2;
    TH2D *h_R_dedx_p, *h_R_dedx_m2;
    TH1D *h_R_s0_dedx;
    TH2D *h_R_s0_dedx_x, *h_R_s0_beta_x;
    TH1D *h_R_s2_pad;
    TH1D *h_R_s2_dedx;
    TH2D *h_R_s2_dedx_x, *h_R_s2_beta_x;
    TH2D *h_R_s2_dedx_pad, *h_R_s2_beta_pad;
    TH1D *h_R_a1_sum, *h_R_a2_sum;
    TH2D *h_R_a1_sum_x, *h_R_a2_sum_x;
    TH2D *h_R_a1_sum_p, *h_R_a2_sum_p;
    TH2D *h_R_a1_sum_m2, *h_R_a2_sum_m2;

    TH1D *h_R_tgt;
    TH2D *h_R_s2pad_tgt;
    TH2D *h_R_p_tgt, *h_R_pathl_tgt, *h_R_tgy_tgt, *h_R_tgth_tgt, *h_R_tgph_tgt;
    TH2D *h_R_x_tgt, *h_R_y_tgt;

//// Coin ////
    TH1D *h_ct;
    TH1D *h_ct_wK, *h_ct_wK_z;
    TH2D* h_ct_Rp;
    TH1D *h_ct_wK_acc, *h_ct_wK_z_acc;
    TH2D *h_Ls2x_ct;
    TH2D *h_Rs2x_ct;
    TH2D *h_a1sum_ct, *h_a2sum_ct;
    TH1D *h_mm, *h_mmall, *h_mmfoil;
    TH1D *h_mmbg, *h_mmallbg, *h_mmfoilbg;
    TH2D *h_Lp_mm, *h_Ll_mm, *h_Ltgy_mm, *h_Ltgth_mm, *h_Ltgph_mm;
    TH2D *h_Lvx_mm, *h_Lvy_mm, *h_Lvz_mm;
    TH2D *h_Lx_mm, *h_Ly_mm, *h_Lth_mm, *h_Lph_mm;
    TH2D *h_Rp_mm, *h_Rl_mm, *h_Rtgy_mm, *h_Rtgth_mm, *h_Rtgph_mm;
    TH2D *h_Rvx_mm, *h_Rvy_mm, *h_Rvz_mm;
    TH2D *h_Rx_mm, *h_Ry_mm, *h_Rth_mm, *h_Rph_mm;
    TH2D *h_Rp_Lp;
    TH1D *h_mm_L;
    TH1D *h_mm_L_ec;
    TH1D *h_mm_nnL;
    TH1D *h_acc_L;
    TH1D *h_acc_nnL;
    TH1D *h_mm_H3L;
    TH1D *h_acc_H3L;
    TH1D *h_peak_H3L;  
    TH1D *h_mm_Al;
    TH1D *h_mm_Al_acc;
    TH1D *h_peak_Al;
    TH1D *h_mm_MgL;
    TH1D *h_mm_MgL_acc;
    TH1D *h_peak_MgL;
    TH1D *h_peak_L;
    TH1D *h_peak_nnL;
    TH1D *h_acc_Al;
    TH1D *h_mm_acc;
    TH1D *h_peak_mm;
    TH1D *h_mm_pi;
    TH1D* h_ct_wK_z_all;
    TH1D* h_ct_acc;
    TH1D* h_mm_Al_bg;
    /// Added by itabashi ///
    TH1D*h_Rz;
    TH1D*h_Rz_cut;
    TH1D*h_Rth;
    TH1D*h_Rph;
    TH1D*h_Rp;
    TH1D*h_Lz;
    TH1D*h_Lth;
    TH1D*h_Lph;
    TH1D*h_Lp;
    TH1D*h_Rz_c;
    TH1D*h_Rth_c;
    TH1D*h_Rph_c;
    TH1D*h_Rp_c;
    TH1D*h_Lz_c;
    TH1D*h_Lth_c;
    TH1D*h_Lph_c;
    TH1D*h_Lp_c;    

    TF1* fAl_R;

//----Survival Ratio-----//
	TEfficiency* pEff1=0;
	TEfficiency* pEff2=0;
	TEfficiency* pEff3=0;
	TEfficiency* pEff4=0;
	TEfficiency* pEff5=0;
	TEfficiency* pEff6=0;
	TEfficiency* pEff7=0;
	TEfficiency* pEff8=0;
	TEfficiency* pEff9=0;
	TEfficiency* pEff10=0;
	TEfficiency* pEff11=0;
	TEfficiency* pEff12=0;
	TEfficiency* pEff13=0;
	TEfficiency* pEff14=0;
	TEfficiency* pEff15=0;
	TEfficiency* pEff111=0;
	TEfficiency* pEff1111=0;

 private:
    double L_s0l_toff    , L_s0r_toff;
    double L_s2l_toff[16], L_s2r_toff[16];
    double R_s0l_toff    , R_s0r_toff;
    double R_s2l_toff[16], R_s2r_toff[16];

    double L_s0l_t    , L_s0r_t    , L_s0_t;
    double L_s2l_t[16], L_s2r_t[16], L_s2_t[16];
    double R_s0l_t    , R_s0r_t    , R_s0_t;
    double R_s2l_t[16], R_s2r_t[16], R_s2_t[16];
    double R_p, L_p, B_p;
public:
  double min_mm,max_mm;
  double min_Lp,max_Lp;
  double mt;
  double mh;
  int bin_mm;
  int bin_Lp=200;
  double coin_offset;
  string param_mt[100];
  bool MT_p[100];
  bool ploss;
  double tdc_time=0.05623;//[ns/ch]
  bool Lp_scale=false;
  bool nnL_flag=false;
  double ac1_off[24],ac1_1pe[24],ac2_off[26],ac2_1pe[26];
  double R_pathl,L_pathl;
  double R_pathl_c, L_pathl_c;
  double ct;
  double coin_shift;
  double R_pz,R_px,R_py;
  double L_pz,L_px,L_py;
  int count=0;
  double ac1_th[24]={5249.78, 5263.03, 5244.2, 5239.62, 5276.17,  5249.84,
		     5274.63, 5291.02, 5259.66, 5305.42, 5262.97, 5333.76, 5275.35, 5190.67, 5303.03, 5414.55, 5224.28, 5320.19, 5253.68, 5242.92, 5219.99, 5296.6, 5345.91, 5357.9};
// 
// //--- Coin Offset -----//
// double pathl_off,s2_offset,coin_offset;
 double pathl_off,s2_offset;
// //----- Cut Parameters ----------//
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
// //-------------------------------//
 bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig,cut_track,cut_s0;
//
//
// //===== Fitting =========//
//
//
// //--- Parameters -----//
//  double bg_min,bg_max;
//  double bgs_min,bgs_max;
//  double Lfom[100][100][2],Sfom[100][100][2];
//  double L0_err[100][100][2],L1_err[100][100][2],L2_err[100][100][2];
//  double S0_err[100][100][2],S1_err[100][100][2],S2_err[100][100][2];
//  double nL_err[100][100][2],nS_err[100][100][2];
//  double bgL_ac1[100][100], bgL_ac2[100][100],bgS_ac1[100][100], bgS_ac2[100][100];
//  double totL_ac1[100][100], totL_ac2[100][100],totS_ac1[100][100], totS_ac2[100][100];
//  double nL[100][100][2],sigL[100][100][2],meanL[100][100][2];
// double nS[100][100][2],sigS[100][100][2],meanS[100][100][2];
// double kmin[100][100][2],kmax[100][100][2];
// double inte_ktot[100][100][2], inte_ksig[100][100][2];
// double p0_acc[100][100][2], p1_acc[100][100][2];
 double n_p[100][100][100],sig_p[100][100][100],mean_p[100][100][100];
 double n_pi[100][100][100],sig_pi[100][100][100],mean_pi[100][100][100];
 double n_k[100][100][100],sig_k[100][100][100],mean_k[100][100][100];
 double n_pi_noAC1,sig_pi_noAC1,mean_pi_noAC1;
 double n_k_noAC1,sig_k_noAC1,mean_k_noAC1;
 double n_p_noAC1,sig_p_noAC1,mean_p_noAC1;
 double n_L_noAC1,sig_L_noAC1,mean_L_noAC1;
 double n_S_noAC1,sig_S_noAC1,mean_S_noAC1;
 double n_pi_noAC2,sig_pi_noAC2,mean_pi_noAC2;
 double n_k_noAC2,sig_k_noAC2,mean_k_noAC2;
 double n_p_noAC2,sig_p_noAC2,mean_p_noAC2;
 double n_L_noAC2,sig_L_noAC2,mean_L_noAC2;
 double n_S_noAC2,sig_S_noAC2,mean_S_noAC2;
 double n_pi_noAC,sig_pi_noAC,mean_pi_noAC;
 double n_k_noAC,sig_k_noAC,mean_k_noAC;
 double n_p_noAC,sig_p_noAC,mean_p_noAC;
 double n_L_noAC,sig_L_noAC,mean_L_noAC;
 double n_S_noAC,sig_S_noAC,mean_S_noAC;
 double n_pi_noZ[100],sig_pi_noZ[100],mean_pi_noZ[100];
 double n_k_noZ[100],sig_k_noZ[100],mean_k_noZ[100];
 double n_p_noZ[100],sig_p_noZ[100],mean_p_noZ[100];
 double n_L_noZ[100],sig_L_noZ[100],mean_L_noZ[100];
 double n_S_noZ[100],sig_S_noZ[100],mean_S_noZ[100];
 double n_L[100][100][100],sig_L[100][100][100],mean_L[100][100][100];
 double n_S[100][100][100],sig_S[100][100][100],mean_S[100][100][100];
// int bin_ac1_adc[100][100],bin_min_ac1,bin_max_ac1,bin_ac2_adc[100][100],bin_max_ac2,bin_min_ac2;
// double sum_k[100][100][2],sum_p[100][100][2],sum_pi[100][100][2]; 
// double sum_k_err[100][100][2],sum_p_err[100][100][2],sum_pi_err[100][100][2]; 
// double inte_acc[100][100][2];
// double th_ac1[100],th_ac2[100];
// int bin_th_ac1[100][100],bin_th_ac2[100][100]; 
// double nk[100][100][100][100][2],npi[100][100][100][100][2],np[100][100][100][100][2];
// double max_nk[100][100][2],max_npi[100][100][2],max_np[100][100][2];
 double n_p_err[100][100][100],n_pi_err[100][100][100],n_k_err[100][100][100];
 double n_p_err_noAC1,n_pi_err_noAC1,n_k_err_noAC1;
 double n_p_err_noAC2,n_pi_err_noAC2,n_k_err_noAC2;
 double n_p_err_noAC,n_pi_err_noAC,n_k_err_noAC;
 double n_L_err[100][100][100],n_S_err[100][100][100];
// double FOM_ac1[100][100],FOM_ac2[100][100];
// double max_fom_ac1,max_fom_ac2;
// int fom_th1,fom_th2;
// double nLam_ac1,nLam_ac2,SNLam_ac1,SNLam_ac2;
// int fom_max_th2,fom_max_th1;
// double FOM_max_ac1[100],FOM_max_ac2[100],FOM_th1[100],FOM_th2[100]; 
// 
 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_sig_k,def_mean_k,def_acc;
 double def_sig_L,def_mean_L,def_mean_S,def_sig_S;
 double def_num_k,def_num_p,def_num_pi,def_acc_k,def_acc_pi,def_acc_p;
 double signal[100][100][100], noise[100][100][100], fom_L[100][100][100];
 double fom_pi1[100], fom_pi2l[100], fom_pi2u[100];
 double fom_k1[100], fom_k2l[100], fom_k2u[100];
 double fom_p1[100], fom_p2l[100], fom_p2u[100];
 double err_fom_pi1[100], err_fom_pi2l[100], err_fom_pi2u[100];
 double err_fom_k1[100], err_fom_k2l[100], err_fom_k2u[100];
 double err_fom_p1[100], err_fom_p2l[100], err_fom_p2u[100];
//
// double def_t1_k[100][100],def_t1_pi[100][100],def_t1_p[100][100],def_t1_acc[100][100];
// double def_t1_k_err[100][100],def_t1_pi_err[100][100],def_t1_p_err[100][100],def_t1_acc_err[100][100];
// double t1sig_k[100][100],t1sig_p[100][100],t1sig_pi[100][100],t1mean_p[100][100],t1mean_k[100][100],t1mean_pi[100][100];
// double t1sum_k[100],t1sum_pi[100],t1sum_p[100];
// double t1sum_k_err[100],t1sum_pi_err[100],t1sum_p_err[100];
// double def_t2_k[100][100],def_t2_pi[100][100],def_t2_p[100][100],def_t2_acc[100][100];
// double t2sig_k[100][100],t2sig_p[100][100],t2sig_pi[100][100],t2mean_p[100][100],t2mean_k[100][100],t2mean_pi[100][100];
// double def_t2_k_err[100][100],def_t2_pi_err[100][100],def_t2_p_err[100][100],def_t2_acc_err[100][100];
// double t2sum_k[100],t2sum_pi[100],t2sum_p[100];
// double t2sum_k_err[100],t2sum_pi_err[100],t2sum_p_err[100];
// double def_t3_k[100][100],def_t3_pi[100][100],def_t3_p[100][100],def_t3_acc[100][100];
// double t3sig_k[100][100],t3sig_p[100][100],t3sig_pi[100][100],t3mean_p[100][100],t3mean_k[100][100],t3mean_pi[100][100];
// double t3sum_k[100][100],t3sum_pi[100][100],t3sum_p[100][100];
// double emp[100];
//
// double rate_k[100][100][2],rate_p[100][100][2],rate_pi[100][100][2];
// double rate_k_err[100][100][2],rate_p_err[100][100][2],rate_pi_err[100][100][2];
// double sum_acc[100][100][2];
// double max_SN_ac1[100],max_SN_ac2[100];
// int SN_ac1[100],SN_ac2[100];
// double bg_0[100][100][2],bg_1[100][100][2],bg_2[100][100][2];
// double bg_s0[100][100][2],bg_s1[100][100][2],bg_s2[100][100][2];
// double L0[100][100][2],L1[100][100][2],L2[100][100][2];
// double S0[100][100][2],S1[100][100][2],S2[100][100][2];
// double sum_k_max=1250.;
// double fom_max=0.0;
// 
//
//
// //====== Tuning ============//
//
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
//
//-- Integral Range Definition --//
double center_pi=0., range_pi=0.;
double center_k=0., range_k=0.;
double center_p=0., range_p=0.;
double center_L=0., range_L=0.;
double center_S=0., range_S=0.;


// //====== Draw ===========//
 TCanvas* c1;
 TCanvas* c2;
 TCanvas* c3;
 TCanvas* c4;
 TCanvas* c5;
 TCanvas* c6;
 TCanvas* c7;
 TCanvas* c8;
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
// TCanvas* c20;
// TCanvas* c21;
// TCanvas* c22;
// TCanvas* c23;
// TCanvas* c24;
// TCanvas* c25; 
// TCanvas* c30;
// TCanvas* c31;
// TCanvas* c32;
// TCanvas* c33;
// TCanvas* c34;
//
//
//   
// 
};




