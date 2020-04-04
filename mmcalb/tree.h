#ifndef tree_h
#define tree_h 1
#include <iostream>
#include <fstream>
#include "Setting.h"
#include "define.h"

//#define max 100
const int Max=100;
class tree{

 public:
  tree();
  ~tree();
 public:

  void SetRun(string ifname);  
  void ChainTree(string ifname);
  void SetBranch();
  void NewBranch(string ofname, bool rarm);

  //----- ChainTree -------//
  TChain* T;
  int ent;
  int ENum;

  //------ SetBranch -------//
  
  double runnum;
  double hallap;
  double DRevtype;
  //----- Scintillation Triger Counters -------// 
  double RF1[Max],LF1[Max];
  double Rs0r_ac[Max],Rs0l_ac[Max],Ls0r_ac[Max],Ls0l_ac[Max];
  double Rs2r_ac[Max],Rs2l_ac[Max],Ls2r_ac[Max],Ls2l_ac[Max];
  double Rs0r_tc[Max],Rs0l_tc[Max],Ls0r_tc[Max],Ls0l_tc[Max];
  double Rs2r_tc[Max],Rs2l_tc[Max],Ls2r_tc[Max],Ls2l_tc[Max];
  double Rtrn,Ltrn;
  //----- PID Detectors --------------------//
  double Ra1t[Max],Ra1a[Max],Ra1a_p[Max],Ra1a_c[Max],Ra1sum;
  double Ra2t[Max],Ra2a[Max],Ra2a_p[Max],Ra2a_c[Max],Ra2sum;
  double La1t[Max],La1a[Max],La1a_p[Max],La1a_c[Max],La1sum;
  double La2t[Max],La2a[Max],La2a_p[Max],La2a_c[Max],La2sum;
  double Rgssum;
  double Lpssum;
  double Lshsum;
  double Lcersum;
  //----- at target ----------//
  double Rp[Max],Rpx[Max],Rpy[Max],Rpz[Max],Lp[Max],Lpx[Max],Lpy[Max],Lpz[Max];
  double Rth[Max],Rph[Max],Rx[Max],Ry[Max],Rz[Max];
  double Lth[Max],Lph[Max],Lx[Max],Ly[Max],Lz[Max];
  double R_Ras_x,L_Ras_x;
  //----- at FP -------------//
  double Rth_fp[Max],Rph_fp[Max],Rx_fp[Max],Ry_fp[Max];
  double Lth_fp[Max],Lph_fp[Max],Lx_fp[Max],Ly_fp[Max];
  
  double Rbeta[Max],Lbeta[Max];
  double rs2pathl[Max],rs0pathl[Max],rtrpathl[Max];
  double ls2pathl[Max],ls0pathl[Max],ltrpathl[Max];
  double trigger[100];
  double Rs2tpads[100],Ls2tpads[100];
  double Rs2trpad[100],Ls2trpad[100];

  //----- SetNewBranch -----//
  
  //  TFile* fnew;
  //  TTree* tnew;


  
};

#endif
