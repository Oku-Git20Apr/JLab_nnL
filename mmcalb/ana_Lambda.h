#ifndef ana_h
#define ana_h 1

using namespace std;

#include "Setting.h"
#include "Tree.h"
#include "ParamMan.h"
#include "define.h"
//#include "momcalib.h"

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

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
class ana : public Tree
{
  public:
    ana();
    ~ana();

  private:
    Setting *set;
    ParamMan *param;
    TFile *ofp;
    //    momcalib* mom;
    int GetMaxEvent()     { return ENumMax; }

  public:
    void ReadParam(string name);
    void Loop();
    void Draw();
    void Swich(bool nnL,bool scale);
    bool Close();
    void MakeHist();
    void SetRoot(string ifname);
    void SetRunList(string ifname);
    void SetMaxEvent( int N )  { ENumMax = N; }
    void Draw_mm();
    void matrix(string mtparam);
    void MTParam_R();
    void MTParam_L();
    void MTParam_G();
    void MTP_mom();
    void GetACParam();
    double  AC_npe(int nac, int seg, double adc);
    void Calib(int rt, int lt);
    double BG_Al(int events);
    void PathCalib(int rhit, int lhit);
    double Eloss(double yp,double z,char* arm);
    void CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
    double CoinCalc_gogami(int RS2_seg, int LS2_seg, int rhit, int lhit);
  //    double CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  double CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit);
 private:
    int ENumMax;
    int ENum;

  private:
// Lines, Textx
    TLine  *line;
    TLatex *text;

    TH2D *h_rbay_rbax, *h_rbby_rbbx;
    TH2D *h_rby_rbx;

    TTree* tree_out;
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
  double tdc_time;
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


};

#endif
