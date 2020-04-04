#ifndef momcalib_h
#define momcalib_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"


bool RHRS;
bool LHRS;
const int MAX=1000;
int nite=0;


extern double s2f1_off(int i,char* ARM,char* MODE, int KINE);
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
//extern double tune(double* pa, int j, int angflag); 
extern  double Calc_ras(double a,double b,double c){return  a *b + c;};  
extern double calcf2t_ang(double* P,double xf, double xpf, double yf, double fpf,double z);
extern double calcf2t_zt(double* P, double xf, double xpf, double yf, double ypf);
extern double calcf2t_mom(double* P, double xf, double xpf, double yf, double ypf, double zt);

extern double calcf2t_mom_RL( double* RP, double Rxf, double Rxpf, double Ryf, double Rypf, double Rzt,
			      double* LP, double Lxf, double Lxpf, double Lyf, double Lypf, double Lzt	     );



//=================================================//
//============= MomCalib Class ====================//
//=================================================//

class momcalib : public tree
{

 public:
  momcalib();
  ~momcalib();

 public:
  //  void SetBranch(string ifname, bool rarm);
  //  void NewBranch(string ofname, bool rarm);
  void EventSelection();
  void SetRoot(string ifname);
  void SingleRoot(string ifname);
  void MakeHist();
  void NewRoot(string ofname);
  void MTParam(string mtparam);
  void MTParam_R();
  void MTParam_L();
  void MTP_mom();
  void ParamCorr();
  void zCorr(bool rarm, bool larm);
  void rasCorr(bool rarm, bool larm );
  void angCorr(bool rth, bool rph, bool lth, bool lph);
  void plossCorr(bool PLoss);
  void momCorr(bool rarm, bool larm);
  double SSx(double z, double th);
  double SSy(double z, double ph);  
  int mode(string arm,char* Target, int F1tdc);
  void MomTuning(string ofname);
  double Eloss(double xp,double z,  char* arm);
  double tune(double* pa, int j, int angflag);   

  void Fill();
  void Close();

  int MODE=5;

   Setting *set = new Setting();


   //------mode -------------//
   char* target;
   int tdc_mode=0;

   //------ MTParam --------//
  string param[10];
  bool MT[10];  
  bool ploss;
  //------ NewROOT -----------//
  TFile* ofp;
  TTree* tnew;




  //----- SSx && SSy ------//
  const double l0 = 100.3; // length target to Sieve slit
  double rssx,rssy,lssx,lssy,rssx_c,rssy_c,lssx_c,lssy_c;
  //----- MakeHist ---------//  
  TH1F* hLz;
  TH1F* hLth;
  TH1F* hLph;
  TH1F* hLz_c;
  TH1F* hLz_rc;
  TH1F* hLth_c;
  TH1F* hLph_c;
  TH1F* hRp;
  TH1F* hRp_c;
  TH1F* hRp_ec;  
  TH1F* hLp;
  TH1F* hLp_c;
  TH1F* hLp_ec;
  TH2F* hLss;
  TH2F* hLss_c;
  TH2F* hRss;
  TH2F* hRss_c;  
  //------------------//

  TH1F* hct; 
  TH1F* hLz_cut;
  TH1F* hRz_cut;  
  TH1F* hct_cut;
  TH1F* hac1_cut;
  TH1F* hac2_cut;
  TH1F* hcer_cut;


  //----- tuning check -----//
  TH1F* hmm_L_t[100];
  TH1F* hmm_S_t[100];
  TH1F* hmm_t[100];  
  TH1F* hmm_L_tc[100];
  TH1F* hmm_S_tc[100];
  TH1F* hmm_tc[100];  

  //------------------------//  
  TH1F* hRz;
  TH1F* hRth;
  TH1F* hRph;
  TH1F* hRz_c;
  TH1F* hRz_rc;  
  TH1F* hRth_c;
  TH1F* hRph_c;    
  TH1F* hmm_select;
  TH1F* hmm_L;
  TH1F* hmm_S;
  TH1F* hmm_cut;
  TH1F* hLam;
  TH1F* hSig;

  double min_z,max_z;
  int bin_z;
  double min_th,max_th;
  int bin_th;
  double min_ph,max_ph;
  int bin_ph;
  double min_lp,max_lp;
  int bin_lp;
  double min_rp,max_rp;
  int bin_rp;
  double min_bp,max_bp;
  int bin_bp;
  double min_ssx,max_ssx,min_ssy,max_ssy;
  int bin_ssx,bin_ssy;
  double min_mm,max_mm;
  int bin_mm;

  double min_ct,max_ct;
  int bin_ct;
  double tdc_time=0.056;//[ns]
  
  //------ MTtuing -----------//


  //------ EventSelection ---//
  div_t d;
  //--- Event Selection mass cut ---// 
  double mmL_range = 0.010; // [GeV/c^2]
  double mmS_range = 0.010; // [GeV/c^2]

  
    int nlam=0;
    int nsig=0;
    
 // Cut Parameters //
 const double coin_cutmin=-248;
 const double coin_cutmax=-244; 
 const double rpathl_cutmin=28.7;
 const double rpathl_cutmax=29.4;
 const double lpathl_cutmin=28.6;
 const double lpathl_cutmax=29.2;
 const double rbeta_cutmin=0.0;
 const double rbeta_cutmax=1.0;
 const double lbeta_cutmin=0.9;
 const double lbeta_cutmax=1.0;
 const double Rvz_cutmin=-0.1;
 const double Rvz_cutmax= 0.1;
 const double Lvz_cutmin=-0.1;
 const double Lvz_cutmax= 0.1;
 const double Rx_cutmin= -0.4;
 const double Rx_cutmax= 0.4;
 const double a1_th=100.0;
 const double a2_th=500.0; 
 // Event Flag //
    bool Rs0_flag;
    bool Ls0_flag;
    bool Rs2_flag;
    bool Ls2_flag;
    bool Scinti_flag;
    bool ac1_flag;
    bool ac2_flag;
    bool gs_flag;
    bool sh_flag;
    bool ps_flag;
    bool trig_flag;
    bool track_flag;
    bool Rpathl_flag;
    bool Lpathl_flag;
    bool coin_flag;
    bool x_flag;
    bool y_flag;
    bool z_flag;
    bool RPID_flag;
    bool LPID_flag;

    
 double mm;
 double mm_L;
 double mm_acc;
 double mm_pi;
 double ct;
 double ctime;
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
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
 //---- MomTuning ----//
 double chi_sq[1000],chi_sq1[100],chi_sq2[100];
 TGraphErrors* gchi_p  = new TGraphErrors();
 TGraphErrors* gchi_Rp = new TGraphErrors();
 TGraphErrors* gchi_Lp = new TGraphErrors();
 //----- ELoss ------//
 

 //----  ParamCorr -----//
 int nPC=0;
 bool PC_flag=false;
 
};

//=====================================//
//======== Momcalib ==================//
//=====================================//

momcalib::momcalib(){

 set->Initialize();
 
};
momcalib::~momcalib(){};


int momcalib::mode(string arm, char* Target, int F1tdc){
  cout<<endl;
  cout<<"======================================="<<endl;
  if(arm=="R"){
    MODE=-1;
  cout<<"====== RHRS momentum tuning ==========="<<endl;

  }else if(arm=="L"){
    MODE=1;
  cout<<"====== LHRS momentum tuning ==========="<<endl;    
  } else if(arm=="C"){
    MODE=0;
  cout<<"====== R & L momentum tuning ==========="<<endl;    
  } else{cout<<"Please Select Tuning Mode !"<<endl; exit(1);}

  cout<<"======================================="<<endl;

    target=Target;

  // tdc_time tdc_time=0.56;//[ns]
    cout<<"F1tdc_mode : "<<typeid(F1tdc).name()<<endl;
  if(F1tdc==1){tdc_time=0.056;// [ns]
    tdc_mode=1;
  }else if(F1tdc==2){tdc_time=0.058; //[ns]
    tdc_mode=2;
  }else{cout<<"Error F1tdc modes "<<F1tdc<<endl;exit(1);}
  
  cout<<"tdc resolution [ns]"<<tdc_time<<endl;
  cout<<"Target : "<<target<<endl;
  
  return MODE;
}



//-----------------------------//
//--------- SetRun  -----------//
//-----------------------------//

void momcalib::SingleRoot(string ifname){

  SetRun(ifname);
  SetBranch();

};



//-----------------------------//
//--------- SetRoot -----------//
//-----------------------------//

void momcalib::SetRoot(string ifname){

  ChainTree(ifname);
  SetBranch();

};

//-----------------------------//
//--------- NewRoot -----------//
//-----------------------------//

void momcalib::NewRoot(string ofname){

  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Momentum matrix tuning");
  tnew =T->CloneTree(0);
  
  tnew->Branch("mm", &mm,"mm/D");  
  tnew->Branch("mm_acc", &mm_acc,"mm_acc/D");  
  tnew->Branch("mm_L", &mm_L,"mm_L/D");  
  tnew->Branch("mm_pi", &mm_pi,"mm_pi/D");
  tnew->Branch("coin", &coin_t,"coin_t/D");  
  tnew->Branch("coin_acc", &ctime,"coin_acc/D");  

};

//=====================================================================//

void momcalib::MakeHist(){

  
  min_z=-0.2;
  max_z=0.2;
  bin_z=400;
  hLz=new TH1F("hLz","",bin_z,min_z,max_z);
  set->SetTH1(hLz,"LHRS z-vertex w/ correction","z [m]","Counts");
  hLz->SetLineColor(4);
  hLz->SetFillColor(4);  
  hLz->SetFillStyle(3002);
  hLz_c=new TH1F("hLz_c","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_c,"LHRS z-vertex z parameter tuing","z [m]","Counts");
  hLz_c->SetLineColor(8);
  hLz_c->SetFillColor(8);  
  hLz_c->SetFillStyle(3002);
  hLz_rc=new TH1F("hLz_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_rc,"LHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hLz_rc->SetLineColor(2);
  hLz_rc->SetFillColor(2);  
  hLz_rc->SetFillStyle(3002);
  
  hRz=new TH1F("hRz","",bin_z,min_z,max_z);
  set->SetTH1(hRz,"RHRS z-vertex ","z [m]","Counts");
  hRz->SetLineColor(4);
  hRz->SetFillColor(4);  
  hRz->SetFillStyle(3002);
  hRz_c=new TH1F("hRz_c","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_c,"RHRS z-vertex z parameter tuing","z [m]","Counts");
  hRz_c->SetLineColor(8);
  hRz_c->SetFillColor(8);  
  hRz_c->SetFillStyle(3002);
  hRz_rc=new TH1F("hRz_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_rc,"RHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hRz_rc->SetLineColor(2);
  hRz_rc->SetFillColor(2);  
  hRz_rc->SetFillStyle(3002);


  
  min_th=-0.1;
  max_th=0.1;
  bin_th=1000;
  
  hLth=new TH1F("hLth","",bin_th,min_th,max_th);
  set->SetTH1(hLth,"LHRS th-vertex ","theta [rad]","Counts");
  hLth->SetLineColor(4);
  hLth->SetFillColor(4);  
  hLth->SetFillStyle(3002);
  hLth_c=new TH1F("hLth_c","",bin_th,min_th,max_th);  
  set->SetTH1(hLth_c,"LHRS theta-vertex w/ parameter tuning","theta [rad]","Counts");
  hLth_c->SetLineColor(2);
  hLth_c->SetFillColor(2);  
  hLth_c->SetFillStyle(3002);
  hRth=new TH1F("hRth","",bin_th,min_th,max_th);
  set->SetTH1(hRth,"RHRS th-vertex ","theta [rad]","Counts");
  hRth->SetLineColor(4);
  hRth->SetFillColor(4);  
  hRth->SetFillStyle(3002);  
  hRth_c=new TH1F("hRth_c","",bin_th,min_th,max_th);  
  set->SetTH1(hRth_c,"RHRS th-vertex w/ paramater tuning","theta [rad]","Counts");  
  hRth_c->SetLineColor(2);
  hRth_c->SetFillColor(2);  
  hRth_c->SetFillStyle(3002);

  
  min_ph=-0.1;
  max_ph=0.1;
  bin_ph=1000;

  hLph=new TH1F("hLph","",bin_ph,min_ph,max_ph);
  set->SetTH1(hLph,"LHRS ph-vertex ","phi [rad]","Counts");
  hLph->SetLineColor(4);
  hLph->SetFillColor(4);  
  hLph->SetFillStyle(3002);  
  hLph_c=new TH1F("hLph_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hLph_c,"LHRS ph-vertex w/ parameter tuning","phi [rad]","Counts");
  hLph_c->SetLineColor(2);
  hLph_c->SetFillColor(2);  
  hLph_c->SetFillStyle(3002);
  hRph=new TH1F("hRph","",bin_ph,min_ph,max_ph);
  set->SetTH1(hRph,"RHRS ph-vertex ","phi [rad]","Counts");
  hRph->SetLineColor(4);
  hRph->SetFillColor(4);  
  hRph->SetFillStyle(3002);  
  hRph_c=new TH1F("hRph_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hRph_c,"RHRS ph-vertex parameter tuning","phi [rad]","Counts");  
  hRph_c->SetLineColor(2);
  hRph_c->SetFillColor(2);  
  hRph_c->SetFillStyle(3002);


  
  min_ssx = -8.0; //[cm]
  max_ssx =  8.0; //[cm]
  min_ssy =  -5.0; //[cm]
  max_ssy =  5.0; //[cm]
  bin_ssx =(int)((max_ssx - min_ssx)/0.05);
  bin_ssy =(int)((max_ssy - min_ssy)/0.05);

  hLss=new TH2F("hLss","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss,"Sieve Slit pattern in LHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hLss_c=new TH2F("hLss_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss_c,"Sieve Slit pattern in LHRS w/  parameter tuning","ssy [cm]","ssx [cm]");

  hRss=new TH2F("hRss","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss,"Sieve Slit pattern in RHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hRss_c=new TH2F("hRss_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss_c,"Sieve Slit pattern in RHRS w/  parameter tuning","ssy [cm]","ssx [cm]");
  
  
  min_lp = 2.0;
  max_lp = 2.3;
  bin_lp =(int)((max_lp - min_lp)/10000.);

  hLp=new TH1F("hLp","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLph,"LHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp->SetLineColor(4);
  hLp->SetFillColor(4);  
  hLp->SetFillStyle(3002);
  hLp_c=new TH1F("hLp_c","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLph_c,"LHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp_c->SetLineColor(8);
  hLp_c->SetFillColor(8);  
  hLp_c->SetFillStyle(3002);
  hLp_ec=new TH1F("hLp_ec","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp_ec,"LHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hLp_ec->SetLineColor(2);
  hLp_ec->SetFillColor(2);  
  hLp_ec->SetFillStyle(3002);
  
  min_rp = 1.7;
  max_rp = 1.9;
  bin_rp =(int)((max_rp - min_rp)/1000);

  hRp=new TH1F("hRp","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRph,"RHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp->SetLineColor(4);
  hRp->SetFillColor(4);  
  hRp->SetFillStyle(3002);
  hRp_c=new TH1F("hRp_c","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRph_c,"RHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp_c->SetLineColor(8);
  hRp_c->SetFillColor(8);  
  hRp_c->SetFillStyle(3002);
  hRp_ec=new TH1F("hRp_ec","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRp_ec,"RHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hRp_ec->SetLineColor(2);
  hRp_ec->SetFillColor(2);  
  hRp_ec->SetFillStyle(3002);  

  
  min_ct=-50.0;
  max_ct=50.0;
  bin_ct=(double)((max_ct-min_ct)/tdc_time);




  
  hct=new TH1F("hct","",bin_ct,min_ct,max_ct);
  set->SetTH1(hct,"Coincidence time","coin-time [ns]","Counts");  
  

  //----- Cut hist --------//
  hLz_cut=new TH1F("hLz_cut","",bin_z,min_z,max_z);
  set->SetTH1(hLz_cut,"LHRS z cut event","z [m]","Counts");
  hLz_cut->SetLineColor(2);
  hLz_cut->SetFillColor(2);  
  hLz_cut->SetFillStyle(3002);

  hRz_cut=new TH1F("hRz_cut","",bin_z,min_z,max_z);
  set->SetTH1(hRz_cut,"RHRS z cut event","z [m]","Counts");
  hRz_cut->SetLineColor(2);
  hRz_cut->SetFillColor(2);  
  hRz_cut->SetFillStyle(3002);    

 
  hct_cut=new TH1F("hct_cut","",bin_ct,min_ct,max_ct);
  set->SetTH1(hct,"Coincidence time w/ cut","coin-time [ns]","Counts");  

  hct_cut->SetLineColor(2);
  hct_cut->SetFillColor(2);  
  hct_cut->SetFillStyle(3002);  


  
  min_mm=1.0; // [GeV]
  max_mm=1.5; // [GeV]
  int mm_width=5; // MeV
  double mmbin_width=(double)mm_width*1.0e-3; //[GeV]
  bin_mm=(int)((max_mm-min_mm)/mmbin_width);

  hmm_select=new TH1F("hmm_select","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_select,"Missing mass Event Selction","mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_cut=new TH1F("hmm_cut","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_cut,"Missing mass  w/ cut hist ","mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  //------ Event Selection hist ----------------//
  hmm_L=new TH1F("hmm_L","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_L,"Lambda missing mass Event Selction","mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));  
  hmm_L->SetLineColor(2);
  hmm_L->SetFillColor(2);  
  hmm_L->SetFillStyle(3002);

  hmm_S=new TH1F("hmm_S","",bin_mm,min_mm,max_mm);    
  set->SetTH1(hmm_S,"Sigma missing mass Event Selction","mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_S->SetLineColor(kBlue);
  hmm_S->SetFillColor(kBlue);  
  hmm_S->SetFillStyle(3002);



 //------- tuning check ----------------------//
  for(int i=0;i<nite;i++){

    hmm_L_t[i]=new TH1F(Form("hmm_L_t_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_t[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_L_tc[i]=new TH1F(Form("hmm_L_tc_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_tc[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_S_t[i]=new TH1F(Form("hmm_S_t_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_S_t[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_S_tc[i]=new TH1F(Form("hmm_S_tc_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_S_tc[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));        
    
    }
  

  
  
}

//==========================================================

void momcalib::MTParam(string mtparam){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }

  if(MODE==0 || MODE==-1) {MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;}
  if(MODE==0 || MODE==1) {MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;}
  MTP_mom();

  for(int i=0;i<10;i++)MT[i]=false;  
  //======= Tuning selection flag =====================//
  //--------- RHRS ------------------------//
  MT[0] = true;  // RHRS z correction
  MT[1] = true;  // RHRS raster correction
  MT[2] = true;  // RHRS theta correction
  MT[3] = true;  // RHRS phi correction
  //--------- LHRS -----------------------//
  MT[4] = true;  // LHRS z correction
  MT[5] = true;  // LHRS raster correction
  // MT[5] = false;  // LHRS raster correction
  MT[6] = true;  // LHRS theta correction
  MT[7] = true;  // LHRS phi correction
  //-------- momentum calibration ---------//
  MT[8] = true; // RHRS momentum correction  
  MT[9] = true; // LHRS momentum correction  
  ploss = true;  // Energy Loss
  //================================================//

  cout<<endl;
  cout<<"======== Correction Parameters ========="<<endl;
  if(MT[0])cout<<" RHRS z      correction "<<endl;
  else     cout<<" RHRS z                    no correction "<<endl;
  if(MT[1])cout<<" RHRS raster correction "<<endl;
  else     cout<<" RHRS raster               no correction "<<endl;
  if(MT[2])cout<<" RHRS theta  correction "<<endl;
  else     cout<<" RHRS theta                no correction "<<endl;
  if(MT[3])cout<<" RHRS phi    correction "<<endl;
  else     cout<<" RHRS phi                  no correction "<<endl;
  if(MT[4])cout<<" LHRS z      correction "<<endl;
  else     cout<<" LHRS z                    no correction "<<endl;
  if(MT[5])cout<<" LHRS raster correction "<<endl;
  else     cout<<" LHRS raster               no correction "<<endl;
  if(MT[6])cout<<" LHRS theta  correction "<<endl;
  else     cout<<" LHRS theta                no correction "<<endl;
  if(MT[7])cout<<" LHRS phi    correction "<<endl;
  else     cout<<" LHRS phi                  no correction "<<endl;
  if(MT[8])cout<<" RHRS mom    correction "<<endl;
  else     cout<<" RHRS mom                  no correction "<<endl;
  if(MT[9])cout<<" LHRS mom    correction "<<endl;
  else     cout<<" LHRS mom                  no correction "<<endl;
  if(ploss)cout<<" Energy Los  correction "<<endl;
  else     cout<<" Energy Los                no correction "<<endl;
  cout<<endl;
}



//=====================================================================//

void momcalib::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//

  
  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);

   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;
    //    cout<<"R Mzt : "<<Pzt[i]<<endl;
   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param[1].c_str()); // optimized
    //    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);

  for (int i=0;i<nParamT_ras;i++){

    Mras >> Pras[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras[i] = 0.0;
  }


  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);

    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param[3].c_str()); // optimized  
  ifstream Mypt(name_Mypt);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
  }
  Mypt.close();    



  
};
//////////////////////////////////////////////////////////////

void momcalib::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L); 
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param[5].c_str()); // optimized
    //    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);

  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> Pras_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras_L[i] = 0.0;
  }
  
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  //  cout<<"LHRS theta parameters file: "<<name_Mxpt_L<<endl;  
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    //   cout<<"LHRS theta : "<<par<<endl;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();

  
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param[7].c_str()); // optimized
  ifstream Mypt_L(name_Mypt_L);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;    
  }
  Mypt_L.close();    


}

//========================================================//


void momcalib::MTP_mom(){

  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, param[8].c_str()); // optimized
    ifstream Mpt(name_Mpt);

   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p >> p;
    Prp[i]=par;
    Opt_par_R[i]=par;
    Opt_par[i]=par;
   }
  Mpt.close();

  //====== LHRS Momentum parameters ========//
    char name_Mpt_L[500];
    sprintf(name_Mpt_L, param[9].c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);

   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p >> p;
    Plp[i]=par;
    Opt_par_L[i]=par;
    Opt_par[i+nParamTp]=par;  // Both momentum paramters
   }
  Mpt_L.close();





}

//===============================================================================//

void momcalib::zCorr(bool rarm, bool larm){

  
  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]   = (Rz[0]-Ztm)/Ztr;
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;    
  //===================================//    

    
    if(rarm) Rz[0]   = calcf2t_zt(Pzt, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0]); // nomalized
    if(larm) Lz[0]   = calcf2t_zt(Pzt_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0]); //nomalized

    

    //========== Scaled =================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;
    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    //===================================//
    
    Rz[0]   = Rz[0] * Ztr + Ztm;     // scaled
    Lz[0]   = Lz[0] * Ztr + Ztm;     // scaled    
    
}

//=================================================================================//

void momcalib::rasCorr(bool rarm, bool larm){

  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]   = (Rz[0]-Ztm)/Ztr;

    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;
  //===================================//    



    //======== Raster Correction ==========================//    

    if(rarm){
    RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled     
    Rz[0] = Rz[0] - RasterCor; // correction
    Rz[0]   = (Rz[0]-Ztm)/Ztr;    // nomalization    
    }
    
    if(larm){
    RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L = RasterCor_L/tan(hrs_ang);
    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled
    Lz[0] = Lz[0]+RasterCor_L;
    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr;    // nomalization
    }
    //====================================================//




    //========== Scaled ==================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    
    //===================================//

    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled
    Lz[0]  = Lz[0]*Ztr +Ztm; // scaled
    
    
}

//===============================================================================//

void momcalib::angCorr(bool rth, bool rph,bool lth, bool lph ){


  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]     = (Rz[0]-Ztm)/Ztr;
    Rth[0]    = (Rth[0] - Xptm)/Xptr;
    Rph[0]    = (Rph[0] - Yptm)/Yptr;    
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]     = (Lz[0]-Ztm)/Ztr;
    Lth[0]    = (Lth[0] - Xptm)/Xptr; 
    Lph[0]    = (Lph[0] - Yptm)/Yptr;     
  //===================================//    


    if(rth)Rth[0]  = calcf2t_ang(Pxpt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]); // nomalized
    if(lph)Rph[0]  = calcf2t_ang(Pypt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]); // nomalized
    
    if(lth)Lth[0]  = calcf2t_ang(Pxpt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]); // nomalized 
    if(lph)Lph[0]  = calcf2t_ang(Pypt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]); // nomalized   


    //========== Scaled ==================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    
    //===================================//

   
    Rph[0]  = Rph[0] * Yptr + Yptm; // scaled
    Rth[0]  = Rth[0] * Xptr + Xptm; // scaled
    Rz[0]   = Rz[0]  *  Ztr + Ztm; // scaled
    Lth[0]  = Lth[0] * Xptr + Xptm;  // scaled    
    Lph[0]  = Lph[0] * Yptr + Yptm;  // scaled    
    Lz[0]   = Lz[0]  * Ztr  + Ztm; // scaled
    
}

//================================================================================//

void momcalib::plossCorr(bool PLoss){

  

    if(PLoss){
    double dpe,dpe_,dpk;
    dpe=0.0; dpe_=0.0; dpk=0.0;
    dpe   = Eloss(0.0,Rz[0],"B");
    dpe_  = Eloss(Lth[0],Lz[0],"L");
    dpk   = Eloss(Rth[0],Rz[0],"R");
    pe    = pe    - dpe;
    Lp[0] = Lp[0] + dpe_;
    Rp[0] = Rp[0] + dpk;
    }// Energy loss

  

}

//===============================================================================//
void momcalib::momCorr(bool rarm, bool larm){


  
  
  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;

    Rz[0]   = (Rz[0]-Ztm)/Ztr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;
    Rp[0]   = (Rp[0]-PRm)/PRr;    
    Lp[0]   = (Lp[0]-PLm)/PLr;    
  //===================================//    
    
    if(rarm)Rp[0] = calcf2t_mom(Opt_par_R, Rx_fp[0], Ry_fp[0], Rth_fp[0], Rph_fp[0],Rz[0]);
    if(larm)Lp[0] = calcf2t_mom(Opt_par_L, Lx_fp[0], Ly_fp[0], Lth_fp[0], Lph_fp[0],Lz[0]);    

    

    //========== Scaled ==================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    
    //===================================//


    Rz[0]  = Rz[0] * Ztr + Ztm; // scaled
    Rp[0]  = Rp[0] * PRr + PRm; // scaled
    Lz[0]  = Lz[0] * Ztr + Ztm; // scaled
    Lp[0]  = Lp[0] * PLr + PLm; // scaled    

    
    
}

//===============================================================================//

double momcalib::SSx(double z,double th){

  double ssx=-100.0;
      for(int j=0 ; j<nfoil ; j++){
	l[j]=0.0;
	l[j] = sqrt(pow(l0,2.0) + pow(fcent_real[j]*100.,2.0) -2.0*l0*fcent_real[j]*100.*cos(hrs_ang));
	dth[j] = asin(l0/l[j]*sin(hrs_ang)) - hrs_ang;
	projectf[j] = cos( dth[j] );
	if(fcent[j]-selection_width<z
	   && z<fcent[j]+selection_width){
	  ssx = - th * l[j] * projectf[j];
	}
      }
      return ssx;
}

//==============================================================================//

double momcalib::SSy(double z,double ph){

  double ssy=-100.0;
      for(int j=0 ; j<nfoil ; j++){
	l[j]=0.0;
	l[j] = sqrt(pow(l0,2.0) + pow(fcent_real[j]*100.,2.0) -2.0*l0*fcent_real[j]*100.*cos(hrs_ang));
	dth[j] = asin(l0/l[j]*sin(hrs_ang)) - hrs_ang;
	projectf[j] = cos( dth[j] );
	
	if(fcent[j]-selection_width<z
	   && z<fcent[j]+selection_width){
	  ssy = - ph * l[j] * projectf[j];
	}
      }
      return ssy;
}


//==============================================================================//

void momcalib::ParamCorr(){


    //===================================================================//
    //============================ Tuning value =========================//
    //===================================================================//



  
 
 // RHRS //

  

    //   
    //    hRth->Fill(Rth[0]);    
    //    hRph->Fill(Rph[0]);

  



  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rth[0]    = (Rth[0] - Xptm)/Xptr;
    Rph[0]    = (Rph[0] - Yptm)/Yptr;    
    Rz[0]   = (Rz[0]-Ztm)/Ztr;
  //===================================//    



    
    if(MT[0])Rz[0] = calcf2t_zt(Pzt, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0]);

    
    if(MT[1]){
    RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled     
    Rz[0] = Rz[0] - RasterCor; // correction    
    hRz_c->Fill(Rz[0]);    
    Rz[0]   = (Rz[0]-Ztm)/Ztr;    // nomalization
    }
    
    
    if(MT[2])Rth[0]  = calcf2t_ang(Pxpt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]);
    if(MT[3])Rph[0] = calcf2t_ang(Pypt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]);
    if(MT[8])Rp[0] = calcf2t_mom(Opt_par_R, Rx_fp[0], Ry_fp[0], Rth_fp[0], Rph_fp[0],Rz[0]);

    
    Rph[0]  = Rph[0]*Yptr +Yptm;
    Rth[0]  = Rth[0]*Xptr +Xptm;
    Rz[0]  = Rz[0]*Ztr +Ztm;


    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;


    //    hRth_c->Fill(Rth[0]);        
    //    hRph_c->Fill(Rph[0]);    


    
 // LHRS //



    //   
    //    hLth->Fill(Lth[0]);    
    //    hLph->Fill(Lph[0]);


    //=========  nomalization  ================//
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr; 
    Lth[0]    = (Lth[0] - Xptm)/Xptr; 
    Lph[0]    = (Lph[0] - Yptm)/Yptr; 
    //========================================//
    
    


    if(MT[4])  Lz[0]   = calcf2t_zt(Pzt_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0]);

    //======== Raster Correction ==========================//
    if(MT[5]){
    RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L = RasterCor_L/tan(hrs_ang);
    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled
    Lz[0] = Lz[0]+RasterCor_L;
    hLz_c->Fill(Lz[0]);    
    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr;    // nomalization

    }
    //====================================================//
    
    if(MT[6])Lth[0]  = calcf2t_ang(Pxpt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]);
    if(MT[7])Lph[0]  = calcf2t_ang(Pypt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]);
    if(MT[9])Lp[0] = calcf2t_mom(Opt_par_L, Lx_fp[0], Ly_fp[0], Lth_fp[0], Lph_fp[0],Lz[0]);

    
    Lth[0]  = Lth[0]*Xptr +Xptm;  // scaled    
    Lph[0]  = Lph[0]*Yptr +Yptm;  // scaled
    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;



    //    hLth_c->Fill(Lth[0]);        
    //    hLph_c->Fill(Lph[0]);    


}


//================================================================================//
void momcalib::EventSelection(){

  cout<<endl;

  cout<<"=========================================="<<endl;
  cout<<"========= Event Selection ================"<<endl;
  cout<<"=========================================="<<endl;  

  int nLam=0;
  int nSig=0;
  
  ENum= T->GetEntries();
  ENum=1000000;
  cout<<"Select Events : "<<ENum<<endl;
  d=div(ENum,100);
  int nd=0;
  
  for (int i=0;i<ENum;i++){

    // ===== Initialization ====== //    
    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rz[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lz[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }
    rssx=-2222.0;
    rssy=-2222.0;
    lssx=-2222.0;
    lssy=-2222.0;    
    rssx_c=-2222.0;
    rssy_c=-2222.0;
    lssx_c=-2222.0;
    lssy_c=-2222.0;    

    mm    =-2222.0;
    mm_L  =-2222.0;
    mm_pi =-2222.0;
    Ra1sum=-2222.0;
    Ra2sum=-2222.0;
    Rgssum=-2222.0;
    Lpssum=-2222.0;    



    Lshsum=-2222.0;
    R_Ras_x=-2222.0;
    L_Ras_x=-2222.0;
    coin_t=-2222.0;
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    RPID_flag=false;
    LPID_flag=false;

    
    T->GetEntry(i);

    pe=hallap/1000.;
    
    hLz  ->Fill(Lz[0]);
    hLth ->Fill(Lth[0]);
    hLph ->Fill(Lph[0]);    
    hLp  ->Fill(Lp[0]);

    lssx=SSx(Lz[0],Lth[0]);    
    lssy=SSy(Lz[0],Lph[0]);        
    hLss->Fill(lssy,lssx);
    rssx=SSx(Rz[0],Rth[0]);    
    rssy=SSy(Rz[0],Rph[0]);        
    hRss->Fill(rssy,rssx);

    
    
    hRz  ->Fill(Rz[0]);
    hRth ->Fill(Rth[0]);
    hRph ->Fill(Rph[0]);    
    hRp  ->Fill(Rp[0]);


    
    
 //===== Parameter Tuning ============//
    zCorr(MT[0],MT[4]);
    hLz_c->Fill(Lz[0]);
    hRz_c->Fill(Rz[0]);
    
    rasCorr(MT[1],MT[5]);
    hLz_rc->Fill(Lz[0]);
    hRz_rc->Fill(Rz[0]);
    
    angCorr(MT[2],MT[3],MT[6],MT[7]);
    hLth_c->Fill(Lth[0]);
    hRth_c->Fill(Rth[0]);
    hLph_c->Fill(Lph[0]);
    hRph_c->Fill(Rph[0]);

    lssx_c=SSx(Lz[0],Lth[0]);    
    lssy_c=SSy(Lz[0],Lph[0]);        
    hLss_c->Fill(lssy_c,lssx_c);
    rssx_c=SSx(Rz[0],Rth[0]);    
    rssy_c=SSy(Rz[0],Rph[0]);        
    hRss_c->Fill(rssy_c,rssx_c);

    momCorr(MT[8],MT[9]);
    hRp_c->Fill(Rp[0]);
    hLp_c->Fill(Lp[0]);    
    //    plossCorr(ploss);
    //    hRp_ec->Fill(Rp[0]);
    //    hLp_ec->Fill(Lp[0]);
    
 //====================================//




    
 //------ Set Physics value --------//
    Ee=sqrt(pow(pe,2)+pow(Me,2));
    Ee_=sqrt(pow(Lp[0],2)+pow(Me,2));
    Ek=sqrt(pow(Rp[0],2)+pow(MK,2));
    Ls2pads=(int)Ls2tpads[0];
    Rs2pads=(int)Rs2tpads[0];
    rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
    lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
    rbeta=Rp[0]/Ek; 
    rpath_corr=rpathl/rbeta/c;
    lbeta=Lp[0]/Ee_; 
    lpath_corr=lpathl/lbeta/c;


 
 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( Lp[0], Lth_fp[0], Lph_fp[0] );
 R_v.SetMagThetaPhi( Rp[0], Rth_fp[0], Rph_fp[0] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );
 mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

 
 //--- Set Coincidence time ---------// 
 // Rs2_off=RS2_off_H1[Rs2pads];
 // Ls2_off=LS2_off_H1[Ls2pads];
 Rs2_off= s2f1_off(Rs2pads,"R",target,tdc_mode);
 Ls2_off= s2f1_off(Rs2pads,"L",target,tdc_mode);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   

 


 
 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true;
  if((Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax) &&  ( Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax)) z_flag=true; 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads) track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  if(Ra1sum<a1_th)      ac1_flag=true;
  if(Ra2sum>a2_th)      ac2_flag=true;

  gs_flag=true;
  ps_flag=true;
  sh_flag=true;
  
  if(ac1_flag==true && ac2_flag==true && gs_flag==true) RPID_flag=true;
  if(Lcersum>2000.)LPID_flag=true;  
  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag)      Scinti_flag=true;

  
 
    // ------------------------------------------- //
    // ------- Event selection for tuning -------- //
    // ------------------------------------------- //

  
  if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){

    if(RPID_flag && LPID_flag){
            hct->Fill(coin_t);
      if(z_flag){
	hLz_cut->Fill(Lz[0]);
	hRz_cut->Fill(Rz[0]);}  
      if(coin_flag)hct_cut->Fill(coin_t);      }
    
      //======== Production Events ============//
    if(coin_flag && RPID_flag && LPID_flag && z_flag){
      //     cout<<"mm: "<<mm<<endl;
      mm_L=mm;
      hmm_select->Fill(mm);
 
      //====== Lambda Selction =========//
      //      if(ML-mmL_range < mm && mm < ML+mmL_range && ntune_event<nmax){
      if(1.12 < mm && mm < 1.15 && ntune_event<nmax){      
      mass_flag[ntune_event]=0;
      mass_ref[ntune_event] = ML;
      MM[ntune_event]=mm;
      beam_p[ntune_event]=pe;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;
      rp[ntune_event] =Rp[0];
      lp[ntune_event] =Lp[0];      
      bp[ntune_event] =pe;
      hmm_L->Fill(mm);
      nLam++;
      ntune_event++;
      
      //      }else if(MS0-mmS_range < mm && mm < MS0+mmS_range && ntune_event<nmax){
      }else if(1.20 < mm && mm < 1.23 && ntune_event<nmax){      
      mass_flag[ntune_event]=1;
      mass_ref[ntune_event] = MS0;
      beam_p[ntune_event]=pe;
      MM[ntune_event]=mm;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;            
      rp[ntune_event] =Rp[0];
      lp[ntune_event] =Lp[0];      
      bp[ntune_event] =pe;
      hmm_S->Fill(mm);
      nSig++;
      ntune_event++;      }



    }
  } // event selection
  
  if(i==d.quot * 10*nd){
    cout<<10*nd<<" % Filled ("<<i<<" / "<<ENum<<")"<<endl; nd++;}

  }//End Fill


  cout<<"Select Events: "<<ntune_event<<endl;
  cout<<"Select Lambda events : "<<nLam<<endl;
  cout<<"Select Sigma  events : "<<nSig<<endl;
  }

//========================================================//




void momcalib::MomTuning(string ofname){


  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started ====================" << endl;
    cout << "======================================================" <<endl;}
    char tempc[500],tempc2[500];
    const  char* new_tempc=ofname.c_str();
    cout<<"new marix file: "<<new_tempc<<endl;
    
    ofstream * ofs1;
    ofstream * ofs2;       
  for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    
    if(MODE==-1){
    cout<<"------- Rp tuning -----"<<endl;
    chi_sq1[i]=0.0;
    chi_sq1[i] = tune(Opt_par_R,i,-1);   // Rp
    }else if(MODE==1){
    cout<<"------- Lp tuning -----"<<endl;    
    chi_sq2[i]=0.0;
    chi_sq2[i] = tune(Opt_par_L,i,1);  // Lp
    }else if(MODE==0){
    cout<<"---------pk & pe tuning ------"<<endl;
    chi_sq[i]=0.0;
    //    for(int s=0;s<40;s++)cout<<"opt_para [ "<<s<<" ] : "<<Opt_par[s]<<endl;
    chi_sq[i] = tune(Opt_par,i,0);
    }

    cout << " Tuning# = " << i << ": chisq = ";

    if(MODE==-1){
      cout  << chi_sq1[i] << endl;
      gchi_Rp->SetPoint(i,i,chi_sq1[i]);          
    }else if(MODE== 1){
      cout  << chi_sq2[i] << endl;
      gchi_Lp->SetPoint(i,i,chi_sq2[i]);     
    }else if(MODE== 0){
      cout  << chi_sq[i]  << endl;
      gchi_p ->SetPoint(i,i,chi_sq[i]);    }

    
    if(MODE==-1 || MODE==0){
    sprintf(tempc,  "%s_Rp_%d_MODE%d.dat",new_tempc,i,MODE); 
    cout<<"new matrix Rp: "<<tempc<<endl;
    ofs1 = new ofstream(tempc);}

    if(MODE==1 ||MODE==0){
    sprintf(tempc2, "%s_Lp_%d_MODE%d.dat",new_tempc,i,MODE);
    cout<<"new matrix Lp: "<<tempc2<<endl;    
    ofs2 = new ofstream(tempc2);}
 
    int nppp = 0;
    for(int i=0 ; i<nnp+1 ; i++){
      for(int e=0 ; e<nnp+1 ; e++){
	for(int d=0 ; d<nnp+1 ; d++){
	  for(int c=0 ; c<nnp+1 ; c++){
	    for(int b=0 ; b<nnp+1 ; b++){
	      for(int a=0 ; a<nnp+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(MODE==-1){		  
		  *ofs1 << Opt_par_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==1){
		  *ofs2 << Opt_par_L[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==0){
		  *ofs1 << Opt_par_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  *ofs2 << Opt_par_L[nppp+nParamTp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  cout<<"Rparam ["<<nppp<<"] "<<Opt_par_R[nppp]<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;		  
		  cout<<"Lparam ["<<nppp<<"] "<<Opt_par_L[nppp]<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<endl;
		  
		  }

		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }

    if(MODE==-1 || MODE==0){
      ofs1->close();
      ofs1->clear();
    }else if(MODE==0 || MODE==1){
      ofs2->close();
      ofs2->clear();}


    
  }

  
  cout<<"========== Tuning is done ============="<<endl;
  
  

}

//========================================================//

double momcalib::Eloss(double xp, double z,char* arm){



  
  double x = - tan(hrs_ang-xp); //xp : theta [rad]
  double ph[3],pl[2];
  double dEloss;
  bool high;
  
  if(z>0.08)high=false;
  else high=true;
  
    //==== thickness 0.400 mm ========//

  if(arm=="R"){
    ph[0] = -1.3175;
    ph[1] = -4.6151;
    ph[2] = 2.0369;
    pl[0] = 3.158e-2;
    pl[1] = 4.058e-1;
  }else if(arm=="L"){
    ph[0] = -1.3576;
    ph[1] = -4.5957;
    ph[2] = 2.0909;
    pl[0] = 6.2341e-3;
    pl[1] = 4.0336e-1;
  }

  double dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];
  double dEloss_l = pl[0]*x +pl[1];

  
  if(high)dEloss = dEloss_h;
  else dEloss = dEloss_l;
  //==== thickness 0.4 mm in beam energy loss ======//
  if(arm=="B")dEloss=0.23; //[MeV/c]
  dEloss=dEloss/1000.; // [GeV/c]
  return dEloss;
  
}

//========================================================//

void momcalib::Fill(){

  cout<<endl;
  cout<<"======================================"<<endl;
  cout<<"============= Event Fill ============="<<endl;
  cout<<"======================================"<<endl;  

  

  ENum= T->GetEntries();
  cout<<"Events :"<<ENum<<endl;
     if(ENum<10000)d=div(ENum,1000);
   else   d=div(ENum,10000);

  for (int i=0;i<ENum;i++){

    // ===== Initialization ====== //    

    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }

    
    mm    =-2222.0;
    mm_L  =-2222.0;
    mm_pi =-2222.0;
    mm_acc =-2222.0;
    Ra1sum=-2222.0;
    Ra2sum=-2222.0;
    Rgssum=-2222.0;
    Lpssum=-2222.0;    
    Lshsum=-2222.0;
    
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    RPID_flag=false;
    LPID_flag=false;

    T->GetEntry(i);
   
    pe=hallap/1000.;
    
 //===== Parameter Tuning ============//
    zCorr(MT[0],MT[4]);
    rasCorr(MT[1],MT[5]);
    angCorr(MT[2],MT[3],MT[6],MT[7]);
    momCorr(MT[8],MT[9]);
    plossCorr(ploss);    
 //====================================//

    
 //------ Set Physics value --------//
 Ee=sqrt(pow(pe,2)+pow(Me,2));
 Ee_=sqrt(pow(Lp[0],2)+pow(Me,2));
 Ek=sqrt(pow(Rp[0],2)+pow(MK,2));
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rbeta=Rp[0]/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=Lp[0]/Ee_; 
 lpath_corr=lpathl/lbeta/c;

    

 
 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( Lp[0], Lth_fp[0], Lph_fp[0] );
 R_v.SetMagThetaPhi( Rp[0], Rth_fp[0], Rph_fp[0] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );

 mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
 
 //--- Set Coincidence time ---------// 
 Rs2_off=RS2_off_H1[Rs2pads];
 Ls2_off=LS2_off_H1[Ls2pads];
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   


 
 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true;
  if(Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax && Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax) z_flag=true; 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads) track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  if(Ra1sum<a1_th)      ac1_flag=true;
  if(Ra2sum>a2_th)      ac2_flag=true;

  gs_flag=true;
  ps_flag=true;
  sh_flag=true;
  
  if(ac1_flag==true && ac2_flag==true && gs_flag==true) RPID_flag=true;
  //  if(ps_flag==true && sh_flag==true)                    LPID_flag=true;
  if(Lcersum>2000.)LPID_flag=true;  
  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag)      Scinti_flag=true;

  
  if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){

    
  //======== Production Events ============//
    if(coin_flag && RPID_flag && LPID_flag && z_flag){
      mm_L=mm;
      hmm_cut->Fill(mm);
    }
  //======== ACCidencel B.G. ==============///
    if(RPID_flag && LPID_flag && z_flag &&
       (-63<coin_t && coin_t <-15) || (15<coin_t && coin_t<63)){
      ct=coin_t;
      mm_acc=mm;
              while(1){
	       if(-3.0<ct && ct<3.0){ctime=ct; break;}
	       else if(ct<-3.0) ct=ct+6;
	       else if(3.0<ct)  ct=ct-6;
	      }
    }
    //=====================================//

    
  }

  
  tnew->Fill();
          if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }


  
}

//========================================================//


//=======================================================//

void momcalib::Close(){

   tnew->Write();
   hLz->Write();
   hLth->Write();
   hLph->Write();
   hLz_c->Write();
   hLz_rc->Write();   
   hLth_c->Write();
   hLph_c->Write();
   hLss->Write();
   hLss_c->Write();
   hLp->Write();
   hLp_c->Write();   
   hLp_ec->Write();
   hRz->Write();
   hRth->Write();
   hRph->Write();
   hRz_c->Write();
   hRz_rc->Write();   
   hRth_c->Write();
   hRph_c->Write();
   hRss->Write();
   hRss_c->Write();
   hRp->Write();
   hRp_c->Write();   
   hRp_ec->Write();
   hLz_cut->Write();
   hRz_cut->Write();
   hct->Write();   
   hct_cut->Write();

   
   hmm_select->Write();
   hmm_L->Write();
   hmm_S->Write();
   hmm_cut->Write();
   
   ofp->Close();

}



//=========================================================//
//================= Function ==============================//
//=========================================================//

// ####################################################
double s2f1_off(int i,char* ARM,char* MODE, int KINE){
// ####################################################

  double RS2_offset[16],LS2_offset[16];
  if(MODE=="H" && KINE==2){
 
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(MODE=="H" && KINE==1){

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
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


// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 0.01;//[GeV/c^2]
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  double mass;
  double Ee,Ee_,Ek;
  double rbeta,lbeta;
  //  double rp,lp;
  double chi2=0.0;
  double rx_c,ry_c,rth_c,rph_c,rz_c,lx_c,ly_c,lth_c,lph_c;

  for(int i=0 ; i<ntune_event ; i++){

    //====== Initialization ========//

    residual =0.0;


    // RHRS //
      
    //======== Scaled paramters ===============//

    rx_fp[i]  = (rx_fp[i]-XFPm)/XFPr;
    rth_fp[i] = (rth_fp[i]-XpFPm)/XpFPr;
    ry_fp[i]  = (ry_fp[i]-YFPm)/YFPr;
    rph_fp[i] = (rph_fp[i]-YpFPm)/YpFPr;
    rth[i]    = (rth[i] - Xptm)/Xptr;
    rph[i]    = (rph[i] - Yptm)/Yptr;    
    rz[i]     = (rz[i]-Ztm)/Ztr;
    rp[i]     = (rp[i] - PRm)/PRr;
									     
    //==== momentum tuning ==========//
    rp[i] = calcf2t_mom(Opt_par_R, rx_fp[i], ry_fp[i], rth_fp[i], rph_fp[i],rz[i]);

    //==== Scaled to Nomal =======//
    
    rth[i]  = rth[i]*Xptr +Xptm;
    rph[i]  = rph[i]*Yptr +Yptm;
    rz[i]   = rz[i]*Ztr +Ztm;
    rp[i]   = rp[i]*PRr +PRm;
    
    rx_fp[i]  = rx_fp[i]  * XFPr + XFPm;
    rth_fp[i] = rth_fp[i] * XpFPr + XpFPm;
    ry_fp[i]  = ry_fp[i]  * YFPr + YFPm;
    rph_fp[i] = rph_fp[i] * YpFPr + YpFPm;
    
    //=============================//
   

    // LHRS //
    
    //======== Scaled paramters ===============//

    lx_fp[i]  = (lx_fp[i]-XFPm)/XFPr;
    lth_fp[i] = (lth_fp[i]-XpFPm)/XpFPr;
    ly_fp[i]  = (ly_fp[i]-YFPm)/YFPr;
    lph_fp[i] = (lph_fp[i]-YpFPm)/YpFPr;
    lth[i]    = (lth[i] - Xptm)/Xptr;
    lph[i]    = (lph[i] - Yptm)/Yptr;    
    lz[i]     = (lz[i]-Ztm)/Ztr;
    lp[i]     = (lp[i] - PLm)/PLr;
    //==== momentum tuning ==========//
    
    lp[i] = calcf2t_mom(Opt_par_L, lx_fp[i], ly_fp[i], lth_fp[i], lph_fp[i], lz[i]);


    //==== Scaled to Nomal =======//
    
    lth[i]  = lth[i]*Xptr +Xptm;
    lph[i]  = lph[i]*Yptr +Yptm;
    lz[i]   = lz[i]*Ztr + Ztm;
    lp[i]   = lp[i]*PLr + PLm;    
    
    lx_fp[i]  = lx_fp[i]  * XFPr + XFPm;
    lth_fp[i] = lth_fp[i] * XpFPr + XpFPm;
    ly_fp[i]  = ly_fp[i]  * YFPr + YFPm;
    lph_fp[i] = lph_fp[i] * YpFPr + YpFPm;

    //=============================//
   
    

    


    
 //------ Set Physics value --------//
 Ee=sqrt(pow(bp[i],2)+pow(Me,2));
 Ee_=sqrt(pow(lp[i],2)+pow(Me,2));
 Ek=sqrt(pow(rp[i],2)+pow(MK,2));


 rbeta=rp[i]/Ek; 
 lbeta=lp[i]/Ee_; 


 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( lp[i], lth_fp[i], lph_fp[i] );
 R_v.SetMagThetaPhi( rp[i], rth_fp[i], rph_fp[i] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );

 mass = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

    
    //======================//
    //====== Residual ======//
    //======================//

    
    residual = mass - mass_ref[i];
    chi2=chi2 + pow(residual/sigma,2.0);
    //    cout<<"mass: "<<mass<<" ref:"<<mass_ref[i]<<" chi2 "<<chi2<<endl;
  }

 
}//end fcn



// #############################################################
double momcalib::tune(double* pa, int j, int MODE) 
// #############################################################
{
  
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamTp;
  int arm=1;
  if(MODE==0){allparam=nParamTp*2;
    arm=2;
  }


  //======= Coincidence analysis ========//
  if(MODE==0){
    for(int i=0;i<nParamTp;i++){
      Opt_par_R[i]=pa[i];
      Opt_par_L[i]=pa[i+nParamTp];
    }
  }
  //======================================//



  TMinuit* minuit= new TMinuit(allparam);
  minuit->SetFCN(fcn); // fcn Chi-square function


  double start[allparam];
  double step[allparam];
  const int nMatT =nnp;  
  const int nXf   =nnp;
  const int nXpf  =nnp;
  const int nYf   =nnp;
  const int nYpf  =nnp;
  const int nZt   =nnp; // The number of order is reduced for test (4-->2)
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;

  for(int f=0;f<arm;f++){
    for (int n=0;n<nMatT+1;n++){
      for(e=0;e<n+1;e++){
	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){ 
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){ 
		if (a+b+c+d+e==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		    start[npar] = pa[npar];
		     step[npar] = 1.0e-3;
		  }
		  else{
		    start[npar] = 0.0;
		    step[npar] = 0.0;
		  }
		  npar++;

		}
	      }
	    }
	  }
	}    
      }
    }
  }

  // ~~~ Chi-square ~~~~
  arglist[0] = 1.0;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  minuit -> SetPrintLevel(-1);
  
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i]; 
    step[i] = 1.0;

    //          LLim[i] = pa[i] - pa[i]*0.8;
    //          ULim[i] = pa[i] + pa[i]*0.8;
    //          LLim[i] = pa[i] - 5.0; // temp
    //          ULim[i] = pa[i] + 5.0; // temp
    
          LLim[i] = pa[i] - 1000.0; // temp
          ULim[i] = pa[i] + 1000.0; // temp
	  minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0; // original
  //arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
 
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  cout<<endl;
  cout<<"amin: "<<amin<<" edm: "<<edm<<" errdef: "<<errdef<<" nvpar: "<<nvpar<<endl;
  cout<<endl;
  
  if(amin>0) chi2=amin;

  
  for(int i=0 ; i<allparam ; i++){
    if(MODE==0){minuit -> GetParameter(i,Opt_par[i],er);
      if(i<allparam/2){minuit -> GetParameter(i,Opt_par_R[i],er);}// RHRS momentum parameters
      else {minuit -> GetParameter(i,Opt_par_L[i],er);}// RHRS momentum parameters
    }else if(MODE==-1){minuit -> GetParameter(i,Opt_par_R[i],er);// RHRS momentum parameters
    }else if(MODE==1)minuit -> GetParameter(i,Opt_par_L[i],er); // LHRS momentum parameters

  }
  
  return chi2;
}


// ###################################################
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
// ###################################################
{
  // ----------------------------------------------------------------- //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp --- //
  // ----------------------------------------------------------------- //

  const int nMatT=nnp;  
  const int nXf=nnp;
  const int nXpf=nnp;
  const int nYf=nnp;
  const int nYpf=nnp;
  const int nZt=nnp;
  const int nXt=nnp;
  const int nXpt=nnp;
  const int nYt=nnp;
  const int nYpt=nnp;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0,f=0,g=0,h=0,i=0;
  for (int n=0;n<nMatT+1;n++){
	    for(e=0;e<n+1;e++){
	      for (d=0;d<n+1;d++){
		for (c=0;c<n+1;c++){ 
		  for (b=0;b<n+1;b++){
		    for (a=0;a<n+1;a++){ 
		      if (a+b+c+d+e+f+g+h+i==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt
		    && f<=nXt && g<=nXpt && h<=nYt && i<=nYpt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
		      }
		    }
		  }
		}
	      }    
	    }
	  }

  return Y; 
  
}



// ###################################################
double calcf2t_ang(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
// ####################################################
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  
  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
	      npar++;
	      }  
	    }
	  }
	}
      }    
    }
  }


  return Y; 
  
}
      	

// #################################################
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ###############################################

  int nnz=3;  

  const int nMatT=nnz; 
  const int nXf=nnz;
  const int nXpf=nnz;
  const int nYf=nnz;
  const int nYpf=nnz;

  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d));
		  }
		  else{
		    x = 0.;
		  }
		  Y += x*P[npar];
		  npar++;
		}
		
	      }
	    }
	  }
  	}
  }
  
  return Y;
}







#endif
