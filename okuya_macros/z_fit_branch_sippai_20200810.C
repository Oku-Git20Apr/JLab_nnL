//---Okuyama---//
//--2020/8/10--//
//
//(1) run analyzer (z_analyzer.cc)
//--> create fout.root and (TFile)file_out and (TTree)tree_out
//
//(2) run fitting macro (z_fit.C, this file)
//--> fit all of the histograms (e.g. cointime, missing mass for every cut condition)
//--> make some efficiency curves
//--> estimate systematic errors

#define MAX 50     // Maximum No. of Tracks
#define RS0 1      // No. of Segments of R-S0
#define RS2 16     // No. of Segments of R-S2
#define RA1 24     // No. of Segments of R-AC1
#define RA2 26     // No. of Segments of R-AC2
#define RCR 10     // No. of Segments of R-GC
#define LCL 10     // No. of Segments of L-GC
#define RPS 48     // No. of Segments of R-Pre-Shower
#define RSH 75     // No. of Segments of R-Shower
#define RF1TDC 64  // No. of ch of R-F1TDC
#define LS0 1      // No. of Segments of L-S0
#define LS2 16     // No. of Segments of L-S2
#define LF1TDC 64  // No. of ch of L-F1TDC



//double CoinCalc_gogami(int RS2_seg, int LS2_seg,int rhit, int lhit);



void z_fit(){


	string ofname = "z_fit.dat";
	string pdfname = "z_fit.pdf";
cout << " output file name is " << ofname << endl;
cout << " output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("fout.root","read");
  TFile *file_new = new TFile("fout_new.root","recreate");
cout<<"aaaaaaaaaaaaaaaaaaa"<<endl;
  TTree *tree_old = (TTree*)file->Get("tree_out");
cout<<"aaaaaaaaaaaaaaaaaaa"<<endl;
if(!tree_old){cout<<"Error"<<endl;return 1;}
cout<<"aaaaaaaaaaaaaaaaaaa"<<endl;
  TTree *tree = tree_old->CloneTree(1000);
cout<<"aaaaaaaaaaaaaaaaaaa"<<endl;
	tree->Write();
//  file->Write();
 


// =================================================== //
// ==== Offset and scaling factors for matrices ====== //
// =================================================== //
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  Ztm = -0.15,Ztr=0.35;
//==== momentum scaled  parameters =====//
const double  PLm = 2.0, PLr=0.22; 
const double  PRm =1.74, PRr=0.2;
const double  PaRm = 25.713, PaRr = 0.4;
const double  PaLm = 25.713, PaLr = 0.4;
//---Physics Constant---//
 
const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_sig_S=0.004;
 const double def_mean_S=MS0-ML;
 const double def_sig_p=0.852;
 const double def_mean_p=-8.0;
 const double def_sig_pi=0.443;
 const double def_mean_pi=3.0;
 const double def_sig_k=0.644;
 const double def_mean_k=0.0;
 const double def_acc=27.7;
 const double min_coin_c=-20.0;
 const double max_coin_c=20.0;
 const double tdc_time=0.056;//ns
 int bin_coin_c=(int)((max_coin_c-min_coin_c)/tdc_time);
 const double min_mm=-0.1;//GeV/c^2
 const double max_mm=0.2;//GeV/c^2
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;
 
//---For Efficiency Result---//
 double center_pi, center_k, center_p, center_L, center_S;
 double range_pi, range_k, range_p, range_L, range_S;
 double n_pi_noZ, n_k_noZ, n_p_noZ, n_L_noZ, n_S_noZ;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_noZ, mean_k_noZ, mean_p_noZ, mean_L_noZ, mean_S_noZ;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_noZ, sig_k_noZ, sig_p_noZ, sig_L_noZ, sig_S_noZ;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];


//---Fitting Function---//
 TF1* fpi_noZ;
 TF1* fk_noZ;
 TF1* fp_noZ;
 TF1* fmmbg_noZ;
 TF1* fL_noZ;
 TF1* fS_noZ;
 TF1* fmm_noZ;
 TF1* fcoin_noZ;
 TF1* facc[100];
 TF1* fpi[100];
 TF1* fk[100];
 TF1* fp[100];
 TF1* fcoin[100];
 TF1* fL[100];
 TF1* fS[100];
 TF1* fmm[100];
 TF1* fmmbg[100];

//---Histograms from analyzer---//
 TH1F *hcoin_k_fom_noZ=(TH1F*)file->Get("hcoin_k_fom_noZ");
 TH1F *hcoin_bg_fom_noZ=(TH1F*)file->Get("hcoin_bg_fom_noZ");
 TH1F *hcoin_wo_bg_fom_noZ=(TH1F*)file->Get("hcoin_wo_bg_fom_noZ");
 TH1F *hcoin_pi_noZ=(TH1F*)file->Get("hcoin_pi_noZ");
 TH1F *hmm_L_fom_noZ=(TH1F*)file->Get("hmm_L_fom_noZ");
 TH1F *hmm_bg_fom_noZ=(TH1F*)file->Get("hmm_bg_fom_noZ");
 TH1F *hmm_wo_bg_fom_noZ=(TH1F*)file->Get("hmm_wo_bg_fom_noZ");
 TH1F *hmm_pi_fom_noZ=(TH1F*)file->Get("hmm_pi_fom_noZ");
 TH1F *hmm_pibg_fom_noZ=(TH1F*)file->Get("hmm_pibg_fom_noZ");
 TH1F *hmm_pi_wobg_fom_noZ=(TH1F*)file->Get("hmm_pi_wobg_fom_noZ");
 TH1F *hcoin_k_fom[100];
 TH1F *hcoin_bg_fom[100];
 TH1F *hcoin_wo_bg_fom[100];
 TH1F *hcoin_pi[100];
 TH1F *hmm_L_fom[100];
 TH1F *hmm_bg_fom[100];
 TH1F *hmm_wo_bg_fom[100];
 TH1F *hmm_pi_fom[100];
 TH1F *hmm_pibg_fom[100];
 TH1F *hmm_pi_wobg_fom[100];
for(int i=0;i<100;i++){
hcoin_k_fom[i]=(TH1F*)file->Get(Form("hcoin_k_fom[%d][0][0]",i));
hcoin_bg_fom[i]=(TH1F*)file->Get(Form("hcoin_bg_fom[%d][0][0]",i));
hcoin_wo_bg_fom[i]=(TH1F*)file->Get(Form("hcoin_wo_bg_fom[%d][0][0]",i));
hcoin_pi[i]=(TH1F*)file->Get(Form("hcoin_pi[%d][0][0]",i));
hmm_L_fom[i]=(TH1F*)file->Get(Form("hmm_L_fom[%d][0][0]",i));
hmm_bg_fom[i]=(TH1F*)file->Get(Form("hmm_bg_fom[%d][0][0]",i));
hmm_wo_bg_fom[i]=(TH1F*)file->Get(Form("hmm_wo_bg_fom[%d][0][0]",i));
hmm_pi_fom[i]=(TH1F*)file->Get(Form("hmm_pi_fom[%d][0][0]",i));
hmm_pibg_fom[i]=(TH1F*)file->Get(Form("hmm_pibg_fom[%d][0][0]",i));
hmm_pi_wobg_fom[i]=(TH1F*)file->Get(Form("hmm_pi_wobg_fom[%d][0][0]",i));
}

//---Efficiency---//
TEfficiency *pEff1;
TEfficiency *pEff2;
TEfficiency *pEff3;
TEfficiency *pEff4;
TEfficiency *pEff5;
//TH1F* h_pisr1;
//TH1F* h_ksr1;
//TH1F* h_psr1;
//TH1F* h_Lsr1;
//TH1F* h_Ssr1;
//TH1F* h_pitot1;
//TH1F* h_ktot1;
//TH1F* h_ptot1; 
//TH1F* h_Ltot1; 
//TH1F* h_Stot1; 
//---Default Result---//
//---Efficiency result from analyzer---//
TH1F* h_pisr1=(TH1F*)file->Get("h_pisr1");
TH1F* h_ksr1=(TH1F*)file->Get("h_ksr1");
TH1F* h_psr1=(TH1F*)file->Get("h_psr1");
TH1F* h_Lsr1=(TH1F*)file->Get("h_Lsr1");
TH1F* h_Ssr1=(TH1F*)file->Get("h_Ssr1");
TH1F* h_pitot1=(TH1F*)file->Get("h_pitot1");
TH1F* h_ktot1=(TH1F*)file->Get("h_ktot1");
TH1F* h_ptot1=(TH1F*)file->Get("h_ktot1"); 
TH1F* h_Ltot1=(TH1F*)file->Get("h_ktot1"); 
TH1F* h_Stot1=(TH1F*)file->Get("h_ktot1"); 
//--------------------------------------------//
 h_pisr1->Sumw2();
 h_pisr1->SetLineColor(kRed);
 h_ksr1->Sumw2();
 h_ksr1->SetLineColor(kRed);
 h_psr1->Sumw2();
 h_psr1->SetLineColor(kRed);
 h_Lsr1->Sumw2();
 h_Lsr1->SetLineColor(kAzure);
 h_Ssr1->Sumw2();
 h_Ssr1->SetLineColor(kCyan);
 h_pitot1->Sumw2();
 h_pitot1->SetLineColor(kBlack);
 h_ktot1->Sumw2();
 h_ktot1->SetLineColor(kBlack);
 h_ptot1->Sumw2();
 h_ptot1->SetLineColor(kBlack);
 h_Ltot1->Sumw2();
 h_Ltot1->SetLineColor(kAzure);
 h_Stot1->Sumw2();
 h_Stot1->SetLineColor(kCyan);

//--Cut Condition--//
 int nth=100;
 double zver[100];
 zver[0]=0.;
 for(int i=1; i<nth; i++) zver[i]=zver[i-1]+0.005;//SUM
 

//---------------------------//
//--Event by Event Analysis--//
//---------------------------//
//-----,if neccessary--------//
//---------------------------//
  int ev = 0;
  int ENum = 0;
  int ENum2 = 0;
    //double ct,mm,ct_den,mm_den,mmbg_den;
    //double ct_eff[100], mm_eff[100], mmbg_eff[100];
    int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
	double ct, ct_eff;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("tr.ntrack_l",1);tree->SetBranchAddress("tr.ntrack_l",&NLtr);
	tree->SetBranchStatus("tr.ntrack_r",1);tree->SetBranchAddress("tr.ntrack_r",&NRtr);
	tree->SetBranchStatus("tr.Ls2_pad[100]",1);tree->SetBranchAddress("tr.Ls2_pad[100]",&Ls2_pad);
	tree->SetBranchStatus("tr.Rs2_pad[100]",1);tree->SetBranchAddress("tr.Rs2_pad[100]",&Rs2_pad);
	//tree->SetBranchStatus("tr.ct_gb",1);tree_new->SetBranchAddress("tr.ct_gb",&ct);
	//tree->SetBranchStatus("tr.ct_eff[100]",1);tree_new->SetBranchAddress("tr.ct_eff[100]",&ct_eff);
	//tree->SetBranchStatus("tr.ct_den",1);tree_new->SetBranchAddress("tr.ct_den",&ct_den);
	//tree->SetBranchStatus("tr.mm_den",1);tree_new->SetBranchAddress("tr.mm_den",&mm_den);
	//tree->SetBranchStatus("tr.mm_eff[100]",1);tree_new->SetBranchAddress("tr.mm_eff[100]",&mm_eff);
	//tree->SetBranchStatus("tr.mmbg_den",1);tree_new->SetBranchAddress("tr.mmbg_den",&mmbg_den);
	//tree->SetBranchStatus("tr.mmbg_eff[100]",1);tree_new->SetBranchAddress("tr.mmbg_eff[100]",&mmbg_eff);
//// Tracking ////
  double R_tr_n;                                                                 // No. of Tracks
  double R_tr_flag[MAX], R_tr_ndof[MAX];                                         // track status, track NDoF
  double R_tr_chi2[MAX];                                                         // track chi2
  double R_tr_beta[MAX];                                                         // beta of track
  double R_tr_d_x[MAX], R_tr_d_y[MAX], R_tr_d_th[MAX], R_tr_d_ph[MAX];           // x, y, theta, phi at Detector
  double R_tr_r_x[MAX], R_tr_r_y[MAX], R_tr_r_th[MAX], R_tr_r_ph[MAX];           // x, y, theta, phi at Rot-Coordinate
  double R_tr_x[MAX], R_tr_y[MAX], R_tr_th[MAX], R_tr_ph[MAX];                   // x, y, theta, phi
  double R_tr_time[MAX];                                                         // theta, time@RF
  double R_tr_p[MAX], R_tr_pathl[MAX], R_tr_px[MAX], R_tr_py[MAX], R_tr_pz[MAX]; // mom(unit 10GeV?), length(TtoP), momx, momy, momz
  double R_tr_tg_dp[MAX], R_tr_tg_y[MAX], R_tr_tg_th[MAX], R_tr_tg_ph[MAX];      // dp, y, theta, phi at target
  double R_tr_vx[MAX], R_tr_vy[MAX], R_tr_vz[MAX]; // vertex X, Y, Z

  //  double R_rtRFtime[6];
  double R_rtRFtime[MAX];
  double L_tr_n;                                                                 // No. of Tracks
  double L_tr_flag[MAX], L_tr_ndof[MAX];                                         // track status, track NDoF
  double L_tr_chi2[MAX];                                                         // track chi2
  double L_tr_beta[MAX];                                                         // beta of track
  double L_tr_d_x[MAX], L_tr_d_y[MAX], L_tr_d_th[MAX], L_tr_d_ph[MAX];           // x, y, theta, phi at Detector
  double L_tr_r_x[MAX], L_tr_r_y[MAX], L_tr_r_th[MAX], L_tr_r_ph[MAX];           // x, y, theta, phi at Rot-Coordinate
  double L_tr_x[MAX], L_tr_y[MAX], L_tr_th[MAX], L_tr_ph[MAX];                   // x, y, theta, phi
  double L_tr_time[MAX];                                                         // theta, time@RF
  double L_tr_p[MAX], L_tr_pathl[MAX], L_tr_px[MAX], L_tr_py[MAX], L_tr_pz[MAX]; // mom, length(TtoP), momx, momy, momz
  double L_tr_tg_dp[MAX], L_tr_tg_y[MAX], L_tr_tg_th[MAX], L_tr_tg_ph[MAX];      // dp, y, theta, phi at target
  double L_tr_vx[MAX], L_tr_vy[MAX], L_tr_vz[MAX]; // vertex X, Y, Z

  double R_s2_trpath[MAX], L_s2_trpath[MAX];
  double RTDC_F1FirstHit[RF1TDC];

  double L_rtRFtime[6];
  tree->SetBranchStatus("R.tr.x"               ,1);  tree->SetBranchAddress("R.tr.x"               , R_tr_x              );
  tree->SetBranchStatus("R.tr.y"               ,1);  tree->SetBranchAddress("R.tr.y"               , R_tr_y              );
  tree->SetBranchStatus("R.tr.th"              ,1);  tree->SetBranchAddress("R.tr.th"              , R_tr_th             );
  tree->SetBranchStatus("R.tr.ph"              ,1);  tree->SetBranchAddress("R.tr.ph"              , R_tr_ph             );
  //  tree->SetBranchStatus("R.tr.time"            ,1);  tree->SetBranchAddress("R.tr.time"            , R_tr_time           );
  tree->SetBranchStatus("R.tr.p"               ,1);  tree->SetBranchAddress("R.tr.p"               , R_tr_p              );
  tree->SetBranchStatus("R.tr.pathl"           ,1);  tree->SetBranchAddress("R.tr.pathl"           , R_tr_pathl          );
  tree->SetBranchStatus("R.tr.px"              ,1);  tree->SetBranchAddress("R.tr.px"              , R_tr_px             );
  tree->SetBranchStatus("R.tr.py"              ,1);  tree->SetBranchAddress("R.tr.py"              , R_tr_py             );
  tree->SetBranchStatus("R.tr.pz"              ,1);  tree->SetBranchAddress("R.tr.pz"              , R_tr_pz             );
  tree->SetBranchStatus("R.tr.tg_dp"           ,1);  tree->SetBranchAddress("R.tr.tg_dp"           , R_tr_tg_dp          );
  tree->SetBranchStatus("R.tr.tg_y"            ,1);  tree->SetBranchAddress("R.tr.tg_y"            , R_tr_tg_y           );
  tree->SetBranchStatus("R.tr.tg_th"           ,1);  tree->SetBranchAddress("R.tr.tg_th"           , R_tr_tg_th          );
  tree->SetBranchStatus("R.tr.tg_ph"           ,1);  tree->SetBranchAddress("R.tr.tg_ph"           , R_tr_tg_ph          );
  tree->SetBranchStatus("R.tr.vx"              ,1);  tree->SetBranchAddress("R.tr.vx"              , R_tr_vx             );
  tree->SetBranchStatus("R.tr.vy"              ,1);  tree->SetBranchAddress("R.tr.vy"              , R_tr_vy             );
  tree->SetBranchStatus("R.tr.vz"              ,1);  tree->SetBranchAddress("R.tr.vz"              , R_tr_vz             );
  tree->SetBranchStatus("L.tr.n"               ,1);  tree->SetBranchAddress("L.tr.n"               ,&L_tr_n              );
    tree->SetBranchStatus("L.tr.flag"            ,1);  tree->SetBranchAddress("L.tr.flag"            , L_tr_flag           );
    tree->SetBranchStatus("L.tr.ndof"            ,1);  tree->SetBranchAddress("L.tr.ndof"            , L_tr_ndof           );
   tree->SetBranchStatus("L.tr.chi2"            ,1);  tree->SetBranchAddress("L.tr.chi2"            , L_tr_chi2           );
     tree->SetBranchStatus("L.tr.beta"            ,1);  tree->SetBranchAddress("L.tr.beta"            , L_tr_beta           );
     tree->SetBranchStatus("L.tr.d_x"             ,1);  tree->SetBranchAddress("L.tr.d_x"             , L_tr_d_x            );
     tree->SetBranchStatus("L.tr.d_y"             ,1);  tree->SetBranchAddress("L.tr.d_y"             , L_tr_d_y            );
     tree->SetBranchStatus("L.tr.d_th"            ,1);  tree->SetBranchAddress("L.tr.d_th"            , L_tr_d_th           );
     tree->SetBranchStatus("L.tr.d_ph"            ,1);  tree->SetBranchAddress("L.tr.d_ph"            , L_tr_d_ph           );
     tree->SetBranchStatus("L.tr.r_x"             ,1);  tree->SetBranchAddress("L.tr.r_x"             , L_tr_r_x            );
     tree->SetBranchStatus("L.tr.r_y"             ,1);  tree->SetBranchAddress("L.tr.r_y"             , L_tr_r_y            );
     tree->SetBranchStatus("L.tr.r_th"            ,1);  tree->SetBranchAddress("L.tr.r_th"            , L_tr_r_th           );
     tree->SetBranchStatus("L.tr.r_ph"            ,1);  tree->SetBranchAddress("L.tr.r_ph"            , L_tr_r_ph           );
  tree->SetBranchStatus("L.tr.x"               ,1);  tree->SetBranchAddress("L.tr.x"               , L_tr_x              );
  tree->SetBranchStatus("L.tr.y"               ,1);  tree->SetBranchAddress("L.tr.y"               , L_tr_y              );
  tree->SetBranchStatus("L.tr.th"              ,1);  tree->SetBranchAddress("L.tr.th"              , L_tr_th             );
  tree->SetBranchStatus("L.tr.ph"              ,1);  tree->SetBranchAddress("L.tr.ph"              , L_tr_ph             );
   tree->SetBranchStatus("L.tr.time"            ,1);  tree->SetBranchAddress("L.tr.time"            , L_tr_time           );
  tree->SetBranchStatus("L.tr.p"               ,1);  tree->SetBranchAddress("L.tr.p"               , L_tr_p              );
  tree->SetBranchStatus("L.tr.pathl"           ,1);  tree->SetBranchAddress("L.tr.pathl"           , L_tr_pathl          );
  tree->SetBranchStatus("L.tr.px"              ,1);  tree->SetBranchAddress("L.tr.px"              , L_tr_px             );
  tree->SetBranchStatus("L.tr.py"              ,1);  tree->SetBranchAddress("L.tr.py"              , L_tr_py             );
  tree->SetBranchStatus("L.tr.pz"              ,1);  tree->SetBranchAddress("L.tr.pz"              , L_tr_pz             );
  tree->SetBranchStatus("L.tr.tg_dp"           ,1);  tree->SetBranchAddress("L.tr.tg_dp"           , L_tr_tg_dp          );
  tree->SetBranchStatus("L.tr.tg_y"            ,1);  tree->SetBranchAddress("L.tr.tg_y"            , L_tr_tg_y           );
  tree->SetBranchStatus("L.tr.tg_th"           ,1);  tree->SetBranchAddress("L.tr.tg_th"           , L_tr_tg_th          );
  tree->SetBranchStatus("L.tr.tg_ph"           ,1);  tree->SetBranchAddress("L.tr.tg_ph"           , L_tr_tg_ph          );
  tree->SetBranchStatus("L.tr.vx"              ,1);  tree->SetBranchAddress("L.tr.vx"              , L_tr_vx             );
  tree->SetBranchStatus("L.tr.vy"              ,1);  tree->SetBranchAddress("L.tr.vy"              , L_tr_vy             );
  tree->SetBranchStatus("L.tr.vz"              ,1);  tree->SetBranchAddress("L.tr.vz"              , L_tr_vz             );
  tree->SetBranchStatus("R.s2.trpath"          ,1);  tree->SetBranchAddress("R.s2.trpath"          , R_s2_trpath         );
  tree->SetBranchStatus("L.s2.trpath"          ,1);  tree->SetBranchAddress("L.s2.trpath"          , L_s2_trpath         );
  tree->SetBranchStatus("RTDC.F1FirstHit"      ,1);  tree->SetBranchAddress("RTDC.F1FirstHit"      , RTDC_F1FirstHit);
  ENum=tree->GetEntries();
  ENum2=ENum;
  cout<<"Events: "<<ENum<<endl; 
	for(int k=0;k<ENum;k++){
		tree->GetEntry(k);

//int RS2_seg[100];
//int LS2_seg[100];
//for(int i=0;i<100;i++){
//	RS2_seg[i]=0;
//	LS2_seg[i]=0;
//}
      for(int lt=0;lt<NLtr;lt++){
      for(int rt=0;rt<NRtr;rt++){
	//int Ls2_seg = Ls2_pad[lt];
	//int Rs2_seg = Rs2_pad[rt];
//  ct=CoinCalc_gogami(Rs2_seg,Ls2_seg,rt,lt);
  int rhit=rt;
  int lhit=lt;
  int LS2_seg = Ls2_pad[lt];
  int RS2_seg = Rs2_pad[rt];
//beta
  double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+Mpi*Mpi);
  //double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  double beta_L  = L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);  

//length
  double LenL  = L_tr_pathl[lhit];// - L_s2_trpath[lhit];
  double LenR  = R_tr_pathl[rhit];// - R_s2_trpath[rhit];
  
//correction timing
  double cor_L = (LenL-3.18)/(beta_L*LightVelocity);
  double cor_R = (LenR- R_s2_trpath[rhit])/(beta_R*LightVelocity);



//  double tref_L  = LTDC_F1FirstHit[40]       * tdc_time;
  double tref_R  = RTDC_F1FirstHit[9]        * tdc_time;
  //  double rf      = LTDC_F1FirstHit[47]       * tdc_time;
  //  double rf_R    = RTDC_F1FirstHit[15]        * tdc_time;


 double timeL_R = RTDC_F1FirstHit[RS2_seg+16] * tdc_time;
 double timeR_R = RTDC_F1FirstHit[RS2_seg+48] * tdc_time; 
// double timeL_L = LTDC_F1FirstHit[LS2_seg] * tdc_time;
// double timeR_L = LTDC_F1FirstHit[LS2_seg+48] * tdc_time; 


 double toffset_R = -364.6-150.; // for H2_1
// double toffset_L = 1762.0;

// double meantime_L = tref_L - (timeL_L+timeR_L)/2.0 + toffset_L + cor_L;
 double meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R + cor_R;

 // meantime_R=100;RS2_seg=7;LS2_seg=7;


 double s2_tzero_R[16]={0,-1.02037,0.0046854,0.0534834,-0.534372,-0.60597,0.343139,-0.293262,0.267898,-0.666823,0.272364,0.0969059,-0.893806,-1.01129,-1.13495,-0.784991};
 double s2_tzero_L[16]={0,2.005,0.654413,1.34976,0.0290891,0.187557,0.00499421,-0.914343,-1.24058,-0.535878,-0.77564,2.22918,0.804909,0.607826,-0.635764,0};

 meantime_R= meantime_R - s2_tzero_R[RS2_seg] -s2_tzero_L[LS2_seg];

 double yfp_cor_R =R_tr_y[rhit]*-0.182869 +R_tr_ph[rhit]*-0.0211276;
 double yfp_cor_L = -10.0* L_tr_y[lhit] -28.0* L_tr_ph[lhit];

 // cout<<"yfp_cor_R "<<yfp_cor_R<<" yfp_corL "<<yfp_cor_L<<endl;
 
 meantime_R= meantime_R + yfp_cor_R + yfp_cor_L;
 meantime_R= meantime_R - cor_L +75.4;
 
// tr.yp_cor=0.0;
// tr.yp_cor= + yfp_cor_R + yfp_cor_L;




  //======= Nomalization ==================//
  R_tr_x[rhit]    = (R_tr_x[rhit]-XFPm)/XFPr;
  R_tr_th[rhit]   = (R_tr_th[rhit]-XpFPm)/XpFPr;
  R_tr_y[rhit]    = (R_tr_y[rhit]-YFPm)/YFPr;
  R_tr_ph[rhit]   = (R_tr_ph[rhit]-YpFPm)/YpFPr;
  R_tr_vz[rhit]   = (R_tr_vz[rhit]-Ztm)/Ztr;


  L_tr_x[lhit]    = (L_tr_x[lhit]-XFPm)/XFPr; 
  L_tr_th[lhit]   = (L_tr_th[lhit]-XpFPm)/XpFPr;
  L_tr_y[lhit]    = (L_tr_y[lhit]-YFPm)/YFPr;
  L_tr_ph[lhit]   = (L_tr_ph[lhit]-YpFPm)/YpFPr;
  L_tr_vz[lhit]   = (L_tr_vz[lhit]-Ztm)/Ztr;




 //double ctimecorR = calcf2t_3rd(PctimeR, R_tr_x[rhit],R_tr_th[rhit],R_tr_y[rhit],R_tr_ph[rhit],R_tr_vz[rhit]);
 //double ctimecorL = calcf2t_3rd(PctimeL, L_tr_x[lhit],L_tr_th[lhit],L_tr_y[lhit],L_tr_ph[lhit],L_tr_vz[lhit]);
 double ctimecorR = 0.;
 double ctimecorL = 0.;


    //========== Scaled at FP ==================//
    R_tr_x[rhit]  = R_tr_x[rhit]  * XFPr + XFPm;
    R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
    R_tr_y[rhit]  = R_tr_y[rhit]  * YFPr + YFPm;
    R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr + YpFPm;

    L_tr_x[lhit]  = L_tr_x[lhit]  * XFPr + XFPm;
    L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
    L_tr_y[lhit]  = L_tr_y[lhit]  * YFPr + YFPm;
    L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr + YpFPm;    



 double ctime = - meantime_R + ctimecorL + ctimecorR;
 double ctime_before = - meantime_R;

// tr.ctimecorR=ctimecorR;
// tr.ctimecorL=ctimecorL;
 // double time_rf = rf - meantime;
 // double time_rf_R = rf - meantime_R -tref_R;
 // double ctime = - meantime_R + mean_time - kcenter;





//  cout<<"======== check coin time ====== "<<endl;
//  cout<<"s2R off "<<s2_tzero_R[RS2_seg]<<" s2L off "<<s2_tzero_L[LS2_seg]<<" meantime_R "<<meantime_R<<" ctime "<<ctime<<endl;
//  cout<<" meantime_L "<<meantime_L<<" meantime_R "<<meantime_R<<" corL "<<cor_L<<" corR "<<cor_R<<" ctime "<<ctime<<endl;
//  cout<<" ctimecorL "<<ctimecorL<<" ctimecorR "<<ctimecorR<<" yfp_cor_R "<<yfp_cor_R<<" yfp_cor_L "<<yfp_cor_L<<endl;

 ctime=ctime -1.4;
 ctime_before=ctime_before -1.4;


// tr.ct_gb=-ctime;

 if(3600.0<ctime && ctime<3665){
   ctime = ctime - 3637.88 - 12.76;
   ctime = ctime - 12.0-3.1;
 }
//else ctime=-9999.;
 if(3600.0<ctime_before && ctime_before<3665){
   ctime_before = ctime_before - 3637.88 - 12.76;
   ctime_before = ctime_before - 12.0-3.1;
 }

 ctime=-ctime;
 ctime_before=-ctime_before;

 ct=ctime;


	if(k==ev*100){
cout << "Event (Fill) : " << k << "/" << ENum << endl;
cout << "ct="<<ct<<endl;
	ev += 1;
	}



			}//NRtr
			}//NLtr
	}//ENum
//---------------------------//
//--Event by Event Analysis--//
//----------END--------------//




//-----No Z Cut-----//
 hcoin_bg_fom_noZ->Scale(40./80.);
 hcoin_wo_bg_fom_noZ->Add(hcoin_k_fom_noZ,hcoin_bg_fom_noZ,1.0,-1.0);
 fp_noZ=new TF1("fp_noZ","gausn(0)",min_coin_c,max_coin_c);
 fp_noZ->SetNpx(2000);
 fpi_noZ =new TF1("fpi_noZ","gausn(0)+gausn(3)",def_mean_pi-3*def_sig_pi,def_mean_pi+6*def_sig_pi);
 fpi_noZ->SetNpx(2000);
 fk_noZ=new TF1("fk_noZ","gausn(0)",min_coin_c,max_coin_c);
 fk_noZ->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noZ->Fit("fp_noZ","Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//n_p_noZ=fp_noZ->GetParameter(0);
 mean_p_noZ=fp_noZ->GetParameter(1);
 sig_p_noZ=fp_noZ->GetParameter(2);
 //center_p=mean_p_noZ;
 //range_p=2*sig_p_noZ;
 //range_p=1.435;
 center_p=def_mean_p;
 range_p=2*def_sig_p;
 n_p_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_p-range_p),hcoin_wo_bg_fom_noZ->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi_noZ->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 //fpi_noZ->SetParameters(1000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom_noZ->Fit("fpi_noZ","Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=fpi_noZ->Integral(def_mean_pi-5*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=n_pi_noZ/((max_coin_c-min_coin_c)/bin_coin_c);
 mean_pi_noZ=fpi_noZ->GetParameter(1);
 sig_pi_noZ=fpi_noZ->GetParameter(2);
 hcoin_pi_noZ->FillRandom("fpi_noZ",n_pi_noZ*1000.);
 hcoin_pi_noZ->Scale(1./1000.);
 hcoin_wo_bg_fom_noZ->Add(hcoin_wo_bg_fom_noZ,hcoin_pi_noZ,1.,-1.);
 //center_pi=mean_pi_noZ;
 //range_pi=2*sig_pi_noZ;
 //range_pi=0.7;
 center_pi=def_mean_pi;
 range_pi=2*def_sig_pi;
 n_pi_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noZ->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",-1,1);
 //n_k_noZ=fk_noZ->GetParameter(0);
 mean_k_noZ=fk_noZ->GetParameter(1);
// mean_k_noZ=0;
 sig_k_noZ=fk_noZ->GetParameter(2);
 //center_k=mean_k_noZ;
 center_k=0.;
 //range_k=2*sig_k_noZ;
 range_k=1.;

 n_k_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_k-range_k),hcoin_wo_bg_fom_noZ->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noZ =new TF1("fcoin_noZ","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noZ->SetNpx(2000);
 fcoin_noZ->SetTitle("Cointime w/o AC cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noZ->SetParameters(n_pi_noZ,mean_pi_noZ,sig_pi_noZ,n_k_noZ,mean_k_noZ,sig_k_noZ,n_p_noZ,mean_p_noZ,sig_p_noZ);
 //hcoin_wo_bg_fom_noZ->Fit("fcoin_woAC","Rq0","0",min_coin_c,max_coin_c);
 //n_pi_noZ=fcoin_noZ->GetParameter(0);//Npi_nocut
 //mean_pi_noZ=fcoin_noZ->GetParameter(1);
 //sig_pi_noZ=fcoin_noZ->GetParameter(2);
 //n_k_noZ=fcoin_noZ->GetParameter(3);//Nk_nocut
 //mean_k_noZ=fcoin_noZ->GetParameter(4);
 //sig_k_noZ=fcoin_noZ->GetParameter(5);
 //n_p_noZ=fcoin_noZ->GetParameter(6);//Np_nocut
 //mean_p_noZ=fcoin_noZ->GetParameter(7);
 //sig_p_noZ=fcoin_noZ->GetParameter(8);
cout<<"n_pi_noZ="<<n_pi_noZ<<"n_k_noZ="<<n_k_noZ<<"n_p_noZ="<<n_p_noZ
<<"mean_pi_noZ="<<mean_pi_noZ<<"sig_pi_noZ="<<sig_pi_noZ<<"mean_k_noZ="<<mean_k_noZ<<"sig_k_noZ="<<sig_k_noZ<<"mean_p_noZ="<<mean_p_noZ<<"sig_p_noZ="<<sig_p_noZ<<endl;

 

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom_noZ->Scale(2./80.);
 hmm_pibg_fom_noZ->Scale(0.7/80.);
 hmm_wo_bg_fom_noZ->Add(hmm_L_fom_noZ,hmm_bg_fom_noZ,1.0,-1.0);
 hmm_pi_wobg_fom_noZ->Add(hmm_pi_fom_noZ,hmm_pibg_fom_noZ,1.0,-1.0);


// fmmbg_noZ=new TF1("fmmbg_noZ","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg_noZ=new TF1("fmmbg_noZ",F_Voigt,min_mm,max_mm,4);
 fmmbg_noZ=new TF1("fmmbg_noZ","pol4",-0.05,0.15);
// fmmbg_noZ->SetParameters(5,0.05,0.05,0.01);
// fmmbg_noZ->SetParLimits(0,0.,100000.);//positive
// fmmbg_noZ->SetParLimits(3,0.,100.);//positive
 fmmbg_noZ->SetNpx(2000);
 //fmmbg_noZ->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg_noZ->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg_noZ->SetParameter(1,0.05);
// fmmbg_noZ->SetParameter(2,0.03);
 fL_noZ=new TF1("fL_noZ","gausn(0)",min_mm,max_mm);
 fL_noZ->SetNpx(2000);
 fL_noZ->SetParLimits(2,0.,0.01);
 fS_noZ=new TF1("fS_noZ","gausn(0)",min_mm,max_mm);
 fS_noZ->SetNpx(2000);
 fS_noZ->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noZ=new TF1("fmm_noZ","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
 fmm_noZ->SetNpx(2000);
 fmm_noZ->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm_noZ->SetParLimits(0,0.,1000000.);//positive
 fmm_noZ->SetParLimits(3,0.,1000000.);//positive
// fmm_noZ->SetParLimits(3,0.,100.);//positive//Voigt
// fmm_noZ->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm_noZ->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm_noZ->SetParameter(1,def_mean_L);
// fmm_noZ->SetParameter(4,def_mean_S);

 hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","0",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
 mean_L_noZ=fL_noZ->GetParameter(1);
 sig_L_noZ=fL_noZ->GetParameter(2);
 center_L=def_mean_L;
 range_L=2*def_sig_L;

 hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","0",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
 mean_S_noZ=fS_noZ->GetParameter(1);
 sig_S_noZ=fS_noZ->GetParameter(2);
 center_S=def_mean_S;
 range_S=2*def_sig_S;


 //------- Fitting ----------//

cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noZ->Fit("fmmbg_noZ","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom_noZ->Fit("fmmbg_noZ","Rq0","",-0.05,0.15);
//double fmmbga = fmmbg_noZ->GetParameter(0);
//double fmmbgb = fmmbg_noZ->GetParameter(1);
//double fmmbgc = fmmbg_noZ->GetParameter(2);
//double fmmbgd = fmmbg_noZ->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm_noZ->FixParameter(1,b);//mean
// fmm_noZ->FixParameter(2,c);//sigma
// fmm_noZ->FixParameter(3,d);//lg
//double e = fmmbg_noZ->GetParameter(4);
//double f = fmmbg_noZ->GetParameter(5);
// fmm_noZ->SetParameter(0,a);
// fmm_noZ->SetParameter(1,b);
// fmm_noZ->SetParameter(2,c);
// fmm_noZ->SetParameter(3,d);
// fmm_noZ->SetParameter(4,e);
// fmm_noZ->SetParameter(5,f);
// fmm_noZ->SetParameter(6,500);
 fmm_noZ->SetParameter(1,def_mean_L);
 fmm_noZ->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm_noZ->SetParameter(2,def_sig_L);
 fmm_noZ->SetParLimits(2,0.,2*def_sig_L);
// fmm_noZ->SetParameters(9,100);
 fmm_noZ->SetParameter(4,def_mean_S);
 fmm_noZ->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm_noZ->SetParameter(5,def_sig_S);
 fmm_noZ->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom_noZ->Fit("fmm_noZ","Rq0","0",-0.05,0.15);
double fmm_noZpar0 = fmm_noZ->GetParameter(0);cout<<"fmm_noZ[0]="<<fmm_noZpar0<<endl;//area(L)
double fmm_noZpar1 = fmm_noZ->GetParameter(1);cout<<"fmm_noZ[1]="<<fmm_noZpar1<<endl;//mean(L)
double fmm_noZpar2 = fmm_noZ->GetParameter(2);cout<<"fmm_noZ[2]="<<fmm_noZpar2<<endl;//sigma(L)
double fmm_noZpar3 = fmm_noZ->GetParameter(3);cout<<"fmm_noZ[3]="<<fmm_noZpar3<<endl;//area(S)
double fmm_noZpar4 = fmm_noZ->GetParameter(4);cout<<"fmm_noZ[4]="<<fmm_noZpar4<<endl;//mean(S)
double fmm_noZpar5 = fmm_noZ->GetParameter(5);cout<<"fmm_noZ[5]="<<fmm_noZpar5<<endl;//sigma(S)
double fmm_noZpar6 = fmm_noZ->GetParameter(6);cout<<"fmm_noZ[6]="<<fmm_noZpar6<<endl;//poly_const
double fmm_noZpar7 = fmm_noZ->GetParameter(7);cout<<"fmm_noZ[7]="<<fmm_noZpar7<<endl;//poly_x
double fmm_noZpar8 = fmm_noZ->GetParameter(8);cout<<"fmm_noZ[8]="<<fmm_noZpar8<<endl;//poly_x^2
double fmm_noZpar9 = fmm_noZ->GetParameter(9);cout<<"fmm_noZ[9]="<<fmm_noZpar9<<endl;//poly_x^3
double fmm_noZpar10 = fmm_noZ->GetParameter(10);cout<<"fmm_noZ[10]="<<fmm_noZpar10<<endl;//poly_x^4
 fmmbg_noZ->SetParameters(fmm_noZpar6,fmm_noZpar7,fmm_noZpar8,fmm_noZpar9,fmm_noZpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noZ=fL_noZ->GetParameter(0);
// mean_L_noZ=fL_noZ->GetParameter(1);
// sig_L_noZ=fL_noZ->GetParameter(2);
 mean_L_noZ=def_mean_L;
 mean_S_noZ=def_mean_S;
 sig_L_noZ=def_sig_L;
 sig_S_noZ=def_sig_S;
 n_L_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_L-range_L),hmm_wo_bg_fom_noZ->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L_noZ<<endl;
 double integralL=fmmbg_noZ->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom_noZ->FindBin(center_L+range_L)-hmm_wo_bg_fom_noZ->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L_noZ=n_L_noZ-integralL;
cout<<"after(L):: "<<n_L_noZ<<endl;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,5)-pow(mean_L_noZ-2*sig_L_noZ,5))*a/5;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,4)-pow(mean_L_noZ-2*sig_L_noZ,4))*b/4;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,3)-pow(mean_L_noZ-2*sig_L_noZ,3))*c/3;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,2)-pow(mean_L_noZ-2*sig_L_noZ,2))*d/2;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,1)-pow(mean_L_noZ-2*sig_L_noZ,1))*e;
//
 n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_S-range_S),hmm_wo_bg_fom_noZ->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S_noZ<<endl;
 double integralS=fmmbg_noZ->Integral(center_S-range_S,center_S+range_S);
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom_noZ->FindBin(center_S+range_S)-hmm_wo_bg_fom_noZ->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
 if(integralS>0)n_S_noZ-=integralS;
cout<<"after(S):: "<<n_S_noZ<<endl;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,5)-pow(mean_S_noZ-2*sig_S_noZ,5))*a/5;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,4)-pow(mean_S_noZ-2*sig_S_noZ,4))*b/4;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,3)-pow(mean_S_noZ-2*sig_S_noZ,3))*c/3;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,2)-pow(mean_S_noZ-2*sig_S_noZ,2))*d/2;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,1)-pow(mean_S_noZ-2*sig_S_noZ,1))*e;

 cout<<"n_L"<<n_L_noZ<<endl;
cout<<"mean_L"<<mean_L_noZ<<endl;
cout<<"sig_L"<<sig_L_noZ<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S_noZ=fS_noZ->GetParameter(0);
// mean_S_noZ=fS_noZ->GetParameter(1);
// mean_S_noZ=def_mean_S;
// sig_S_noZ=def_sig_S;
// sig_S_noZ=fS_noZ->GetParameter(2);
// n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ-2*sig_S_noZ),hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ+2*sig_S_noZ));
 cout<<"n_S"<<n_S_noZ<<endl;
cout<<"mean_S"<<mean_S_noZ<<endl;
cout<<"sig_S"<<sig_S_noZ<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//



 //----------------------------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//

	for(int i=0;i<nth;i++){

			
//-----Background subtraction-----//
//------------cointime------------//
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom[i]->Scale(40./80.);
 hcoin_wo_bg_fom[i]->Add(hcoin_k_fom[i],hcoin_bg_fom[i],1.0,-1.0);

 fp[i]=new TF1(Form("fp[%d]",i),"gausn(0)",min_coin_c,max_coin_c);
 fp[i]->SetNpx(2000);
 fpi[i] =new TF1(Form("fpi[%d]",i),"gausn(0)+gausn(3)",min_coin_c,max_coin_c);
 fpi[i]->SetNpx(2000);
 fk[i]=new TF1(Form("fk[%d]",i),"gausn(0)",min_coin_c,max_coin_c);
 fk[i]->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom[i]->Fit(Form("fp[%d]",i),"Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 //n_p[i]=fp[i]->GetParameter(0);
 mean_p[i]=fp[i]->GetParameter(1);
 sig_p[i]=fp[i]->GetParameter(2);
 n_p[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_p-range_p),hcoin_wo_bg_fom[i]->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi[i]->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom[i]->Fit(Form("fpi[%d]",i),"Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i]=fpi[i]->Integral(min_coin_c,max_coin_c);
 n_pi[i]=n_pi[i]/((max_coin_c-min_coin_c)/bin_coin_c);
 //n_pi[i]=n_pi[i]/(2*range_pi/(hcoin_wo_bg_fom[i]->FindBin(center_pi+range_pi)-hcoin_wo_bg_fom[i]->FindBin(center_pi-range_pi)));
 mean_pi[i]=fpi[i]->GetParameter(1);
cout << mean_pi[i] << endl;
 sig_pi[i]=fpi[i]->GetParameter(2);
 hcoin_pi[i]->FillRandom(Form("fpi[%d]",i),n_pi[i]);
 //hcoin_wo_bg_fom[i]->Add(hcoin_wo_bg_fom[i],hcoin_pi[i],1.,-1.);
n_pi[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_pi-range_pi),hcoin_wo_bg_fom[i]->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
// hcoin_wo_bg_fom[i]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom[i]->Fit(Form("fk[%d]",i),"Rq0","0",-1,1);
// n_k[i]=fk[i]->GetParameter(0);
 mean_k[i]=fk[i]->GetParameter(1);
 sig_k[i]=fk[i]->GetParameter(2);
 n_k[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_k-range_k),hcoin_wo_bg_fom[i]->FindBin(center_k+range_k));

 // n_k[i][th1][1]=hcoin_k_fom[i][th1]->Integral(hcoin_k_fom[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_k_fom[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 fcoin[i] =new TF1(Form("fcoin[%d]",i),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i]->SetNpx(2000);
 fcoin[i]->SetParameters(n_pi[i],mean_pi[i],sig_pi[i],n_k[i],mean_k[i],sig_k[i],n_p[i],mean_p[i],sig_p[i]);
 //hcoin_wo_bg_fom[i]->Fit(Form("fcoin[%d][%d][%d]",i,j,l),"Rq0","0",min_coin_c,max_coin_c);
 //n_pi[i]=fcoin[i]->GetParameter(0);//Npi_nocut
 //mean_pi[i]=fcoin[i]->GetParameter(1);
 //sig_pi[i]=fcoin[i]->GetParameter(2);
 //n_k[i]=fcoin[i]->GetParameter(3);//Nk_nocut
 //mean_k[i]=fcoin[i]->GetParameter(4);
 //sig_k[i]=fcoin[i]->GetParameter(5);
 //n_p[i]=fcoin[i]->GetParameter(6);//Np_nocut
 //mean_p[i]=fcoin[i]->GetParameter(7);
 //sig_p[i]=fcoin[i]->GetParameter(8);

cout<<"n_pi["<<i<<"]="<<n_pi[i]<<endl;
cout<<"n_k["<<i<<"]="<<n_k[i]<<endl;
cout<<"n_p["<<i<<"]="<<n_p[i]<<endl;

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom[i]->Scale(2./80.);
 hmm_pibg_fom[i]->Scale(0.7/80.);
 hmm_wo_bg_fom[i]->Add(hmm_L_fom[i],hmm_bg_fom[i],1.0,-1.0);
 hmm_pi_wobg_fom[i]->Add(hmm_pi_fom[i],hmm_pibg_fom[i],1.0,-1.0);


// fmmbg[i]=new TF1("fmmbg[i]","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg[i]=new TF1("fmmbg[i]",F_Voigt,min_mm,max_mm,4);
 fmmbg[i]=new TF1("fmmbg[i]","pol4",-0.05,0.15);
// fmmbg[i]->SetParameters(5,0.05,0.05,0.01);
// fmmbg[i]->SetParLimits(0,0.,100000.);//positive
// fmmbg[i]->SetParLimits(3,0.,100.);//positive
 fmmbg[i]->SetNpx(2000);
 //fmmbg[i]->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg[i]->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg[i]->SetParameter(1,0.05);
// fmmbg[i]->SetParameter(2,0.03);
 fL[i]=new TF1("fL[i]","gausn(0)",min_mm,max_mm);
 fL[i]->SetNpx(2000);
 fL[i]->SetParLimits(2,0.,0.01);
 fS[i]=new TF1("fS[i]","gausn(0)",min_mm,max_mm);
 fS[i]->SetNpx(2000);
 fS[i]->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm[i]=new TF1("fmm[i]","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
 fmm[i]->SetNpx(2000);
 fmm[i]->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm[i]->SetParLimits(0,0.,1000000.);//positive
 fmm[i]->SetParLimits(3,0.,1000000.);//positive
// fmm[i]->SetParLimits(3,0.,100.);//positive//Voigt
// fmm[i]->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm[i]->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm[i]->SetParameter(1,def_mean_L);
// fmm[i]->SetParameter(4,def_mean_S);



 //------- Fitting ----------//

//cout<<"fmmbg fit start"<<endl;
// mmbg_eff[i]->Fit("fmmbg[i]","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom[i]->Fit("fmmbg[i]","Rq0","",-0.05,0.15);
//double fmmbga = fmmbg[i]->GetParameter(0);
//double fmmbgb = fmmbg[i]->GetParameter(1);
//double fmmbgc = fmmbg[i]->GetParameter(2);
//double fmmbgd = fmmbg[i]->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm[i]->FixParameter(1,b);//mean
// fmm[i]->FixParameter(2,c);//sigma
// fmm[i]->FixParameter(3,d);//lg
//double e = fmmbg[i]->GetParameter(4);
//double f = fmmbg[i]->GetParameter(5);
// fmm[i]->SetParameter(0,a);
// fmm[i]->SetParameter(1,b);
// fmm[i]->SetParameter(2,c);
// fmm[i]->SetParameter(3,d);
// fmm[i]->SetParameter(4,e);
// fmm[i]->SetParameter(5,f);
// fmm[i]->SetParameter(6,500);
 fmm[i]->SetParameter(1,def_mean_L);
 fmm[i]->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm[i]->SetParameter(2,def_sig_L);
 fmm[i]->SetParLimits(2,0.,2*def_sig_L);
// fmm[i]->SetParameters(9,100);
 fmm[i]->SetParameter(4,def_mean_S);
 fmm[i]->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm[i]->SetParameter(5,def_sig_S);
 fmm[i]->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom[i]->Fit("fmm[i]","Rq0","0",-0.05,0.15);
//double fmmpar0 = fmm[i]->GetParameter(0);//cout<<"fmm[0]="<<fmmpar0<<endl;//area(L)
//double fmmpar1 = fmm[i]->GetParameter(1);//cout<<"fmm[1]="<<fmmpar1<<endl;//mean(L)
//double fmmpar2 = fmm[i]->GetParameter(2);//cout<<"fmm[2]="<<fmmpar2<<endl;//sigma(L)
//double fmmpar3 = fmm[i]->GetParameter(3);//cout<<"fmm[3]="<<fmmpar3<<endl;//area(S)
//double fmmpar4 = fmm[i]->GetParameter(4);//cout<<"fmm[4]="<<fmmpar4<<endl;//mean(S)
//double fmmpar5 = fmm[i]->GetParameter(5);//cout<<"fmm[5]="<<fmmpar5<<endl;//sigma(S)
mean_L[i]=fmm[i]->GetParameter(1);
sig_L[i]=fmm[i]->GetParameter(2);
mean_S[i]=fmm[i]->GetParameter(4);
sig_S[i]=fmm[i]->GetParameter(5);
double fmmpar6 = fmm[i]->GetParameter(6);//cout<<"fmm[6]="<<fmmpar6<<endl;//poly_const
double fmmpar7 = fmm[i]->GetParameter(7);//cout<<"fmm[7]="<<fmmpar7<<endl;//poly_x
double fmmpar8 = fmm[i]->GetParameter(8);//cout<<"fmm[8]="<<fmmpar8<<endl;//poly_x^2
double fmmpar9 = fmm[i]->GetParameter(9);//cout<<"fmm[9]="<<fmmpar9<<endl;//poly_x^3
double fmmpar10 = fmm[i]->GetParameter(10);//cout<<"fmm[10]="<<fmmpar10<<endl;//poly_x^4
 fmmbg[i]->SetParameters(fmmpar6,fmmpar7,fmmpar8,fmmpar9,fmmpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom[i]->Fit("fL[i]","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L[i]=fL[i]->GetParameter(0);
// mean_L[i]=fL[i]->GetParameter(1);
// sig_L[i]=fL[i]->GetParameter(2);
 n_L[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(center_L-range_L),hmm_wo_bg_fom[i]->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L[i]<<endl;
 double integralL=fmmbg[i]->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom[i]->FindBin(center_L+range_L)-hmm_wo_bg_fom[i]->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L[i]=n_L[i]-integralL;
else cout<<"negative BG: ("<<i<<")"<<endl;
cout<<"after(L):: "<<n_L[i]<<endl;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],5)-pow(mean_L[i]-2*sig_L[i],5))*a/5;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],4)-pow(mean_L[i]-2*sig_L[i],4))*b/4;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],3)-pow(mean_L[i]-2*sig_L[i],3))*c/3;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],2)-pow(mean_L[i]-2*sig_L[i],2))*d/2;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],1)-pow(mean_L[i]-2*sig_L[i],1))*e;
//
 n_S[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(center_S-range_S),hmm_wo_bg_fom[i]->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S[i]<<endl;
 double integralS=fmmbg[i]->Integral(center_S-range_S,center_S+range_S);
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom[i]->FindBin(center_S+range_S)-hmm_wo_bg_fom[i]->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
 if(integralS>0)n_S[i]-=integralS;
else cout<<"negative BG: ("<<i<<")"<<endl;
cout<<"after(S):: "<<n_S[i]<<endl;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],5)-pow(mean_S[i]-2*sig_S[i],5))*a/5;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],4)-pow(mean_S[i]-2*sig_S[i],4))*b/4;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],3)-pow(mean_S[i]-2*sig_S[i],3))*c/3;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],2)-pow(mean_S[i]-2*sig_S[i],2))*d/2;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],1)-pow(mean_S[i]-2*sig_S[i],1))*e;

 cout<<"n_L"<<n_L[i]<<endl;
cout<<"mean_L"<<mean_L[i]<<endl;
cout<<"sig_L"<<sig_L[i]<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom[i]->Fit("fS[i]","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S[i]=fS[i]->GetParameter(0);
// mean_S[i]=fS[i]->GetParameter(1);
// mean_S[i]=def_mean_S;
// sig_S[i]=def_sig_S;
// sig_S[i]=fS[i]->GetParameter(2);
// n_S[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(mean_S[i]-2*sig_S[i]),hmm_wo_bg_fom[i]->FindBin(mean_S[i]+2*sig_S[i]));
 cout<<"n_S"<<n_S[i]<<endl;
cout<<"mean_S"<<mean_S[i]<<endl;
cout<<"sig_S"<<sig_S[i]<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//

	
	
	if(n_pi[i]>0.){}else{n_pi[i]=1.;}
	if(n_k[i]>0.){}else{n_k[i]=1.;}
	if(n_p[i]>0.){}else{n_p[i]=1.;}
	if(n_L[i]>0.){}else{n_L[i]=1.;}
	if(n_S[i]>0.){}else{n_S[i]=1.;}




	if(n_pi[i]>n_pi_noZ){
		for(int fill=0;fill<n_pi_noZ;fill++){h_pisr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_pi[i];fill++){h_pisr1->Fill(zver[i]+0.0025);}}

	if(n_k[i]>n_k_noZ){
		for(int fill=0;fill<n_k_noZ;fill++){h_ksr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_k[i];fill++){h_ksr1->Fill(zver[i]+0.0025);}}
	if(n_p[i]>n_p_noZ){
		for(int fill=0;fill<n_p_noZ;fill++){h_psr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_p[i];fill++){h_psr1->Fill(zver[i]+0.0025);}}

	if(n_L[i]>n_L_noZ){
		for(int fill=0;fill<n_L_noZ;fill++){h_Lsr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_L[i];fill++){h_Lsr1->Fill(zver[i]+0.0025);}}
	if(n_S[i]>n_S_noZ){
		for(int fill=0;fill<n_S_noZ;fill++){h_Ssr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_S[i];fill++){h_Ssr1->Fill(zver[i]+0.0025);}}

	for(int fill=0;fill<n_pi_noZ;fill++) h_pitot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_k_noZ;fill++) h_ktot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_p_noZ;fill++) h_ptot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_L_noZ;fill++) h_Ltot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_S_noZ;fill++) h_Stot1->Fill(zver[i]+0.0025);





	

	}//for i

//ofstream fout(Form("SR_z_%d.dat",(int)zver[0]));
ofstream fout(pdfname.c_str());
cout<<"writing in "<<pdfname<<"..."<<endl;

fout<<n_pi_noZ<<" "<<n_k_noZ<<" "<<n_p_noZ<<" "<<n_L_noZ<<" "<<n_S_noZ<<endl;
for(int i=0;i<nth;i++){
		fout<<n_pi[i]<<" "<<n_k[i]<<" "<<n_p[i]<<" "<<n_L[i]<<" "<<n_S[i]<<"|Ave(z)|<"<<zver[i] <<endl;
}
		


cout << "TEfficiency START!" << endl;
cout<<"pEff1:"<<endl;
if(TEfficiency::CheckConsistency(*h_pisr1,*h_pitot1,"w")){
pEff1 = new TEfficiency(*h_pisr1,*h_pitot1);
}
cout<<"pEff2:"<<endl;
if(TEfficiency::CheckConsistency(*h_ksr1,*h_ktot1,"w")){
pEff2 = new TEfficiency(*h_ksr1,*h_ktot1);
}
cout<<"pEff3:"<<endl;
if(TEfficiency::CheckConsistency(*h_psr1,*h_ptot1,"w")){
pEff3 = new TEfficiency(*h_psr1,*h_ptot1);
}

//-----------Lambda-------------------//
cout<<"pEff4:"<<endl;
if(TEfficiency::CheckConsistency(*h_Lsr1,*h_Ltot1,"w")){
pEff4 = new TEfficiency(*h_Lsr1,*h_Ltot1);
}

//-----------Sigma-------------------//
cout<<"pEff5:"<<endl;
if(TEfficiency::CheckConsistency(*h_Ssr1,*h_Stot1,"w")){
pEff5 = new TEfficiency(*h_Ssr1,*h_Stot1);
}





//------------DRAW--------------//
  cout<<"DRAW START"<<endl;
TCanvas *c1 = new TCanvas("c1","c1",800.,800.);
TCanvas *c2 = new TCanvas("c2","c2",800.,800.);
TCanvas *c3 = new TCanvas("c3","c3",800.,800.);
TCanvas *c4 = new TCanvas("c4","c4",800.,800.);
TCanvas *c5 = new TCanvas("c5","c5",800.,800.);
TCanvas *c6 = new TCanvas("c6","c6",800.,800.);
TCanvas *c7 = new TCanvas("c7","c7",800.,800.);
TCanvas *c8 = new TCanvas("c8","c8",800.,800.);
TCanvas *c9 = new TCanvas("c9","c9",800.,800.);
TCanvas *c10 = new TCanvas("c10","c10",800.,800.);

c1->cd();
 pEff1->Draw("");
c2->cd();
 pEff2->Draw("");
c3->cd();
 pEff3->Draw("");
c4->cd();
 pEff4->Draw("");
c5->cd();
 pEff5->Draw("");
	




//------------PRINT--------------//
  cout<<"Print START"<<endl;
  cout<<"pdf name: "<<pdfname<<endl;
 c1->Print(Form("%s[",pdfname.c_str()));
 c1->Print(Form("%s",pdfname.c_str()));
 c2->Print(Form("%s",pdfname.c_str()));
 c3->Print(Form("%s",pdfname.c_str()));
 c4->Print(Form("%s",pdfname.c_str()));
 c5->Print(Form("%s",pdfname.c_str()));
 c6->Print(Form("%s",pdfname.c_str()));
 c7->Print(Form("%s",pdfname.c_str()));
 c8->Print(Form("%s",pdfname.c_str()));
 c9->Print(Form("%s",pdfname.c_str()));
 c10->Print(Form("%s",pdfname.c_str()));
 c10->Print(Form("%s]",pdfname.c_str()));




  cout<<"Well done!"<<endl;
 return 0;

}


//double CoinCalc_gogami(int RS2_seg, int LS2_seg,int rhit, int lhit){
//
////beta
//  double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+Mpi*Mpi);
//  //double beta_R  = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
//  double beta_L  = L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);  
//
////length
//  double LenL  = L_tr_pathl[lhit];// - L_s2_trpath[lhit];
//  double LenR  = R_tr_pathl[rhit];// - R_s2_trpath[rhit];
//  
////correction timing
//  double cor_L = (LenL-3.18)/(beta_L*LightVelocity);
//  double cor_R = (LenR- R_s2_trpath[rhit])/(beta_R*LightVelocity);
//
//
//
////  double tref_L  = LTDC_F1FirstHit[40]       * tdc_time;
//  double tref_R  = RTDC_F1FirstHit[9]        * tdc_time;
//  //  double rf      = LTDC_F1FirstHit[47]       * tdc_time;
//  //  double rf_R    = RTDC_F1FirstHit[15]        * tdc_time;
//
//
// double timeL_R = RTDC_F1FirstHit[RS2_seg+16] * tdc_time;
// double timeR_R = RTDC_F1FirstHit[RS2_seg+48] * tdc_time; 
//// double timeL_L = LTDC_F1FirstHit[LS2_seg] * tdc_time;
//// double timeR_L = LTDC_F1FirstHit[LS2_seg+48] * tdc_time; 
//
//
// double toffset_R = -364.6-150.; // for H2_1
//// double toffset_L = 1762.0;
//
//// double meantime_L = tref_L - (timeL_L+timeR_L)/2.0 + toffset_L + cor_L;
// double meantime_R = tref_R - (timeL_R+timeR_R)/2.0 + toffset_R + cor_R;
//
// // meantime_R=100;RS2_seg=7;LS2_seg=7;
//
//
// double s2_tzero_R[16]={0,-1.02037,0.0046854,0.0534834,-0.534372,-0.60597,0.343139,-0.293262,0.267898,-0.666823,0.272364,0.0969059,-0.893806,-1.01129,-1.13495,-0.784991};
// double s2_tzero_L[16]={0,2.005,0.654413,1.34976,0.0290891,0.187557,0.00499421,-0.914343,-1.24058,-0.535878,-0.77564,2.22918,0.804909,0.607826,-0.635764,0};
//
// meantime_R= meantime_R - s2_tzero_R[RS2_seg] -s2_tzero_L[LS2_seg];
//
// double yfp_cor_R =R_tr_y[rhit]*-0.182869 +R_tr_ph[rhit]*-0.0211276;
// double yfp_cor_L = -10.0* L_tr_y[lhit] -28.0* L_tr_ph[lhit];
//
// // cout<<"yfp_cor_R "<<yfp_cor_R<<" yfp_corL "<<yfp_cor_L<<endl;
// 
// meantime_R= meantime_R + yfp_cor_R + yfp_cor_L;
// meantime_R= meantime_R - cor_L +75.4;
// 
// tr.yp_cor=0.0;
// tr.yp_cor= + yfp_cor_R + yfp_cor_L;
//
//
//
//
//  //======= Nomalization ==================//
//  R_tr_x[rhit]    = (R_tr_x[rhit]-XFPm)/XFPr;
//  R_tr_th[rhit]   = (R_tr_th[rhit]-XpFPm)/XpFPr;
//  R_tr_y[rhit]    = (R_tr_y[rhit]-YFPm)/YFPr;
//  R_tr_ph[rhit]   = (R_tr_ph[rhit]-YpFPm)/YpFPr;
//  R_tr_vz[rhit]   = (R_tr_vz[rhit]-Ztm)/Ztr;
//
//
//  L_tr_x[lhit]    = (L_tr_x[lhit]-XFPm)/XFPr; 
//  L_tr_th[lhit]   = (L_tr_th[lhit]-XpFPm)/XpFPr;
//  L_tr_y[lhit]    = (L_tr_y[lhit]-YFPm)/YFPr;
//  L_tr_ph[lhit]   = (L_tr_ph[lhit]-YpFPm)/YpFPr;
//  L_tr_vz[lhit]   = (L_tr_vz[lhit]-Ztm)/Ztr;
//
//
//
//
// //double ctimecorR = calcf2t_3rd(PctimeR, R_tr_x[rhit],R_tr_th[rhit],R_tr_y[rhit],R_tr_ph[rhit],R_tr_vz[rhit]);
// //double ctimecorL = calcf2t_3rd(PctimeL, L_tr_x[lhit],L_tr_th[lhit],L_tr_y[lhit],L_tr_ph[lhit],L_tr_vz[lhit]);
// double ctimecorR = 0.;
// double ctimecorL = 0.;
//
//
//    //========== Scaled at FP ==================//
//    R_tr_x[rhit]  = R_tr_x[rhit]  * XFPr + XFPm;
//    R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
//    R_tr_y[rhit]  = R_tr_y[rhit]  * YFPr + YFPm;
//    R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr + YpFPm;
//
//    L_tr_x[lhit]  = L_tr_x[lhit]  * XFPr + XFPm;
//    L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
//    L_tr_y[lhit]  = L_tr_y[lhit]  * YFPr + YFPm;
//    L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr + YpFPm;    
//
//
//
// double ctime = - meantime_R + ctimecorL + ctimecorR;
// double ctime_before = - meantime_R;
//
// tr.ctimecorR=ctimecorR;
// tr.ctimecorL=ctimecorL;
// // double time_rf = rf - meantime;
// // double time_rf_R = rf - meantime_R -tref_R;
// // double ctime = - meantime_R + mean_time - kcenter;
//
//
//
//
//
////  cout<<"======== check coin time ====== "<<endl;
////  cout<<"s2R off "<<s2_tzero_R[RS2_seg]<<" s2L off "<<s2_tzero_L[LS2_seg]<<" meantime_R "<<meantime_R<<" ctime "<<ctime<<endl;
////  cout<<" meantime_L "<<meantime_L<<" meantime_R "<<meantime_R<<" corL "<<cor_L<<" corR "<<cor_R<<" ctime "<<ctime<<endl;
////  cout<<" ctimecorL "<<ctimecorL<<" ctimecorR "<<ctimecorR<<" yfp_cor_R "<<yfp_cor_R<<" yfp_cor_L "<<yfp_cor_L<<endl;
//
// ctime=ctime -1.4;
// ctime_before=ctime_before -1.4;
//
//
// tr.ct_gb=-ctime;
//
// if(3600.0<ctime && ctime<3665){
//   ctime = ctime - 3637.88 - 12.76;
//   ctime = ctime - 12.0-3.1;
// }
////else ctime=-9999.;
// if(3600.0<ctime_before && ctime_before<3665){
//   ctime_before = ctime_before - 3637.88 - 12.76;
//   ctime_before = ctime_before - 12.0-3.1;
// }
//
// ctime=-ctime;
// ctime_before=-ctime_before;
////kazuki
//	//  if(bestcut){
//	//	hct_test3->Fill(ct_test);
//	//	hct_test2->Fill(ctime);
//	//	hct_test->Fill(ctime_before);
//	//	}
//	//	h_ctct->Fill(ctime,ctime_before);
//	//	h_gbetaR->Fill(beta_R);
//	//	h_gbetaL->Fill(beta_L);
//	//	h_gLenR->Fill(LenR);
//	//	h_gLenL->Fill(LenL);
//	//	h_gcorR->Fill(cor_R);
//	//	h_gcorL->Fill(cor_L);
//	//	h_gcorLR->Fill(cor_L,cor_R);
//	//	//time reference?
//	//	h_gtref_R->Fill(tref_R);
//	//	h_gtimeR_R->Fill(timeR_R);
//	//	h_gtimeL_R->Fill(timeL_R);
//	//	h_gmeantime->Fill(meantime_R);
//	//	h_gctcorR->Fill(ctimecorR);
//	//	h_gctcorL->Fill(ctimecorL);
//
// return ctime;
//
//
//}
