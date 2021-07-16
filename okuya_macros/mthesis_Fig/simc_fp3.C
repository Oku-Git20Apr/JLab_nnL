//--  Data vs SIMC (FP)  --//
//
//K. Okuyama (Feb. 21, 2021)
//K. Okuyama (Jul. 10, 2021)
//	-I forgot what I wanted to do with simc_fp_??? series.
//	-I rebuilt "simc_fp3.C" and "simc_fp_momcut.C"
//
//This is taken over from simc_fp.C

double FP_cut1(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double flag = par[0];//L or R
	return 0.17*x0/100.+0.025;
}

double FP_cut2(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double flag = par[0];//L or R
	return 0.17*x0/100.-0.035;
}

double FP_cut3(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double flag = par[0];//L or R
	return 0.40*x0/100.+0.130;
}
double FP_cut4(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double flag = par[0];//L or R
	return 0.40*x0/100.-0.130;
}

void simc_fp3(){
  
/*-- Input file --*/
  TFile *file = new TFile("../h2all_2020Nov.root","read");
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS.root","read");// L:S0=3:1 (0.9M vs 0.3M) 2020/12/10
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_Alupdate.root","read");// L:S0=3:1 (0.9M vs 0.3M) 
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_update_momminus10.root","read");// L:S0=3:1 (0.9M vs 0.3M) 
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");

/*-- Mixed Event Analysis --*/
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_2020Nov.root","read");
  double nbunch = 6000.;//effetive bunches (6 bunches x 5 mixtures)

/*-- Output file --*/
  TTree *tree = (TTree*)file->Get("tree_out");

/*-- Global Settings --*/
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);



//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_n_L=250.; 
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_n_S=70.; 
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
 int bin_mm=(max_mm-min_mm)/0.001; //Counts/2 MeV
 bin_mm=(int)bin_mm;
 //const double fit_min_mm=-0.006;
 const double fit_min_mm=-0.01;
 const double fit_max_mm=0.095;
 const int fit_bin_mm = (fit_max_mm-fit_min_mm)/0.001;
 const double fit_bin_width = (fit_max_mm-fit_min_mm)/fit_bin_mm;


//---------------------------------------//
//               Branch                  //
//---------------------------------------//

 int nrun, NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;

	double L_tr_chi2;
	double L_tr_x, L_tr_y, L_tr_th, L_tr_ph;
	double L_tr_p;
	double L_tr_tg_th, L_tr_tg_ph;
	double L_tr_vz, L_tr_vz2, L_tr_vz3;
	double L_tr_vz_saved;
	double R_tr_chi2;
	double R_tr_x, R_tr_y, R_tr_th, R_tr_ph;
	double R_tr_p;
	double R_tr_tg_th, R_tr_tg_ph;
	double R_tr_vz, R_tr_vz2, R_tr_vz3;
	double L_mom, R_mom, B_mom; 
	double L_ene, R_ene, B_ene; 
	double ac1sum, ac2sum;//NPE SUM

	//chain_dummy->SetBranchStatus("*",0);
  	//chain_dummy->SetBranchStatus("L.tr.vz",1);  chain_dummy->SetBranchAddress("L.tr.vz", &L_tr_vz2);
  	//chain_dummy->SetBranchStatus("R.tr.vz",1);  chain_dummy->SetBranchAddress("R.tr.vz", &R_tr_vz2);
	//chain_true->SetBranchStatus("*",0);
  	//chain_true->SetBranchStatus("L.tr.vz",1);  chain_true->SetBranchAddress("L.tr.vz", &L_tr_vz3);
  	//chain_true->SetBranchStatus("R.tr.vz",1);  chain_true->SetBranchAddress("R.tr.vz", &R_tr_vz3);

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("nrun",1);tree->SetBranchAddress("nrun",&nrun);
	tree->SetBranchStatus("tr.ntrack_l",1);tree->SetBranchAddress("tr.ntrack_l",&NLtr);
	tree->SetBranchStatus("tr.ntrack_r",1);tree->SetBranchAddress("tr.ntrack_r",&NRtr);
	
  	tree->SetBranchStatus("ac1_npe_sum",1);  tree->SetBranchAddress("ac1_npe_sum", &ac1sum);
  	tree->SetBranchStatus("ac2_npe_sum",1);  tree->SetBranchAddress("ac2_npe_sum", &ac2sum);
  	tree->SetBranchStatus("Lp_c",1);  tree->SetBranchAddress("Lp_c", &L_mom);
  	tree->SetBranchStatus("Rp_c",1);  tree->SetBranchAddress("Rp_c", &R_mom);
  	tree->SetBranchStatus("Bp_c",1);  tree->SetBranchAddress("Bp_c", &B_mom);
  	tree->SetBranchStatus("ct_orig",1);  tree->SetBranchAddress("ct_orig", &ct);

  	tree->SetBranchStatus("L.tr.chi2",1);  tree->SetBranchAddress("L.tr.chi2", &L_tr_chi2);
  	tree->SetBranchStatus("L.tr.x",1);  tree->SetBranchAddress("L.tr.x", &L_tr_x);
  	tree->SetBranchStatus("L.tr.y",1);  tree->SetBranchAddress("L.tr.y", &L_tr_y);
  	tree->SetBranchStatus("L.tr.th",1);  tree->SetBranchAddress("L.tr.th", &L_tr_th);
  	tree->SetBranchStatus("L.tr.ph",1);  tree->SetBranchAddress("L.tr.ph", &L_tr_ph);
  	tree->SetBranchStatus("L.tr.p",1);  tree->SetBranchAddress("L.tr.p", &L_tr_p);
  	tree->SetBranchStatus("L.tr.tg_th",1);  tree->SetBranchAddress("L.tr.tg_th", &L_tr_tg_th );
  	tree->SetBranchStatus("L.tr.tg_ph",1);  tree->SetBranchAddress("L.tr.tg_ph", &L_tr_tg_ph );
  	tree->SetBranchStatus("L.tr.vz",1);  tree->SetBranchAddress("L.tr.vz", &L_tr_vz);

  	tree->SetBranchStatus("R.tr.chi2",1);  tree->SetBranchAddress("R.tr.chi2", &R_tr_chi2);
	tree->SetBranchStatus("R.tr.x" ,1);  tree->SetBranchAddress("R.tr.x" , &R_tr_x );
  	tree->SetBranchStatus("R.tr.y" ,1);  tree->SetBranchAddress("R.tr.y" , &R_tr_y );
  	tree->SetBranchStatus("R.tr.th",1);  tree->SetBranchAddress("R.tr.th", &R_tr_th);
  	tree->SetBranchStatus("R.tr.ph",1);  tree->SetBranchAddress("R.tr.ph", &R_tr_ph);
  	tree->SetBranchStatus("R.tr.p",1);  tree->SetBranchAddress("R.tr.p", &R_tr_p);
  	tree->SetBranchStatus("R.tr.tg_th",1);  tree->SetBranchAddress("R.tr.tg_th", &R_tr_tg_th);
  	tree->SetBranchStatus("R.tr.tg_ph",1);  tree->SetBranchAddress("R.tr.tg_ph", &R_tr_tg_ph);
  	tree->SetBranchStatus("R.tr.vz",1);  tree->SetBranchAddress("R.tr.vz", &R_tr_vz);




  TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  h1->GetXaxis()->SetTitle("coin time (ns)");
  h1->GetYaxis()->SetTitle("Counts / 100 ps");
  h1->GetXaxis()->SetRangeUser(-14.0,17.);
  double xmin = -0.1, xmax = 0.2; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_L_fom_best  = new TH1F("hmm_L_fom_best","hmm_L_fom_best",xbin,xmin,xmax);
  hmm_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_L_fom_best->GetYaxis()->SetTitle("Counts / MeV");
  hmm_L_fom_best->SetLineColor(1);
  TH1F* hcs_L_fom_best  = new TH1F("hcs_L_fom_best","hcs_L_fom_best",xbin,xmin,xmax);
  hcs_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} [GeV/c^{2}]");
  hcs_L_fom_best->GetYaxis()->SetTitle("d#sigma/d#Omega (C.M.F.) [nb/sr]");
  hcs_L_fom_best->SetLineColor(1);
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin/2,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");
  TH1F* hmm_L_strict  = new TH1F("hmm_L_strict","",xbin,xmin*1000.,xmax*1000.);
  hmm_L_strict->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  hmm_L_strict->GetYaxis()->SetTitle("Counts/(MeV/c^{2})");
  hmm_L_strict->SetLineColor(kBlack);
  TH1F* hmm_wobg_strict  = new TH1F("hmm_wobg_strict","",xbin,xmin*1000.,xmax*1000.);
  hmm_wobg_strict->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  hmm_wobg_strict->GetYaxis()->SetTitle("Counts/(MeV/c^{2})");
  hmm_wobg_strict->SetLineColor(kBlack);
  TH1F* hmm_wobg_fom_best  = new TH1F("hmm_wobg_fom_best","",xbin,xmin,xmax);
  hmm_wobg_fom_best->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [GeV/c^{2}]");
  hmm_wobg_fom_best->GetYaxis()->SetTitle("Counts/(GeV/c^{2})");
  hmm_wobg_fom_best->SetLineColor(kBlack);
  TH1F* h_ct  = new TH1F("h_ct","h_ct",1000/0.056,-5000,5000.);
  TH1F* h_ct2  = new TH1F("h_ct2","h_ct2",1000/0.056,-5000,5000.);
  TH1F* hmm_ctout_only  = new TH1F("hmm_ctout_only","hmm_ctout_only",300,-0.1,0.2);
  TH1F* hmm_ctout  = new TH1F("hmm_ctout","hmm_ctout",300,-0.1,0.2);
  TH1F* hmm_ctout2  = new TH1F("hmm_ctout2","hmm_ctout2",300,-0.1,0.2);
  TH1F* hmm_ctout3  = new TH1F("hmm_ctout3","hmm_ctout3",300,-0.1,0.2);
  TH1F* hmm_ctout4  = new TH1F("hmm_ctout4","hmm_ctout4",300,-0.1,0.2);
  TH2F* h_ctmm = new TH2F("h_ctmm","Cointime vs MM",100/0.056,-20.,20.,100,-0.1,0.2);
  TH2F* h_ctmm2 = new TH2F("h_ctmm2","Cointime vs MM (strict)",100/0.056,-20.,20.,100,-0.1,0.2);
  TH2F* gklab_gkcm = new TH2F("gklab_gkcm","#theta_{gk}^{lab} vs #theta_{gk}^{CM}",100,0.,0.15,100,0.,0.4);
  TH2F* gklab_eklab = new TH2F("gklab_eklab","#theta_{gk}^{lab} vs #theta_{ek}^{lab}",100,0.,0.15,100,0.15,0.3);
  TH2F* eelab_eklab = new TH2F("eelab_eklab","#theta_{ee}^{lab} vs #theta_{ek}^{lab}",100,0.15,0.3,100,0.15,0.3);
  TH2F* eklab_gkcm = new TH2F("eklab_gkcm","#theta_{gk}^{CM} vs #theta_{ek}^{lab}",100,0.,0.4,100,0.15,0.3);
  TH2F* ekcm_gkcm = new TH2F("ekcm_gkcm","#theta_{gk}^{CM} vs #theta_{ek}^{CM}",100,0.,0.25,100,0.15,0.4);
  TH2F* cos_gklab_gkcm = new TH2F("cos_gklab_gkcm","#theta_{gk}^{lab} vs #theta_{gk}^{CM}",100,0.98,1.0,100,0.92,1.0);
  TH2F* cos_gklab_eklab = new TH2F("cos_gklab_eklab","#theta_{gk}^{lab} vs #theta_{ek}^{lab}",100,0.98,1.0,100,0.95,1.0);
  TH2F* cos_eelab_eklab = new TH2F("cos_eelab_eklab","#theta_{ee}^{lab} vs #theta_{ek}^{lab}",100,0.95,1.0,100,0.95,1.0);
  TH2F* cos_eklab_gkcm = new TH2F("cos_eklab_gkcm","#theta_{gk}^{CM} vs #theta_{ek}^{lab}",100,0.95,1.0,100,0.95,1.0);
  TH2F* cos_ekcm_gkcm = new TH2F("cos_ekcm_gkcm","#theta_{gk}^{CM} vs #theta_{ek}^{CM}",100,0.9,1.0,100,0.9,1.0);

  TH1D* h_theta_ee = new TH1D("h_theta_ee", "theta_ee",1000,0.1,0.35);
  TH1D* h_phi_ee = new TH1D("h_phi_ee", "phi_ee",1000,0.,PI);
  TH1D* h_theta_ek = new TH1D("h_theta_ek", "theta_ek",1000,0.1,0.35);
  TH1D* h_phi_ek = new TH1D("h_phi_ek", "phi_ek",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_g = new TH1D("h_theta_g", "theta_g",1000,0.1,0.35);
  TH1D* h_phi_g = new TH1D("h_phi_g", "phi_g",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_gk_lab = new TH1D("h_theta_gk_lab", "theta_gk_lab",1000,0.,0.2);
  TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "theta_gk_cm",1000,0.,0.3);
  TH1D* h_cos_gk_lab = new TH1D("h_cos_gk_lab", "cos_gk_lab",1000,0.97,1.0);
  TH1D* h_cos_gk_cm = new TH1D("h_cos_gk_cm", "cos_gk_cm",1000,0.8,1.0);
  TH1D* h_mom_g = new TH1D("h_mom_g", "mom_g",1000,1.8,2.5);
  TH1D* h_qsq = new TH1D("h_qsq", "Q^2",1000,0.,0.8);
  TH1D* h_w = new TH1D("h_w", "W",1000,0.,0.8);
  TH2D* h_thph_ee = new TH2D("h_thph_ee", "theta_ee:phi_ee" ,1000,0.1,0.35,1000,PI/2-1.,PI/2+1.);
  TH2D* h_thph_ek = new TH2D("h_thph_ek", "theta_ek:phi_ek" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
  TH2D* h_thph_g = new TH2D("h_thph_g", "theta_g:phi_g" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
  
  TH1F* h_nltrack  = new TH1F("h_nltrack","NLtr",10,-2,8);
  TH1F* h_nrtrack  = new TH1F("h_nrtrack","NRtr",10,-2,8);
  //TH2F* h_pepk  = new TH2F("h_pepk","h_pepk (tight)",40,1.73,1.93,40,1.95,2.25);
  TH2D* h_pepk = new TH2D("h_pepk", "p_{K} vs p_{e'}" ,50,1720.,1940.,50,1980.,2220.);
  h_pepk->GetXaxis()->SetTitle("p_{K} [MeV/c]");
  h_pepk->GetYaxis()->SetTitle("p_{e'} [MeV/c]");
  h_pepk->SetNdivisions(505);
  //h_pepk->GetZaxis()->SetLabelOffset(-0.005);
  TH1D* h_pe = new TH1D("h_pe", "p_{e'}" ,200,1980.,2220.);
  h_pe->GetXaxis()->SetTitle("Momentum [MeV/c]");
  h_pe->GetYaxis()->SetTitle("Counts");
  TH1D* h_pk = new TH1D("h_pk", "p_{K}" ,200,1720.,1940.);
  h_pk->GetXaxis()->SetTitle("Momentum [MeV/c]");
  h_pk->GetYaxis()->SetTitle("Counts");
  TH2D* h_pepk_simc = new TH2D("h_pepk_simc", "p_{K} vs p_{e'}" ,50,1720.,1940.,50,1980.,2220.);
  h_pepk_simc->GetXaxis()->SetTitle("p_{K} [MeV/c]");
  h_pepk_simc->GetYaxis()->SetTitle("p_{e'} [MeV/c]");
  h_pepk_simc->SetNdivisions(505);
  //h_pepk_simc->GetZaxis()->SetLabelOffset(-0.005);
  tree_simc->Project("h_pepk_simc","Lp_rec:Rp_rec","");
  TH1D* h_pe_simc = new TH1D("h_pe_simc", "p_{e'}" ,200,1980.,2220.);
  h_pe_simc->GetXaxis()->SetTitle("Momentum [MeV/c]");
  h_pe_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_pe_simc","Lp_rec","");
  TH1D* h_pk_simc = new TH1D("h_pk_simc", "p_{K}" ,200,1720.,1940.);
  h_pk_simc->GetXaxis()->SetTitle("Momentum [MeV/c]");
  h_pk_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_pk_simc","Rp_rec","");
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",100,-0.15,0.15,100,-0.15,0.15);
 
  
//Focal Plane (Feb. 21, 2021)
  TH1D* h_R_y_data       = new TH1D("h_R_y_data"      ,"h_R_y_data"      ,80,   -6.,  6.);
  h_R_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  h_R_y_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_y_simc       = new TH1D("h_R_y_simc"      ,"h_R_y_simc"      ,80,   -6.,  6.);
  h_R_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  h_R_y_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_y_simc","h_yfp+0.516","");
  TH1D* h_R_x_data       = new TH1D("h_R_x_data"      ,"h_R_x_data"      ,80,   -80.,  80.);
  h_R_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  h_R_x_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_x_simc       = new TH1D("h_R_x_simc"      ,"h_R_x_simc"      ,80,   -80.,  80.);
  h_R_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  h_R_x_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_x_simc","h_xfp","");
  TH1D* h_R_th_data       = new TH1D("h_R_th_data"      ,"h_R_th_data"      ,80,   -0.2,  0.2);
  h_R_th_data->SetNdivisions(505);
  h_R_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  h_R_th_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_th_simc       = new TH1D("h_R_th_simc"      ,"h_R_th_simc"      ,80,   -0.2,  0.2);
  h_R_th_simc->SetNdivisions(505);
  h_R_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  h_R_th_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_th_simc","h_xpfp","");
  TH1D* h_R_ph_data       = new TH1D("h_R_ph_data"      ,"h_R_ph_data"      ,80,   -0.05,  0.05);
  h_R_ph_data->SetNdivisions(505);
  h_R_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_ph_simc       = new TH1D("h_R_ph_simc"      ,"h_R_ph_simc"      ,80,   -0.05,  0.05);
  h_R_ph_simc->SetNdivisions(505);
  h_R_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_ph_simc","h_ypfp","");
  TH1D* h_R_tg_th_data       = new TH1D("h_R_tg_th_data"      ,"h_R_tg_th_data"      ,80,   -0.1,  0.1);
  h_R_tg_th_data->SetNdivisions(505);
  h_R_tg_th_data->GetXaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_tg_th_simc       = new TH1D("h_R_tg_th_simc"      ,"h_R_tg_th_simc"      ,80,   -0.1,  0.1);
  h_R_tg_th_simc->SetNdivisions(505);
  h_R_tg_th_simc->GetXaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_tg_th_simc","h_xptar","");
  TH1D* h_R_tg_ph_data       = new TH1D("h_R_tg_ph_data"      ,"h_R_tg_ph_data"      ,80,   -0.05,  0.05);
  h_R_tg_ph_data->SetNdivisions(505);
  h_R_tg_ph_data->GetXaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_R_tg_ph_simc       = new TH1D("h_R_tg_ph_simc"      ,"h_R_tg_ph_simc"      ,80,   -0.05,  0.05);
  h_R_tg_ph_simc->SetNdivisions(505);
  h_R_tg_ph_simc->GetXaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_R_tg_ph_simc","h_yptar","");

  TH1D* h_L_y_data       = new TH1D("h_L_y_data"      ,"h_L_y_data"      ,80,   -6.,  6.);
  h_L_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  h_L_y_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_y_simc       = new TH1D("h_L_y_simc"      ,"h_L_y_simc"      ,80,   -6.,  6.);
  h_L_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  h_L_y_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_y_simc","e_yfp+0.807","");
  TH1D* h_L_x_data       = new TH1D("h_L_x_data"      ,"h_L_x_data"      ,80,   -80.,  80.);
  h_L_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  h_L_x_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_x_simc       = new TH1D("h_L_x_simc"      ,"h_L_x_simc"      ,80,   -80.,  80.);
  h_L_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  h_L_x_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_x_simc","e_xfp","");
  TH1D* h_L_th_data       = new TH1D("h_L_th_data"      ,"h_L_th_data"      ,80,   -0.2,  0.2);
  h_L_th_data->SetNdivisions(505);
  h_L_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  h_L_th_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_th_simc       = new TH1D("h_L_th_simc"      ,"h_L_th_simc"      ,80,   -0.2,  0.2);
  h_L_th_simc->SetNdivisions(505);
  h_L_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  h_L_th_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_th_simc","e_xpfp","");
  TH1D* h_L_ph_data       = new TH1D("h_L_ph_data"      ,"h_L_ph_data"      ,80,   -0.05,  0.05);
  h_L_ph_data->SetNdivisions(505);
  h_L_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_ph_simc       = new TH1D("h_L_ph_simc"      ,"h_L_ph_simc"      ,80,   -0.05,  0.05);
  h_L_ph_simc->SetNdivisions(505);
  h_L_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_ph_simc","e_ypfp","");
  TH1D* h_L_tg_th_data       = new TH1D("h_L_tg_th_data"      ,"h_L_tg_th_data"      ,80,   -0.1,  0.1);
  h_L_tg_th_data->SetNdivisions(505);
  h_L_tg_th_data->GetXaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_tg_th_simc       = new TH1D("h_L_tg_th_simc"      ,"h_L_tg_th_simc"      ,80,   -0.1,  0.1);
  h_L_tg_th_simc->SetNdivisions(505);
  h_L_tg_th_simc->GetXaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_tg_th_simc","e_xptar","");
  TH1D* h_L_tg_ph_data       = new TH1D("h_L_tg_ph_data"      ,"h_L_tg_ph_data"      ,80,   -0.05,  0.05);
  h_L_tg_ph_data->SetNdivisions(505);
  h_L_tg_ph_data->GetXaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_data->GetYaxis()->SetTitle("Counts");
  TH1D* h_L_tg_ph_simc       = new TH1D("h_L_tg_ph_simc"      ,"h_L_tg_ph_simc"      ,80,   -0.05,  0.05);
  h_L_tg_ph_simc->SetNdivisions(505);
  h_L_tg_ph_simc->GetXaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_simc->GetYaxis()->SetTitle("Counts");
  tree_simc->Project("h_L_tg_ph_simc","e_yptar","");

//== FP x FP ==//
//X-Y, X-theta, Y-phi, theta-phi 
//
//X-Y
  TH2D* h_L_x_y_data       = new TH2D("h_L_x_y_data"      ,"h_L_x_y_data"   ,80,  -6., 6.   ,80,   -80.,  80.);
  h_L_x_y_data->GetYaxis()->SetTitle("X(FP) [cm]");
  h_L_x_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_L_x_y_simc       = new TH2D("h_L_x_y_simc"      ,"h_L_x_y_simc"   ,80,  -6., 6.   ,80,   -80.,  80.);
  h_L_x_y_simc->GetYaxis()->SetTitle("X(FP) [cm]");
  h_L_x_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_L_x_y_simc","e_xfp:e_yfp+0.807","");
  TH2D* h_R_x_y_data       = new TH2D("h_R_x_y_data"      ,"h_R_x_y_data"   ,80,  -6., 6.   ,80,   -80.,  80.);
  h_R_x_y_data->GetYaxis()->SetTitle("X(FP) [cm]");
  h_R_x_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_R_x_y_simc       = new TH2D("h_R_x_y_simc"      ,"h_R_x_y_simc"   ,80,  -6., 6.   ,80,   -80.,  80.);
  h_R_x_y_simc->GetYaxis()->SetTitle("X(FP) [cm]");
  h_R_x_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_R_x_y_simc","h_xfp:h_yfp+0.516","");
//X-theta
  TH2D* h_L_th_x_data       = new TH2D("h_L_th_x_data"      ,"h_L_th_x_data"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_L_th_x_data->GetYaxis()->SetTitle("#theta(FP) [rad]");
  h_L_th_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  TH2D* h_L_th_x_simc       = new TH2D("h_L_th_x_simc"      ,"h_L_th_x_simc"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_L_th_x_simc->GetYaxis()->SetTitle("#theta(FP) [rad]");
  h_L_th_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  tree_simc->Project("h_L_th_x_simc","e_xpfp:e_xfp","");
  //tree_simc->Project("h_L_th_x_simc","e_xpfp:e_xfp","e_xpfp<0.17*e_xfp/100.+0.025&&e_xpfp>0.17*e_xfp/100.-0.035&&e_xpfp<0.40*e_xfp/100.+0.130");
  TH2D* h_R_th_x_data       = new TH2D("h_R_th_x_data"      ,"h_R_th_x_data"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_R_th_x_data->GetYaxis()->SetTitle("#theta(FP) [rad]");
  h_R_th_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  TH2D* h_R_th_x_simc       = new TH2D("h_R_th_x_simc"      ,"h_R_th_x_simc"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_R_th_x_simc->GetYaxis()->SetTitle("#theta(FP) [rad]");
  h_R_th_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  tree_simc->Project("h_R_th_x_simc","h_xpfp:h_xfp","");
  //tree_simc->Project("h_R_th_x_simc","h_xpfp:h_xfp","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130");
//Y-phi
  TH2D* h_L_ph_y_data       = new TH2D("h_L_ph_y_data"      ,"h_L_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_L_ph_y_data->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_L_ph_y_simc       = new TH2D("h_L_ph_y_simc"      ,"h_L_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_L_ph_y_simc->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_L_ph_y_simc","e_ypfp:e_yfp+0.807","");
  TH2D* h_R_ph_y_data       = new TH2D("h_R_ph_y_data"      ,"h_R_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_R_ph_y_data->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_R_ph_y_simc       = new TH2D("h_R_ph_y_simc"      ,"h_R_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_R_ph_y_simc->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_R_ph_y_simc","h_ypfp:h_yfp+0.516","");
//theta-phi
  TH2D* h_L_ph_th_data       = new TH2D("h_L_ph_th_data"      ,"h_L_ph_th_data"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_L_ph_th_data->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  TH2D* h_L_ph_th_simc       = new TH2D("h_L_ph_th_simc"      ,"h_L_ph_th_simc"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_L_ph_th_simc->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_L_ph_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  tree_simc->Project("h_L_ph_th_simc","e_ypfp:e_xpfp","");
  TH2D* h_R_ph_th_data       = new TH2D("h_R_ph_th_data"      ,"h_R_ph_th_data"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_R_ph_th_data->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  TH2D* h_R_ph_th_simc       = new TH2D("h_R_ph_th_simc"      ,"h_R_ph_th_simc"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_R_ph_th_simc->GetYaxis()->SetTitle("#phi(FP) [rad]");
  h_R_ph_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  tree_simc->Project("h_R_ph_th_simc","h_ypfp:h_xpfp","");

//== tar x tar ==//
//theta-phi, phi-Z
//
//theta-phi
  TH2D* h_L_tg_ph_tg_th_data       = new TH2D("h_L_tg_ph_tg_th_data"      ,"h_L_tg_ph_tg_th_data"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_L_tg_ph_tg_th_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_tg_th_data->GetXaxis()->SetTitle("#theta(tar) [rad]");
  TH2D* h_L_tg_ph_tg_th_simc       = new TH2D("h_L_tg_ph_tg_th_simc"      ,"h_L_tg_ph_tg_th_simc"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_L_tg_ph_tg_th_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_tg_th_simc->GetXaxis()->SetTitle("#theta(tar) [rad]");
  tree_simc->Project("h_L_tg_ph_tg_th_simc","e_yptar:e_xptar","");
  TH2D* h_R_tg_ph_tg_th_data       = new TH2D("h_R_tg_ph_tg_th_data"      ,"h_R_tg_ph_tg_th_data"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_R_tg_ph_tg_th_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_tg_th_data->GetXaxis()->SetTitle("#theta(tar) [rad]");
  TH2D* h_R_tg_ph_tg_th_simc       = new TH2D("h_R_tg_ph_tg_th_simc"      ,"h_R_tg_ph_tg_th_simc"   ,80,    -0.12,  0.12   ,80,   -0.05,  0.05);
  h_R_tg_ph_tg_th_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_tg_th_simc->GetXaxis()->SetTitle("#theta(tar) [rad]");
  tree_simc->Project("h_R_tg_ph_tg_th_simc","h_yptar:h_xptar","");
//phi-Z
  TH2D* h_L_z_tg_ph_data       = new TH2D("h_L_z_tg_ph_data"      ,"h_L_z_tg_ph_data"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_L_z_tg_ph_data->SetNdivisions(505);
  h_L_z_tg_ph_data->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_tg_ph_data->GetXaxis()->SetTitle("#phi(tar) [rad]");
  TH2D* h_L_z_tg_ph_simc       = new TH2D("h_L_z_tg_ph_simc"      ,"h_L_z_tg_ph_simc"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_L_z_tg_ph_simc->SetNdivisions(505);
  h_L_z_tg_ph_simc->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_tg_ph_simc->GetXaxis()->SetTitle("#phi(tar) [rad]");
  tree_simc->Project("h_L_z_tg_ph_simc","zposi:e_yptar","");                                                              
  TH2D* h_R_z_tg_ph_data       = new TH2D("h_R_z_tg_ph_data"      ,"h_R_z_tg_ph_data"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_R_z_tg_ph_data->SetNdivisions(505);
  h_R_z_tg_ph_data->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_tg_ph_data->GetXaxis()->SetTitle("#phi(tar) [rad]");
  TH2D* h_R_z_tg_ph_simc       = new TH2D("h_R_z_tg_ph_simc"      ,"h_R_z_tg_ph_simc"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_R_z_tg_ph_simc->SetNdivisions(505);
  h_R_z_tg_ph_simc->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_tg_ph_simc->GetXaxis()->SetTitle("#phi(tar) [rad]");
  tree_simc->Project("h_R_z_tg_ph_simc","zposi:h_yptar","");

//== FP x tar ==//
//X-theta, Y-phi, theta-theta, phi-phi, Y-Z, phi-Z
//X-theta
  TH2D* h_L_tg_th_x_data       = new TH2D("h_L_tg_th_x_data"      ,"h_L_tg_th_x_data"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_L_tg_th_x_data->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  TH2D* h_L_tg_th_x_simc       = new TH2D("h_L_tg_th_x_simc"      ,"h_L_tg_th_x_simc"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_L_tg_th_x_simc->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  tree_simc->Project("h_L_tg_th_x_simc","e_xptar:e_xfp","");
  //tree_simc->Project("h_L_tg_th_x_simc","e_xpfp:e_xfp","e_xpfp<0.17*e_xfp/100.+0.025&&e_xpfp>0.17*e_xfp/100.-0.035&&e_xpfp<0.40*e_xfp/100.+0.130");
  TH2D* h_R_tg_th_x_data       = new TH2D("h_R_tg_th_x_data"      ,"h_R_tg_th_x_data"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_R_tg_th_x_data->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_x_data->GetXaxis()->SetTitle("X(FP) [cm]");
  TH2D* h_R_tg_th_x_simc       = new TH2D("h_R_tg_th_x_simc"      ,"h_R_tg_th_x_simc"   ,80,  -80., 80.   ,80,   -0.12,  0.12);
  h_R_tg_th_x_simc->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_x_simc->GetXaxis()->SetTitle("X(FP) [cm]");
  tree_simc->Project("h_R_tg_th_x_simc","h_xptar:h_xfp","");
  //tree_simc->Project("h_R_th_x_simc","h_xpfp:h_xfp","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130");
//Y-phi
  TH2D* h_L_tg_ph_y_data       = new TH2D("h_L_tg_ph_y_data"      ,"h_L_tg_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_L_tg_ph_y_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_L_tg_ph_y_simc       = new TH2D("h_L_tg_ph_y_simc"      ,"h_L_tg_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_L_tg_ph_y_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_L_tg_ph_y_simc","e_ypfp:e_yfp+0.807","");
  TH2D* h_R_tg_ph_y_data       = new TH2D("h_R_tg_ph_y_data"      ,"h_R_tg_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_R_tg_ph_y_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_R_tg_ph_y_simc       = new TH2D("h_R_tg_ph_y_simc"      ,"h_R_tg_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  h_R_tg_ph_y_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_R_tg_ph_y_simc","h_ypfp:h_yfp+0.516","");
//theta-theta
  TH2D* h_L_tg_th_th_data       = new TH2D("h_L_tg_th_th_data"      ,"h_L_tg_th_th_data"   ,80,  -0.12,  0.12  ,80,   -0.12,  0.12);
  h_L_tg_th_th_data->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  TH2D* h_L_tg_th_th_simc       = new TH2D("h_L_tg_th_th_simc"      ,"h_L_tg_th_th_simc"   ,80,  -0.12,  0.12  ,80,   -0.12,  0.12);
  h_L_tg_th_th_simc->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_L_tg_th_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  tree_simc->Project("h_L_tg_th_th_simc","e_xptar:e_xpfp","");
  //tree_simc->Project("h_L_tg_th_th_simc","e_thpfp:e_thfp","e_thpfp<0.17*e_thfp/100.+0.025&&e_thpfp>0.17*e_thfp/100.-0.035&&e_thpfp<0.40*e_thfp/100.+0.130");
  TH2D* h_R_tg_th_th_data       = new TH2D("h_R_tg_th_th_data"      ,"h_R_tg_th_th_data"   ,80,  -0.12,  0.12  ,80,   -0.12,  0.12);
  h_R_tg_th_th_data->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_th_data->GetXaxis()->SetTitle("#theta(FP) [rad]");
  TH2D* h_R_tg_th_th_simc       = new TH2D("h_R_tg_th_th_simc"      ,"h_R_tg_th_th_simc"   ,80,  -0.12,  0.12  ,80,   -0.12,  0.12);
  tree_simc->Project("h_R_tg_th_th_simc","h_xptar:h_xpfp","");
  h_R_tg_th_th_simc->GetYaxis()->SetTitle("#theta(tar) [rad]");
  h_R_tg_th_th_simc->GetXaxis()->SetTitle("#theta(FP) [rad]");
  //tree_simc->Project("h_R_th_x_simc","h_xpfp:h_xfp","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130");
//phi-phi
  TH2D* h_L_tg_ph_ph_data       = new TH2D("h_L_tg_ph_ph_data"      ,"h_L_tg_ph_ph_data"   ,80,   -0.05,  0.05  ,80,   -0.05,  0.05);
  h_L_tg_ph_ph_data->SetNdivisions(505);
  h_L_tg_ph_ph_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  TH2D* h_L_tg_ph_ph_simc       = new TH2D("h_L_tg_ph_ph_simc"      ,"h_L_tg_ph_ph_simc"   ,80,   -0.05,  0.05  ,80,   -0.05,  0.05);
  h_L_tg_ph_ph_simc->SetNdivisions(505);
  h_L_tg_ph_ph_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_L_tg_ph_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  tree_simc->Project("h_L_tg_ph_ph_simc","e_yptar:e_ypfp","");                                               
  TH2D* h_R_tg_ph_ph_data       = new TH2D("h_R_tg_ph_ph_data"      ,"h_R_tg_ph_ph_data"   ,80,   -0.05,  0.05  ,80,   -0.05,  0.05);
  h_R_tg_ph_ph_data->SetNdivisions(505);
  h_R_tg_ph_ph_data->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  TH2D* h_R_tg_ph_ph_simc       = new TH2D("h_R_tg_ph_ph_simc"      ,"h_R_tg_ph_ph_simc"   ,80,   -0.05,  0.05  ,80,   -0.05,  0.05);
  h_R_tg_ph_ph_simc->SetNdivisions(505);
  h_R_tg_ph_ph_simc->GetYaxis()->SetTitle("#phi(tar) [rad]");
  h_R_tg_ph_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  tree_simc->Project("h_R_tg_ph_ph_simc","h_yptar:h_ypfp","");
//Y-Z
  TH2D* h_L_z_y_data       = new TH2D("h_L_z_y_data"      ,"h_L_z_y_data"   ,80,    -6., 6.   ,80,   -15., 15.  );
  h_L_z_y_data->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_L_z_y_simc       = new TH2D("h_L_z_y_simc"      ,"h_L_z_y_simc"   ,80,    -6., 6.   ,80,   -15., 15.  );
  h_L_z_y_simc->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_L_z_y_simc","zposi:e_yfp+0.807","");                                                          
  TH2D* h_R_z_y_data       = new TH2D("h_R_z_y_data"      ,"h_R_z_y_data"   ,80,    -6., 6.   ,80,   -15., 15.  );
  h_R_z_y_data->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_y_data->GetXaxis()->SetTitle("Y(FP) [cm]");
  TH2D* h_R_z_y_simc       = new TH2D("h_R_z_y_simc"      ,"h_R_z_y_simc"   ,80,    -6., 6.   ,80,   -15., 15.  );
  h_R_z_y_simc->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_y_simc->GetXaxis()->SetTitle("Y(FP) [cm]");
  tree_simc->Project("h_R_z_y_simc","zposi:h_yfp+0.516","");
//phi-Z
  TH2D* h_L_z_ph_data       = new TH2D("h_L_z_ph_data"      ,"h_L_z_ph_data"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_L_z_ph_data->SetNdivisions(505);
  h_L_z_ph_data->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  TH2D* h_L_z_ph_simc       = new TH2D("h_L_z_ph_simc"      ,"h_L_z_ph_simc"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_L_z_ph_simc->SetNdivisions(505);
  h_L_z_ph_simc->GetYaxis()->SetTitle("Z [cm]");
  h_L_z_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  tree_simc->Project("h_L_z_ph_simc","zposi:e_ypfp","");                                                              
  TH2D* h_R_z_ph_data       = new TH2D("h_R_z_ph_data"      ,"h_R_z_ph_data"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_R_z_ph_data->SetNdivisions(505);
  h_R_z_ph_data->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_ph_data->GetXaxis()->SetTitle("#phi(FP) [rad]");
  TH2D* h_R_z_ph_simc       = new TH2D("h_R_z_ph_simc"      ,"h_R_z_ph_simc"   ,80,   -0.05,  0.05   ,80,   -15., 15.  );
  h_R_z_ph_simc->SetNdivisions(505);
  h_R_z_ph_simc->GetYaxis()->SetTitle("Z [cm]");
  h_R_z_ph_simc->GetXaxis()->SetTitle("#phi(FP) [rad]");
  tree_simc->Project("h_R_z_ph_simc","zposi:h_ypfp","");


  TH2D* h_R_y_vz_data       = new TH2D("h_R_y_vz_data"      ,"h_R_y_vz_data"   ,80,   -6.,  6.,80,  -15., 15.   );
  TH2D* h_R_y_vz_simc       = new TH2D("h_R_y_vz_simc"      ,"h_R_y_vz_simc"   ,80,   -6.,  6.,80,  -15., 15.   );
  tree_simc->Project("h_R_y_vz_simc","zposi:h_yfp+0.516","");                                                         
  TH2D* h_L_y_vz_data       = new TH2D("h_L_y_vz_data"      ,"h_L_y_vz_data"   ,80,   -6.,  6.,80,  -15., 15.   );
  TH2D* h_L_y_vz_simc       = new TH2D("h_L_y_vz_simc"      ,"h_L_y_vz_simc"   ,80,   -6.,  6.,80,  -15., 15.   );
  tree_simc->Project("h_L_y_vz_simc","zposi:e_yfp+0.807","");
  TH2D* h_R_ph_vz_data       = new TH2D("h_R_ph_vz_data"      ,"h_R_ph_vz_data"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  TH2D* h_R_ph_vz_simc       = new TH2D("h_R_ph_vz_simc"      ,"h_R_ph_vz_simc"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  tree_simc->Project("h_R_ph_vz_simc","zposi:h_ypfp","");                                                         
  TH2D* h_L_ph_vz_data       = new TH2D("h_L_ph_vz_data"      ,"h_L_ph_vz_data"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  TH2D* h_L_ph_vz_simc       = new TH2D("h_L_ph_vz_simc"      ,"h_L_ph_vz_simc"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  tree_simc->Project("h_L_ph_vz_simc","zposi:e_ypfp","");
  TH2D* h_R_p_x_data       = new TH2D("h_R_p_x_data"      ,"h_R_p_x_data"   ,80,   1700.,  1950.,80,  -80., 80.  );
  h_R_p_x_data->GetYaxis()->SetTitle("X(FP) [cm]");
  h_R_p_x_data->GetXaxis()->SetTitle("Momentum [MeV/c]");
  TH2D* h_R_p_x_simc       = new TH2D("h_R_p_x_simc"      ,"h_R_p_x_simc"   ,80,   1700.,  1950.,80,  -80., 80.  );
  h_R_p_x_simc->GetYaxis()->SetTitle("X(FP) [cm]");
  h_R_p_x_simc->GetXaxis()->SetTitle("Momentum [MeV/c]");
  tree_simc->Project("h_R_p_x_simc","h_xfp:Rp_rec","");
  TH2D* h_L_p_x_data       = new TH2D("h_L_p_x_data"      ,"h_L_p_x_data"   ,80,   1950.,  2250.,80,  -80., 80.  );
  h_L_p_x_data->GetYaxis()->SetTitle("X(FP) [cm]");
  h_L_p_x_data->GetXaxis()->SetTitle("Momentum [MeV/c]");
  TH2D* h_L_p_x_simc       = new TH2D("h_L_p_x_simc"      ,"h_L_p_x_simc"   ,80,   1950.,  2250.,80,  -80., 80.  );
  h_L_p_x_simc->GetYaxis()->SetTitle("X(FP) [cm]");
  h_L_p_x_simc->GetXaxis()->SetTitle("Momentum [MeV/c]");
  tree_simc->Project("h_L_p_x_simc","e_xfp:Lp_rec","");

// Data => Line Color is Red
  h_R_y_data->SetLineColor(kRed);
  h_R_x_data->SetLineColor(kRed);
  h_R_th_data->SetLineColor(kRed);
  h_R_ph_data->SetLineColor(kRed);
  h_R_tg_th_data->SetLineColor(kRed);
  h_R_tg_ph_data->SetLineColor(kRed);
  h_L_y_data->SetLineColor(kRed);
  h_L_x_data->SetLineColor(kRed);
  h_L_th_data->SetLineColor(kRed);
  h_L_ph_data->SetLineColor(kRed);
  h_L_tg_th_data->SetLineColor(kRed);
  h_L_tg_ph_data->SetLineColor(kRed);

// SIMC => Line Color is Azure 
  h_R_y_simc->SetLineColor(kAzure);
  h_R_x_simc->SetLineColor(kAzure);
  h_R_th_simc->SetLineColor(kAzure);
  h_R_ph_simc->SetLineColor(kAzure);
  h_R_tg_th_simc->SetLineColor(kAzure);
  h_R_tg_ph_simc->SetLineColor(kAzure);
  h_L_y_simc->SetLineColor(kAzure);
  h_L_x_simc->SetLineColor(kAzure);
  h_L_th_simc->SetLineColor(kAzure);
  h_L_ph_simc->SetLineColor(kAzure);
  h_L_tg_th_simc->SetLineColor(kAzure);
  h_L_tg_ph_simc->SetLineColor(kAzure);


  TH1D* h_pe_simc_FPcut = new TH1D("h_pe_simc_FPcut", "p_{e'} (FP cut)" ,200,1980.,2220.);
  tree_simc->Project("h_pe_simc_FPcut","Lp_rec","e_xpfp<0.17*e_xfp/100.+0.025&&e_xpfp>0.17*e_xfp/100.-0.035&&e_xpfp<0.40*e_xfp/100.+0.130&&e_xpfp>0.40*e_xfp/100.-0.130");
  TH1D* h_pk_simc_FPcut = new TH1D("h_pk_simc_FPcut", "p_{K} (FP cut)" ,200,1720.,1940.);
  tree_simc->Project("h_pk_simc_FPcut","Rp_rec","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130&&h_xpfp>0.40*h_xfp/100.-0.130");

  TH1F* h_zave  = new TH1F("h_zave","Z-vertex (Ave.)",1000,-0.25,0.25);
  TH1F* h_zave_dummy  = new TH1F("h_zave_dummy","Z-vertex (Ave.)",1000,-0.15,0.15);
  TH1F* h_zave_true  = new TH1F("h_zave_true","Z-vertex (Ave.)",1000,-0.15,0.15);
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool ct_cut_ctout = false;
  bool event_selection = false;
  bool event_selection_ctout = false;
  bool event_selection_ctout2 = false;
  bool event_selection_ctout3 = false;
  bool event_selection_ctout4 = false;
  bool event_selection_nocut = false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum

/***********************************/
/**	What parameter do you change? **/
/***********************************/
  bool zcut_flag  = true;
  bool accut_flag = false;
  bool fpcut_flag = false;
  bool trcut_flag = false;
/***********************************/
/***********************************/
/***********************************/
  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  //tree->Draw(">>elist" , "fabs(ct_orig)<3.");//ctsum (does NOT dintinguish #track)
  //TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  //int ENum = elist->GetN(); 
  int ENum = tree->GetEntries(); 
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(int i=0;i<ENum;i++){
	//tree->GetEntry(elist->GetEntry(i));
	tree->GetEntry(i);

    if(i%time_div==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (i+1) - diff;
      cout<<i<<" / "<<ENum<<" ("<<i*100/ENum<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }

      
        L_Tr = L_FP = false;
        if( L_tr_chi2<0.01 ) L_Tr = true;
        if( L_tr_th<0.17*L_tr_x+0.025
         && L_tr_th>0.17*L_tr_x-0.035
         && L_tr_th<0.40*L_tr_x+0.130 ) L_FP = true;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2<0.01 ) R_Tr = true;
        if( R_tr_th<0.17*R_tr_x+0.025
         && R_tr_th>0.17*R_tr_x-0.035
         && R_tr_th<0.40*R_tr_x+0.130 ) R_FP = true;


	

		double pi_pos = 3.18;//ns

		if(fabs(ct)<1.006)ct_cut=true;
		else ct_cut=false;
		//if(fabs(ct+pi_pos+147.455)<1.006||fabs(ct+pi_pos+185.081)<1.006||fabs(ct+pi_pos+348.277)<1.006||fabs(ct+pi_pos+687.262)<1.006||fabs(ct+pi_pos+1151.74)<1.006||fabs(ct+pi_pos+2181.05)<1.006||fabs(ct+pi_pos+2017.86)<1.006||fabs(ct+pi_pos+2520.02)<1.006||fabs(ct+pi_pos+2984.54)<1.006||fabs(ct+pi_pos+3673.55)<1.006)ct_cut_ctout=true;
		if(fabs(ct+pi_pos+3650.07)<10.)ct_cut_ctout=true;
		else ct_cut_ctout=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum>3.&&ac2sum>5.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		//else event_selection=false;
		if(ct<-80.)event_selection_ctout=true;
		else event_selection_ctout=false;
		if(ct<-80.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_ctout2=true;
		else event_selection_ctout2=false;
		if(ct>150.)event_selection_ctout3=true;
		else event_selection_ctout3=false;
		if(ct>150.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_ctout4=true;
		else event_selection_ctout4=false;

		event_selection_nocut=false;
		if(zcut_flag){
			if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(accut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(fpcut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&L_Tr)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(trcut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_FP&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}


	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
	    double R_px = R_pz * (R_tr_tg_th );
	    double R_py = R_pz * ( R_tr_tg_ph );

	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(( L_tr_tg_th ), 2.0) + pow(( L_tr_tg_ph),2.0));
	    double L_px = L_pz * ( L_tr_tg_th );
	    double L_py = L_pz * ( L_tr_tg_ph );

	    double B_E =sqrt(B_mom*B_mom + Me*Me);
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
	    double L_E =sqrt(L_mom*L_mom + Me*Me);

		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);

	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );

        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
	    mm=mass - mh;//shift by ML
		
		h_ctmm->Fill(ct,mm);
		h_ct->Fill(ct);
		if(event_selection)h_ctmm2->Fill(ct,mm);
		if(event_selection)h_ct2->Fill(ct);
		if(event_selection&&ct_cut_ctout)hmm_ctout_only->Fill(mm);
		if(event_selection_ctout)hmm_ctout->Fill(mm);
		if(event_selection_ctout2)hmm_ctout2->Fill(mm);
		if(event_selection_ctout3)hmm_ctout3->Fill(mm);
		if(event_selection_ctout4)hmm_ctout4->Fill(mm);
		if(event_selection&&ct_cut)hmm_L_fom_best->Fill(mm);
		if(event_selection&&ct_cut)hmm_L_strict->Fill(mm*1000.);
		if(event_selection_nocut&&ct_cut)hmm_L_fom_nocut->Fill(mm);
		double theta_ee = L_4vec.Theta();
		//test double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double phi_L = atan((phi*cos(phi0)+sin(phi0))/theta);//LHRS frame
		double phi_ee = L_4vec.Phi();//original frame
		double phi_ek = R_4vec.Phi()+2*PI;//original frame

		G_4vec = B_4vec - L_4vec;
		//double mom_g = G_4vec.Rho();
		double mom_g=sqrt(G_4vec.Px()*G_4vec.Px()+G_4vec.Py()*G_4vec.Py()+G_4vec.Pz()*G_4vec.Pz());
		double Qsq = G_4vec.M()*G_4vec.M();
		double phi_g = G_4vec.Phi()+2*PI;
		double theta_g = G_4vec.Theta();
		double theta_gk_lab = G_4vec.Angle(R_4vec.Vect());
		double omega=G_4vec.E();
		double beta=mom_g/(omega+Mp);
	
		TVector3 boost;
		TLorentzVector GT_4vec;
		GT_4vec=G_4vec+T_4vec;
		boost=GT_4vec.BoostVector();
		R_4vec.Boost(-boost);
		L_4vec.Boost(-boost);
		B_4vec.Boost(-boost);
		double theta_gk_cm = G_4vec.Angle(R_4vec.Vect());
		double theta_ek_cm = R_4vec.Theta();
		double pR_cm=sqrt(R_4vec.Px()*R_4vec.Px()+R_4vec.Py()*R_4vec.Py()+R_4vec.Pz()*R_4vec.Pz());
		double pL_cm=sqrt(L_4vec.Px()*L_4vec.Px()+L_4vec.Py()*L_4vec.Py()+L_4vec.Pz()*L_4vec.Pz());
		double pB_cm=sqrt(B_4vec.Px()*B_4vec.Px()+B_4vec.Py()*B_4vec.Py()+B_4vec.Pz()*B_4vec.Pz());

		double n = MK/ML;
		double p_cm=sqrt(GT_4vec.Px()*GT_4vec.Px()+GT_4vec.Py()*GT_4vec.Py()+GT_4vec.Pz()*GT_4vec.Pz());
		double E_cm = GT_4vec.E();
//beta=2.3/(2.2+Mp);
//pR_cm=0.65;
//theta_gk_cm=0.12;
		double gamma=1./sqrt(1-beta*beta);
		double ER_cm=sqrt(pR_cm*pR_cm+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm = (gamma*pR_cm*pR_cm*(pR_cm*cos(theta_gk_cm)+beta*ER_cm))/(pow(sqrt(pR_cm*pR_cm*sin(theta_gk_cm)*sin(theta_gk_cm)+gamma*gamma*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)),3.));
//cout<<"labtocm="<<labtocm<<endl;
		double tan_lab1 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+beta*sqrt(MK*MK+pR_cm*pR_cm)/pR_cm));
		double tan_lab2 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+(omega*Mp-Qsq*Qsq)/(omega*Mp+Mp*Mp)));
		//if(tan_lab1!=tan_lab2)cout<<"tan1="<<atan(tan_lab1)<<", tan2="<<atan(tan_lab2)<<"theta_gk_lab="<<theta_gk_lab<<endl;


		if(event_selection&&ct_cut){
			gklab_gkcm->Fill(theta_gk_lab,theta_gk_cm);
			gklab_eklab->Fill(theta_gk_lab,theta_ek);
			eelab_eklab->Fill(theta_ee,theta_ek);
			eklab_gkcm->Fill(theta_gk_cm,theta_ek);
			ekcm_gkcm->Fill(theta_gk_cm,theta_ek_cm);
			cos_gklab_gkcm->Fill(cos(theta_gk_lab),cos(theta_gk_cm));
			cos_gklab_eklab->Fill(cos(theta_gk_lab),cos(theta_ek));
			cos_eelab_eklab->Fill(cos(theta_ee),cos(theta_ek));
			cos_eklab_gkcm->Fill(cos(theta_gk_cm),cos(theta_ek));
			cos_ekcm_gkcm->Fill(cos(theta_gk_cm),cos(theta_ek_cm));
			h_pepk->Fill(R_mom*1000.,L_mom*1000.);
			h_pe->Fill(L_mom*1000.);
			h_pk->Fill(R_mom*1000.);
			h_R_x_data->Fill(R_tr_x*100.);
			h_R_y_data->Fill(R_tr_y*100.);
			h_R_th_data->Fill(R_tr_th);
			h_R_ph_data->Fill(R_tr_ph);
			h_R_tg_th_data->Fill(R_tr_tg_th);
			h_R_tg_ph_data->Fill(R_tr_ph);
			h_L_x_data->Fill(L_tr_x*100.);
			h_L_y_data->Fill(L_tr_y*100.);
//== FP x FP ==//
//X-Y, X-theta, Y-phi, theta-phi 
			h_L_x_y_data->Fill(L_tr_y*100.,L_tr_x*100.);
			h_R_x_y_data->Fill(R_tr_y*100.,R_tr_x*100.);
			h_L_th_x_data->Fill(L_tr_x*100.,L_tr_th);
			h_R_th_x_data->Fill(R_tr_x*100.,R_tr_th);
			h_L_ph_y_data->Fill(L_tr_y*100.,L_tr_ph);
			h_R_ph_y_data->Fill(R_tr_y*100.,R_tr_ph);
			h_L_ph_th_data->Fill(L_tr_th,L_tr_ph);
			h_R_ph_th_data->Fill(R_tr_th,R_tr_ph);

//== tar x tar ==//
//theta-phi, phi-Z
			h_L_tg_ph_tg_th_data->Fill(L_tr_tg_th,L_tr_tg_ph);
			h_R_tg_ph_tg_th_data->Fill(R_tr_tg_th,R_tr_tg_ph);
			h_L_z_tg_ph_data->Fill(L_tr_tg_ph,L_tr_vz*100.);
			h_R_z_tg_ph_data->Fill(R_tr_tg_ph,R_tr_vz*100.);

//== FP x tar ==//
//X-theta, Y-phi, theta-theta, phi-phi, Y-Z, phi-Z
			h_L_tg_th_x_data->Fill(L_tr_x*100.,L_tr_tg_th);
			h_R_tg_th_x_data->Fill(R_tr_x*100.,R_tr_tg_th);
			h_L_tg_th_th_data->Fill(L_tr_th,L_tr_tg_th);
			h_R_tg_th_th_data->Fill(R_tr_th,R_tr_tg_th);
			h_L_tg_ph_y_data->Fill(L_tr_y*100.,L_tr_tg_ph);
			h_R_tg_ph_y_data->Fill(R_tr_y*100.,R_tr_tg_ph);
			h_L_tg_ph_ph_data->Fill(L_tr_ph,L_tr_tg_ph);
			h_R_tg_ph_ph_data->Fill(R_tr_ph,R_tr_tg_ph);
			h_L_z_y_data->Fill(L_tr_y*100.,L_tr_vz*100.);
			h_R_z_y_data->Fill(R_tr_y*100.,R_tr_vz*100.);
			h_L_z_ph_data->Fill(L_tr_ph,L_tr_vz*100.);
			h_R_z_ph_data->Fill(R_tr_ph,R_tr_vz*100.);

			//h_L_y_vz_data->Fill(L_tr_y*100.,L_tr_vz*100.);
			//h_R_y_vz_data->Fill(R_tr_y*100.,R_tr_vz*100.);
			//h_L_ph_vz_data->Fill(L_tr_ph,L_tr_vz*100.);
			//h_R_ph_vz_data->Fill(R_tr_ph,R_tr_vz*100.);
			h_L_p_x_data->Fill(L_mom*1000.,L_tr_x*100.);
			h_R_p_x_data->Fill(R_mom*1000.,R_tr_x*100.);
			h_L_th_data->Fill(L_tr_th);
			h_L_ph_data->Fill(L_tr_ph);
			h_L_tg_th_data->Fill(L_tr_tg_th);
			h_L_tg_ph_data->Fill(L_tr_ph);
			//h_L_tg_ph_y_data->Fill(-1.*L_tr_y*100.,L_tr_tg_ph);
			//h_R_tg_ph_y_data->Fill(-1.*R_tr_y*100.,R_tr_tg_ph);
			//h_L_tg_th_x_data->Fill(L_tr_x*100.,L_tr_th);
			//h_R_tg_th_x_data->Fill(R_tr_x*100.,R_tr_th);
			//if(R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)h_pe->Fill(L_mom*1000.);
			//if(R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)h_pk->Fill(R_mom*1000.);
		}

		//if(abs(R_tr_vz-L_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		if(abs(R_tr_vz-L_tr_vz)<0.025)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		h_nltrack->Fill(NLtr);
		h_nrtrack->Fill(NRtr);

}//ENum

cout<<"Integral"<<endl;
cout<<"h_pe(mom_cut)="<<h_pe->Integral(h_pe->FindBin(2010.),h_pe->FindBin(2160.))<<endl;
cout<<"h_pe="<<h_pe->Integral()<<endl;
cout<<"h_pk(mom_cut)="<<h_pk->Integral(h_pk->FindBin(1760.),h_pk->FindBin(1900.))<<endl;
cout<<"h_pk="<<h_pk->Integral()<<endl;
cout<<"h_pe_simc(mom_cut)="<<h_pe_simc->Integral(h_pe_simc->FindBin(2010.),h_pe_simc->FindBin(2160.))<<endl;
cout<<"h_pe_simc="<<h_pe_simc->Integral()<<endl;
cout<<"h_pk_simc(mom_cut)="<<h_pk_simc->Integral(h_pk_simc->FindBin(1760.),h_pk_simc->FindBin(1900.))<<endl;
cout<<"h_pk_simc="<<h_pk_simc->Integral()<<endl;
cout<<"GetEntries"<<endl;
cout<<"h_pe="<<h_pe->GetEntries()<<endl;
cout<<"h_pk="<<h_pk->GetEntries()<<endl;
cout<<"h_pe_simc="<<h_pe_simc->GetEntries()<<endl;
cout<<"h_pk_simc="<<h_pk_simc->GetEntries()<<endl;

///-----********-----///
// c10: pe, pk, simc
// c20: pe, pk, data
// c30: pe, pk, same
// c40: pe, pk, chi-square
//
// c50: FP_HRS-R 1D
// c51: FP_HRS-L 1D
//
// c60: 
// c61: 
// c64: 
///-----********-----///


	double dtemp, stemp;
	TCanvas* c10 = new TCanvas("c10","c10",900.,900.);
	c10->Divide(2,2);
	c10->cd(1);
	h_pe_simc->SetNdivisions(505);
	h_pe_simc->Draw("");
	c10->cd(2);
	h_pk_simc->SetNdivisions(505);
	h_pk_simc->Draw("");
	c10->cd(3);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.15);
	gPad->SetBottomMargin(0.15);
	h_pepk_simc->Draw("colz");
	TCanvas* c20 = new TCanvas("c20","c20",900.,900.);
	c20->Divide(2,2);
	c20->cd(1);
	h_pe->Draw("");
	c20->cd(2);
	h_pk->Draw("");
	c20->cd(3);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.15);
	gPad->SetBottomMargin(0.15);
	h_pepk->Draw("colz");
	TCanvas* c30 = new TCanvas("c30","c30",900.,900.);
	h_pe->SetLineColor(kAzure);
	h_pk->SetLineColor(kAzure);
	h_pe_simc->SetLineColor(kRed);
	h_pk_simc->SetLineColor(kRed);
	c30->Divide(2,2);
	c30->cd(1);
	h_pe->Draw("");
	//h_pe_simc->Scale(1863./922483.);//orig
	//h_pe_simc->Scale(2421./1.2e+6);
	//h_pe_simc->Scale(2081./923474.);//rec
	stemp=h_pe_simc->Integral();
	dtemp=h_pe->Integral();
	cout<<"stemp="<<stemp<<endl;
	cout<<"dtemp="<<dtemp<<endl;
	h_pe_simc->Scale(dtemp/stemp);
	h_pe_simc->Draw("same");
	c30->cd(2);
	h_pk->SetNdivisions(505);
	h_pk->Draw("");
	//h_pk_simc->Scale(2226./1041430.);//rec
	//h_pk_simc->Scale(2286./1040300.);//orig
	//h_pk_simc->Scale(1863./1040300.);
	//h_pk_simc->Scale(2421./1.2e+6);
	stemp=h_pk_simc->Integral();
	dtemp=h_pk->Integral();
	h_pk_simc->Scale(dtemp/stemp);
	h_pk_simc->Draw("same");
	c30->cd(3);
	h_pe->Draw("e");
	h_pe_simc->Draw("same");
	//h_pe_simc_FPcut->Scale(2081./923474.);//rec
	stemp=h_pe_simc_FPcut->Integral();
	dtemp=h_pe->Integral();
	h_pe_simc_FPcut->Scale(dtemp/stemp);
	h_pe_simc_FPcut->SetLineColor(kGreen);
	h_pe_simc_FPcut->Draw("same");
	c30->cd(4);
	h_pk->Draw("e");
	h_pk_simc->Draw("same");
	//h_pk_simc_FPcut->Scale(2226./1041430.);//rec
	stemp=h_pk_simc_FPcut->Integral();
	dtemp=h_pk->Integral();
	h_pk_simc_FPcut->Scale(dtemp/stemp);
	h_pk_simc_FPcut->SetLineColor(kGreen);
	h_pk_simc_FPcut->Draw("same");
	
	TCanvas* c40 = new TCanvas("c40","c40",900.,900.);
    TH1D* h_pe_chisq = new TH1D("h_pe_chisq", "#chi_{e'}^{2}" ,200,1980.,2220.);
    TH1D* h_pk_chisq = new TH1D("h_pk_chisq", "#chi_{K}^{2}" ,200,1720.,1940.);
	double e_data, e_expt, e_temp, e_chisq=0.;
	double k_data, k_expt, k_temp, k_chisq=0.;
	int e_bin=0, k_bin=0;
	for(int i=0;i<200;i++){
		e_data = h_pe->GetBinContent(i+1);
		e_expt = h_pe_simc->GetBinContent(i+1);
		if(e_expt==0.)e_temp=0.01;
		else e_temp = (e_data - e_expt)*(e_data - e_expt)/e_expt;
		//e_temp = (e_data - e_expt)*(e_data - e_expt)/e_expt;
		if(h_pe->GetBinCenter(i+1)>2010.&&h_pe->GetBinCenter(i+1)<2160.){e_chisq += e_temp;e_bin++;}
		//e_chisq += e_temp;e_bin++;
		h_pe_chisq->SetBinContent(i+1,e_temp);
		k_data = h_pk->GetBinContent(i+1);
		k_expt = h_pk_simc->GetBinContent(i+1);
		if(k_expt==0.)k_temp=0.01;
		else k_temp = (k_data - k_expt)*(k_data - k_expt)/k_expt;
		//k_temp = (k_data - k_expt)*(k_data - k_expt)/k_expt;
		if(h_pk->GetBinCenter(i+1)>1760.&&h_pk->GetBinCenter(i+1)<1900.){k_chisq += k_temp;k_bin++;}
		//k_chisq += k_temp;k_bin++;
		h_pk_chisq->SetBinContent(i+1,k_temp);
	}
	c40->Divide(2,2);
	c40->cd(1);
	h_pe->Draw("e");
	h_pe_simc->Draw("same");
	c40->cd(2);
	h_pk->Draw("e");
	h_pk_simc->Draw("same");
	c40->cd(3);
	h_pe_chisq->Draw("");
	c40->cd(4);
	h_pk_chisq->SetNdivisions(505);
	h_pk_chisq->Draw("");
	
	cout<<"e chi-square = "<<e_chisq<<endl;
	cout<<"e bin = "<<e_bin<<endl;
	cout<<"k chi-square = "<<k_chisq<<endl;
	cout<<"k bin = "<<k_bin<<endl;
//	c3->Print("./pdf/mm_tight.pdf");
//	c4->Print("./pdf/mm_wobg_tight.pdf");

cout << "h_R_y_simc = "<<h_R_y_simc->Integral()<<endl;
cout << "h_R_y_data = "<<h_R_y_data->Integral()<<endl;
cout << "h_R_x_simc = "<<h_R_x_simc->Integral()<<endl;
cout << "h_R_x_data = "<<h_R_x_data->Integral()<<endl;
cout << "h_R_th_simc = "<<h_R_th_simc->Integral()<<endl;
cout << "h_R_ph_data = "<<h_R_ph_data->Integral()<<endl;
cout << "h_R_tg_th_simc = "<<h_R_tg_th_simc->Integral()<<endl;
cout << "h_R_tg_ph_data = "<<h_R_tg_ph_data->Integral()<<endl;
cout << "h_L_y_simc = "<<h_L_y_simc->Integral()<<endl;
cout << "h_L_y_data = "<<h_L_y_data->Integral()<<endl;
cout << "h_L_x_simc = "<<h_L_x_simc->Integral()<<endl;
cout << "h_L_x_data = "<<h_L_x_data->Integral()<<endl;
cout << "h_L_th_simc = "<<h_L_th_simc->Integral()<<endl;
cout << "h_L_ph_data = "<<h_L_ph_data->Integral()<<endl;
cout << "h_L_tg_th_simc = "<<h_L_tg_th_simc->Integral()<<endl;
cout << "h_L_tg_ph_data = "<<h_L_tg_ph_data->Integral()<<endl;

	TCanvas* c50 = new TCanvas("c50","c50",900.,900.);
	c50->Divide(3,2);
	c50->cd(1);
	stemp=h_R_y_simc->Integral();
	dtemp=h_R_y_data->Integral();
	cout<<"stemp="<<stemp<<endl;
	cout<<"dtemp="<<dtemp<<endl;
	h_R_y_simc->Scale(dtemp/stemp);
	h_R_y_data->Draw("e");
	//h_R_y_data->Fit("gausn");
	h_R_y_simc->Draw("histsame");
	//h_R_y_simc->Fit("gausn");
	c50->cd(2);
	stemp=h_R_x_simc->Integral();
	dtemp=h_R_x_data->Integral();
	h_R_x_simc->Scale(dtemp/stemp);
	h_R_x_data->Draw("e");
	h_R_x_simc->Draw("histsame");
	c50->cd(3);
	stemp=h_R_th_simc->Integral();
	dtemp=h_R_th_data->Integral();
	h_R_th_simc->Scale(dtemp/stemp);
	h_R_th_data->Draw("e");
	h_R_th_simc->Draw("histsame");
	c50->cd(4);
	stemp=h_R_ph_simc->Integral();
	dtemp=h_R_ph_data->Integral();
	h_R_ph_simc->Scale(dtemp/stemp);
	h_R_ph_data->Draw("e");
	h_R_ph_simc->Draw("histsame");
	c50->cd(5);
	stemp=h_R_tg_th_simc->Integral();
	dtemp=h_R_tg_th_data->Integral();
	h_R_tg_th_simc->Scale(dtemp/stemp);
	h_R_tg_th_data->Draw("e");
	h_R_tg_th_simc->Draw("histsame");
	c50->cd(6);
	stemp=h_R_tg_ph_simc->Integral();
	dtemp=h_R_tg_ph_data->Integral();
	h_R_tg_ph_simc->Scale(dtemp/stemp);
	h_R_tg_ph_data->Draw("e");
	h_R_tg_ph_simc->Draw("histsame");
	
	TCanvas* c51 = new TCanvas("c51","c51",900.,900.);
	c51->Divide(3,2);
	c51->cd(1);
	stemp=h_L_y_simc->Integral();
	dtemp=h_L_y_data->Integral();
	h_L_y_simc->Scale(dtemp/stemp);
	h_L_y_data->Draw("e");
	//h_L_y_data->Fit("gausn");
	h_L_y_simc->Draw("histsame");
	//h_L_y_simc->Fit("gausn");
	c51->cd(2);
	stemp=h_L_x_simc->Integral();
	dtemp=h_L_x_data->Integral();
	h_L_x_simc->Scale(dtemp/stemp);
	h_L_x_data->Draw("e");
	h_L_x_simc->Draw("histsame");
	c51->cd(3);
	stemp=h_L_th_simc->Integral();
	dtemp=h_L_th_data->Integral();
	h_L_th_simc->Scale(dtemp/stemp);
	h_L_th_data->Draw("e");
	h_L_th_simc->Draw("histsame");
	c51->cd(4);
	stemp=h_L_ph_simc->Integral();
	dtemp=h_L_ph_data->Integral();
	h_L_ph_simc->Scale(dtemp/stemp);
	h_L_ph_data->Draw("e");
	h_L_ph_simc->Draw("histsame");
	c51->cd(5);
	stemp=h_L_tg_th_simc->Integral();
	dtemp=h_L_tg_th_data->Integral();
	h_L_tg_th_simc->Scale(dtemp/stemp);
	h_L_tg_th_data->Draw("e");
	h_L_tg_th_simc->Draw("histsame");
	c51->cd(6);
	stemp=h_L_tg_ph_simc->Integral();
	dtemp=h_L_tg_ph_data->Integral();
	h_L_tg_ph_simc->Scale(dtemp/stemp);
	h_L_tg_ph_data->Draw("e");
	h_L_tg_ph_simc->Draw("histsame");

cout << "h_R_y_simc = "<<h_R_y_simc->Integral()<<endl;
cout << "h_R_y_data = "<<h_R_y_data->Integral()<<endl;
cout << "h_R_x_simc = "<<h_R_x_simc->Integral()<<endl;
cout << "h_R_x_data = "<<h_R_x_data->Integral()<<endl;
cout << "h_R_th_simc = "<<h_R_th_simc->Integral()<<endl;
cout << "h_R_ph_data = "<<h_R_ph_data->Integral()<<endl;
cout << "h_R_tg_th_simc = "<<h_R_tg_th_simc->Integral()<<endl;
cout << "h_R_tg_ph_data = "<<h_R_tg_ph_data->Integral()<<endl;
cout << "h_L_y_simc = "<<h_L_y_simc->Integral()<<endl;
cout << "h_L_y_data = "<<h_L_y_data->Integral()<<endl;
cout << "h_L_x_simc = "<<h_L_x_simc->Integral()<<endl;
cout << "h_L_x_data = "<<h_L_x_data->Integral()<<endl;
cout << "h_L_th_simc = "<<h_L_th_simc->Integral()<<endl;
cout << "h_L_ph_data = "<<h_L_ph_data->Integral()<<endl;
cout << "h_L_tg_th_simc = "<<h_L_tg_th_simc->Integral()<<endl;
cout << "h_L_tg_ph_data = "<<h_L_tg_ph_data->Integral()<<endl;

	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadTopMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

//	TCanvas* c60 = new TCanvas("c60","c60",900.,900.);
//	c60->Divide(2,2);
//	c60->cd(1);
//	h_L_tg_ph_y_simc->Draw("colz");
//	c60->cd(2);
//	h_L_tg_ph_y_data->Draw("colz");
//	c60->cd(3);
//	h_L_tg_th_x_simc->Draw("colz");
//	c60->cd(4);
//	h_L_tg_th_x_data->Draw("colz");
//	TCanvas* c61 = new TCanvas("c61","c61",900.,900.);
//	c61->Divide(2,2);
//	c61->cd(1);
//	h_R_tg_ph_y_simc->Draw("colz");
//	c61->cd(2);
//	h_R_tg_ph_y_data->Draw("colz");
//	c61->cd(3);
//	h_R_tg_th_x_simc->Draw("colz");
//	c61->cd(4);
//	h_R_tg_th_x_data->Draw("colz");
//	
//	TCanvas* c62 = new TCanvas("c62","c62",900.,900.);
//	c62->Divide(2,2);
//	c62->cd(1);
//  	h_R_y_vz_simc->Draw("colz");    
//	c62->cd(2);
//	h_R_y_vz_data->Draw("colz");    
//	c62->cd(3);
//  	h_L_y_vz_simc->Draw("colz");
//	c62->cd(4);
//  	h_L_y_vz_data->Draw("colz");
//
//	TCanvas* c63 = new TCanvas("c63","c63",900.,900.);
//	c63->Divide(2,2);
//	c63->cd(1);
//  	h_R_ph_vz_simc->Draw("colz");    
//	c63->cd(2);
//	h_R_ph_vz_data->Draw("colz");    
//	c63->cd(3);
//  	h_L_ph_vz_simc->Draw("colz");
//	c63->cd(4);
//  	h_L_ph_vz_data->Draw("colz");




//Each Pad Size are changed.
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadTopMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);

	TCanvas* c64 = new TCanvas("c64","c64",900.,900.);
	c64->Divide(2,2);
	c64->cd(1);
  	h_R_p_x_simc->Draw("colz");    
	c64->cd(2);
  	h_R_p_x_data->Draw("colz");    
	c64->cd(3);
  	h_L_p_x_simc->Draw("colz");    
	c64->cd(4);
  	h_L_p_x_data->Draw("colz");    

	TCanvas* c70 = new TCanvas("c70","c70 (FP)",900.,900.);
	c70->Divide(2,2);
	c70->cd(1);
  	h_R_x_y_simc->Draw("colz");    
	c70->cd(2);
  	h_R_x_y_data->Draw("colz");    
	c70->cd(3);
  	h_L_x_y_simc->Draw("colz");    
	c70->cd(4);
  	h_L_x_y_data->Draw("colz");    

	TCanvas* c71 = new TCanvas("c71","c71 (FP)",900.,900.);
	c71->Divide(2,2);
	TF1* func_fp1 = new TF1("func_fp1",FP_cut1,-80.,80.,1);
	func_fp1->SetNpx(600);
	func_fp1->SetParameter(1,0);
	func_fp1->SetLineColor(kRed);
	func_fp1->SetLineWidth(4);
	TF1* func_fp2 = new TF1("func_fp2",FP_cut2,-80.,80.,1);
	func_fp2->SetNpx(600);
	func_fp2->SetParameter(1,0);
	func_fp2->SetLineColor(kRed);
	func_fp2->SetLineWidth(4);
	TF1* func_fp3 = new TF1("func_fp3",FP_cut3,-80.,80.,1);
	func_fp3->SetNpx(600);
	func_fp3->SetParameter(1,0);
	func_fp3->SetLineColor(kRed);
	func_fp3->SetLineWidth(4);
	TF1* func_fp4 = new TF1("func_fp4",FP_cut4,-80.,80.,1);
	func_fp4->SetNpx(600);
	func_fp4->SetParameter(1,0);
	func_fp4->SetLineColor(kViolet);
	func_fp4->SetLineWidth(4);
	c71->cd(1);
  	h_R_th_x_simc->Draw("colz");    
	func_fp1->Draw("same");
	func_fp2->Draw("same");
	func_fp3->Draw("same");
	func_fp4->Draw("same");
	c71->cd(2);
  	h_R_th_x_data->Draw("colz");    
	func_fp1->Draw("same");
	func_fp2->Draw("same");
	func_fp3->Draw("same");
	func_fp4->Draw("same");
	c71->cd(3);
  	h_L_th_x_simc->Draw("colz");    
	func_fp1->Draw("same");
	func_fp2->Draw("same");
	func_fp3->Draw("same");
	func_fp4->Draw("same");
	c71->cd(4);
  	h_L_th_x_data->Draw("colz");    
	func_fp1->Draw("same");
	func_fp2->Draw("same");
	func_fp3->Draw("same");
	func_fp4->Draw("same");

	TCanvas* c72 = new TCanvas("c72","c72 (FP)",900.,900.);
	c72->Divide(2,2);
	c72->cd(1);
  	h_R_ph_y_simc->Draw("colz");    
	c72->cd(2);
  	h_R_ph_y_data->Draw("colz");    
	c72->cd(3);
  	h_L_ph_y_simc->Draw("colz");    
	c72->cd(4);
  	h_L_ph_y_data->Draw("colz");    

	TCanvas* c73 = new TCanvas("c73","c73 (FP)",900.,900.);
	c73->Divide(2,2);
	c73->cd(1);
  	h_R_ph_th_simc->Draw("colz");    
	c73->cd(2);
  	h_R_ph_th_data->Draw("colz");    
	c73->cd(3);
  	h_L_ph_th_simc->Draw("colz");    
	c73->cd(4);
  	h_L_ph_th_data->Draw("colz");    

	TCanvas* c80 = new TCanvas("c80","c80 (tar)",900.,900.);
	c80->Divide(2,2);
	c80->cd(1);
  	h_R_tg_ph_tg_th_simc->Draw("colz");    
	c80->cd(2);
  	h_R_tg_ph_tg_th_data->Draw("colz");    
	c80->cd(3);
  	h_L_tg_ph_tg_th_simc->Draw("colz");    
	c80->cd(4);
  	h_L_tg_ph_tg_th_data->Draw("colz");    

	TCanvas* c81 = new TCanvas("c81","c81 (tar)",900.,900.);
	c81->Divide(2,2);
	c81->cd(1);
  	h_R_z_tg_ph_simc->Draw("colz");    
	c81->cd(2);
  	h_R_z_tg_ph_data->Draw("colz");    
	c81->cd(3);
  	h_L_z_tg_ph_simc->Draw("colz");    
	c81->cd(4);
  	h_L_z_tg_ph_data->Draw("colz");    

	TCanvas* c90 = new TCanvas("c90","c90 (FP x tar)",900.,900.);
	c90->Divide(2,2);
	c90->cd(1);
	h_R_tg_th_x_simc->Draw("colz");
	c90->cd(2);
	h_R_tg_th_x_data->Draw("colz");
	c90->cd(3);
	h_L_tg_th_x_simc->Draw("colz");
	c90->cd(4);
	h_L_tg_th_x_data->Draw("colz");

	TCanvas* c91 = new TCanvas("c91","c91 (FP x tar)",900.,900.);
	c91->Divide(2,2);
	c91->cd(1);
	h_R_tg_th_th_simc->Draw("colz");
	c91->cd(2);
	h_R_tg_th_th_data->Draw("colz");
	c91->cd(3);
	h_L_tg_th_th_simc->Draw("colz");
	c91->cd(4);
	h_L_tg_th_th_data->Draw("colz");

	TCanvas* c92 = new TCanvas("c92","c92 (FP x tar)",900.,900.);
	c92->Divide(2,2);
	c92->cd(1);
	h_R_tg_ph_y_simc->Draw("colz");
	c92->cd(2);
	h_R_tg_ph_y_data->Draw("colz");
	c92->cd(3);
	h_L_tg_ph_y_simc->Draw("colz");
	c92->cd(4);
	h_L_tg_ph_y_data->Draw("colz");

	TCanvas* c93 = new TCanvas("c93","c93 (FP x tar)",900.,900.);
	c93->Divide(2,2);
	c93->cd(1);
	h_R_tg_ph_ph_simc->Draw("colz");
	c93->cd(2);
	h_R_tg_ph_ph_data->Draw("colz");
	c93->cd(3);
	h_L_tg_ph_ph_simc->Draw("colz");
	c93->cd(4);
	h_L_tg_ph_ph_data->Draw("colz");

	TCanvas* c94 = new TCanvas("c94","c94 (tar)",900.,900.);
	c94->Divide(2,2);
	c94->cd(1);
  	h_R_z_y_simc->Draw("colz");    
	c94->cd(2);
  	h_R_z_y_data->Draw("colz");    
	c94->cd(3);
  	h_L_z_y_simc->Draw("colz");    
	c94->cd(4);
  	h_L_z_y_data->Draw("colz");    

	TCanvas* c95 = new TCanvas("c95","c95 (tar)",900.,900.);
	c95->Divide(2,2);
	c95->cd(1);
  	h_R_z_ph_simc->Draw("colz");    
	c95->cd(2);
  	h_R_z_ph_data->Draw("colz");    
	c95->cd(3);
  	h_L_z_ph_simc->Draw("colz");    
	c95->cd(4);
  	h_L_z_ph_data->Draw("colz");    

/*--- Print ---*/
cout << "Print is starting" << endl;
	//"SIMC_vs_DATA.pdf
	//"SIMC_vs_DATA_newZ.pdf
	//"SIMC_vs_DATA_cell_x10.pdf
	//"SIMC_vs_noCalib.pdf
	//"SIMC_FP_Rebuild_all.pdf //2021/7/11 simc_fp3.C
	//"SIMC_FP_Rebuild_mom.pdf //2021/7/11 simc_fp_momcut.C
	
	c10->Print("SIMC_FP_Rebuilt_all.pdf[");
	c10->Print("SIMC_FP_Rebuilt_all.pdf[");//pepk simc
	c20->Print("SIMC_FP_Rebuilt_all.pdf");//pepk data
	c30->Print("SIMC_FP_Rebuilt_all.pdf");//pepk same
	c50->Print("SIMC_FP_Rebuilt_all.pdf");//FP 1D
	c51->Print("SIMC_FP_Rebuilt_all.pdf");//FP 1D
	//c60->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	//c61->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c64->Print("SIMC_FP_Rebuilt_all.pdf");//FP vs Mom. 2D
	c70->Print("SIMC_FP_Rebuilt_all.pdf");//FP 2D
	c71->Print("SIMC_FP_Rebuilt_all.pdf");//FP 2D
	c72->Print("SIMC_FP_Rebuilt_all.pdf");//FP 2D
	c73->Print("SIMC_FP_Rebuilt_all.pdf");//FP 2D
	c80->Print("SIMC_FP_Rebuilt_all.pdf");//tar 2D
	c81->Print("SIMC_FP_Rebuilt_all.pdf");//tar 2D
	c90->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c91->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c92->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c93->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c94->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c95->Print("SIMC_FP_Rebuilt_all.pdf");//FPxtar 2D
	c95->Print("SIMC_FP_Rebuilt_all.pdf]");//FPxtar 2D
	
cout << "Well done!" << endl;
}//fit
