//--------------------------------//
//--  Fitting w/ Response func. --//
//--------------------------------//
//
//K. Okuyama (Aug. 22, 2020)
//
//This is taken over from fit_many.C
//No array branch mode 

double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

double F_2Gauss( double *x, double *par )
  {
	double val = par[0] * (TMath::Gaus(x[0],par[1],par[2]) + par[3] * TMath::Gaus(x[0],par[4],par[5]));
    return val;
  }

double F_1Gauss( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
    return val;
  }

double FMM_Voigt( double *x, double *par )
  {
    // par[6] : area
    // par[7] : location
    // par[8] : gaussian sigma
    // par[9] : lorentz fwhm
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    val += par[6] * TMath::Voigt(x[0]-par[7],par[8],par[9],4);
    return val;
  }

double FMM_1Gauss( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
	val += par[6] * TMath::Gaus(x[0],par[7],par[8]);
    return val;
  }

double FMM_2Gauss( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
	val += par[6] * (TMath::Gaus(x[0],par[7],par[8]) + par[9] * TMath::Gaus(x[0],par[10],par[11]));
    return val;
  }

double FMM_2Poly( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    val += par[6] + par[7]*x[0] + par[8]*x[0]*x[0];
    return val;
  }

double FMM_3Poly( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    val += par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0];
    return val;
  }

double FMM_4Poly( double *x, double *par )
  {
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    val += par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0] + par[10]*TMath::Power(x[0],4.);
    return val;
  }

double FMM_4Poly_wRes( double *x, double *par )
  {
	double Napier = 2.7182818;
	double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);//Lambda Gaussian
	//double val = par[0] * TMath::Landau(x[0],par[1],par[2]);//Lambda Landau 
	val += par[3] * TMath::Gaus(x[0],par[4],par[5]);//Sigma Gaussian
    val += par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0] + par[10]*TMath::Power(x[0],4.);//4Poly BG
	if((x[0]-par[1])>0)val += (par[11]) * TMath::Power(Napier,-par[12]*(x[0]-par[13]))*TMath::Gaus(x[0],par[13],par[16]);//Lambda Radiative tail
	//if((x[0]-par[1])>0)val += (par[11]/4.) * TMath::Power(Napier,-par[17]*(x[0]-par[18]))*TMath::Gaus(x[0],par[13],par[16]);//Lambda Radiative tail
	//if((x[0]-par[4]-2*par[5])>0)val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
	if((x[0]-par[4]>0))val += par[14] * TMath::Power(Napier,-par[12]*(x[0]-par[15]))*TMath::Gaus(x[0],par[15],par[16]);//Simga Radiative tail//relevant to Lambda Radiative tail
	//val += par[11] * TMath::Power(Napier,-par[12]*(x[0]-par[13]));//Lambda Radiative tail
	//val += par[14] * TMath::Power(Napier,-par[15]*(x[0]-par[16]));//Simga Radiative tail
    return val;
  }

double expgaus2(double *x, double *par, int num) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;
// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[num+2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[num+3];
     fland = TMath::Gaus(xx,x[0],par[num+2]);
     sum += fland * TMath::Exp(-xx/par[num+1]);
     xx = xupp - (i-.5) * step - par[num+3];
     fland = TMath::Gaus(xx,x[0],par[num+2]);
     sum += fland * TMath::Exp(-xx/par[num+1]);
  }
  //val = par[2] * step * sum * invsq2pi / par[3];
  val = par[num] * step * sum * invsq2pi / (par[num+2]*par[num+1]*exp(-par[num+3]/par[num+1]));
  return val;
}

double FMM_Lambda_Sigma( double *x, double *par , int num)
  {
  double val = par[num] * TMath::Gaus(x[0],par[num+1],par[num+2]);//Lambda Gaussian
  val += par[num+3] * TMath::Gaus(x[0],par[num+4],par[num+5]);//Sigma Gaussian
  return val;
}
double FMM_Res( double *x, double *par ){

	return FMM_Lambda_Sigma(x,par,0)+expgaus2(x,par,6)+expgaus2(x,par,10)+expgaus2(x,par,14);

}
double FMM_Res_Lambdaonly( double *x, double *par ){

	return par[0]*TMath::Gaus(x[0],par[1],par[2])+expgaus2(x,par,3)+expgaus2(x,par,7);

}



void fit_resp(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("h2all.root","read");//input file of all H2 run(default: h2all4.root)
	//ACCBGの引き算はmea_hist.ccから
  TFile *file_mea = new TFile("bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  double nbunch = 600.;//effetive bunches (6 bunches x 5 mixtures)
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
//	tree->Write();

    
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
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;


//---For Efficiency Result---//
 double center_pi, center_k, center_p, center_L, center_S;
 double range_pi, range_k, range_p, range_L, range_S;
 double n_pi_nocut, n_k_nocut, n_p_nocut, n_L_nocut, n_S_nocut;
 double const_pi_nocut, const_k_nocut, const_p_nocut, const_L_nocut, const_S_nocut;
 double n_pi_best, n_k_best, n_p_best, n_L_best, n_S_best;
 double const_pi_best, const_k_best, const_p_best, const_L_best, const_S_best;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_nocut, mean_k_nocut, mean_p_nocut, mean_L_nocut, mean_S_nocut;
 double mean_pi_best, mean_k_best, mean_p_best, mean_L_best, mean_S_best;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_nocut, sig_k_nocut, sig_p_nocut, sig_L_nocut, sig_S_nocut;
 double sig_pi_best, sig_k_best, sig_p_best, sig_L_best, sig_S_best;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];


//---Fitting Function---//
 TF1* fmmbg_nocut;
 TF1* fmm_pi_nocut_Voigt;
 TF1* fmm_pi_nocut_4Poly;
 TF1* fmm_pi_nocut_3Poly;
 TF1* fmm_pi_nocut_2Poly;
 TF1* fmm_pi_nocut_2Gauss;
 TF1* fmm_pi_nocut_1Gauss;
 TF1* fL_nocut;
 TF1* fS_nocut;
 TF1* fmm_nocut;
 TF1* fmm_nocut_Voigt;
 TF1* fmm_nocut_4Poly;
 TF1* fmm_nocut_3Poly;
 TF1* fmm_nocut_2Poly;
 TF1* fmm_nocut_2Gauss;
 TF1* fmm_nocut_1Gauss;
 TF1* fmmbg_nocut_Voigt;
 TF1* fmmbg_nocut_4Poly;
 TF1* fmmbg_nocut_3Poly;
 TF1* fmmbg_nocut_2Poly;
 TF1* fmmbg_nocut_2Gauss;
 TF1* fmmbg_nocut_1Gauss;
 TF1* fmmbg_best;
 TF1* fmm_pi_best_Voigt;
 TF1* fmm_pi_best_4Poly;
 TF1* fmm_pi_best_3Poly;
 TF1* fmm_pi_best_2Poly;
 TF1* fmm_pi_best_2Gauss;
 TF1* fmm_pi_best_1Gauss;
 TF1* fL_best;
 TF1* fS_best;
 TF1* fmm_best;
 TF1* fmm_best_Voigt;
 TF1* fmm_best_4Poly;
 TF1* fmm_best_3Poly;
 TF1* fmm_best_2Poly;
 TF1* fmm_best_2Gauss;
 TF1* fmm_best_1Gauss;
 TF1* fmmbg_best_Voigt;
 TF1* fmmbg_best_4Poly;
 TF1* fmmbg_best_3Poly;
 TF1* fmmbg_best_2Poly;
 TF1* fmmbg_best_2Gauss;
 TF1* fmmbg_best_1Gauss;
 TF1* facc[100];
 TF1* fpi[100];
 TF1* fk[100];
 TF1* fp[100];
 TF1* fcoin[100];
 TF1* fL[100];
 TF1* fS[100];
 TF1* fmm[100];
 TF1* fmmbg[100];



//---------------------------------------//
//               Branch                  //
//---------------------------------------//

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;

	double L_tr_chi2;
	double L_tr_x, L_tr_y, L_tr_th, L_tr_ph;
	double L_tr_p;
	double L_tr_tg_th, L_tr_tg_ph;
	double L_tr_vz;
	double L_tr_vz_saved;
	double R_tr_chi2;
	double R_tr_x, R_tr_y, R_tr_th, R_tr_ph;
	double R_tr_p;
	double R_tr_tg_th, R_tr_tg_ph;
	double R_tr_vz;
	double L_mom, R_mom, B_mom; 
	double L_ene, R_ene, B_ene; 
	double ac1sum, ac2sum;//NPE SUM

	tree->SetBranchStatus("*",0);
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
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");

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
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
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
  tree->Draw(">>elist" , "fabs(ct_orig)<1.0");//ctsum (does NOT dintinguish #track)
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  int ENum = elist->GetN(); 
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(int i=0;i<ENum;i++){
	tree->GetEntry(elist->GetEntry(i));

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


	


		if(fabs(ct)<1)ct_cut=true;
		else ct_cut=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;

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
		if(event_selection&&ct_cut)hmm_L_fom_best->Fill(mm);
		if(event_selection_nocut&&ct_cut)hmm_L_fom_nocut->Fill(mm);



}//ENum

	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2");
	TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	TH1F* hmm_Al_fom_nocut=(TH1F*)file->Get("hmm_Al_fom_best");
	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	//TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	//TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
	TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut");
	TH1F* hmm_Albg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_forAl");
	//TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_nocut_result");
    int fitmin = hmm_L_fom_best->FindBin(0.10);
    int fitmax = hmm_L_fom_best->FindBin(0.15);
    double num1 = hmm_L_fom_best->Integral(fitmin,fitmax);
    double num2 = hmm_bg_fom_best->Integral(fitmin,fitmax);
    double mixscale = num1/num2;
	cout<<"hmm_L integral ="<<num1<<endl;
	cout<<"hmm_bg integral ="<<num2<<endl;
	cout<<"mixscale(mixed/original)="<<1/mixscale<<endl;
	hmm_bg_fom_best->Sumw2();
	hmm_bg_fom_nocut->Sumw2();
	hmm_Albg_fom_nocut->Sumw2();
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_nocut->Scale(1./nbunch);
	hmm_Albg_fom_nocut->Scale(1./nbunch);
	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
	hmm_wo_bg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
	hmm_wo_bg_fom_nocut->Add(hmm_L_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
	//hmm_pi_wobg_fom_nocut->Add(hmm_pi_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	hmm_pi_wobg_fom_nocut->Add(hmm_Al_fom_nocut,hmm_Albg_fom_nocut,1.0,-1.0);


	 double constL=0.;
	 double meanL=0.;
	 double sigL=0.;
	 double constS=0.;
	 double meanS=0.;
	 double sigS=0.;
	 double par1 = 0.;
	 double par2 = 0.;
	 double par3 = 0.;
	 double par4 = 0.;
	 double par5 = 0.;
	 double par6 = 0.;
	 double integralL_nocut = 0.;
	 double integralS_nocut = 0.;
	 n_L_nocut=0.;
	 n_S_nocut=0.;
	 double fmin = -0.05;
	 double fmax =  0.15;
/****************************************/
/************BEST CUT********************/
/****************************************/
cout<<"BEST CUT START"<<endl;

	fmmbg_best=new TF1("fmmbg_best","pol4",min_mm,max_mm);
	fmmbg_best->SetNpx(2000);
	 fL_best=new TF1("fL_best","gaus(0)",min_mm,max_mm);
	 fL_best->SetNpx(2000);
	 fL_best->SetParameters(def_n_L,def_mean_L,def_sig_L);
	 fL_best->SetParLimits(0,0.,100000);
	 fL_best->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fL_best->SetParLimits(2,0.,0.01);
	 fS_best=new TF1("fS_best","gaus(0)",min_mm,max_mm);
	 fS_best->SetNpx(2000);
	 fS_best->SetParameters(def_n_S,def_mean_S,def_sig_S);
	 fS_best->SetParLimits(0,0.,100000);
	 fS_best->SetParLimits(1,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fS_best->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
	 
	 fmm_best=new TF1("fmm_best","gaus(0)+gaus(3)+pol4(6)",min_mm,max_mm);
	 fmm_best->SetNpx(2000);
	 fmm_best->SetTitle("Missing Mass (best cut)");
	 fmm_best->SetParLimits(0,0.,1000000.);//positive
	 fmm_best->SetParLimits(3,0.,1000000.);//positive
		
	 hmm_wo_bg_fom_best->Fit("fL_best","N","",def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 const_L_best=fL_best->GetParameter(0);
	 mean_L_best=fL_best->GetParameter(1);
	 sig_L_best=fL_best->GetParameter(2);
	
	 hmm_wo_bg_fom_best->Fit("fS_best","N","",def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 const_S_best=fS_best->GetParameter(0);
	 mean_S_best=fS_best->GetParameter(1);
	 sig_S_best=fS_best->GetParameter(2);


	 constL=0.;
	 meanL=0.;
	 sigL=0.;
	 constS=0.;
	 meanS=0.;
	 sigS=0.;
	 par1 = 0.;
	 par2 = 0.;
	 par3 = 0.;
	 par4 = 0.;
	 par5 = 0.;
	 par6 = 0.;
	 double integralL_best = 0.;
	 double integralS_best = 0.;
	 n_L_best=0.;
	 n_S_best=0.;

/*%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%    4th Polynomial	%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%*/
	//--- w/ 4th Polynomial func.
	 cout<<"4Poly MODE START"<<endl;
	 fmm_best_4Poly=new TF1("fmm_best_4Poly",FMM_Res,min_mm,max_mm,18);
	 fmmbg_best_4Poly=new TF1("fmmbg_best_4Poly","pol4",min_mm,max_mm);
	 fmm_best_4Poly->SetNpx(200);
	 fmm_best_4Poly->SetTitle("Missing Mass (best)");
	 fmm_best_4Poly->SetParLimits(0,0.,1000.);//positive
	 fmm_best_4Poly->SetParLimits(3,0.,300.);//positive
	 fmm_best_4Poly->SetParameter(0,const_L_best*0.85);
	 fmm_best_4Poly->SetParameter(1,mean_L_best);
	 fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L*0.4,def_mean_L+def_sig_L*0.4);
	 fmm_best_4Poly->SetParameter(2,sig_L_best);
	 fmm_best_4Poly->SetParLimits(2,0.,0.01);
	 fmm_best_4Poly->SetParameter(3,const_S_best*0.85);
	 fmm_best_4Poly->SetParameter(4,mean_S_best);
	 fmm_best_4Poly->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_best_4Poly->SetParameter(5,sig_S_best);
	 fmm_best_4Poly->SetParLimits(5,0.,0.003);
	 fmm_best_4Poly->SetParameter(6,0.7);//scale
	 fmm_best_4Poly->SetParLimits(6,0.,2.);
	 fmm_best_4Poly->SetParameter(7,0.04);//att.
	 fmm_best_4Poly->SetParLimits(7,0.01,0.1);
	 fmm_best_4Poly->SetParameter(8,0.01);//sigma
	 fmm_best_4Poly->SetParLimits(8,0.,0.1);
	 fmm_best_4Poly->SetParameter(9,-0.0001);//peak pos.
	 fmm_best_4Poly->SetParLimits(9,-0.002,0.002);
	 fmm_best_4Poly->SetParameter(10,10.);
	 fmm_best_4Poly->SetParLimits(10,0.,20.);
	 fmm_best_4Poly->SetParameter(11,0.04);
	 fmm_best_4Poly->SetParLimits(11,0.01,0.1);
	 fmm_best_4Poly->SetParameter(12,0.002);
	 fmm_best_4Poly->SetParLimits(12,0.,0.01);
	 fmm_best_4Poly->SetParameter(13,-0.077);
	 fmm_best_4Poly->SetParLimits(13,-0.085,-0.065);
	 fmm_best_4Poly->SetParameter(14,0.3);
	 fmm_best_4Poly->SetParLimits(14,0.,1.0);
	 fmm_best_4Poly->SetParameter(15,0.004);
	 fmm_best_4Poly->SetParLimits(15,0.001,0.01);
	 fmm_best_4Poly->SetParameter(16,0.002);
	 fmm_best_4Poly->SetParLimits(16,0.,0.01);
	 fmm_best_4Poly->SetParameter(17,0.0);
	 fmm_best_4Poly->SetParLimits(17,-0.002,0.002);
	 hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","L","",-0.05,0.1);//Total fitting w/ 4Poly BG
	 double chisq = fmm_best_4Poly->GetChisquare();
	 double dof  = fmm_best_4Poly->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;
	 TF1 *fmm_Lambdaonly_4Poly=new TF1("fmm_Lambdaonly_4Poly",FMM_Res_Lambdaonly,min_mm,max_mm,11);
	 double p0=fmm_best_4Poly->GetParameter(0);
	 double p1=fmm_best_4Poly->GetParameter(1);
	 double p2=fmm_best_4Poly->GetParameter(2);
	 double p6=fmm_best_4Poly->GetParameter(6);
	 double p7=fmm_best_4Poly->GetParameter(7);
	 double p8=fmm_best_4Poly->GetParameter(8);
	 double p9=fmm_best_4Poly->GetParameter(9);
	 double p14=fmm_best_4Poly->GetParameter(14);
	 double p15=fmm_best_4Poly->GetParameter(15);
	 double p16=fmm_best_4Poly->GetParameter(16);
	 double p17=fmm_best_4Poly->GetParameter(17);
	 fmm_Lambdaonly_4Poly->SetParameter(0,p0);
	 fmm_Lambdaonly_4Poly->SetParameter(1,p1);
	 fmm_Lambdaonly_4Poly->SetParameter(2,p2);
	 fmm_Lambdaonly_4Poly->SetParameter(3,p6);
	 fmm_Lambdaonly_4Poly->SetParameter(4,p7);
	 fmm_Lambdaonly_4Poly->SetParameter(5,p8);
	 fmm_Lambdaonly_4Poly->SetParameter(6,p9);
	 fmm_Lambdaonly_4Poly->SetParameter(7,p14);
	 fmm_Lambdaonly_4Poly->SetParameter(8,p15);
	 fmm_Lambdaonly_4Poly->SetParameter(9,p16);
	 fmm_Lambdaonly_4Poly->SetParameter(10,p17);
	 double nLambda = fmm_Lambdaonly_4Poly->Integral(-0.05,0.10);
	 nLambda = nLambda/0.001;
	 cout<<"nLambda="<<nLambda<<endl;
//	 fmm_best_4Poly->FixParameter(6,0.);
//	 fmm_best_4Poly->FixParameter(7,0.);
//	 fmm_best_4Poly->FixParameter(8,0.);
//	 fmm_best_4Poly->FixParameter(9,0.);
//	 fmm_best_4Poly->FixParameter(10,0.);
//	 fmm_best_4Poly->SetParameter(11,40.);//Resp. func. scale
//	 fmm_best_4Poly->SetParLimits(11,0.,100.);
//	 fmm_best_4Poly->SetParameter(12,2.);//Resp. func. att.
//	 fmm_best_4Poly->SetParameter(13,def_mean_L+def_sig_L);//Resp. func. peak
//	 fmm_best_4Poly->SetParameter(14,10.);//Resp. func. scale 
//	 fmm_best_4Poly->SetParLimits(14,0.,100.);
//	 //fmm_best_4Poly->SetParameter(15,2.);//Resp. func. att.
//	 fmm_best_4Poly->SetParameter(15,def_mean_S+def_sig_S);//Resp. func. peak 
//	 fmm_best_4Poly->SetParameter(16,def_sig_L);//Resp. func. peak 
////	 fmm_best_4Poly->SetParameter(17,def_sig_S);//Resp. func. peak 
//	// fmm_best_4Poly->SetParameter(17,16.);//Resp. func. att. 
//	// fmm_best_4Poly->SetParameter(18,def_mean_L);//Resp. func. peak 

//	 fmm_best_4Poly=new TF1("fmm_best_4Poly",FMM_4Poly,min_mm,max_mm,11);
//	 fmmbg_best_4Poly=new TF1("fmmbg_best_4Poly","pol4",min_mm,max_mm);
//	 fmm_best_4Poly->SetNpx(2000);
//	 fmm_best_4Poly->SetTitle("Missing Mass (best)");
//	 fmm_best_4Poly->SetParLimits(0,0.,1000000.);//positive
//	 fmm_best_4Poly->SetParLimits(3,0.,1000000.);//positive
//	 fmm_best_4Poly->SetParameter(0,const_L_best);
//	 fmm_best_4Poly->SetParameter(1,mean_L_best);
//	 fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 fmm_best_4Poly->SetParameter(2,sig_L_best);
//	 fmm_best_4Poly->SetParLimits(2,0.,0.01);
//	 fmm_best_4Poly->SetParameter(3,const_S_best);
//	 fmm_best_4Poly->SetParameter(4,mean_S_best);
//	 fmm_best_4Poly->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fmm_best_4Poly->SetParameter(5,sig_S_best);
//	 fmm_best_4Poly->SetParLimits(5,0.,0.01);
	 constL=fmm_best_4Poly->GetParameter(0);
	 meanL =fmm_best_4Poly->GetParameter(1);
	 sigL  =fmm_best_4Poly->GetParameter(2);
	 constS=fmm_best_4Poly->GetParameter(3);
	 meanS =fmm_best_4Poly->GetParameter(4);
	 sigS  =fmm_best_4Poly->GetParameter(5);
	 par1  =fmm_best_4Poly->GetParameter(6);
	 par2  =fmm_best_4Poly->GetParameter(7);
	 par3  =fmm_best_4Poly->GetParameter(8);
	 par4  =fmm_best_4Poly->GetParameter(9);
	 par5  =fmm_best_4Poly->GetParameter(10);
	 fmmbg_best_4Poly->SetParameters(par1, par2, par3, par4, par5);
	
	 n_L_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_L-range_L),hmm_wo_bg_fom_best->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_best<<endl;
	 integralL_best=fmmbg_best_4Poly->Integral(center_L-range_L,center_L+range_L);
	 integralL_best=integralL_best/(2*range_L/(hmm_wo_bg_fom_best->FindBin(center_L+range_L)-hmm_wo_bg_fom_best->FindBin(center_L-range_L)));
	cout<<"integralL_best="<<integralL_best<<endl;
	 if(integralL_best>0)n_L_best=n_L_best-integralL_best;
	cout<<"after(L):: "<<n_L_best<<endl;
	 n_S_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_S-range_S),hmm_wo_bg_fom_best->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_best<<endl;
	 integralS_best=fmmbg_best_4Poly->Integral(center_S-range_S,center_S+range_S);
	 integralS_best=integralS_best/(2*range_S/(hmm_wo_bg_fom_best->FindBin(center_S+range_S)-hmm_wo_bg_fom_best->FindBin(center_S-range_S)));
	cout<<"integralS_best="<<integralS_best<<endl;
	 if(integralS_best>0)n_S_best=n_S_best-integralS_best;
	cout<<"after(S):: "<<n_S_best<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_best=0.; integralS_best=0.;
	 n_L_best=0.; n_S_best=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

	
/*%%%%%%%%%%%%%%%%*/
/*%%    Voigt	%%*/
/*%%%%%%%%%%%%%%%%*/
	//--- w/ Voigt func.
	 cout<<"Voigt MODE START"<<endl;
	 fmmbg_best_Voigt=new TF1("fmmbg_best_Voigt",F_Voigt,min_mm,max_mm,4);
	 fmm_pi_best_Voigt=new TF1("fmm_pi_best_Voigt",F_Voigt,min_mm,max_mm,4);
	 fmm_pi_best_Voigt->SetNpx(2000);
	 fmm_pi_best_Voigt->SetParameters(30,0.05,0.05,0.00006);
	 fmm_pi_best_Voigt->SetParLimits(0,0.,1000000);//positive
	 hmm_pi_wobg_fom_best->Fit("fmm_pi_best_Voigt","N","",fmin,fmax);//1st Fit
	 double fmmpibest1=fmm_pi_best_Voigt->GetParameter(1);
	 double fmmpibest2=fmm_pi_best_Voigt->GetParameter(2);
	 double fmmpibest3=fmm_pi_best_Voigt->GetParameter(3);
	 fmm_pi_best_Voigt->FixParameter(1,fmmpibest1);
	 fmm_pi_best_Voigt->FixParameter(2,fmmpibest2);
	 fmm_pi_best_Voigt->FixParameter(3,fmmpibest3);
	 fmm_pi_best_Voigt->SetParLimits(0,0.,1000000);//positive
	 //hmm_wo_bg_fom_best->Fit("fmm_pi_best_Voigt","N","",-0.1,-0.02);//2nd Fit
	 fmm_best_Voigt=new TF1("fmm_best_Voigt",FMM_Voigt,min_mm,max_mm,10);
	 fmm_best_Voigt->SetNpx(2000);
	 fmm_best_Voigt->SetTitle("Missing Mass (best)");
	 fmm_best_Voigt->SetParLimits(0,0.,1000000.);//positive
	 fmm_best_Voigt->SetParLimits(3,0.,1000000.);//positive
	 fmm_best_Voigt->SetParLimits(6,0.,1000000.);//positive
	 fmm_best_Voigt->SetParameter(0,const_L_best);
	 fmm_best_Voigt->SetParameter(1,mean_L_best);
	 fmm_best_Voigt->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_best_Voigt->SetParameter(2,sig_L_best);
	 fmm_best_Voigt->SetParLimits(2,0.,0.01);
	 fmm_best_Voigt->SetParameter(3,const_S_best);
	 fmm_best_Voigt->SetParameter(4,mean_S_best);
	 fmm_best_Voigt->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_best_Voigt->SetParameter(5,sig_S_best);
	 fmm_best_Voigt->SetParLimits(5,0.,0.01);
	 fmm_best_Voigt->SetParameter(6,(fmm_pi_best_Voigt->GetParameter(0))/40.);
	 fmm_best_Voigt->FixParameter(7,fmmpibest1);
	 fmm_best_Voigt->FixParameter(8,fmmpibest2);
	 fmm_best_Voigt->FixParameter(9,fmmpibest3);
	 hmm_wo_bg_fom_best->Fit("fmm_best_Voigt","N","",-0.02,0.08);//Total fitting w/ Voigt BG
	 constL=fmm_best_Voigt->GetParameter(0);
	 meanL =fmm_best_Voigt->GetParameter(1);
	 sigL  =fmm_best_Voigt->GetParameter(2);
	 constS=fmm_best_Voigt->GetParameter(3);
	 meanS =fmm_best_Voigt->GetParameter(4);
	 sigS  =fmm_best_Voigt->GetParameter(5);
	 par1  =fmm_best_Voigt->GetParameter(6);
	 par2  =fmm_best_Voigt->GetParameter(7);
	 par3  =fmm_best_Voigt->GetParameter(8);
	 par4  =fmm_best_Voigt->GetParameter(9);
	 fmmbg_best_Voigt->SetParameters(par1, par2, par3, par4);
	
	 n_L_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_L-range_L),hmm_wo_bg_fom_best->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_best<<endl;
	 integralL_best=fmmbg_best_Voigt->Integral(center_L-range_L,center_L+range_L);
	 integralL_best=integralL_best/(2*range_L/(hmm_wo_bg_fom_best->FindBin(center_L+range_L)-hmm_wo_bg_fom_best->FindBin(center_L-range_L)));
	cout<<"integralL_best="<<integralL_best<<endl;
	 if(integralL_best>0)n_L_best=n_L_best-integralL_best;
	cout<<"after(L):: "<<n_L_best<<endl;
	 n_S_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_S-range_S),hmm_wo_bg_fom_best->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_best<<endl;
	 integralS_best=fmmbg_best_Voigt->Integral(center_S-range_S,center_S+range_S);
	 integralS_best=integralS_best/(2*range_S/(hmm_wo_bg_fom_best->FindBin(center_S+range_S)-hmm_wo_bg_fom_best->FindBin(center_S-range_S)));
	cout<<"integralS_best="<<integralS_best<<endl;
	 if(integralS_best>0)n_S_best=n_S_best-integralS_best;
	cout<<"after(S):: "<<n_S_best<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_best=0.; integralS_best=0.;
	 n_L_best=0.; n_S_best=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

/*%%%%%%%%%%%%%%%%%%%%%%*/
/*%%    2 Gaussian    %%*/
/*%%%%%%%%%%%%%%%%%%%%%%*/
	//--- w/ 2Gauss func.
	 cout<<"2Gauss MODE START"<<endl;
	 fmm_pi_best_2Gauss=new TF1("fmm_pi_best_2Gauss",F_2Gauss,min_mm,max_mm,6);
	 fmmbg_best_2Gauss=new TF1("fmmbg_best_2Gauss",F_2Gauss,min_mm,max_mm,6);
	 fmm_pi_best_2Gauss->SetNpx(2000);
	 fmm_pi_best_2Gauss->SetParameters(50,0.05,0.05,0.5,0.1,0.05);
	 hmm_pi_wobg_fom_best->Fit("fmm_pi_best_2Gauss","N","",fmin,fmax);//1st Fit
	 double fmmpibest2Gauss1=fmm_pi_best_2Gauss->GetParameter(1);
	 double fmmpibest2Gauss2=fmm_pi_best_2Gauss->GetParameter(2);
	 double fmmpibest2Gauss3=fmm_pi_best_2Gauss->GetParameter(3);
	 double fmmpibest2Gauss4=fmm_pi_best_2Gauss->GetParameter(4);
	 double fmmpibest2Gauss5=fmm_pi_best_2Gauss->GetParameter(5);
	 fmm_pi_best_2Gauss->FixParameter(1,fmmpibest2Gauss1);
	 fmm_pi_best_2Gauss->FixParameter(2,fmmpibest2Gauss2);
	 fmm_pi_best_2Gauss->FixParameter(4,fmmpibest2Gauss4);
	 fmm_pi_best_2Gauss->FixParameter(5,fmmpibest2Gauss5);
	 fmm_pi_best_2Gauss->SetParLimits(0,0.,1000000);//positive
	 fmm_pi_best_2Gauss->SetParLimits(3,0.,1000000);//positive
	 //hmm_wo_bg_fom_best->Fit("fmm_pi_best_2Gauss","N","",-0.1,-0.02);//2nd Fit
	 fmm_best_2Gauss=new TF1("fmm_best_2Gauss",FMM_2Gauss,min_mm,max_mm,12);
	 fmm_best_2Gauss->SetNpx(2000);
	 fmm_best_2Gauss->SetTitle("Missing Mass (best)");
	 fmm_best_2Gauss->SetParLimits(0,0.,1000000.);//positive
	 fmm_best_2Gauss->SetParLimits(3,0.,1000000.);//positive
	 fmm_best_2Gauss->SetParLimits(6,0.,1000000.);//positive
	 fmm_best_2Gauss->SetParLimits(9,0.,1000000.);//positive
	 fmm_best_2Gauss->SetParameter(0,const_L_best);
	 fmm_best_2Gauss->SetParameter(1,mean_L_best);
	 fmm_best_2Gauss->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_best_2Gauss->SetParameter(2,sig_L_best);
	 fmm_best_2Gauss->SetParLimits(2,0.,0.01);
	 fmm_best_2Gauss->SetParameter(3,const_S_best);
	 fmm_best_2Gauss->SetParameter(4,mean_S_best);
	 fmm_best_2Gauss->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_best_2Gauss->SetParameter(5,sig_S_best);
	 fmm_best_2Gauss->SetParLimits(5,0.,0.01);
	 fmm_best_2Gauss->SetParameter(6,1.);
	 fmm_best_2Gauss->FixParameter(7,fmmpibest2Gauss1);
	 fmm_best_2Gauss->FixParameter(8,fmmpibest2Gauss2);
	 fmm_best_2Gauss->FixParameter(9,fmmpibest2Gauss3);
	 fmm_best_2Gauss->FixParameter(10,fmmpibest2Gauss4);
	 fmm_best_2Gauss->FixParameter(11,fmmpibest2Gauss5);
	 hmm_wo_bg_fom_best->Fit("fmm_best_2Gauss","N","",-0.02,0.08);//Total fitting w/ 2Gauss BG
	 constL=fmm_best_2Gauss->GetParameter(0);
	 meanL =fmm_best_2Gauss->GetParameter(1);
	 sigL  =fmm_best_2Gauss->GetParameter(2);
	 constS=fmm_best_2Gauss->GetParameter(3);
	 meanS =fmm_best_2Gauss->GetParameter(4);
	 sigS  =fmm_best_2Gauss->GetParameter(5);
	 par1  =fmm_best_2Gauss->GetParameter(6);
	 par2  =fmm_best_2Gauss->GetParameter(7);
	 par3  =fmm_best_2Gauss->GetParameter(8);
	 par4  =fmm_best_2Gauss->GetParameter(9);
	 par5  =fmm_best_2Gauss->GetParameter(10);
	 par6  =fmm_best_2Gauss->GetParameter(11);
	 fmmbg_best_2Gauss->SetParameters(par1, par2, par3, par4, par5, par6);
	
	 n_L_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_L-range_L),hmm_wo_bg_fom_best->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_best<<endl;
	 integralL_best=fmmbg_best_2Gauss->Integral(center_L-range_L,center_L+range_L);
	 integralL_best=integralL_best/(2*range_L/(hmm_wo_bg_fom_best->FindBin(center_L+range_L)-hmm_wo_bg_fom_best->FindBin(center_L-range_L)));
	cout<<"integralL_best="<<integralL_best<<endl;
	 if(integralL_best>0)n_L_best=n_L_best-integralL_best;
	cout<<"after(L):: "<<n_L_best<<endl;
	 n_S_best=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_S-range_S),hmm_wo_bg_fom_best->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_best<<endl;
	 integralS_best=fmmbg_best_2Gauss->Integral(center_S-range_S,center_S+range_S);
	 integralS_best=integralS_best/(2*range_S/(hmm_wo_bg_fom_best->FindBin(center_S+range_S)-hmm_wo_bg_fom_best->FindBin(center_S-range_S)));
	cout<<"integralS_best="<<integralS_best<<endl;
	 if(integralS_best>0)n_S_best=n_S_best-integralS_best;
	cout<<"after(S):: "<<n_S_best<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_best=0.; integralS_best=0.;
	 n_L_best=0.; n_S_best=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

/****************************************/
/************NO CUT**********************/
/****************************************/
cout<<"NO CUT START"<<endl;
	fmmbg_nocut=new TF1("fmmbg_nocut","pol4",min_mm,max_mm);
	fmmbg_nocut->SetNpx(2000);
	 fL_nocut=new TF1("fL_nocut","gaus(0)",min_mm,max_mm);
	 fL_nocut->SetNpx(2000);
	 fL_nocut->SetParameters(def_n_L,def_mean_L,def_sig_L);
	 fL_nocut->SetParLimits(0,0.,100000);
	 fL_nocut->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fL_nocut->SetParLimits(2,0.,0.01);
	 fS_nocut=new TF1("fS_nocut","gaus(0)",min_mm,max_mm);
	 fS_nocut->SetNpx(2000);
	 fS_nocut->SetParameters(def_n_S,def_mean_S,def_sig_S);
	 fS_nocut->SetParLimits(0,0.,100000);
	 fS_nocut->SetParLimits(1,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fS_nocut->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
	 
	 fmm_nocut=new TF1("fmm_nocut","gaus(0)+gaus(3)+pol4(6)",min_mm,max_mm);
	 fmm_nocut->SetNpx(2000);
	 fmm_nocut->SetTitle("Missing Mass (nocut cut)");
	 fmm_nocut->SetParLimits(0,0.,1000000.);//positive
	 fmm_nocut->SetParLimits(3,0.,1000000.);//positive
		
	 hmm_wo_bg_fom_nocut->Fit("fL_nocut","N","",def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 const_L_nocut=fL_nocut->GetParameter(0);
	 mean_L_nocut=fL_nocut->GetParameter(1);
	 sig_L_nocut=fL_nocut->GetParameter(2);
	 center_L=def_mean_L;
	 range_L=2.*def_sig_L;
	
	 hmm_wo_bg_fom_nocut->Fit("fS_nocut","N","",def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 const_S_nocut=fS_nocut->GetParameter(0);
	 mean_S_nocut=fS_nocut->GetParameter(1);
	 sig_S_nocut=fS_nocut->GetParameter(2);
	 center_S=def_mean_S;
	 range_S=2.*def_sig_S;


/*%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%    4th Polynomial	%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%*/
	//--- w/ 4th Polynomial func.
	 cout<<"4Poly MODE START"<<endl;
	 fmm_nocut_4Poly=new TF1("fmm_nocut_4Poly",FMM_4Poly,min_mm,max_mm,11);
	 fmmbg_nocut_4Poly=new TF1("fmmbg_nocut_4Poly","pol4",min_mm,max_mm);
	 fmm_nocut_4Poly->SetNpx(2000);
	 fmm_nocut_4Poly->SetTitle("Missing Mass (nocut)");
	 fmm_nocut_4Poly->SetParLimits(0,0.,1000000.);//positive
	 fmm_nocut_4Poly->SetParLimits(3,0.,1000000.);//positive
	 fmm_nocut_4Poly->SetParameter(0,const_L_nocut);
	 fmm_nocut_4Poly->SetParameter(1,mean_L_nocut);
	 fmm_nocut_4Poly->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_nocut_4Poly->SetParameter(2,sig_L_nocut);
	 fmm_nocut_4Poly->SetParLimits(2,0.,0.01);
	 fmm_nocut_4Poly->SetParameter(3,const_S_nocut);
	 fmm_nocut_4Poly->SetParameter(4,mean_S_nocut);
	 fmm_nocut_4Poly->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_nocut_4Poly->SetParameter(5,sig_S_nocut);
	 fmm_nocut_4Poly->SetParLimits(5,0.,0.01);
	 hmm_wo_bg_fom_nocut->Fit("fmm_nocut_4Poly","N","",fmin,fmax);//Total fitting w/ 4Poly BG
	 constL=fmm_nocut_4Poly->GetParameter(0);
	 meanL =fmm_nocut_4Poly->GetParameter(1);
	 sigL  =fmm_nocut_4Poly->GetParameter(2);
	 constS=fmm_nocut_4Poly->GetParameter(3);
	 meanS =fmm_nocut_4Poly->GetParameter(4);
	 sigS  =fmm_nocut_4Poly->GetParameter(5);
	 par1  =fmm_nocut_4Poly->GetParameter(6);
	 par2  =fmm_nocut_4Poly->GetParameter(7);
	 par3  =fmm_nocut_4Poly->GetParameter(8);
	 par4  =fmm_nocut_4Poly->GetParameter(9);
	 par5  =fmm_nocut_4Poly->GetParameter(10);
	 fmmbg_nocut_4Poly->SetParameters(par1, par2, par3, par4, par5);
	
	 n_L_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_L-range_L),hmm_wo_bg_fom_nocut->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_nocut<<endl;
	 integralL_nocut=fmmbg_nocut_4Poly->Integral(center_L-range_L,center_L+range_L);
	 integralL_nocut=integralL_nocut/(2*range_L/(hmm_wo_bg_fom_nocut->FindBin(center_L+range_L)-hmm_wo_bg_fom_nocut->FindBin(center_L-range_L)));
	cout<<"integralL_nocut="<<integralL_nocut<<endl;
	 if(integralL_nocut>0)n_L_nocut=n_L_nocut-integralL_nocut;
	cout<<"after(L):: "<<n_L_nocut<<endl;
	 n_S_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_S-range_S),hmm_wo_bg_fom_nocut->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_nocut<<endl;
	 integralS_nocut=fmmbg_nocut_4Poly->Integral(center_S-range_S,center_S+range_S);
	 integralS_nocut=integralS_nocut/(2*range_S/(hmm_wo_bg_fom_nocut->FindBin(center_S+range_S)-hmm_wo_bg_fom_nocut->FindBin(center_S-range_S)));
	cout<<"integralS_nocut="<<integralS_nocut<<endl;
	 if(integralS_nocut>0)n_S_nocut=n_S_nocut-integralS_nocut;
	cout<<"after(S):: "<<n_S_nocut<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_nocut=0.; integralS_nocut=0.;
	 n_L_nocut=0.; n_S_nocut=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

	
/*%%%%%%%%%%%%%%%%*/
/*%%    Voigt	%%*/
/*%%%%%%%%%%%%%%%%*/
	//--- w/ Voigt func.
	 cout<<"Voigt MODE START"<<endl;
	 fmmbg_nocut_Voigt=new TF1("fmmbg_nocut_Voigt",F_Voigt,min_mm,max_mm,4);
	 fmm_pi_nocut_Voigt=new TF1("fmm_pi_nocut_Voigt",F_Voigt,min_mm,max_mm,4);
	 fmm_pi_nocut_Voigt->SetNpx(2000);
	 fmm_pi_nocut_Voigt->SetParameters(3.7,0.051,0.045,0.00001);
	 fmm_pi_nocut_Voigt->SetParLimits(0,0.,1000000);//positive
	 hmm_pi_wobg_fom_nocut->Fit("fmm_pi_nocut_Voigt","N","",-0.04,0.10);//1st Fit
	 double fmmpinocut1=fmm_pi_nocut_Voigt->GetParameter(1);
	 double fmmpinocut2=fmm_pi_nocut_Voigt->GetParameter(2);
	 double fmmpinocut3=fmm_pi_nocut_Voigt->GetParameter(3);
	 fmm_pi_nocut_Voigt->FixParameter(1,fmmpinocut1);
	 fmm_pi_nocut_Voigt->FixParameter(2,fmmpinocut2);
	 fmm_pi_nocut_Voigt->FixParameter(3,fmmpinocut3);
	 fmm_pi_nocut_Voigt->SetParLimits(0,0.,1000000);//positive
	 //hmm_wo_bg_fom_nocut->Fit("fmm_pi_nocut_Voigt","N","",-0.1,-0.02);//2nd Fit
	 fmm_nocut_Voigt=new TF1("fmm_nocut_Voigt",FMM_Voigt,min_mm,max_mm,10);
	 fmm_nocut_Voigt->SetNpx(2000);
	 fmm_nocut_Voigt->SetTitle("Missing Mass (nocut)");
	 fmm_nocut_Voigt->SetParLimits(0,0.,1000000.);//positive
	 fmm_nocut_Voigt->SetParLimits(3,0.,1000000.);//positive
	 fmm_nocut_Voigt->SetParLimits(6,0.,1000000.);//positive
	 fmm_nocut_Voigt->SetParameter(0,const_L_nocut);
	 fmm_nocut_Voigt->SetParameter(1,mean_L_nocut);
	 fmm_nocut_Voigt->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_nocut_Voigt->SetParameter(2,sig_L_nocut);
	 fmm_nocut_Voigt->SetParLimits(2,0.,0.01);
	 fmm_nocut_Voigt->SetParameter(3,const_S_nocut);
	 fmm_nocut_Voigt->SetParameter(4,mean_S_nocut);
	 fmm_nocut_Voigt->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_nocut_Voigt->SetParameter(5,sig_S_nocut);
	 fmm_nocut_Voigt->SetParLimits(5,0.,0.01);
	 fmm_nocut_Voigt->SetParameter(6,(fmm_pi_nocut_Voigt->GetParameter(0))/30.);
	 fmm_nocut_Voigt->FixParameter(7,fmmpinocut1);
	 fmm_nocut_Voigt->FixParameter(8,fmmpinocut2);
	 fmm_nocut_Voigt->FixParameter(9,fmmpinocut3);
	 hmm_wo_bg_fom_nocut->Fit("fmm_nocut_Voigt","N","",fmin,fmax);//Total fitting w/ Voigt BG
	 constL=fmm_nocut_Voigt->GetParameter(0);
	 meanL =fmm_nocut_Voigt->GetParameter(1);
	 sigL  =fmm_nocut_Voigt->GetParameter(2);
	 constS=fmm_nocut_Voigt->GetParameter(3);
	 meanS =fmm_nocut_Voigt->GetParameter(4);
	 sigS  =fmm_nocut_Voigt->GetParameter(5);
	 par1  =fmm_nocut_Voigt->GetParameter(6);
	 par2  =fmm_nocut_Voigt->GetParameter(7);
	 par3  =fmm_nocut_Voigt->GetParameter(8);
	 par4  =fmm_nocut_Voigt->GetParameter(9);
	 fmmbg_nocut_Voigt->SetParameters(par1, par2, par3, par4);
	
	 n_L_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_L-range_L),hmm_wo_bg_fom_nocut->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_nocut<<endl;
	 integralL_nocut=fmmbg_nocut_Voigt->Integral(center_L-range_L,center_L+range_L);
	 integralL_nocut=integralL_nocut/(2*range_L/(hmm_wo_bg_fom_nocut->FindBin(center_L+range_L)-hmm_wo_bg_fom_nocut->FindBin(center_L-range_L)));
	cout<<"integralL_nocut="<<integralL_nocut<<endl;
	 if(integralL_nocut>0)n_L_nocut=n_L_nocut-integralL_nocut;
	cout<<"after(L):: "<<n_L_nocut<<endl;
	 n_S_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_S-range_S),hmm_wo_bg_fom_nocut->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_nocut<<endl;
	 integralS_nocut=fmmbg_nocut_Voigt->Integral(center_S-range_S,center_S+range_S);
	 integralS_nocut=integralS_nocut/(2*range_S/(hmm_wo_bg_fom_nocut->FindBin(center_S+range_S)-hmm_wo_bg_fom_nocut->FindBin(center_S-range_S)));
	cout<<"integralS_nocut="<<integralS_nocut<<endl;
	 if(integralS_nocut>0)n_S_nocut=n_S_nocut-integralS_nocut;
	cout<<"after(S):: "<<n_S_nocut<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_nocut=0.; integralS_nocut=0.;
	 n_L_nocut=0.; n_S_nocut=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

/*%%%%%%%%%%%%%%%%%%%%%%*/
/*%%    2 Gaussian    %%*/
/*%%%%%%%%%%%%%%%%%%%%%%*/
	//--- w/ 2Gauss func.
	 cout<<"2Gauss MODE START"<<endl;
	 fmm_pi_nocut_2Gauss=new TF1("fmm_pi_nocut_2Gauss",F_2Gauss,min_mm,max_mm,6);
	 fmmbg_nocut_2Gauss=new TF1("fmmbg_nocut_2Gauss",F_2Gauss,min_mm,max_mm,6);
	 fmm_pi_nocut_2Gauss->SetNpx(2000);
	 fmm_pi_nocut_2Gauss->SetParameters(30,0.05,0.04,0.1,0.08,0.04);
	 fmm_pi_nocut_2Gauss->SetParLimits(0,0.,10000000.);//positive
	 fmm_pi_nocut_2Gauss->SetParLimits(3,0.,10000000.);//positive
	 fmm_pi_nocut_2Gauss->SetParLimits(2,0.01,1.);//Not too narrow
	 fmm_pi_nocut_2Gauss->SetParLimits(5,0.01,1.);//Not too narrow
	 hmm_pi_wobg_fom_nocut->Fit("fmm_pi_nocut_2Gauss","N","",-0.02,0.08);//1st Fit
	 double fmmpinocut2Gauss1=fmm_pi_nocut_2Gauss->GetParameter(1);
	 double fmmpinocut2Gauss2=fmm_pi_nocut_2Gauss->GetParameter(2);
	 double fmmpinocut2Gauss3=fmm_pi_nocut_2Gauss->GetParameter(3);
	 double fmmpinocut2Gauss4=fmm_pi_nocut_2Gauss->GetParameter(4);
	 double fmmpinocut2Gauss5=fmm_pi_nocut_2Gauss->GetParameter(5);
	 fmm_pi_nocut_2Gauss->FixParameter(1,fmmpinocut2Gauss1);
	 fmm_pi_nocut_2Gauss->FixParameter(2,fmmpinocut2Gauss2);
	 fmm_pi_nocut_2Gauss->FixParameter(4,fmmpinocut2Gauss4);
	 fmm_pi_nocut_2Gauss->FixParameter(5,fmmpinocut2Gauss5);
	 fmm_pi_nocut_2Gauss->SetParLimits(0,0.,1000000);//positive
	 fmm_pi_nocut_2Gauss->SetParLimits(3,0.,1000000);//positive
	 //hmm_wo_bg_fom_nocut->Fit("fmm_pi_nocut_2Gauss","N","",-0.1,-0.02);//2nd Fit
	 fmm_nocut_2Gauss=new TF1("fmm_nocut_2Gauss",FMM_2Gauss,min_mm,max_mm,12);
	 fmm_nocut_2Gauss->SetNpx(2000);
	 fmm_nocut_2Gauss->SetTitle("Missing Mass (nocut)");
	 fmm_nocut_2Gauss->SetParLimits(0,0.,1000000.);//positive
	 fmm_nocut_2Gauss->SetParLimits(3,0.,1000000.);//positive
	 fmm_nocut_2Gauss->SetParLimits(6,0.,1000000.);//positive
	 fmm_nocut_2Gauss->SetParLimits(9,0.,1000000.);//positive
	 fmm_nocut_2Gauss->SetParameter(0,const_L_nocut);
	 fmm_nocut_2Gauss->SetParameter(1,mean_L_nocut);
	 fmm_nocut_2Gauss->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_nocut_2Gauss->SetParameter(2,sig_L_nocut);
	 fmm_nocut_2Gauss->SetParLimits(2,0.,0.01);
	 fmm_nocut_2Gauss->SetParameter(3,const_S_nocut);
	 fmm_nocut_2Gauss->SetParameter(4,mean_S_nocut);
	 fmm_nocut_2Gauss->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_nocut_2Gauss->SetParameter(5,sig_S_nocut);
	 fmm_nocut_2Gauss->SetParLimits(5,0.,0.01);
	 fmm_nocut_2Gauss->SetParameter(6,20.);
	 fmm_nocut_2Gauss->FixParameter(7,fmmpinocut2Gauss1);
	 fmm_nocut_2Gauss->FixParameter(8,fmmpinocut2Gauss2);
	 fmm_nocut_2Gauss->FixParameter(9,fmmpinocut2Gauss3);
	 fmm_nocut_2Gauss->FixParameter(10,fmmpinocut2Gauss4);
	 fmm_nocut_2Gauss->FixParameter(11,fmmpinocut2Gauss5);
	 hmm_wo_bg_fom_nocut->Fit("fmm_nocut_2Gauss","N","",fmin,fmax);//Total fitting w/ 2Gauss BG
	 constL=fmm_nocut_2Gauss->GetParameter(0);
	 meanL =fmm_nocut_2Gauss->GetParameter(1);
	 sigL  =fmm_nocut_2Gauss->GetParameter(2);
	 constS=fmm_nocut_2Gauss->GetParameter(3);
	 meanS =fmm_nocut_2Gauss->GetParameter(4);
	 sigS  =fmm_nocut_2Gauss->GetParameter(5);
	 par1  =fmm_nocut_2Gauss->GetParameter(6);
	 par2  =fmm_nocut_2Gauss->GetParameter(7);
	 par3  =fmm_nocut_2Gauss->GetParameter(8);
	 par4  =fmm_nocut_2Gauss->GetParameter(9);
	 par5  =fmm_nocut_2Gauss->GetParameter(10);
	 par6  =fmm_nocut_2Gauss->GetParameter(11);
	 fmmbg_nocut_2Gauss->SetParameters(par1, par2, par3, par4, par5, par6);
	
	 n_L_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_L-range_L),hmm_wo_bg_fom_nocut->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_nocut<<endl;
	 integralL_nocut=fmmbg_nocut_2Gauss->Integral(center_L-range_L,center_L+range_L);
	 integralL_nocut=integralL_nocut/(2*range_L/(hmm_wo_bg_fom_nocut->FindBin(center_L+range_L)-hmm_wo_bg_fom_nocut->FindBin(center_L-range_L)));
	cout<<"integralL_nocut="<<integralL_nocut<<endl;
	 if(integralL_nocut>0)n_L_nocut=n_L_nocut-integralL_nocut;
	cout<<"after(L):: "<<n_L_nocut<<endl;
	 n_S_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_S-range_S),hmm_wo_bg_fom_nocut->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_nocut<<endl;
	 integralS_nocut=fmmbg_nocut_2Gauss->Integral(center_S-range_S,center_S+range_S);
	 integralS_nocut=integralS_nocut/(2*range_S/(hmm_wo_bg_fom_nocut->FindBin(center_S+range_S)-hmm_wo_bg_fom_nocut->FindBin(center_S-range_S)));
	cout<<"integralS_nocut="<<integralS_nocut<<endl;
	 if(integralS_nocut>0)n_S_nocut=n_S_nocut-integralS_nocut;
	cout<<"after(S):: "<<n_S_nocut<<endl;
	 cout<<"constL"<<constL<<endl;
	cout<<"meanL"<<meanL<<endl;
	cout<<"sigL"<<sigL<<endl;
	 cout<<"constS"<<constS<<endl;
	cout<<"meanS"<<meanS<<endl;
	cout<<"sigS"<<sigS<<endl;

     par1=0.; par2=0.; par3=0.; par4=0.; par5=0.; par6=0.;
	 integralL_nocut=0.; integralS_nocut=0.;
	 n_L_nocut=0.; n_S_nocut=0.;
	 constL=0.;meanL=0.;sigL=0.;constS=0.;meanS=0.;sigS=0.;

//	 fmm_nocut->SetParameter(0,const_L_nocut);
//	 fmm_nocut->SetParameter(1,mean_L_nocut);
//	 fmm_nocut->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 fmm_nocut->SetParameter(2,sig_L_nocut);
//	 fmm_nocut->SetParLimits(2,0.,2*def_sig_L);
//	// fmm_nocut->SetParameters(9,100);
//	 fmm_nocut->SetParameter(3,const_S_nocut);
//	 fmm_nocut->SetParameter(4,mean_S_nocut);
//	 fmm_nocut->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fmm_nocut->SetParameter(5,sig_S_nocut);
//	 fmm_nocut->SetParLimits(5,0.,2*def_sig_S);
//	 hmm_wo_bg_fom_nocut->Fit("fmm_nocut","N","0",-0.05,0.15);
//	double fmm_nocutpar0 = fmm_nocut->GetParameter(0);cout<<"fmm_nocut[0]="<<fmm_nocutpar0<<endl;//area(L)
//	double fmm_nocutpar1 = fmm_nocut->GetParameter(1);cout<<"fmm_nocut[1]="<<fmm_nocutpar1<<endl;//mean(L)
//	double fmm_nocutpar2 = fmm_nocut->GetParameter(2);cout<<"fmm_nocut[2]="<<fmm_nocutpar2<<endl;//sigma(L)
//	double fmm_nocutpar3 = fmm_nocut->GetParameter(3);cout<<"fmm_nocut[3]="<<fmm_nocutpar3<<endl;//area(S)
//	double fmm_nocutpar4 = fmm_nocut->GetParameter(4);cout<<"fmm_nocut[4]="<<fmm_nocutpar4<<endl;//mean(S)
//	double fmm_nocutpar5 = fmm_nocut->GetParameter(5);cout<<"fmm_nocut[5]="<<fmm_nocutpar5<<endl;//sigma(S)
//	double fmm_nocutpar6 = fmm_nocut->GetParameter(6);cout<<"fmm_nocut[6]="<<fmm_nocutpar6<<endl;//poly_const
//	double fmm_nocutpar7 = fmm_nocut->GetParameter(7);cout<<"fmm_nocut[7]="<<fmm_nocutpar7<<endl;//poly_x
//	double fmm_nocutpar8 = fmm_nocut->GetParameter(8);cout<<"fmm_nocut[8]="<<fmm_nocutpar8<<endl;//poly_x^2
//	double fmm_nocutpar9 = fmm_nocut->GetParameter(9);cout<<"fmm_nocut[9]="<<fmm_nocutpar9<<endl;//poly_x^3
//	double fmm_nocutpar10 = fmm_nocut->GetParameter(10);cout<<"fmm_nocut[10]="<<fmm_nocutpar10<<endl;//poly_x^4
//	 fmmbg_nocut->SetParameters(fmm_nocutpar6,fmm_nocutpar7,fmm_nocutpar8,fmm_nocutpar9,fmm_nocutpar10);
//	
//	 n_L_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_L-range_L),hmm_wo_bg_fom_nocut->FindBin(center_L+range_L));
//	cout<<"before(L):: "<<n_L_nocut<<endl;
//	 integralL_nocut=fmmbg_nocut->Integral(center_L-range_L,center_L+range_L);
//	 integralL_nocut=integralL_nocut/(2*range_L/(hmm_wo_bg_fom_nocut->FindBin(center_L+range_L)-hmm_wo_bg_fom_nocut->FindBin(center_L-range_L)));
//	cout<<"integralL_nocut="<<integralL_nocut<<endl;
//	 if(integralL_nocut>0)n_L_nocut=n_L_nocut-integralL_nocut;
//	cout<<"after(L):: "<<n_L_nocut<<endl;
//	 n_S_nocut=hmm_wo_bg_fom_nocut->Integral(hmm_wo_bg_fom_nocut->FindBin(center_S-range_S),hmm_wo_bg_fom_nocut->FindBin(center_S+range_S));
//	cout<<"before(S):: "<<n_S_nocut<<endl;
//	 integralS_nocut=fmmbg_nocut->Integral(center_S-range_S,center_S+range_S);
//	 integralS_nocut=integralS_nocut/(2*range_S/(hmm_wo_bg_fom_nocut->FindBin(center_S+range_S)-hmm_wo_bg_fom_nocut->FindBin(center_S-range_S)));
//	cout<<"integralS_nocut="<<integralS_nocut<<endl;
//	 if(integralS_nocut>0)n_S_nocut-=integralS_nocut;
//	cout<<"after(S):: "<<n_S_nocut<<endl;
//	 cout<<"n_L"<<n_L_nocut<<endl;
//	cout<<"mean_L"<<mean_L_nocut<<endl;
//	cout<<"sig_L"<<sig_L_nocut<<endl;
//	 cout<<"n_S"<<n_S_nocut<<endl;
//	cout<<"mean_S"<<mean_S_nocut<<endl;
//	cout<<"sig_S"<<sig_S_nocut<<endl;





	hmm_wo_bg_fom_best->Draw("");

/****************************************/
	hmm_L_fom_best->SetLineColor(kAzure);
	hmm_L_fom_nocut->SetLineColor(kAzure);
	hmm_wo_bg_fom_best->SetLineColor(kAzure);
	hmm_wo_bg_fom_nocut->SetLineColor(kAzure);
	hmm_bg_fom_best->SetLineColor(kRed);
	hmm_bg_fom_nocut->SetLineColor(kRed);
	hmm_pi_wobg_fom_best->SetLineColor(kAzure);
	hmm_pi_wobg_fom_nocut->SetLineColor(kAzure);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TCanvas* c3 = new TCanvas("c3","BG (Mixed Event Analysis)");
	c3->Divide(2,2);
	c3->cd(1);
	hmm_L_fom_best->Draw();
	c3->cd(2);
	hmm_L_fom_best->Draw();
	hmm_bg_fom_best->SetStats(0);
	hmm_bg_fom_best->Draw("same");
	c3->cd(3);
	hmm_L_fom_nocut->Draw();
	c3->cd(4);
	hmm_L_fom_nocut->Draw();
	hmm_bg_fom_nocut->SetStats(0);
	hmm_bg_fom_nocut->Draw("same");
	
	TCanvas* c4 = new TCanvas("c4","4th Polynomial");
	c4->Divide(2,2);
	c4->cd(1);
	fmm_best_4Poly->SetLineColor(kGreen);
	hmm_wo_bg_fom_best->Draw("");
	fmm_best_4Poly->Draw("same");
	c4->cd(2);
	fmmbg_best_4Poly->SetLineColor(kGreen);
	fmmbg_best_4Poly->SetFillColor(kGreen);
	fmmbg_best_4Poly->SetFillStyle(3018);
	hmm_wo_bg_fom_best->Draw("");
	fmmbg_best_4Poly->Draw("same");
	c4->cd(3);
	fmm_nocut_4Poly->SetLineColor(kGreen);
	hmm_wo_bg_fom_nocut->Draw("");
	fmm_nocut_4Poly->Draw("same");
	c4->cd(4);
	fmmbg_nocut_4Poly->SetLineColor(kGreen);
	fmmbg_nocut_4Poly->SetFillColor(kGreen);
	fmmbg_nocut_4Poly->SetFillStyle(3018);
	hmm_wo_bg_fom_nocut->Draw("");
	fmmbg_nocut_4Poly->Draw("same");

	TCanvas* c5 = new TCanvas("c5","Voigt func.");
	c5->Divide(2,2);
	c5->cd(1);
	fmm_best_Voigt->SetLineColor(kOrange);
	hmm_wo_bg_fom_best->Draw("");
	fmm_best_Voigt->Draw("same");
	c5->cd(2);
	fmmbg_best_Voigt->SetLineColor(kOrange);
	fmmbg_best_Voigt->SetFillColor(kOrange);
	fmmbg_best_Voigt->SetFillStyle(3018);
	hmm_wo_bg_fom_best->Draw("");
	fmmbg_best_Voigt->Draw("same");
	c5->cd(3);
	fmm_nocut_Voigt->SetLineColor(kOrange);
	hmm_wo_bg_fom_nocut->Draw("");
	fmm_nocut_Voigt->Draw("same");
	c5->cd(4);
	fmmbg_nocut_Voigt->SetLineColor(kOrange);
	fmmbg_nocut_Voigt->SetFillColor(kOrange);
	fmmbg_nocut_Voigt->SetFillStyle(3018);
	hmm_wo_bg_fom_nocut->Draw("");
	fmmbg_nocut_Voigt->Draw("same");

	TCanvas* c6 = new TCanvas("c6","Double Gaussian");
	c6->Divide(2,2);
	c6->cd(1);
	fmm_best_2Gauss->SetLineColor(kRed);
	hmm_wo_bg_fom_best->Draw("");
	fmm_best_2Gauss->Draw("same");
	c6->cd(2);
	fmmbg_best_2Gauss->SetLineColor(kRed);
	fmmbg_best_2Gauss->SetFillColor(kRed);
	fmmbg_best_2Gauss->SetFillStyle(3018);
	hmm_wo_bg_fom_best->Draw("");
	fmmbg_best_2Gauss->Draw("same");
	c6->cd(3);
	fmm_nocut_2Gauss->SetLineColor(kRed);
	hmm_wo_bg_fom_nocut->Draw("");
	fmm_nocut_2Gauss->Draw("same");
	c6->cd(4);
	fmmbg_nocut_2Gauss->SetLineColor(kRed);
	fmmbg_nocut_2Gauss->SetFillColor(kRed);
	fmmbg_nocut_2Gauss->SetFillStyle(3018);
	hmm_wo_bg_fom_nocut->Draw("");
	fmmbg_nocut_2Gauss->Draw("same");
//		TCanvas* c4 = new TCanvas("c4","c4");
//	hmm_wo_bg_fom_best->Draw("");
//	fmmbg_best->Draw("same");
//	fL_best->Draw("same");
//	fS_best->Draw("same");
//		TCanvas* c5 = new TCanvas("c5","c5");
//	hmm_wo_bg_fom_best->Draw("");
//	fmm_best->Draw("same");
//		TCanvas* c6 = new TCanvas("c6","c6");
//	//hmm_pi_wobg_fom_best->Draw("");
//	hmm_wo_bg_fom_best->Draw("");
//	//fmm_pi_best->SetLineColor(kOrange);
//	//fmm_pi_best->Draw("same");
		TCanvas* c7 = new TCanvas("c7","Pion (best)");
	c7->Divide(2,2);
	c7->cd(1);
	hmm_pi_wobg_fom_best->Draw("");
	//fmm_pi_best_4Poly->SetLineColor(kGreen);
	//fmm_pi_best_4Poly->Draw("same");
	c7->cd(2);
	hmm_pi_wobg_fom_best->Draw("");
	fmm_pi_best_Voigt->SetLineColor(kOrange);
	fmm_pi_best_Voigt->Draw("same");
	c7->cd(3);
	hmm_pi_wobg_fom_best->Draw("");
	fmm_pi_best_2Gauss->SetLineColor(kRed);
	fmm_pi_best_2Gauss->Draw("same");
		TCanvas* c8 = new TCanvas("c8","Pion (nocut)");
	c8->Divide(2,2);
	c8->cd(1);
	hmm_pi_wobg_fom_nocut->Draw("");
	//fmm_pi_nocut_4Poly->SetLineColor(kGreen);
	//fmm_pi_nocut_4Poly->Draw("same");
	c8->cd(2);
	hmm_pi_wobg_fom_nocut->Draw("");
	fmm_pi_nocut_Voigt->SetLineColor(kOrange);
	fmm_pi_nocut_Voigt->Draw("same");
	c8->cd(3);
	hmm_pi_wobg_fom_nocut->Draw("");
	fmm_pi_nocut_2Gauss->SetLineColor(kRed);
	fmm_pi_nocut_2Gauss->Draw("same");
	c8->cd(4);
	//hmm_pi_fom_nocut->Draw("");
	//hmm_pi_fom_best->Draw("same");
		TCanvas* c9 = new TCanvas("c9","Total(nocut)");
	hmm_wo_bg_fom_nocut->Draw("");
	fmm_nocut_4Poly->Draw("same");
	fmm_nocut_Voigt->Draw("same");
	fmm_nocut_2Gauss->Draw("same");
		TCanvas* c10 = new TCanvas("c10","Total(best)");
	hmm_wo_bg_fom_best->Draw("");
	fmm_best_4Poly->Draw("same");
	fmm_best_Voigt->Draw("same");
	fmm_best_2Gauss->Draw("same");
		TCanvas* c11 = new TCanvas("c11","Lambda only");
	fmm_Lambdaonly_4Poly->Draw();


/*--- Print ---*/
cout << "Print is starting" << endl;
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


cout << "Well done!" << endl;
}//fit
