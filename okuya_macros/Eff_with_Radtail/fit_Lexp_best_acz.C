//--------------------------------//
//--  Fitting w/ Response func. --//
//--------------------------------//
//
//K. Okuyama (Sep. 14, 2020)
//
//This is taken over from fit_landau.C
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
  val = par[num] * step * sum * invsq2pi / (par[num+2]*par[num+1]*exp(par[num+3]/par[num+1]));
  return val;
}

double landaugaus2(double *x, double *par, int num) {
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val;

// MP shift correction
  mpc = par[num+1] - mpshift * par[num];
// Range of convolution integral
  xlow = x[0] - sc * par[num+3];
  xupp = x[0] + sc * par[num+3];
  step = (xupp-xlow) / np;
// Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[num]) / par[num];
     sum += fland * TMath::Gaus(x[0],xx,par[num+3]);

     xx = xupp - (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[num]) / par[num];
     sum += fland * TMath::Gaus(x[0],xx,par[num+3]);
  }
  val = par[num+2] * step * sum * invsq2pi / par[num+3];

  return val;
}


double FMM_Lambda_Sigma( double *x, double *par , int num)
  {
  double val = par[num] * TMath::Gaus(x[0],par[num+1],par[num+2]);//Lambda Gaussian
  val += par[num+3] * TMath::Gaus(x[0],par[num+4],par[num+5]);//Sigma Gaussian
  return val;
}

double FMM_Response( double *x, double *par ){

   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val1, val2;

// MP shift correction
  mpc = par[1] - mpshift * par[0];
// Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
// Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);

     xx = xupp - (i-.5) * step;
     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
     sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  val1 = step * sum * invsq2pi / par[3];

/*------Landau * Gauss convluted------*/

// Range of convolution integral
  sum  = 0.;
  xlow = 0.;
  xupp = x[0] + 1.6 * sc * par[3];
  step = (xupp-xlow) / np;
  if(step<0.)step = 0.;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[5];
     fland = TMath::Gaus(xx,x[0],par[3]);
     sum += fland * TMath::Exp(-xx/par[4]);
     xx = xupp - (i-.5) * step - par[5];
     fland = TMath::Gaus(xx,x[0],par[3]);
     sum += fland * TMath::Exp(-xx/par[4]);
  }
  //val = par[2] * step * sum * invsq2pi / par[3];
  val2 =  step * sum * invsq2pi / (par[3]*par[4]*exp(par[5]/par[4]));
  //val2 =  step * sum * invsq2pi / (par[3]*par[4]);
/*------Exp * Gauss convluted------*/

  return par[2]*(val1+par[6]*val2);//N x (Landau*Gauss) + N' x (Exp*Gauss)

}

double FMM_Res( double *x, double *par ){

	return FMM_Response(x,par)+FMM_Response(x,&par[7]);

}

double FMM_Res_nocut( double *x, double *par ){

	return FMM_Response(x,par)+FMM_Response(x,&par[7])+par[21]*par[14] * TMath::Voigt(x[0]-par[15],par[16],par[17],4)+(1.-par[21])*par[14] * TMath::Voigt(x[0]-par[18],par[19],par[20],4);

}

double FMM_2BG( double *x, double *par ){

	return par[7]*par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4)+(1.-par[7])*par[0] * TMath::Voigt(x[0]-par[4],par[5],par[6],4);

}

void fit_Lexp_best_acz(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("../h2all.root","read");//input file of all H2 run(default: h2all4.root)
	//ACCBGの引き算はmea_hist.ccから
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_llccrr_new_new.root","read");//input file of BG(MEA) histo.(default: bgmea6.root)
  double nbunch = 6000.;//effetive bunches (6 bunches x 1000 mixtures)
  
//Systematic study
//TFile *file_mea = new TFile("./temp/bgmea_rrr.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
//double nbunch = 3000.;//effetive bunches (3 bunches x 1000 mixtures)


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
 int bin_mm=(max_mm-min_mm)/0.001; //Counts/2 MeV
 bin_mm=(int)bin_mm;
 //const double fit_min_mm=-0.006;
 const double fit_min_mm=-0.02;
 const double fit_max_mm=0.12;
 const int fit_bin_mm = (fit_max_mm-fit_min_mm)/0.001;
 const double fit_bin_width = (fit_max_mm-fit_min_mm)/fit_bin_mm;


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
 TF1* fL_noZ;
 TF1* fS_noZ;
 TF1* fmm_noZ_Lexp;
 TF1* fmm_noAC_Lexp;
 TF1* fmm_strict_Lexp;
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
 TF1* fmm_best_Lexp;
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
  TH1F* hmm_L_fom_strict  = new TH1F("hmm_L_fom_strict","hmm_L_fom_strict",xbin,xmin,xmax);
  TH1F* hmm_L_fom_noZ  = new TH1F("hmm_L_fom_noZ","hmm_L_fom_noZ",xbin,xmin,xmax);
  TH1F* hmm_L_fom_noZ_new  = new TH1F("hmm_L_fom_noZ_new","hmm_L_fom_noZ_new",xbin,xmin,xmax);
  TH1F* hmm_L_fom_noAC  = new TH1F("hmm_L_fom_noAC","hmm_L_fom_noAC",xbin,xmin,xmax);
  TH1F* hmm_L_fom_zdiff  = new TH1F("hmm_L_fom_zdiff","hmm_L_fom_zdiff",xbin,xmin,xmax);
  TH1F* hmm_Al_fom_noZ  = new TH1F("hmm_Al_fom_noZ","hmm_Al_fom_noZ",xbin,xmin,xmax);
  TH1F* hmm_Al_fom_noZ_new  = new TH1F("hmm_Al_fom_noZ_new","hmm_Al_fom_noZ_new",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_strict  = new TH1F("hmm_wo_bg_fom_strict","hmm_wo_bg_fom_strict",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_noZ  = new TH1F("hmm_wo_bg_fom_noZ","hmm_wo_bg_fom_noZ",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_noZ_new  = new TH1F("hmm_wo_bg_fom_noZ_new","hmm_wo_bg_fom_noZ_new",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_noAC  = new TH1F("hmm_wo_bg_fom_noAC","hmm_wo_bg_fom_noAC",xbin,xmin,xmax);
  TH1F* hmm_Al_wobg_fom_noZ  = new TH1F("hmm_Al_wobg_fom_noZ","hmm_Al_wobg_fom",xbin,xmin,xmax);
  TH1F* hmm_Al_wobg_fom_noZ_new  = new TH1F("hmm_Al_wobg_fom_noZ_new","hmm_Al_wobg_fom_new",xbin,xmin,xmax);
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
  bool event_selection_noZ = false;
  bool event_selection_noZ_new=false;
  bool event_selection_strict=false;
  bool event_selection_noAC=false;
  bool event_selection_zdiff=false;
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
  tree->Draw(">>elist" , "fabs(ct_orig)<1.006");//ctsum (does NOT dintinguish #track)
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


	


		if(fabs(ct)<1.006)ct_cut=true;
		else ct_cut=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_strict=true;
		else event_selection_strict=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_noZ_new=true;
		else event_selection_noZ_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_noAC=true;
		else event_selection_noAC=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_zdiff=true;
		else event_selection_zdiff=false;
		event_selection_noZ=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_noZ=true;
		else event_selection_noZ=false;


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
		if(event_selection_noZ&&ct_cut)hmm_L_fom_noZ->Fill(mm);
		if(event_selection_noZ_new&&ct_cut)hmm_L_fom_noZ_new->Fill(mm);
		if(event_selection_strict)hmm_L_fom_strict->Fill(mm);
		if(event_selection_noAC)hmm_L_fom_noAC->Fill(mm);
		if(event_selection_zdiff)hmm_L_fom_zdiff->Fill(mm);
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01))hmm_Al_fom_noZ->Fill(mm);//Al selection
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01))hmm_Al_fom_noZ_new->Fill(mm);//Al selection



}//ENum

	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2");
	TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	TH1F* hmm_Al_fom_nocut=(TH1F*)file->Get("hmm_Al_fom_best");
	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");//best cut pion
	//TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	//TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");//best cut
	TH1F* hmm_bg_fom_noZ=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut");//noZ
	TH1F* hmm_bg_fom_noZ_new=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_new");//strict AC, w/o Z
	TH1F* hmm_bg_fom_strict=(TH1F*)file_mea->Get("hmm_mixacc_result_new");//strict cut 
	TH1F* hmm_bg_fom_noAC=(TH1F*)file_mea->Get("hmm_mixacc_result_woac");// w/o AC, w/Z
	TH1F* hmm_bg_fom_zdiff=(TH1F*)file_mea->Get("hmm_mixacc_result_zdiff");// w/Zdiff
	TH1F* hmm_Albg_fom_noZ=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_forAl");
	TH1F* hmm_Albg_fom_noZ_new=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_new_forAl");//strict AC
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
	hmm_bg_fom_strict->Sumw2();
	hmm_bg_fom_noZ->Sumw2();
	hmm_bg_fom_noZ_new->Sumw2();
	hmm_bg_fom_noAC->Sumw2();
	hmm_Albg_fom_noZ->Sumw2();
	hmm_Albg_fom_noZ_new->Sumw2();
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_strict->Scale(1./nbunch);
	hmm_bg_fom_noZ->Scale(1./nbunch);
	hmm_bg_fom_noZ_new->Scale(1./nbunch);
	hmm_bg_fom_noAC->Scale(1./nbunch);
	hmm_Albg_fom_noZ->Scale(1./nbunch);
	hmm_Albg_fom_noZ_new->Scale(1./nbunch);
	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
	hmm_wo_bg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
	hmm_wo_bg_fom_strict->Add(hmm_L_fom_strict,hmm_bg_fom_strict,1.0,-1.0);
	hmm_wo_bg_fom_noZ->Add(hmm_L_fom_noZ,hmm_bg_fom_noZ,1.0,-1.0);
	hmm_wo_bg_fom_noZ_new->Add(hmm_L_fom_noZ_new,hmm_bg_fom_noZ_new,1.0,-1.0);
	hmm_wo_bg_fom_noAC->Add(hmm_L_fom_noAC,hmm_bg_fom_noAC,1.0,-1.0);
	hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
	//hmm_pi_wobg_fom_nocut->Add(hmm_pi_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	hmm_Al_wobg_fom_noZ->Add(hmm_Al_fom_noZ,hmm_Albg_fom_noZ,1.0,-1.0);
	hmm_Al_wobg_fom_noZ_new->Add(hmm_Al_fom_noZ_new,hmm_Albg_fom_noZ_new,1.0,-1.0);


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

	 TF1 *fAl=new TF1("fAl",F_Voigt,fit_min_mm,fit_max_mm,4);
	 fAl->SetNpx(2000);
	 fAl->SetTitle("Al selected");
	 fAl->SetParameters(3.,0.05,0.04,0.001);
	 fAl->SetParLimits(0,0.,10000.);
	 fAl->SetLineColor(kRed);
	 hmm_Al_wobg_fom_noZ->Fit("fAl","N","",0.,0.1);
	 double Al_par0 = fAl->GetParameter(0);
	 double Al_par1 = fAl->GetParameter(1);
	 double Al_par2 = fAl->GetParameter(2);
	 double Al_par3 = fAl->GetParameter(3);

	 TF1 *fAl_new=new TF1("fAl_new",F_Voigt,fit_min_mm,fit_max_mm,4);
	 fAl_new->SetNpx(2000);
	 fAl_new->SetTitle("Al selected (new)");
	 fAl_new->SetParameters(2.,0.05,0.04,0.001);
	 fAl_new->SetParLimits(0,0.,10000.);
	 fAl_new->SetLineColor(kRed);
	 hmm_Al_wobg_fom_noZ_new->Fit("fAl_new","N","",0.,0.1);
	 double Al_new_par0 = fAl_new->GetParameter(0);
	 double Al_new_par1 = fAl_new->GetParameter(1);
	 double Al_new_par2 = fAl_new->GetParameter(2);
	 double Al_new_par3 = fAl_new->GetParameter(3);

	 TF1 *fpion=new TF1("fpion",F_Voigt,fit_min_mm,fit_max_mm,4);
	 fpion->SetNpx(2000);
	 fpion->SetTitle("Pion selected");
	 fpion->SetParameters(3.,0.05,0.04,0.001);
	 fpion->SetParLimits(0,0.,10000.);
	 fpion->SetLineColor(kRed);
	 hmm_pi_wobg_fom_best->Fit("fpion","N","",0.,0.1);
	 double pion_par0 = fpion->GetParameter(0);
	 double pion_par1 = fpion->GetParameter(1);
	 double pion_par2 = fpion->GetParameter(2);
	 double pion_par3 = fpion->GetParameter(3);

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

/*%%%%%%%%%%%%%%%%*/
/*%%    Lexp	%%*/
/*%%%%%%%%%%%%%%%%*/
	 cout<<"Best Cut START"<<endl;
	 cout<<"(Landau+Exp)*(Gauss) START"<<endl;
	 fmm_best_Lexp=new TF1("fmm_best_Lexp",FMM_Res,fit_min_mm,fit_max_mm,14);
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
   //
   //
   //
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
	 fmm_best_Lexp->SetNpx(20000);
	 fmm_best_Lexp->SetTitle("Missing Mass (best)");
	 fmm_best_Lexp->SetParLimits(2,0.,1000.);//positive
	 fmm_best_Lexp->SetParLimits(9,0.,300.);//positive
	 fmm_best_Lexp->SetParameter(0,0.0007);//Landau width
	 fmm_best_Lexp->SetParameter(1,mean_L_best);
	 fmm_best_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_best_Lexp->SetParameter(2,1.5);//total scale
	 fmm_best_Lexp->SetParameter(3,0.001);//sigma
	 fmm_best_Lexp->SetParLimits(3,0.,0.01);
	 fmm_best_Lexp->SetParameter(4,0.05);//att.
	 fmm_best_Lexp->SetParLimits(4,0.005,0.08);
	 fmm_best_Lexp->SetParameter(5,-0.004);//peak pos.
	 fmm_best_Lexp->SetParLimits(5,-0.05,0.05);
	 fmm_best_Lexp->SetParameter(6,0.6);//relative strength
	 fmm_best_Lexp->SetParLimits(6,0.,1.5);//relative strength

	 fmm_best_Lexp->SetParameter(7,0.0003);//Landau width
	 fmm_best_Lexp->SetParameter(8,mean_S_best);//MPV
	 fmm_best_Lexp->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_best_Lexp->SetParameter(9,0.4);//total scale
	 fmm_best_Lexp->SetParameter(10,0.0015);//sigma
	 fmm_best_Lexp->SetParLimits(10,0.,0.01);
	 fmm_best_Lexp->SetParameter(11,0.05);//att
	 fmm_best_Lexp->SetParLimits(11,0.03,0.12);
	 fmm_best_Lexp->SetParameter(12,0.080);//peak pos.
	 //fmm_best_Lexp->SetParLimits(15,-0.085,-0.055);
	 fmm_best_Lexp->SetParameter(13,0.6);
	 fmm_best_Lexp->SetParLimits(13,0.,1.5);//relative strength

	 hmm_wo_bg_fom_best->Fit("fmm_best_Lexp","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
	 double chisq_best = fmm_best_Lexp->GetChisquare();
	 double dof_best  = fmm_best_Lexp->GetNDF();
	 cout<<"chisq_best="<<chisq_best<<endl;
	 cout<<"dof="<<dof_best<<endl;
	 cout<<"Reduced chi-square = "<<chisq_best/dof_best<<endl;


	 TF1* fmm_Lambda_only_best = new TF1("fmm_Lambda_only_best",FMM_Response,fit_min_mm,fit_max_mm,7);
	 TF1* fmm_Sigma_only_best  = new TF1("fmm_Sigma_only_best" ,FMM_Response, fit_min_mm,fit_max_mm,7);
//Lambda_only_best
	 fmm_Lambda_only_best->SetNpx(20000);
	 fmm_Lambda_only_best->SetParameter(0,fmm_best_Lexp->GetParameter(0));
	 fmm_Lambda_only_best->SetParameter(1,fmm_best_Lexp->GetParameter(1));
	 fmm_Lambda_only_best->SetParameter(2,fmm_best_Lexp->GetParameter(2));
	 fmm_Lambda_only_best->SetParameter(3,fmm_best_Lexp->GetParameter(3));
	 fmm_Lambda_only_best->SetParameter(4,fmm_best_Lexp->GetParameter(4));
	 fmm_Lambda_only_best->SetParameter(5,fmm_best_Lexp->GetParameter(5));
	 fmm_Lambda_only_best->SetParameter(6,fmm_best_Lexp->GetParameter(6));
//Sigma_only_best
	 fmm_Sigma_only_best->SetNpx(20000);
	 fmm_Sigma_only_best->SetParameter(0,fmm_best_Lexp->GetParameter(7));
	 fmm_Sigma_only_best->SetParameter(1,fmm_best_Lexp->GetParameter(8));
	 fmm_Sigma_only_best->SetParameter(2,fmm_best_Lexp->GetParameter(9));
	 fmm_Sigma_only_best->SetParameter(3,fmm_best_Lexp->GetParameter(10));
	 fmm_Sigma_only_best->SetParameter(4,fmm_best_Lexp->GetParameter(11));
	 fmm_Sigma_only_best->SetParameter(5,fmm_best_Lexp->GetParameter(12));
	 fmm_Sigma_only_best->SetParameter(6,fmm_best_Lexp->GetParameter(13));

	double nofL_best = fmm_Lambda_only_best->Integral(fit_min_mm,fit_max_mm);
	double nofL_old_best = fmm_Lambda_only_best->Integral(-0.006,0.006);
	nofL_best = nofL_best/fit_bin_width;
	nofL_old_best = nofL_old_best/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL_best<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old_best<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(-0.006),hmm_wo_bg_fom_best->FindBin(0.006))<<endl;

	double nofS_best = fmm_Sigma_only_best->Integral(fit_min_mm,fit_max_mm);
	double nofS_old_best = fmm_Sigma_only_best->Integral(def_mean_S-0.008,def_mean_S+0.008);
	nofS_best = nofS_best/fit_bin_width;
	nofS_old_best = nofS_old_best/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS_best<<endl;
	cout<<"Number of Sigma w/o radiative tail (TF1 Integral) = "<<nofS_old_best<<endl;
	cout<<"Number of Sigma w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(def_mean_S-0.008),hmm_wo_bg_fom_best->FindBin(def_mean_S+0.008))<<endl;

	 hmm_wo_bg_fom_best->Draw();
	 fmm_best_Lexp->SetLineColor(kRed);
	 fmm_best_Lexp->Draw("same");
	 fmm_Lambda_only_best->SetLineColor(kAzure);
	 fmm_Sigma_only_best->SetLineColor(kCyan);
	 fmm_Lambda_only_best->Draw("same");
	 fmm_Sigma_only_best->Draw("same");
	 //fL_best->SetLineColor(kGreen);
	 //fS_best->SetLineColor(kGreen);
	 //fL_best->Draw("same");
	 //fS_best->Draw("same");


/*%%%%%%%%%%%%%%%%%%%%*/
/*%%    No Z cut	%%*/
/*%%%%%%%%%%%%%%%%%%%%*/
	TCanvas* c3 = new TCanvas("c3","c3");
	 cout<<"noZ START"<<endl;
	 cout<<"(Landau+Exp)*(Gauss) MODE START"<<endl;
	 fmm_noZ_Lexp=new TF1("fmm_noZ_Lexp",FMM_Res_nocut,fit_min_mm,fit_max_mm,22);
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
   //
   //
   //
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
	 fmm_noZ_Lexp->SetNpx(20000);
	 fmm_noZ_Lexp->SetTitle("Missing Mass (noZ)");
	 fmm_noZ_Lexp->FixParameter(0,fmm_best_Lexp->GetParameter(0));
	 fmm_noZ_Lexp->FixParameter(1,fmm_best_Lexp->GetParameter(1));
	 fmm_noZ_Lexp->SetParameter(2,fmm_best_Lexp->GetParameter(2));//L scale
	 fmm_noZ_Lexp->FixParameter(3,fmm_best_Lexp->GetParameter(3));
	 fmm_noZ_Lexp->FixParameter(4,fmm_best_Lexp->GetParameter(4));
	 fmm_noZ_Lexp->FixParameter(5,fmm_best_Lexp->GetParameter(5));
	 fmm_noZ_Lexp->FixParameter(6,fmm_best_Lexp->GetParameter(6));
	 fmm_noZ_Lexp->FixParameter(7,fmm_best_Lexp->GetParameter(7));
	 fmm_noZ_Lexp->FixParameter(8,fmm_best_Lexp->GetParameter(8));
	 fmm_noZ_Lexp->SetParameter(9,fmm_best_Lexp->GetParameter(9));//S scale
	 fmm_noZ_Lexp->FixParameter(10,fmm_best_Lexp->GetParameter(10));
	 fmm_noZ_Lexp->FixParameter(11,fmm_best_Lexp->GetParameter(11));
	 fmm_noZ_Lexp->FixParameter(12,fmm_best_Lexp->GetParameter(12));
	 fmm_noZ_Lexp->FixParameter(13,fmm_best_Lexp->GetParameter(13));
	 //fmm_noZ_Lexp->SetParameter(14,500.);//scale
	 //fmm_noZ_Lexp->SetParLimits(14,0.,1000000.);//scale
	 //fmm_noZ_Lexp->SetParameter(15,0.05);//mean
	 //fmm_noZ_Lexp->SetParameter(16,0.04);//Gsigma
	 //fmm_noZ_Lexp->SetParameter(17,0.01);//Lfwhm
	 fmm_noZ_Lexp->SetParameter(14,6.);//scale
	 fmm_noZ_Lexp->SetParLimits(14,0.,1000000.);//scale
	 fmm_noZ_Lexp->FixParameter(15,Al_par1);//mean
	 fmm_noZ_Lexp->FixParameter(16,Al_par2);//Gsigma
	 fmm_noZ_Lexp->FixParameter(17,Al_par3);//Lfwhm
	 fmm_noZ_Lexp->FixParameter(18,pion_par1);//mean
	 fmm_noZ_Lexp->FixParameter(19,pion_par2);//Gsigma
	 fmm_noZ_Lexp->FixParameter(20,pion_par3);//Lfwhm
	 fmm_noZ_Lexp->FixParameter(21,0.9);//Al vs Pi
	 fmm_noZ_Lexp->SetParLimits(21,0.8,1.);
	 //fmm_noZ_Lexp->SetParLimits(2,0.,1000.);//positive
	 //fmm_noZ_Lexp->SetParLimits(9,0.,300.);//positive
	 //fmm_noZ_Lexp->SetParameter(0,0.0007);//Landau width
	 //fmm_noZ_Lexp->SetParameter(1,mean_L_noZ);
	 //fmm_noZ_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_noZ_Lexp->SetParameter(2,1.5);//total scale
	 //fmm_noZ_Lexp->SetParameter(3,0.001);//sigma
	 //fmm_noZ_Lexp->SetParLimits(3,0.,0.01);
	 //fmm_noZ_Lexp->SetParameter(4,0.05);//att.
	 //fmm_noZ_Lexp->SetParLimits(4,0.005,0.08);
	 //fmm_noZ_Lexp->SetParameter(5,-0.004);//peak pos.
	 //fmm_noZ_Lexp->SetParLimits(5,-0.05,0.05);
	 //fmm_noZ_Lexp->SetParameter(6,0.6);//relative strength
	 //fmm_noZ_Lexp->SetParLimits(6,0.,1.5);//relative strength

	 //fmm_noZ_Lexp->SetParameter(7,0.0003);//Landau width
	 //fmm_noZ_Lexp->SetParameter(8,mean_S_noZ);//MPV
	 //fmm_noZ_Lexp->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 //fmm_noZ_Lexp->SetParameter(9,0.4);//total scale
	 //fmm_noZ_Lexp->SetParameter(10,0.0015);//sigma
	 //fmm_noZ_Lexp->SetParLimits(10,0.,0.01);
	 //fmm_noZ_Lexp->SetParameter(11,0.05);//att
	 //fmm_noZ_Lexp->SetParLimits(11,0.03,0.12);
	 //fmm_noZ_Lexp->SetParameter(12,0.080);//peak pos.
	 ////fmm_noZ_Lexp->SetParLimits(15,-0.085,-0.055);
	 //fmm_noZ_Lexp->SetParameter(13,0.6);
	 //fmm_noZ_Lexp->SetParLimits(13,0.,1.5);//relative strength

	 hmm_wo_bg_fom_noZ->Fit("fmm_noZ_Lexp","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
	 double chisq_noZ = fmm_noZ_Lexp->GetChisquare();
	 double dof_noZ  = fmm_noZ_Lexp->GetNDF();
	 cout<<"chisq_noZ="<<chisq_noZ<<endl;
	 cout<<"dof="<<dof_noZ<<endl;
	 cout<<"Reduced chi-square = "<<chisq_noZ/dof_noZ<<endl;


	 TF1* fmm_Lambda_only_noZ = new TF1("fmm_Lambda_only_noZ",FMM_Response,fit_min_mm,fit_max_mm,7);
	 TF1* fmm_Sigma_only_noZ  = new TF1("fmm_Sigma_only_noZ" ,FMM_Response, fit_min_mm,fit_max_mm,7);
	 TF1* fmm_bg_only_noZ  = new TF1("fmm_bg_only_noZ" ,FMM_2BG, fit_min_mm,fit_max_mm,8);
//Lambda_only_noZ
	 fmm_Lambda_only_noZ->SetNpx(20000);
	 fmm_Lambda_only_noZ->SetParameter(0,fmm_noZ_Lexp->GetParameter(0));
	 fmm_Lambda_only_noZ->SetParameter(1,fmm_noZ_Lexp->GetParameter(1));
	 fmm_Lambda_only_noZ->SetParameter(2,fmm_noZ_Lexp->GetParameter(2));
	 fmm_Lambda_only_noZ->SetParameter(3,fmm_noZ_Lexp->GetParameter(3));
	 fmm_Lambda_only_noZ->SetParameter(4,fmm_noZ_Lexp->GetParameter(4));
	 fmm_Lambda_only_noZ->SetParameter(5,fmm_noZ_Lexp->GetParameter(5));
	 fmm_Lambda_only_noZ->SetParameter(6,fmm_noZ_Lexp->GetParameter(6));
//Sigma_only_noZ
	 fmm_Sigma_only_noZ->SetNpx(20000);
	 fmm_Sigma_only_noZ->SetParameter(0,fmm_noZ_Lexp->GetParameter(7));
	 fmm_Sigma_only_noZ->SetParameter(1,fmm_noZ_Lexp->GetParameter(8));
	 fmm_Sigma_only_noZ->SetParameter(2,fmm_noZ_Lexp->GetParameter(9));
	 fmm_Sigma_only_noZ->SetParameter(3,fmm_noZ_Lexp->GetParameter(10));
	 fmm_Sigma_only_noZ->SetParameter(4,fmm_noZ_Lexp->GetParameter(11));
	 fmm_Sigma_only_noZ->SetParameter(5,fmm_noZ_Lexp->GetParameter(12));
	 fmm_Sigma_only_noZ->SetParameter(6,fmm_noZ_Lexp->GetParameter(13));
//bg_only_noZ
	 fmm_bg_only_noZ->SetNpx(20000);
	 fmm_bg_only_noZ->SetParameter(0,fmm_noZ_Lexp->GetParameter(14));
	 fmm_bg_only_noZ->SetParameter(1,fmm_noZ_Lexp->GetParameter(15));
	 fmm_bg_only_noZ->SetParameter(2,fmm_noZ_Lexp->GetParameter(16));
	 fmm_bg_only_noZ->SetParameter(3,fmm_noZ_Lexp->GetParameter(17));
	 fmm_bg_only_noZ->SetParameter(4,fmm_noZ_Lexp->GetParameter(18));
	 fmm_bg_only_noZ->SetParameter(5,fmm_noZ_Lexp->GetParameter(19));
	 fmm_bg_only_noZ->SetParameter(6,fmm_noZ_Lexp->GetParameter(20));
	 fmm_bg_only_noZ->SetParameter(7,fmm_noZ_Lexp->GetParameter(21));

	double nofL_noZ = fmm_Lambda_only_noZ->Integral(fit_min_mm,fit_max_mm);
	double nofL_old_noZ = fmm_Lambda_only_noZ->Integral(-0.006,0.006);
	nofL_noZ = nofL_noZ/fit_bin_width;
	nofL_old_noZ = nofL_old_noZ/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL_noZ<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old_noZ<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(-0.006),hmm_wo_bg_fom_noZ->FindBin(0.006))<<endl;

	double nofS_noZ = fmm_Sigma_only_noZ->Integral(fit_min_mm,fit_max_mm);
	double nofS_old_noZ = fmm_Sigma_only_noZ->Integral(def_mean_S-0.008,def_mean_S+0.008);
	nofS_noZ = nofS_noZ/fit_bin_width;
	nofS_old_noZ = nofS_old_noZ/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS_noZ<<endl;
	cout<<"Number of Sigma w/o radiative tail (TF1 Integral) = "<<nofS_old_noZ<<endl;
	cout<<"Number of Sigma w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(def_mean_S-0.008),hmm_wo_bg_fom_noZ->FindBin(def_mean_S+0.008))<<endl;

	 hmm_wo_bg_fom_noZ->Draw();
	 fmm_noZ_Lexp->SetLineColor(kRed);
	 fmm_noZ_Lexp->Draw("same");
	 fmm_Lambda_only_noZ->SetLineColor(kAzure);
	 fmm_Sigma_only_noZ->SetLineColor(kCyan);
	 fmm_bg_only_noZ->SetLineColor(kOrange);
	 fmm_Lambda_only_noZ->Draw("same");
	 fmm_Sigma_only_noZ->Draw("same");
	 fmm_bg_only_noZ->Draw("same");
	 //fL_noZ->SetLineColor(kGreen);
	 //fS_noZ->SetLineColor(kGreen);
	 //fL_noZ->Draw("same");
	 //fS_noZ->Draw("same");
	TCanvas* c4 = new TCanvas("c4","c4");
	hmm_Al_wobg_fom_noZ->Draw("");
	fAl->Draw("same");

/*%%%%%%%%%%%%%%%%%%%%*/
/*%%    Strict cut	%%*/
/*%%%%%%%%%%%%%%%%%%%%*/
	TCanvas* c5 = new TCanvas("c5","c5");
	 cout<<"strict cut START"<<endl;
	 cout<<"(Landau+Exp)*(Gauss) MODE START"<<endl;
	 fmm_strict_Lexp=new TF1("fmm_strict_Lexp",FMM_Res_nocut,fit_min_mm,fit_max_mm,22);
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
   //
   //
   //
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
	 fmm_strict_Lexp->SetNpx(20000);
	 fmm_strict_Lexp->SetTitle("Missing Mass (strict)");
	 fmm_strict_Lexp->FixParameter(0,fmm_best_Lexp->GetParameter(0));
	 fmm_strict_Lexp->FixParameter(1,fmm_best_Lexp->GetParameter(1));
	 fmm_strict_Lexp->SetParameter(2,fmm_best_Lexp->GetParameter(2));//L scale
	 fmm_strict_Lexp->SetParLimits(2,0.,100000.);//L scale
	 fmm_strict_Lexp->FixParameter(3,fmm_best_Lexp->GetParameter(3));
	 fmm_strict_Lexp->FixParameter(4,fmm_best_Lexp->GetParameter(4));
	 fmm_strict_Lexp->FixParameter(5,fmm_best_Lexp->GetParameter(5));
	 fmm_strict_Lexp->FixParameter(6,fmm_best_Lexp->GetParameter(6));
	 fmm_strict_Lexp->FixParameter(7,fmm_best_Lexp->GetParameter(7));
	 fmm_strict_Lexp->FixParameter(8,fmm_best_Lexp->GetParameter(8));
	 fmm_strict_Lexp->SetParameter(9,fmm_best_Lexp->GetParameter(9));//S scale
	 fmm_strict_Lexp->SetParLimits(9,0.,100000.);//L scale
	 fmm_strict_Lexp->FixParameter(10,fmm_best_Lexp->GetParameter(10));
	 fmm_strict_Lexp->FixParameter(11,fmm_best_Lexp->GetParameter(11));
	 fmm_strict_Lexp->FixParameter(12,fmm_best_Lexp->GetParameter(12));
	 fmm_strict_Lexp->FixParameter(13,fmm_best_Lexp->GetParameter(13));
	 //fmm_strict_Lexp->SetParameter(14,500.);//scale
	 //fmm_strict_Lexp->SetParLimits(14,0.,1000000.);//scale
	 //fmm_strict_Lexp->SetParameter(15,0.05);//mean
	 //fmm_strict_Lexp->SetParameter(16,0.04);//Gsigma
	 //fmm_strict_Lexp->SetParameter(17,0.01);//Lfwhm
	 fmm_strict_Lexp->SetParameter(14,6.);//scale
	 fmm_strict_Lexp->SetParLimits(14,0.,1000000.);//scale
	 fmm_strict_Lexp->FixParameter(15,Al_par1);//mean
	 fmm_strict_Lexp->FixParameter(16,Al_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(17,Al_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(18,pion_par1);//mean
	 fmm_strict_Lexp->FixParameter(19,pion_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(20,pion_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(21,0.1);//Al vs Pi
	 fmm_strict_Lexp->SetParLimits(21,0.0,0.2);
	 //fmm_strict_Lexp->SetParLimits(2,0.,1000.);//positive
	 //fmm_strict_Lexp->SetParLimits(9,0.,300.);//positive
	 //fmm_strict_Lexp->SetParameter(0,0.0007);//Landau width
	 //fmm_strict_Lexp->SetParameter(1,mean_L_strict);
	 //fmm_strict_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_strict_Lexp->SetParameter(2,1.5);//total scale
	 //fmm_strict_Lexp->SetParameter(3,0.001);//sigma
	 //fmm_strict_Lexp->SetParLimits(3,0.,0.01);
	 //fmm_strict_Lexp->SetParameter(4,0.05);//att.
	 //fmm_strict_Lexp->SetParLimits(4,0.005,0.08);
	 //fmm_strict_Lexp->SetParameter(5,-0.004);//peak pos.
	 //fmm_strict_Lexp->SetParLimits(5,-0.05,0.05);
	 //fmm_strict_Lexp->SetParameter(6,0.6);//relative strength
	 //fmm_strict_Lexp->SetParLimits(6,0.,1.5);//relative strength

	 //fmm_strict_Lexp->SetParameter(7,0.0003);//Landau width
	 //fmm_strict_Lexp->SetParameter(8,mean_S_strict);//MPV
	 //fmm_strict_Lexp->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 //fmm_strict_Lexp->SetParameter(9,0.4);//total scale
	 //fmm_strict_Lexp->SetParameter(10,0.0015);//sigma
	 //fmm_strict_Lexp->SetParLimits(10,0.,0.01);
	 //fmm_strict_Lexp->SetParameter(11,0.05);//att
	 //fmm_strict_Lexp->SetParLimits(11,0.03,0.12);
	 //fmm_strict_Lexp->SetParameter(12,0.080);//peak pos.
	 ////fmm_strict_Lexp->SetParLimits(15,-0.085,-0.055);
	 //fmm_strict_Lexp->SetParameter(13,0.6);
	 //fmm_strict_Lexp->SetParLimits(13,0.,1.5);//relative strength

	 hmm_wo_bg_fom_strict->Fit("fmm_strict_Lexp","","",fit_min_mm,0.097);//Total fitting w/ 4Poly BG
	 double chisq_strict = fmm_strict_Lexp->GetChisquare();
	 double dof_strict  = fmm_strict_Lexp->GetNDF();
	 cout<<"chisq_strict="<<chisq_strict<<endl;
	 cout<<"dof="<<dof_strict<<endl;
	 cout<<"Reduced chi-square = "<<chisq_strict/dof_strict<<endl;


	 TF1* fmm_Lambda_only_strict = new TF1("fmm_Lambda_only_strict",FMM_Response,fit_min_mm,fit_max_mm,7);
	 TF1* fmm_Sigma_only_strict  = new TF1("fmm_Sigma_only_strict" ,FMM_Response, fit_min_mm,fit_max_mm,7);
	 TF1* fmm_bg_only_strict  = new TF1("fmm_bg_only_strict" ,FMM_2BG, fit_min_mm,fit_max_mm,8);
//Lambda_only_strict
	 fmm_Lambda_only_strict->SetNpx(20000);
	 fmm_Lambda_only_strict->SetParameter(0,fmm_strict_Lexp->GetParameter(0));
	 fmm_Lambda_only_strict->SetParameter(1,fmm_strict_Lexp->GetParameter(1));
	 fmm_Lambda_only_strict->SetParameter(2,fmm_strict_Lexp->GetParameter(2));
	 fmm_Lambda_only_strict->SetParameter(3,fmm_strict_Lexp->GetParameter(3));
	 fmm_Lambda_only_strict->SetParameter(4,fmm_strict_Lexp->GetParameter(4));
	 fmm_Lambda_only_strict->SetParameter(5,fmm_strict_Lexp->GetParameter(5));
	 fmm_Lambda_only_strict->SetParameter(6,fmm_strict_Lexp->GetParameter(6));
//Sigma_only_strict
	 fmm_Sigma_only_strict->SetNpx(20000);
	 fmm_Sigma_only_strict->SetParameter(0,fmm_strict_Lexp->GetParameter(7));
	 fmm_Sigma_only_strict->SetParameter(1,fmm_strict_Lexp->GetParameter(8));
	 fmm_Sigma_only_strict->SetParameter(2,fmm_strict_Lexp->GetParameter(9));
	 fmm_Sigma_only_strict->SetParameter(3,fmm_strict_Lexp->GetParameter(10));
	 fmm_Sigma_only_strict->SetParameter(4,fmm_strict_Lexp->GetParameter(11));
	 fmm_Sigma_only_strict->SetParameter(5,fmm_strict_Lexp->GetParameter(12));
	 fmm_Sigma_only_strict->SetParameter(6,fmm_strict_Lexp->GetParameter(13));
//bg_only_strict
	 fmm_bg_only_strict->SetNpx(20000);
	 fmm_bg_only_strict->SetParameter(0,fmm_strict_Lexp->GetParameter(14));
	 fmm_bg_only_strict->SetParameter(1,fmm_strict_Lexp->GetParameter(15));
	 fmm_bg_only_strict->SetParameter(2,fmm_strict_Lexp->GetParameter(16));
	 fmm_bg_only_strict->SetParameter(3,fmm_strict_Lexp->GetParameter(17));
	 fmm_bg_only_strict->SetParameter(4,fmm_strict_Lexp->GetParameter(18));
	 fmm_bg_only_strict->SetParameter(5,fmm_strict_Lexp->GetParameter(19));
	 fmm_bg_only_strict->SetParameter(6,fmm_strict_Lexp->GetParameter(20));
	 fmm_bg_only_strict->SetParameter(7,fmm_strict_Lexp->GetParameter(21));

	double nofL_strict = fmm_Lambda_only_strict->Integral(fit_min_mm,fit_max_mm);
	double nofL_old_strict = fmm_Lambda_only_strict->Integral(-0.006,0.006);
	nofL_strict = nofL_strict/fit_bin_width;
	nofL_old_strict = nofL_old_strict/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL_strict<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old_strict<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_strict->Integral(hmm_wo_bg_fom_strict->FindBin(-0.006),hmm_wo_bg_fom_strict->FindBin(0.006))<<endl;

	double nofS_strict = fmm_Sigma_only_strict->Integral(fit_min_mm,fit_max_mm);
	double nofS_old_strict = fmm_Sigma_only_strict->Integral(def_mean_S-0.008,def_mean_S+0.008);
	nofS_strict = nofS_strict/fit_bin_width;
	nofS_old_strict = nofS_old_strict/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS_strict<<endl;
	cout<<"Number of Sigma w/o radiative tail (TF1 Integral) = "<<nofS_old_strict<<endl;
	cout<<"Number of Sigma w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_strict->Integral(hmm_wo_bg_fom_strict->FindBin(def_mean_S-0.008),hmm_wo_bg_fom_strict->FindBin(def_mean_S+0.008))<<endl;

	 hmm_wo_bg_fom_strict->Draw();
	 fmm_strict_Lexp->SetLineColor(kRed);
	 fmm_strict_Lexp->Draw("same");
	 fmm_Lambda_only_strict->SetLineColor(kAzure);
	 fmm_Sigma_only_strict->SetLineColor(kCyan);
	 fmm_bg_only_strict->SetLineColor(kOrange);
	 fmm_Lambda_only_strict->Draw("same");
	 fmm_Sigma_only_strict->Draw("same");
	 fmm_bg_only_strict->Draw("same");
	 //fL_strict->SetLineColor(kGreen);
	 //fS_strict->SetLineColor(kGreen);
	 //fL_strict->Draw("same");
	 //fS_strict->Draw("same");

/*%%%%%%%%%%%%%%%%%%%%*/
/*%%    No AC cut	%%*/
/*%%%%%%%%%%%%%%%%%%%%*/
	TCanvas* c6 = new TCanvas("c6","c6");
	 cout<<"noAC START"<<endl;
	 cout<<"(Landau+Exp)*(Gauss) MODE START"<<endl;
	 fmm_noAC_Lexp=new TF1("fmm_noAC_Lexp",FMM_Res_nocut,fit_min_mm,fit_max_mm,22);
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=tau of exp function
   //par[5]=Shift of Function Peak
   //par[6]=Relative Strength
   //
   //
   //
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
	 fmm_noAC_Lexp->SetNpx(20000);
	 fmm_noAC_Lexp->SetTitle("Missing Mass (noAC)");
	 fmm_noAC_Lexp->FixParameter(0,fmm_best_Lexp->GetParameter(0));
	 fmm_noAC_Lexp->FixParameter(1,fmm_best_Lexp->GetParameter(1));
	 fmm_noAC_Lexp->SetParameter(2,fmm_best_Lexp->GetParameter(2));//L scale
	 fmm_noAC_Lexp->FixParameter(3,fmm_best_Lexp->GetParameter(3));
	 fmm_noAC_Lexp->FixParameter(4,fmm_best_Lexp->GetParameter(4));
	 fmm_noAC_Lexp->FixParameter(5,fmm_best_Lexp->GetParameter(5));
	 fmm_noAC_Lexp->FixParameter(6,fmm_best_Lexp->GetParameter(6));
	 fmm_noAC_Lexp->FixParameter(7,fmm_best_Lexp->GetParameter(7));
	 fmm_noAC_Lexp->FixParameter(8,fmm_best_Lexp->GetParameter(8));
	 fmm_noAC_Lexp->SetParameter(9,fmm_best_Lexp->GetParameter(9));//S scale
	 fmm_noAC_Lexp->FixParameter(10,fmm_best_Lexp->GetParameter(10));
	 fmm_noAC_Lexp->FixParameter(11,fmm_best_Lexp->GetParameter(11));
	 fmm_noAC_Lexp->FixParameter(12,fmm_best_Lexp->GetParameter(12));
	 fmm_noAC_Lexp->FixParameter(13,fmm_best_Lexp->GetParameter(13));
	 //fmm_noAC_Lexp->SetParameter(14,500.);//scale
	 //fmm_noAC_Lexp->SetParLimits(14,0.,1000000.);//scale
	 //fmm_noAC_Lexp->SetParameter(15,0.05);//mean
	 //fmm_noAC_Lexp->SetParameter(16,0.04);//Gsigma
	 //fmm_noAC_Lexp->SetParameter(17,0.01);//Lfwhm
	 fmm_noAC_Lexp->SetParameter(14,6.);//scale
	 fmm_noAC_Lexp->SetParLimits(14,0.,1000000.);//scale
	 fmm_noAC_Lexp->FixParameter(15,Al_par1);//mean
	 fmm_noAC_Lexp->FixParameter(16,Al_par2);//Gsigma
	 fmm_noAC_Lexp->FixParameter(17,Al_par3);//Lfwhm
	 fmm_noAC_Lexp->FixParameter(18,pion_par1);//mean
	 fmm_noAC_Lexp->FixParameter(19,pion_par2);//Gsigma
	 fmm_noAC_Lexp->FixParameter(20,pion_par3);//Lfwhm
	 fmm_noAC_Lexp->FixParameter(21,0.1);//Al vs Pi
	 fmm_noAC_Lexp->SetParLimits(21,0.,0.2);
	 //fmm_noAC_Lexp->SetParLimits(2,0.,1000.);//positive
	 //fmm_noAC_Lexp->SetParLimits(9,0.,300.);//positive
	 //fmm_noAC_Lexp->SetParameter(0,0.0007);//Landau width
	 //fmm_noAC_Lexp->SetParameter(1,mean_L_noAC);
	 //fmm_noAC_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_noAC_Lexp->SetParameter(2,1.5);//total scale
	 //fmm_noAC_Lexp->SetParameter(3,0.001);//sigma
	 //fmm_noAC_Lexp->SetParLimits(3,0.,0.01);
	 //fmm_noAC_Lexp->SetParameter(4,0.05);//att.
	 //fmm_noAC_Lexp->SetParLimits(4,0.005,0.08);
	 //fmm_noAC_Lexp->SetParameter(5,-0.004);//peak pos.
	 //fmm_noAC_Lexp->SetParLimits(5,-0.05,0.05);
	 //fmm_noAC_Lexp->SetParameter(6,0.6);//relative strength
	 //fmm_noAC_Lexp->SetParLimits(6,0.,1.5);//relative strength

	 //fmm_noAC_Lexp->SetParameter(7,0.0003);//Landau width
	 //fmm_noAC_Lexp->SetParameter(8,mean_S_noAC);//MPV
	 //fmm_noAC_Lexp->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 //fmm_noAC_Lexp->SetParameter(9,0.4);//total scale
	 //fmm_noAC_Lexp->SetParameter(10,0.0015);//sigma
	 //fmm_noAC_Lexp->SetParLimits(10,0.,0.01);
	 //fmm_noAC_Lexp->SetParameter(11,0.05);//att
	 //fmm_noAC_Lexp->SetParLimits(11,0.03,0.12);
	 //fmm_noAC_Lexp->SetParameter(12,0.080);//peak pos.
	 ////fmm_noAC_Lexp->SetParLimits(15,-0.085,-0.055);
	 //fmm_noAC_Lexp->SetParameter(13,0.6);
	 //fmm_noAC_Lexp->SetParLimits(13,0.,1.5);//relative strength

	 hmm_wo_bg_fom_noAC->Fit("fmm_noAC_Lexp","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
	 double chisq_noAC = fmm_noAC_Lexp->GetChisquare();
	 double dof_noAC  = fmm_noAC_Lexp->GetNDF();
	 cout<<"chisq_noAC="<<chisq_noAC<<endl;
	 cout<<"dof="<<dof_noAC<<endl;
	 cout<<"Reduced chi-square = "<<chisq_noAC/dof_noAC<<endl;


	 TF1* fmm_Lambda_only_noAC = new TF1("fmm_Lambda_only_noAC",FMM_Response,fit_min_mm,fit_max_mm,7);
	 TF1* fmm_Sigma_only_noAC  = new TF1("fmm_Sigma_only_noAC" ,FMM_Response, fit_min_mm,fit_max_mm,7);
	 TF1* fmm_bg_only_noAC  = new TF1("fmm_bg_only_noAC" ,FMM_2BG, fit_min_mm,fit_max_mm,8);
//Lambda_only_noAC
	 fmm_Lambda_only_noAC->SetNpx(20000);
	 fmm_Lambda_only_noAC->SetParameter(0,fmm_noAC_Lexp->GetParameter(0));
	 fmm_Lambda_only_noAC->SetParameter(1,fmm_noAC_Lexp->GetParameter(1));
	 fmm_Lambda_only_noAC->SetParameter(2,fmm_noAC_Lexp->GetParameter(2));
	 fmm_Lambda_only_noAC->SetParameter(3,fmm_noAC_Lexp->GetParameter(3));
	 fmm_Lambda_only_noAC->SetParameter(4,fmm_noAC_Lexp->GetParameter(4));
	 fmm_Lambda_only_noAC->SetParameter(5,fmm_noAC_Lexp->GetParameter(5));
	 fmm_Lambda_only_noAC->SetParameter(6,fmm_noAC_Lexp->GetParameter(6));
//Sigma_only_noAC
	 fmm_Sigma_only_noAC->SetNpx(20000);
	 fmm_Sigma_only_noAC->SetParameter(0,fmm_noAC_Lexp->GetParameter(7));
	 fmm_Sigma_only_noAC->SetParameter(1,fmm_noAC_Lexp->GetParameter(8));
	 fmm_Sigma_only_noAC->SetParameter(2,fmm_noAC_Lexp->GetParameter(9));
	 fmm_Sigma_only_noAC->SetParameter(3,fmm_noAC_Lexp->GetParameter(10));
	 fmm_Sigma_only_noAC->SetParameter(4,fmm_noAC_Lexp->GetParameter(11));
	 fmm_Sigma_only_noAC->SetParameter(5,fmm_noAC_Lexp->GetParameter(12));
	 fmm_Sigma_only_noAC->SetParameter(6,fmm_noAC_Lexp->GetParameter(13));
//bg_only_noAC
	 fmm_bg_only_noAC->SetNpx(20000);
	 fmm_bg_only_noAC->SetParameter(0,fmm_noAC_Lexp->GetParameter(14));
	 fmm_bg_only_noAC->SetParameter(1,fmm_noAC_Lexp->GetParameter(15));
	 fmm_bg_only_noAC->SetParameter(2,fmm_noAC_Lexp->GetParameter(16));
	 fmm_bg_only_noAC->SetParameter(3,fmm_noAC_Lexp->GetParameter(17));
	 fmm_bg_only_noAC->SetParameter(4,fmm_noAC_Lexp->GetParameter(18));
	 fmm_bg_only_noAC->SetParameter(5,fmm_noAC_Lexp->GetParameter(19));
	 fmm_bg_only_noAC->SetParameter(6,fmm_noAC_Lexp->GetParameter(20));
	 fmm_bg_only_noAC->SetParameter(7,fmm_noAC_Lexp->GetParameter(21));

	double nofL_noAC = fmm_Lambda_only_noAC->Integral(fit_min_mm,fit_max_mm);
	double nofL_old_noAC = fmm_Lambda_only_noAC->Integral(-0.006,0.006);
	nofL_noAC = nofL_noAC/fit_bin_width;
	nofL_old_noAC = nofL_old_noAC/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL_noAC<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old_noAC<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_noAC->Integral(hmm_wo_bg_fom_noAC->FindBin(-0.006),hmm_wo_bg_fom_noAC->FindBin(0.006))<<endl;

	double nofS_noAC = fmm_Sigma_only_noAC->Integral(fit_min_mm,fit_max_mm);
	double nofS_old_noAC = fmm_Sigma_only_noAC->Integral(def_mean_S-0.008,def_mean_S+0.008);
	nofS_noAC = nofS_noAC/fit_bin_width;
	nofS_old_noAC = nofS_old_noAC/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS_noAC<<endl;
	cout<<"Number of Sigma w/o radiative tail (TF1 Integral) = "<<nofS_old_noAC<<endl;
	cout<<"Number of Sigma w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_noAC->Integral(hmm_wo_bg_fom_noAC->FindBin(def_mean_S-0.008),hmm_wo_bg_fom_noAC->FindBin(def_mean_S+0.008))<<endl;

	 hmm_wo_bg_fom_noAC->Draw();
	 fmm_noAC_Lexp->SetLineColor(kRed);
	 fmm_noAC_Lexp->Draw("same");
	 fmm_Lambda_only_noAC->SetLineColor(kAzure);
	 fmm_Sigma_only_noAC->SetLineColor(kCyan);
	 fmm_bg_only_noAC->SetLineColor(kOrange);
	 fmm_Lambda_only_noAC->Draw("same");
	 fmm_Sigma_only_noAC->Draw("same");
	 fmm_bg_only_noAC->Draw("same");
	 //fL_noAC->SetLineColor(kGreen);
	 //fS_noAC->SetLineColor(kGreen);
	 //fL_noAC->Draw("same");
	 //fS_noAC->Draw("same");
	 
	TCanvas* c7 = new TCanvas("c7","c7");
	 fAl->SetLineColor(kRed);
	 hmm_Al_wobg_fom_noZ->Draw("");
	 fAl->Draw("same");

	TCanvas* c8 = new TCanvas("c8","c8");
	 fAl_new->SetLineColor(kRed);
	 hmm_Al_wobg_fom_noZ_new->Draw("");
	 fAl_new->Draw("same");


	TCanvas* c9 = new TCanvas("c9","c9");
	 fpion->SetLineColor(kRed);
	 hmm_pi_wobg_fom_best->Draw("");
	 fpion->Draw("same");
	
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
	c9->Print(Form("%s]",pdfname.c_str()));
	 

cout << "Well done!" << endl;
}//fit
