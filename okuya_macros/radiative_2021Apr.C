//-- Radiative tail  --//
//comparison
//%Data
//%SIMC
//%Geant4
//
//K. Okuyama (April 27, 2021)
//
//taken over from radiative3.C
//-- after Aluminum cell modeling (April 9, 2021)
//
double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

double pol2gaus(double *x, double *par) {
   //par[0]=x^2
   //par[1]=x^1
   //par[2]=x^0
   //par[3]=Norm
   //par[4]=Width (sigma) of convoluted Gaussian function
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  double np = 500.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, mpc, fland, sum = 0.0, xlow,xupp, step, i;
  double val;

// Range of convolution integral
  xlow = x[0] - sc * par[4];
  xupp = x[0] + sc * par[4];
  step = (xupp-xlow) / np;
  for(i=1.0; i<=np/2; i++) {
     xx = xlow + (i-.5) * step;
     fland = TMath::Gaus(xx,x[0],par[4]);
if(xlow>-0.125&&xupp<0.125){
     sum += fland * (par[0]*x[0]*x[0]+par[1]*x[0]+par[2])*(1./(2.*par[0]*pow(0.125,3.)/3.+2.*par[2]*0.125));
}else sum += fland;

     xx = xupp - (i-.5) * step;
     fland = TMath::Gaus(xx,x[0],par[4]);
if(xlow>-0.125&&xupp<0.125){
     sum += fland * (par[0]*x[0]*x[0]+par[1]*x[0]+par[2])*(1./(2.*par[0]*pow(0.125,3.)/3.+2.*par[2]*0.125));
}else sum += fland;
  }
  val = par[3] * step * sum * invsq2pi / par[4];

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

double F_VZ( double *x, double *par )
{
  return pol2gaus(x,par)+par[5]*TMath::Gaus(x[0],par[6],par[7],1)+par[8]*TMath::Gaus(x[0],par[9],par[10],1)+par[11]*TMath::Gaus(x[0],par[12],par[13],1)+par[14]*TMath::Gaus(x[0],par[15],par[16],1);
}

double F_VZ_cell( double *x, double *par)
{
double a1, a2;
double b1, b2;
a1=par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
a2=par[0]*par[3]*exp(-0.5*pow((x[0]-par[1]+par[4])/par[5],2.));
b1=par[6]*exp(-0.5*pow((x[0]-par[7])/par[2],2.));
b2=par[6]*par[3]*exp(-0.5*pow((x[0]-par[7]+par[4])/par[5],2.));
double a = a1+a2;
double b = b1+b2;
return a+b;
}

double F_VZ2( double *x, double *par)
{
//par[0]: Al front scale
//par[1]: Al front pos.
//par[2]: Al Gauss sigma
//par[3]: Al second gauss strength
//par[4]: Al(front) second gauss pos. 
//par[5]: Al second gauss sigma 
//par[6]: Al rear scale
//par[7]: Al(rear) second gauss pos. 
//par[8]: pol2 coeff.1
//par[9]: pol2 coeff.2
//par[10]: total scale 

double a1, a2;
double b1, b2;
a1=par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
a2=par[0]*par[3]*exp(-0.5*pow((x[0]-par[1]+par[4])/par[5],2.));
b1=par[6]*exp(-0.5*pow((x[0]-par[7])/par[2],2.));
b2=par[6]*par[3]*exp(-0.5*pow((x[0]-par[7]+par[4])/par[5],2.));
double a = a1+a2;
double b = b1+b2;

double c = 0.;
int np = 2000;
for(int i=0;i<np;i++){
 double d = par[8]*pow((x[0]-par[9]),2.)+1.;
 double step = -1.+(double)i/1000.;
 if(step<-0.125) d = 0.;
 if(step> 0.125) d = 0.;

 double aa;
 aa = exp(-0.5*pow((x[0]-step)/par[2],2.))+par[3]*exp(-0.5*pow((x[0]+par[4]-step)/par[5],2.));
 c = c + d*aa;
}
c = par[10]*c;
return a+b+c;
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

void radiative_2021Apr(){
	//BOTH_LS.root === cell thickness ~ 400 um
	//BOTH_LS_cell.root === cell / sin(13.2deg);
	//BOTH_LS_cell_x10.root === cell * 10 / sin(13.2deg);
	
	//rad_Alupdate.root === Al cell modeling update
	//radL_Alupdate.root
	//radS_Alupdate.root
  
  TFile *file = new TFile("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/h2all_2020Nov.root","read");//w/o internal radiation 2020/12/08
  TFile *file_mea = new TFile("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/MixedEventAnalysis/bgmea_2021Jan.root","read");//w/o internal radiation 2020/12/08
  //TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/Suzuki_20201208/G4_temp/data/H2_500um_woInRad.root","read");//w/o internal radiation 2020/12/08
  TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/Suzuki_20201208/G4_temp/data/tree_H2_500um_wInRad.root","read");//w/ internal radiation 2020/12/08

//SIMC//
  TFile *file_simcL = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/radL_Alupdate_double.root","read");
  TFile *file_simcS = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/radS_Alupdate_double.root","read");
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_Alupdate.root","read");
  TTree *tree = (TTree*)file->Get("tree_out");//input
  TTree *SNT = (TTree*)file_simc->Get("SNT");//SIMC all
  TTree *SNTL = (TTree*)file_simcL->Get("SNT");//SIMC Lambda
  TTree *SNTS = (TTree*)file_simcS->Get("SNT");//SIMC Sigma0
  TTree *tree_G4L = (TTree*)file_G4->Get("tree0_12_0");//Lambda
  TTree *tree_G4S = (TTree*)file_G4->Get("tree0_12_1");//Sigma0


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




//-- DATA --//
  TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  h1->GetXaxis()->SetTitle("coin time (ns)");
  h1->GetYaxis()->SetTitle("Counts / 100 ps");
  h1->GetXaxis()->SetRangeUser(-14.0,17.);
  double xmin = -100., xmax = 200.; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_L_fom_best  = new TH1F("hmm_L_fom_best","hmm_L_fom_best",xbin,xmin,xmax);
  hmm_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_L_fom_best->GetYaxis()->SetTitle("Counts / MeV");
  hmm_L_fom_best->SetLineColor(1);
  TH1F* hcs_L_fom_best  = new TH1F("hcs_L_fom_best","hcs_L_fom_best",xbin,xmin,xmax);
  hcs_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} [GeV/c^{2}]");
  hcs_L_fom_best->GetYaxis()->SetTitle("d#sigma/d#Omega (C.M.F.) [nb/sr]");
  hcs_L_fom_best->SetLineColor(1);
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wobg_fom_best  = new TH1F("hmm_wobg_fom_best","hmm_wobg_fom_best",xbin,xmin,xmax);
  hmm_wobg_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (MeV/c^{2})");
  hmm_wobg_fom_best->GetYaxis()->SetTitle("Counts/(MeV/c^{2})");
  hmm_wobg_fom_best->SetLineColor(kBlack);



//-- SIMC --//
  TH1F* hmm_simcL  = new TH1F("hmm_simcL","hmm_simcL",xbin,xmin,xmax);
  TH1F* hmm_simcS  = new TH1F("hmm_simcS","hmm_simcS",xbin,xmin,xmax);
  TH1F* hmm_simc  = new TH1F("hmm_simc","hmm_simc",xbin,xmin,xmax);
	float mm_simcL, mm_simcS;
	float L_momL, L_momR;
	float S_momL, S_momR;
	float L_thL, L_thR;
	float L_phL, L_phR;
	float S_thL, S_thR;
	float S_phL, S_phR;
	SNTL->SetBranchStatus("*",0);
	SNTL->SetBranchStatus("missmass",1);SNTL->SetBranchAddress("missmass",&mm_simcL);
	//SNTL->SetBranchStatus("Lp_rec",1);SNTL->SetBranchAddress("Lp_rec",&L_momL);
	//SNTL->SetBranchStatus("e_xptar",1);SNTL->SetBranchAddress("e_xptar",&L_thL);
    //SNTL->SetBranchStatus("e_yptar"    ,1);SNTL->SetBranchAddress("e_yptar"    ,&L_phL     );
	//SNTL->SetBranchStatus("Rp_rec",1);SNTL->SetBranchAddress("Rp_rec",&L_momR);
	//SNTL->SetBranchStatus("h_xptar",1);SNTL->SetBranchAddress("h_xptar",&L_thR);
    //SNTL->SetBranchStatus("h_yptar"    ,1);SNTL->SetBranchAddress("h_yptar"    ,&L_phR     );
	SNTL->SetBranchStatus("Lp_orig",1);SNTL->SetBranchAddress("Lp_orig",&L_momL);
	SNTL->SetBranchStatus("Lth_gen",1);SNTL->SetBranchAddress("Lth_gen",&L_thL);
    SNTL->SetBranchStatus("Lph_gen"    ,1);SNTL->SetBranchAddress("Lph_gen"    ,&L_phL     );
	SNTL->SetBranchStatus("Rp_orig",1);SNTL->SetBranchAddress("Rp_orig",&L_momR);
	SNTL->SetBranchStatus("Rth_gen",1);SNTL->SetBranchAddress("Rth_gen",&L_thR);
    SNTL->SetBranchStatus("Rph_gen"    ,1);SNTL->SetBranchAddress("Rph_gen"    ,&L_phR     );
	SNTS->SetBranchStatus("*",0);
	SNTS->SetBranchStatus("missmass",1);SNTS->SetBranchAddress("missmass",&mm_simcS);
	SNTS->SetBranchStatus("Lp_orig",1);SNTS->SetBranchAddress("Lp_orig",&S_momL);
	SNTS->SetBranchStatus("Lth_gen",1);SNTS->SetBranchAddress("Lth_gen",&S_thL);
    SNTS->SetBranchStatus("Lph_gen"    ,1);SNTS->SetBranchAddress("Lph_gen"    ,&S_phL     );
	SNTS->SetBranchStatus("Rp_orig",1);SNTS->SetBranchAddress("Rp_orig",&S_momR);
	SNTS->SetBranchStatus("Rth_gen",1);SNTS->SetBranchAddress("Rth_gen",&S_thR);
    SNTS->SetBranchStatus("Rph_gen"    ,1);SNTS->SetBranchAddress("Rph_gen"    ,&S_phR     );


//SIMC Event Loop//
  int ENum_simcL = SNTL->GetEntries(); 
  int ENum_simcS = SNTS->GetEntries(); 

cout<<"Entries(SIMC Lambda): "<<ENum_simcL<<endl;
  TRandom3 ranL_simc;
  for(int i=0;i<ENum_simcL;i++){
	SNTL->GetEntry(i);
	double ran = ranL_simc.Gaus((mm_simcL-ML),0.001);

//KINEMATICS CALCULATION
		L_momL/=1000.;//MeV-->GeV
  		L_momR/=1000.;//MeV-->GeV
  		double L_momB=4.318;//GeV
  		double mh = ML;//hypernuclei
  		double mt = Mp;//target mass

	    double R_pz = L_momR/sqrt(1.0*1.0 + pow(L_thR, 2.0) + pow( L_phR,2.0));
	    double R_px = R_pz * ( L_thR );
	    double R_py = R_pz * ( L_phR );

	    double L_pz = L_momL/sqrt(1.0*1.0 + pow(L_thL, 2.0) + pow(L_phL,2.0));
	    double L_px = L_pz * ( L_thL );
	    double L_py = L_pz * ( L_phL );


	    double B_E =sqrt(L_momB*L_momB + Me*Me);
	    double R_E =sqrt(L_momR*L_momR + MK*MK);
	    double L_E =sqrt(L_momL*L_momL + Me*Me);


		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,L_momB, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);




	    double pL    = L_momL;//MeV
	    double pR    = L_momR;//MeV
		double theta = L_thL;
		double theta_R = L_thR;
		double phi = L_phL;
		double phi_R = L_phR;
		double phi0=13.2*PI/180;//rad
		double phi_L = L_4vec.Phi();//LHRS frame
		double phi_RHRS = R_4vec.Phi();//RHRS frame
	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );
        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML

		double theta_ee = L_4vec.Theta();
		//double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double phi_ee = L_4vec.Phi();//original frame
//cout<<"phi_ee="<<phi_ee<<endl;
		double phi_ek = R_4vec.Phi()+2*PI;//original frame
//cout<<"phi_ek="<<phi_ek<<endl;

		G_4vec = B_4vec - L_4vec;
		double mom_g=sqrt(G_4vec.Px()*G_4vec.Px()+G_4vec.Py()*G_4vec.Py()+G_4vec.Pz()*G_4vec.Pz());
		double Qsq = G_4vec.M()*G_4vec.M();
		double phi_g = G_4vec.Phi()+2*PI;
		double theta_g = G_4vec.Theta();
//cout<<"theta_gk="<<(theta_ek-theta_g)*180./PI<<endl;
		double theta_gk_lab = G_4vec.Angle(R_4vec.Vect());
//cout<<"theta_gk_lab(TLorentz)="<<theta_gk_lab*180./PI<<endl;
		double omega=G_4vec.E();
		double pY = sqrt((omega+Mp-sqrt(pR*pR+MK*MK))*(omega+Mp-sqrt(pR*pR+MK*MK))-ML*ML);
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		//double theta_gk_lab_test = acos((mom_g*mom_g+pR*pR-pY*pY)/(2.*pR*mom_g));
		double theta_eg_lab = acos((B_E-L_E*(1.-Qsq/2./B_E/L_E))/mom_g);
//cout<<"theta_eg_lab="<<theta_eg_lab*180./PI<<endl;
		double theta_gk_lab_test = 13.2*PI/180.-theta_eg_lab;
//cout<<"theta_gk_lab="<<theta_gk_lab_test*180./PI<<endl;
		double beta=mom_g/(omega+Mp);
		double ER=sqrt(pR*pR+MK*MK);
		double gamma=1./sqrt(1-beta*beta);
		double theta_gk_cm_test = atan((pR*sin(theta_gk_lab_test))/(-1.*gamma*beta*ER+gamma*pR*cos(theta_gk_lab_test)));
//cout<<"theta_gk_cm="<<theta_gk_cm_test*180./PI<<endl;
	
		TVector3 boost;
		TLorentzVector GT_4vec;
		GT_4vec=G_4vec+T_4vec;
		boost=GT_4vec.BoostVector();
		R_4vec.Boost(-boost);
		L_4vec.Boost(-boost);
		B_4vec.Boost(-boost);
		double theta_gk_cm = G_4vec.Angle(R_4vec.Vect());


//FILL into histograms//

	//if(L_momR>1760.&&L_momR<1900.&&L_momL>2092.&&L_momL<2160.)hmm_simcL->Fill(ran*1000.);
	if(L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160&&ran<0.15)hmm_simcL->Fill(ran*1000.);//with mom cut

	//cout<<"theta_gk_cm="<<theta_gk_cm*180./PI<<endl;
	//change
	//if(theta_gk_cm*180./PI>=8.&&L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160)hmm_simcL->Fill(ran*1000.);
	//if(Qsq>=0.5&&L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160)hmm_simcL->Fill(ran*1000.);
	}



cout<<"Entries(SIMC Sigma0): "<<ENum_simcS<<endl;
  TRandom3 ranS_simc;
  for(int i=0;i<ENum_simcS;i++){
	SNTS->GetEntry(i);
	double ran = ranS_simc.Gaus((mm_simcS-ML-0.001),0.001);

//KINEMATICS CALCULATION
		S_momL/=1000.;//MeV-->GeV
  		S_momR/=1000.;//MeV-->GeV
  		double S_momB=4.318;//GeV
  		double mh = ML;//hypernuclei
  		double mt = Mp;//target mass
	    double R_pz = S_momR/sqrt(1.0*1.0 + pow(S_thR, 2.0) + pow( S_phR,2.0));
	    double R_px = R_pz * ( S_thR );
	    double R_py = R_pz * ( S_phR );

	    double L_pz = S_momL/sqrt(1.0*1.0 + pow(S_thL, 2.0) + pow(S_phL,2.0));
	    double L_px = L_pz * ( S_thL );
	    double L_py = L_pz * ( S_phL );


	    double B_E =sqrt(S_momB*S_momB + Me*Me);
	    double R_E =sqrt(S_momR*S_momR + MK*MK);
	    double L_E =sqrt(S_momL*S_momL + Me*Me);


		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,S_momB, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);




	    double pL    = S_momL;//GeV
	    double pR    = S_momR;//GeV
		double theta = S_thL;
		double theta_R = S_thR;
		double phi = S_phL;
		double phi_R = S_phR;
		double phi0=13.2*PI/180;//rad
		double phi_L = L_4vec.Phi();//LHRS frame
		double phi_RHRS = R_4vec.Phi();//RHRS frame
	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );
        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML

		double theta_ee = L_4vec.Theta();
		//double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double phi_ee = L_4vec.Phi();//original frame
//cout<<"phi_ee="<<phi_ee<<endl;
		double phi_ek = R_4vec.Phi()+2*PI;//original frame
//cout<<"phi_ek="<<phi_ek<<endl;

		G_4vec = B_4vec - L_4vec;
		double mom_g=sqrt(G_4vec.Px()*G_4vec.Px()+G_4vec.Py()*G_4vec.Py()+G_4vec.Pz()*G_4vec.Pz());
		double Qsq = G_4vec.M()*G_4vec.M();
		double phi_g = G_4vec.Phi()+2*PI;
		double theta_g = G_4vec.Theta();
//cout<<"theta_gk="<<(theta_ek-theta_g)*180./PI<<endl;
		double theta_gk_lab = G_4vec.Angle(R_4vec.Vect());
//cout<<"theta_gk_lab(TLorentz)="<<theta_gk_lab*180./PI<<endl;
		double omega=G_4vec.E();
		double pY = sqrt((omega+Mp-sqrt(pR*pR+MK*MK))*(omega+Mp-sqrt(pR*pR+MK*MK))-ML*ML);
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		//double theta_gk_lab_test = acos((mom_g*mom_g+pR*pR-pY*pY)/(2.*pR*mom_g));
		double theta_eg_lab = acos((B_E-L_E*(1.-Qsq/2./B_E/L_E))/mom_g);
//cout<<"theta_eg_lab="<<theta_eg_lab*180./PI<<endl;
		double theta_gk_lab_test = 13.2*PI/180.-theta_eg_lab;
//cout<<"theta_gk_lab="<<theta_gk_lab_test*180./PI<<endl;
		double beta=mom_g/(omega+Mp);
		double ER=sqrt(pR*pR+MK*MK);
		double gamma=1./sqrt(1-beta*beta);
		double theta_gk_cm_test = atan((pR*sin(theta_gk_lab_test))/(-1.*gamma*beta*ER+gamma*pR*cos(theta_gk_lab_test)));
//cout<<"theta_gk_cm="<<theta_gk_cm_test*180./PI<<endl;
	
		TVector3 boost;
		TLorentzVector GT_4vec;
		GT_4vec=G_4vec+T_4vec;
		boost=GT_4vec.BoostVector();
		R_4vec.Boost(-boost);
		L_4vec.Boost(-boost);
		B_4vec.Boost(-boost);
		double theta_gk_cm = G_4vec.Angle(R_4vec.Vect());

//FILL into histograms//

	//if(S_momR>1760.&&S_momR<1900.&&S_momL>2010.&&S_momL<2108.)hmm_simcS->Fill(ran*1000.);
	if(S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160&&ran<0.15)hmm_simcS->Fill(ran*1000.);//with mom cut

	//change
	//if(theta_gk_cm*180./PI>=8.&&S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160)hmm_simcS->Fill(ran*1000.);
	//if(Qsq>=0.5&&S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160)hmm_simcS->Fill(ran*1000.);
	}

//  TH1F* hmm_simc  = new TH1F("hmm_simc","hmm_simc",xbin,xmin,xmax);
//    //char condi[1000];
//	//SNT->Project("hmm_simc", "missmass*1000.");
//	//SNT->Project("hmm_simc", "missmass");
//	float mm_simc;
//	SNT->SetBranchStatus("*",0);
//	SNT->SetBranchStatus("missmass",1);SNT->SetBranchAddress("missmass",&mm_simc);
//  int ENum_simc = SNT->GetEntries(); 
//cout<<"Entries(SIMC): "<<ENum_simc<<endl;
//  TRandom3 ran_simc;
//  for(int l=0;l<ENum_simc;l++){
//	SNT->GetEntry(l);
//	double ran = ran_simc.Gaus((mm_simc-ML),0.001);
//	hmm_simc->Fill(ran*1000.);
//	}

//-- Geant4 --//
  TH1F* hmm_G4L  = new TH1F("hmm_G4L","hmm_G4L",xbin,xmin,xmax);
  TH1F* hmm_G4S  = new TH1F("hmm_G4S","hmm_G4S",xbin,xmin,xmax);
  TH1F* hmm_G4  = new TH1F("hmm_G4","hmm_G4",xbin,xmin,xmax);
	double mm_G4L, mm_G4S;
	tree_G4L->SetBranchStatus("*",0);
	tree_G4L->SetBranchStatus("MM1",1);tree_G4L->SetBranchAddress("MM1",&mm_G4L);
	tree_G4S->SetBranchStatus("*",0);
	tree_G4S->SetBranchStatus("MM1",1);tree_G4S->SetBranchAddress("MM1",&mm_G4S);
  int ENum_G4L = tree_G4L->GetEntries(); 
  int ENum_G4S = tree_G4S->GetEntries(); 
cout<<"Entries(G4 Lambda): "<<ENum_G4L<<endl;
  TRandom3 ranL_G4;
  for(int i=0;i<ENum_G4L;i++){
	tree_G4L->GetEntry(i);
	double ran = ranL_G4.Gaus((mm_G4L-ML-0.0015),0.001);
	hmm_G4L->Fill(ran*1000.);
	}
cout<<"Entries(G4 Sigma0): "<<ENum_G4S<<endl;
  TRandom3 ranS_G4;
  for(int i=0;i<ENum_G4S;i++){
	tree_G4S->GetEntry(i);
	double ran = ranS_G4.Gaus((mm_G4S-ML-0.0015),0.001);
	hmm_G4S->Fill(ran*1000.);
	}


  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");
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
  TH2F* h_pepk  = new TH2F("h_pepk","h_pepk (tight)",40,1.73,1.93,40,1.95,2.25);
  TH1F* h_ac1_sum  = new TH1F("h_ac1_sum","",1000,-1.,40.);
  h_ac1_sum->GetXaxis()->SetTitle("NPE (AC1)");
  h_ac1_sum->GetYaxis()->SetTitle("Counts");
  h_ac1_sum->SetLineColor(kAzure);
  TH1F* h_ac2_sum  = new TH1F("h_ac2_sum","",1000,-1.,80.);
  h_ac2_sum->GetXaxis()->SetTitle("NPE (AC2)");
  h_ac2_sum->GetYaxis()->SetTitle("Counts");
  h_ac2_sum->SetLineColor(kAzure);
  TH1F* h_coin_nocut  = new TH1F("h_coin_nocut","",40.*1000./56.,-20.,20.);
  h_coin_nocut->GetXaxis()->SetTitle("Coincidence Time [ns]");
  h_coin_nocut->GetYaxis()->SetTitle("Counts");
  h_coin_nocut->SetLineColor(kAzure);
  TH1F* h_coin_strict  = new TH1F("h_coin_strict","",40.*1000./56.,-20.,20.);
  h_coin_strict->GetXaxis()->SetTitle("Coincidence Time [ns]");
  h_coin_strict->GetYaxis()->SetTitle("Counts");
  h_coin_strict->SetLineColor(kAzure);
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",400,-25.*0.01,25.*0.01,400,-25.*0.01,25.*0.01);
  TH2F* h_zz    = new TH2F("h_zz"   ,""   , 1000,-25.,25.,1000,  -25., 25.); 
  h_zz->GetXaxis()->SetTitle("Z-vertex (HRS-R) [cm]");
  h_zz->GetYaxis()->SetTitle("Z-vertex (HRS-L) [cm]");
  h_zz->GetXaxis()->SetTitleColor(kGreen+2);
  h_zz->GetYaxis()->SetTitleColor(kRed);
  h_zz->GetZaxis()->SetLabelOffset(-0.005);
  TH1F* h_zL  = new TH1F("h_zL","",1000.,-15.,15.);
  h_zL->GetXaxis()->SetTitle("Z-vertex (HRS-L) [cm]");
  h_zL->GetXaxis()->SetTitleColor(kRed);
  h_zL->GetXaxis()->SetTitleSize(0.07);
  h_zL->GetXaxis()->SetTitleOffset(0.6);
  h_zL->GetYaxis()->SetTitle("Counts");
  h_zL->GetYaxis()->SetTitleSize(0.07);
  h_zL->GetYaxis()->SetTitleOffset(0.3);
  h_zL->SetLineColor(kAzure);
  TH1F* h_zR  = new TH1F("h_zR","",1000.,-15.,15.);
  h_zR->GetXaxis()->SetTitle("Z-vertex (HRS-R) [cm]");
  h_zR->GetXaxis()->SetTitleColor(kGreen+2);
  h_zR->GetXaxis()->SetTitleSize(0.07);
  h_zR->GetXaxis()->SetTitleOffset(0.6);
  h_zR->GetYaxis()->SetTitle("Counts");
  h_zR->GetYaxis()->SetTitleSize(0.07);
  h_zR->GetYaxis()->SetTitleOffset(0.3);
  h_zR->SetLineColor(kAzure);

  TH1F* h_zave  = new TH1F("h_zave","Z-vertex (Ave.)",1000,-0.25,0.25);
  TH1F* h_zave_dummy  = new TH1F("h_zave_dummy","Z-vertex (Ave.)",1000,-0.15,0.15);
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
		double NA = 6.023*pow(10.,23.);
		double density_h2=70.8;//[mg/cm^2]
		double ntar_h2 = (0.001*density_h2/2.)*NA*2.;//[/cm^2]
		double effAC = 0.750;
		double effZ  = 0.789;
		double effFP = 1.000;
		double effch2= 0.999;
		double effct = 0.966;
		double effDAQ= 0.960;
		double efftr = 0.810;
		double effK	 = 0.170;
		double efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		double RHRS  = 0.006;
		//double LHRS  = 0.006;
		double Charge= 4.6486;//[C]
		double ee	 = 1.602*pow(10.,-19);
		double Ne	 = Charge/ee;//Num of e
		//double Ng	 = 2.48*pow(10.,-6)*Ne;//Geant4
		double Ng	 = 1.94*pow(10.,-6)*Ne;//SIMC
		double cs	 = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
cout<<"Efficiency="<<efficiency<<endl;
cout<<"Ntar(H2)="<<ntar_h2<<endl;
cout<<"Ne="<<Ne<<endl;
cout<<"Ng="<<Ng<<endl;
cout<<"cs="<<cs<<endl;
	//double csL[xbin];




  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  //tree->Draw(">>elist" , "fabs(ct_orig)<1.006");//ctsum (does NOT dintinguish #track)
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


	


		if(fabs(ct)<1.006)ct_cut=true;
		else ct_cut=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.760&&R_mom<1.900&&L_mom>2.010&&L_mom<2.160)event_selection=true;
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

		//change
		if(event_selection&&ct_cut)hmm_L_fom_best->Fill(mm*1000.);
		//if(event_selection&&ct_cut&&theta_gk_cm*180./PI>=8.)hmm_L_fom_best->Fill(mm*1000.);
		//if(event_selection&&ct_cut&&Qsq>=0.5)hmm_L_fom_best->Fill(mm*1000.);
		if(event_selection_nocut&&ct_cut)hmm_L_fom_nocut->Fill(mm);

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
			if(theta_gk_cm*180./PI>=8.)h_pepk->Fill(R_mom,L_mom);
		}
		if(R_Tr&&R_FP&&L_Tr&&L_FP){
		h_zz->Fill(R_tr_vz*100.,L_tr_vz*100.);
		h_zL->Fill(L_tr_vz*100.);
		h_zR->Fill(R_tr_vz*100.);
		}
		//if(R_Tr&&R_FP&&L_Tr&&L_FP&&abs(R_tr_vz-L_tr_vz)<0.025&&abs(R_tr_vz+L_tr_vz)<0.1)
		if(R_Tr&&R_FP&&L_Tr&&L_FP){
		h_coin_nocut->Fill(ct);
		}
		if(event_selection)h_coin_strict->Fill(ct);
		h_ac1_sum->Fill(ac1sum);
		h_ac2_sum->Fill(ac2sum);
		//if(abs(R_tr_vz-L_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		if(abs(R_tr_vz-L_tr_vz)<0.025)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		h_nltrack->Fill(NLtr);
		h_nrtrack->Fill(NRtr);

}//ENum
	//THStack *hs = (THStack*)file_G4->Get("new_mm1stack0_12");
	TH1F* hmm_bg_temp = (TH1F*)file_mea->Get("hmm_mixacc_result_new");
	//change
	//TH1F* hmm_bg_temp = (TH1F*)file_mea->Get("hmm_mixacc_result_new_cm2_2");
	//TH1F* hmm_bg_temp = (TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq2_2");

	for(int i=0;i<300;i++){
		double temp = hmm_bg_temp->GetBinContent(i+1);
		hmm_bg_fom_best->SetBinContent(i+1,temp);
	}
	hmm_bg_fom_best->Scale(1./6000.);
	hmm_wobg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.0,-1.0);

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	//hmm_G4L->Scale(809./3108.);
	//hmm_G4S->Scale(258./442.);
	hmm_G4L->Scale(1140.53/4270.);
	hmm_G4S->Scale(349.678/565.43);
	hmm_G4->Add(hmm_G4L,hmm_G4S,1.0,1.0);
	hmm_G4->SetLineColor(kGreen);
	hmm_wobg_fom_best->Draw("");
	hmm_G4->Draw("same");
cout<<"Data vs Geant4 (Lambda)"<<endl;
cout<<"hmm_G4: "<<hmm_G4->Integral(hmm_G4->FindBin(-6),hmm_G4->FindBin(6.))<<endl;
cout<<"hmm_L(data): "<<hmm_wobg_fom_best->Integral(hmm_wobg_fom_best->FindBin(-6),hmm_wobg_fom_best->FindBin(6.))<<endl;
cout<<"Data vs Geant4 (Sigma0)"<<endl;
cout<<"hmm_G4: "<<hmm_G4->Integral(hmm_G4->FindBin(def_mean_S*1000.-6),hmm_G4->FindBin(def_mean_S*1000.+6.))<<endl;
cout<<"hmm_L(data): "<<hmm_wobg_fom_best->Integral(hmm_wobg_fom_best->FindBin(def_mean_S*1000.-6),hmm_wobg_fom_best->FindBin(def_mean_S*1000.+6.))<<endl;
	////hs->GetHists()->At(0)->Draw("hist")->Integral();
	//TH1F* h_L = (TH1F*)hs->GetHists()->At(0);
	//double NL_G4 = h_L->Integral();
	//h_L->Draw();
	//TCanvas *c2 = new TCanvas("c2","c2",800,800);
	////hs->GetHists()->At(1)->Draw("hist");
	//TH1F* h_S = (TH1F*)hs->GetHists()->At(1);
	//double NS_G4 = h_S->Integral();
	//h_S->Draw();
	//cout<<"NL_G4="<<NL_G4<<endl;
	//cout<<"NS_G4="<<NS_G4<<endl;
	//h_L->Scale(220./500.);
	//h_S->Scale(220./500.);
	////h_L->Scale(1988./NL_G4);
	////h_S->Scale(739./NL_G4);
	//TCanvas* c3 = new TCanvas("c3","c3");
	////hs->Draw("hist");
	//h_L->Draw("");
	//h_S->Draw("same");
	//TCanvas* c4 = new TCanvas("c4","c4");
	//h_L->Draw("");
	//h_S->Draw("same");
	//hmm_wobg_fom_best->SetLineColor(kRed);
	//hmm_wobg_fom_best->SetFillColor(kRed);
	//hmm_wobg_fom_best->SetFillStyle(3004);
	//hmm_wobg_fom_best->Draw("same");
			//h_coin_nocut->Draw("");
			//cout<<"Coin_all_ENum="<<h_coin_nocut->Integral(h_coin_nocut->FindBin(-20.),h_coin_nocut->FindBin(20.))<<endl;
			//TCanvas* c2 = new TCanvas("c2","c2");
			//h_coin_strict->Draw("");
			//cout<<"Coin_kaon_ENum="<<h_coin_strict->Integral(h_coin_strict->FindBin(-1.006),h_coin_strict->FindBin(1.006))<<endl;
			//
			//c1->Print("./pdf/cointime_nocut.pdf");
			//c2->Print("./pdf/cointime_strict.pdf");
	TCanvas* c5 = new TCanvas("c5","c5",800,800);

	//hmm_simcL->Scale(833.828*833.828/357870./833.851);//2021Apr.//Full
	//hmm_simcS->Scale(283.282*283.282/205538./294.881);//2021Apr.//Full
	//hmm_simcL->Scale(833.828/351142.);//2021Apr.//MAX
	//hmm_simcS->Scale(283.282*283.282/204959./295.921);//2021Apr.//MAX
	hmm_simcL->Scale(833.828/330870.);//2021Apr.//MAX
	hmm_simcS->Scale(283.282*283.282/202332./299.470);//2021Apr.//MAX

	hmm_simc->Add(hmm_simcL,hmm_simcS,1.0,1.0);
	hmm_simc->SetLineColor(kRed);

	hmm_wobg_fom_best->Draw("");
	hmm_simc->SetLineWidth(4);
	hmm_simc->Draw("Psame");
cout<<"Data vs SIMC (Lambda)"<<endl;
cout<<"hmm_simc: "<<hmm_simc->Integral(hmm_simc->FindBin(-6),hmm_simc->FindBin(6.))<<endl;
cout<<"hmm_L(data): "<<hmm_wobg_fom_best->Integral(hmm_wobg_fom_best->FindBin(-6),hmm_wobg_fom_best->FindBin(6.))<<endl;
cout<<"Data vs SIMC (Sigma0)"<<endl;
cout<<"hmm_simc: "<<hmm_simc->Integral(hmm_simc->FindBin(def_mean_S*1000.-6),hmm_simc->FindBin(def_mean_S*1000.+6.))<<endl;
cout<<"hmm_L(data): "<<hmm_wobg_fom_best->Integral(hmm_wobg_fom_best->FindBin(def_mean_S*1000.-6),hmm_wobg_fom_best->FindBin(def_mean_S*1000.+6.))<<endl;

cout<<"Total"<<endl;
cout<<"hmm_G4L: "<<hmm_G4L->Integral()<<endl;
cout<<"hmm_G4S: "<<hmm_G4S->Integral()<<endl;
cout<<"hmm_simcL: "<<hmm_simcL->Integral()<<endl;
cout<<"hmm_simcS: "<<hmm_simcS->Integral()<<endl;
cout<<"hmm_data: "<<hmm_wobg_fom_best->Integral()<<endl;
cout<<"hmm_G4: "<<hmm_G4->Integral()<<endl;
cout<<"hmm_simc: "<<hmm_simc->Integral()<<endl;
cout << "Well done!" << endl;
}
