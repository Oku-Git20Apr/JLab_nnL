//-- Radiative tail  --//
//comparison
//%Data
//%SIMC
//%Geant4
//
//K. Okuyama (Jan. 14, 2024)
//taken over from radiative_2023.C
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

void radiative_Q2(){
	//BOTH_LS.root === cell thickness ~ 400 um
	//BOTH_LS_cell.root === cell / sin(13.2deg);
	//BOTH_LS_cell_x10.root === cell * 10 / sin(13.2deg);
	
	//rad_Alupdate.root === Al cell modeling update
	//radL_Alupdate.root
	//radS_Alupdate.root
	
	//rad_2023.root === SIMC 4313MeV
	//radL_2023.root
	//radS_2023.root
  
  TFile *file = new TFile("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/h2all_2020Nov.root","read");
  //TFile *file_mea = new TFile("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/MixedEventAnalysis/bgmea_llccrr_new_2022effK.root","read");//2022/6/11
  TFile *file_mea = new TFile("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/MixedEventAnalysis/bgmea_llccrr_new_2023.root","read");//2023/10/6
  //TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/Suzuki_20201208/G4_temp/data/H2_500um_woInRad.root","read");//w/o internal radiation 2020/12/08
  TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/Suzuki_20201208/G4_temp/data/tree_H2_500um_wInRad.root","read");//w/ internal radiation 2020/12/08

//SIMC//
  TFile *file_simcL = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/radL_2023.root","read");
  TFile *file_simcS = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/radS_2023.root","read");
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_2023.root","read");
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
  TH2F* hmm_Q2_simcL  = new TH2F("hmm_Q2_simcL","hmm_Q2_simcL",xbin,xmin,xmax,50,0.3,0.7);
  hmm_Q2_simcL->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (MeV/c^{2})");
  hmm_Q2_simcL->GetYaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
  TH2F* hmm_Q2_simcS  = new TH2F("hmm_Q2_simcS","hmm_Q2_simcS",xbin,xmin,xmax,50,0.3,0.7);
  hmm_Q2_simcS->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (MeV/c^{2})");
  hmm_Q2_simcS->GetYaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
  TH2F* hmm_Q2_simc  = new TH2F("hmm_Q2_simc","hmm_Q2_simc",xbin,xmin,xmax,50,0.3,0.7);
  hmm_Q2_simc->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (MeV/c^{2})");
  hmm_Q2_simc->GetYaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	float mm_simcL, mm_simcS;
	float L_momL, L_momR;
	float S_momL, S_momR;
	float L_thL, L_thR;
	float L_phL, L_phR;
	float S_thL, S_thR;
	float S_phL, S_phR;
	float L_L_xfp, L_L_xpfp;
	float L_R_xfp, L_R_xpfp;
	float S_L_xfp, S_L_xpfp;
	float S_R_xfp, S_R_xpfp;

	SNTL->SetBranchStatus("*",0);
	SNTL->SetBranchStatus("missmass",1);SNTL->SetBranchAddress("missmass",&mm_simcL);
	SNTL->SetBranchStatus("Lp_orig",1);SNTL->SetBranchAddress("Lp_orig",&L_momL);
	SNTL->SetBranchStatus("Lth_orig",1);SNTL->SetBranchAddress("Lth_orig",&L_thL);
    SNTL->SetBranchStatus("Lph_orig"    ,1);SNTL->SetBranchAddress("Lph_orig"    ,&L_phL     );
	SNTL->SetBranchStatus("Rp_orig",1);SNTL->SetBranchAddress("Rp_orig",&L_momR);
	SNTL->SetBranchStatus("Rth_orig",1);SNTL->SetBranchAddress("Rth_orig",&L_thR);
    SNTL->SetBranchStatus("Rph_orig"    ,1);SNTL->SetBranchAddress("Rph_orig"    ,&L_phR     );
    SNTL->SetBranchStatus("h_xfp"  ,1);SNTL->SetBranchAddress("h_xfp"  ,&L_R_xfp );
    SNTL->SetBranchStatus("h_xpfp" ,1);SNTL->SetBranchAddress("h_xpfp" ,&L_R_xpfp);
    SNTL->SetBranchStatus("e_xfp"  ,1);SNTL->SetBranchAddress("e_xfp"  ,&L_L_xfp );
    SNTL->SetBranchStatus("e_xpfp" ,1);SNTL->SetBranchAddress("e_xpfp" ,&L_L_xpfp);
	SNTS->SetBranchStatus("*",0);
	SNTS->SetBranchStatus("missmass",1);SNTS->SetBranchAddress("missmass",&mm_simcS);
	SNTS->SetBranchStatus("Lp_orig",1);SNTS->SetBranchAddress("Lp_orig",&S_momL);
	SNTS->SetBranchStatus("Lth_orig",1);SNTS->SetBranchAddress("Lth_orig",&S_thL);
    SNTS->SetBranchStatus("Lph_orig"    ,1);SNTS->SetBranchAddress("Lph_orig"    ,&S_phL     );
	SNTS->SetBranchStatus("Rp_orig",1);SNTS->SetBranchAddress("Rp_orig",&S_momR);
	SNTS->SetBranchStatus("Rth_orig",1);SNTS->SetBranchAddress("Rth_orig",&S_thR);
    SNTS->SetBranchStatus("Rph_orig"    ,1);SNTS->SetBranchAddress("Rph_orig"    ,&S_phR     );
    SNTS->SetBranchStatus("h_xfp"  ,1);SNTS->SetBranchAddress("h_xfp"  ,&S_R_xfp );
    SNTS->SetBranchStatus("h_xpfp" ,1);SNTS->SetBranchAddress("h_xpfp" ,&S_R_xpfp);
    SNTS->SetBranchStatus("e_xfp"  ,1);SNTS->SetBranchAddress("e_xfp"  ,&S_L_xfp );
    SNTS->SetBranchStatus("e_xpfp" ,1);SNTS->SetBranchAddress("e_xpfp" ,&S_L_xpfp);


//SIMC Event Loop//
  int ENum_simcL = SNTL->GetEntries(); 
  int ENum_simcS = SNTS->GetEntries(); 

cout<<"Entries(SIMC Lambda): "<<ENum_simcL<<endl;
  TRandom3 ranL_simc;
  int incl = 0;
  for(int i=0;i<ENum_simcL;i++){
	SNTL->GetEntry(i);
	double ran = ranL_simc.Gaus((mm_simcL-ML),0.001);

//KINEMATICS CALCULATION
		L_momL/=1000.;//MeV-->GeV
  		L_momR/=1000.;//MeV-->GeV
  		double L_momB=4.313;//GeV
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
	//if(L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160&&ran<0.15)hmm_simcL->Fill(ran*1000.);//with mom cut

	//cout<<"theta_gk_cm="<<theta_gk_cm*180./PI<<endl;
	if(L_L_xpfp<0.17*L_L_xfp/100.+0.025&&L_L_xpfp>0.17*L_L_xfp/100.-0.035&&L_L_xpfp<0.40*L_L_xfp/100.+0.130&&L_L_xpfp>0.40*L_L_xfp/100.-0.130&&L_R_xpfp<0.17*L_R_xfp/100.+0.025&&L_R_xpfp>0.17*L_R_xfp/100.-0.035&&L_R_xpfp<0.40*L_R_xfp/100.+0.130&&L_R_xpfp>0.40*L_R_xfp/100.-0.130){
	//change
		if(L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160){
			hmm_simcL->Fill(ran*1000.);
			hmm_Q2_simcL->Fill(ran*1000.,Qsq);
			hmm_Q2_simc->Fill(ran*1000.,Qsq);
		//if(Qsq>=0.5&&L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160)hmm_simcL->Fill(ran*1000.);
		//if(theta_gk_cm*180./PI>=8.&&L_momR>1.760&&L_momR<1.900&&L_momL>2.010&&L_momL<2.160)hmm_simcL->Fill(ran*1000.);
	}//FP cut
}
}


cout<<"Entries(SIMC Sigma0): "<<ENum_simcS<<endl;
  TRandom3 ranS_simc;
  for(int i=0;i<ENum_simcS;i++){
	SNTS->GetEntry(i);
	double ran = ranS_simc.Gaus((mm_simcS-ML-0.001),0.001);

//KINEMATICS CALCULATION
		S_momL/=1000.;//MeV-->GeV
  		S_momR/=1000.;//MeV-->GeV
  		double S_momB=4.313;//GeV
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
	//if(S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160&&ran<0.15)hmm_simcS->Fill(ran*1000.);//with mom cut

	if(S_L_xpfp<0.17*S_L_xfp/100.+0.025&&S_L_xpfp>0.17*S_L_xfp/100.-0.035&&S_L_xpfp<0.40*S_L_xfp/100.+0.130&&S_L_xpfp>0.40*S_L_xfp/100.-0.130&&S_R_xpfp<0.17*S_R_xfp/100.+0.025&&S_R_xpfp>0.17*S_R_xfp/100.-0.035&&S_R_xpfp<0.40*S_R_xfp/100.+0.130&&S_R_xpfp>0.40*S_R_xfp/100.-0.130){
	//change
		if(S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160){
			hmm_simcS->Fill(ran*1000.);
			hmm_Q2_simcS->Fill(ran*1000.,Qsq);
			if(incl%3==1)hmm_Q2_simc->Fill(ran*1000.,Qsq);
			incl++;
		//if(Qsq>=0.5&&S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160)hmm_simcS->Fill(ran*1000.);
		//if(theta_gk_cm*180./PI>=8.&&S_momR>1.760&&S_momR<1.900&&S_momL>2.010&&S_momL<2.160)hmm_simcS->Fill(ran*1000.);
	}//FP cut
}
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
//-- Geant4 --//
TCanvas *cmm = new TCanvas("cmm","",800,800);
	hmm_simcL->Draw("");
	hmm_simcS->Draw("same");

TCanvas *cmmQ2L = new TCanvas("cmmQ2L","",800,800);
	hmm_Q2_simcL->Draw("colz");
TCanvas *cmmQ2S = new TCanvas("cmmQ2S","",800,800);
	hmm_Q2_simcS->Draw("colz");
TCanvas *cmmQ2 = new TCanvas("cmmQ2","",800,800);
	cmmQ2->SetLogz(1);
	hmm_Q2_simc->Draw("colz");
	 cmmQ2->SetLeftMargin(0.13);
	 cmmQ2->SetRightMargin(0.13);
	 cmmQ2->SetTopMargin(0.13);
	 cmmQ2->SetBottomMargin(0.13);
	 cmmQ2->Modified();
	 cmmQ2->Update();
	 gPad->Modified();
	 gPad->Update();
cmmQ2->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/dthesis_Fig/pdf/hmm_rad_Q2.pdf");

//c5->Print("SIMC_rad20220118.pdf");
//	 c5->SetLeftMargin(0.11);
//	 c5->SetRightMargin(0.11);
//	 c5->SetTopMargin(0.11);
//	 c5->SetBottomMargin(0.11);
//	 c5->Modified();
//	 c5->Update();
//	 gPad->Modified();
//	 gPad->Update();
//c5->Print("/data/41a/ELS/okuyama/JLab_nnL/okuya_macros/dthesis_Fig/pdf/hmm_rad_simc.pdf");

cout << "Well done!" << endl;
}
