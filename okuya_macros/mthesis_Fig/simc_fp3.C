//--  Data vs SIMC (FP)  --//
//
//K. Okuyama (Feb. 21, 2021)
//
//This is taken over from data_vs_simc.C
//No array branch mode 
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

void simc_fp2(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_cell.root","read");// L:S0=3:1 (0.9M vs 0.3M) 2020/12/10
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/RHRS_grnd.root","read");// L:S0=3:1 (0.9M vs 0.3M) 2020/12/10
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_datafit.root","read");// L:S0=3:1 (0.9M vs 0.3M), pe'=2102.5MeV/c, pK=1825MeV/c,  2020/12/29
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/NONE_LS.root","read");// L:S0=3:1 (0.9M vs 0.3M)  2020/1/5
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_momminus20.root","read");// L:S0=3:1 (0.9M vs 0.3M), pe'=2080MeV/c, pK=1800MeV/c,  2020/12/29
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_500um.root","read");// L:S0=3:1 (0.9M vs 0.3M), Al thickiness = 500 um, 2020/12/29
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_ElossCor.root","read");// L:S0=3:1 (0.9M vs 0.3M), w/ Eloss Correction, 2020/12/10
  //TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/BOTH_LS_Rad3.root","read");// L:S0=3:1 (0.9M vs 0.3M), w/ extrad_flag = 3(Friedrich approximation), 2020/12/10
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");
	//ACCBGの引き算はmea_hist.ccから
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_2020Nov.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TFile *file_mea_mthesis = new TFile("../MixedEventAnalysis/bgmea_mthesis_2020Nov.root","read");//MeV
  double nbunch = 6000.;//effetive bunches (6 bunches x 5 mixtures)
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");


//---  DAQ Efficiency ---//
//H2 run (run111157~111222 & run111480~542)
	string daq_file = "../information/daq.dat";//DAQ Efficiency from ELOG
	int runnum;
	double daq_eff;
	double daq_table[600];
	for(int nnn=0;nnn<600;nnn++){
		daq_table[nnn]=0.0;
	}
	//for(int nnn=157;nnn<221;nnn++){
	//	daq_table[nnn]=0.2;
	//}
	//for(int nnn=480;nnn<543;nnn++){
	//	daq_table[nnn]=0.2;
	//}
	string buf;

	ifstream ifp(daq_file.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << daq_file.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> runnum >> daq_eff;
		//cout << runnum << ", " << daq_eff <<endl;

		daq_table[runnum] = daq_eff;
	}

//----------------HRS-R Acceptance-----------------//

	string AcceptanceR_table = "../information/RHRS_SIMC.dat";//Acceptance Table (SIMC)
	int RHRS_bin;
	double RHRS_SIMC;
	double RHRS_table[150];//1.5<pk[GeV/c]<2.1, 150 partition --> 1bin=4MeV/c
	double RHRS_total=0.;
	int RHRS_total_bin=0;
	string buf2;

/*----- HRS-R Acceptance Table -----*/
	ifstream ifp2(AcceptanceR_table.c_str(),ios::in);
	if (ifp2.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table.c_str() << endl;
	while(1){
		getline(ifp2,buf2);
		if(buf2[0]=='#'){continue;}
		if(ifp2.eof())break;
		stringstream sbuf2(buf2);
		sbuf2 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (average)="<<RHRS_total/(double)RHRS_total_bin<<endl;


//----------------Mistake-----------------//
//VP Flux should be calculated separately.
//It is no use to make tables
//
//	string vpflux_table = "./vpflux_SIMC.dat";//VP Flux Table (SIMC)
//	int vp_bin;
//	double vpflux_SIMC;
//	double Ng_table[150];//1.8<pe[GeV/c]<2.4, 150 partition --> 1bin=4MeV/c
//	double Ng_total=0.;
//	string buf;
//
///*----- VP Flux Table -----*/
//	ifstream ifp(vpflux_table.c_str(),ios::in);
//	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
//cout << "Param file : " << vpflux_table.c_str() << endl;
//	while(1){
//		getline(ifp,buf);
//		if(buf[0]=='#'){continue;}
//		if(ifp.eof())break;
//		stringstream sbuf(buf);
//		sbuf >> vp_bin >> vpflux_SIMC;
//		cout << vp_bin << ", " << vpflux_SIMC <<endl;
//
//		Ng_table[vp_bin-1] = vpflux_SIMC;
//		Ng_total+=vpflux_SIMC;
//	}
//	cout<<"Ng_total="<<Ng_total<<endl;


    
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
  h_pepk->SetNdivisions(505);
  h_pepk->GetZaxis()->SetLabelOffset(-0.005);
  TH1D* h_pe = new TH1D("h_pe", "p_{e'}" ,200,1980.,2220.);
  TH1D* h_pk = new TH1D("h_pk", "p_{K}" ,200,1720.,1940.);
  TH2D* h_pepk_simc = new TH2D("h_pepk_simc", "p_{K} vs p_{e'}" ,50,1720.,1940.,50,1980.,2220.);
  h_pepk_simc->SetNdivisions(505);
  //h_pepk_simc->GetZaxis()->SetLabelOffset(-0.005);
 
  
//Focal Plane (Feb. 21, 2021)
  TH1D* h_R_y_data       = new TH1D("h_R_y_data"      ,"h_R_y_data [cm]"      ,80,   -6.,  6.);
  TH1D* h_R_y_simc       = new TH1D("h_R_y_simc"      ,"h_R_y_simc [cm]"      ,80,   -6.,  6.);
  tree_simc->Project("h_R_y_simc","h_yfp+0.516","");
  TH1D* h_R_x_data       = new TH1D("h_R_x_data"      ,"h_R_x_data [cm]"      ,80,   -80.,  80.);
  TH1D* h_R_x_simc       = new TH1D("h_R_x_simc"      ,"h_R_x_simc [cm]"      ,80,   -80.,  80.);
  tree_simc->Project("h_R_x_simc","h_xfp","");
  TH1D* h_R_th_data       = new TH1D("h_R_th_data"      ,"h_R_th_data"      ,80,   -0.2,  0.2);
  TH1D* h_R_th_simc       = new TH1D("h_R_th_simc"      ,"h_R_th_simc"      ,80,   -0.2,  0.2);
  tree_simc->Project("h_R_th_simc","h_xpfp","");
  TH1D* h_R_ph_data       = new TH1D("h_R_ph_data"      ,"h_R_ph_data"      ,80,   -0.05,  0.05);
  TH1D* h_R_ph_simc       = new TH1D("h_R_ph_simc"      ,"h_R_ph_simc"      ,80,   -0.05,  0.05);
  tree_simc->Project("h_R_ph_simc","h_ypfp","");
  TH1D* h_R_tg_th_data       = new TH1D("h_R_tg_th_data"      ,"h_R_tg_th_data"      ,80,   -0.1,  0.1);
  TH1D* h_R_tg_th_simc       = new TH1D("h_R_tg_th_simc"      ,"h_R_tg_th_simc"      ,80,   -0.1,  0.1);
  tree_simc->Project("h_R_tg_th_simc","h_xptar","");
  TH1D* h_R_tg_ph_data       = new TH1D("h_R_tg_ph_data"      ,"h_R_tg_ph_data"      ,80,   -0.05,  0.05);
  TH1D* h_R_tg_ph_simc       = new TH1D("h_R_tg_ph_simc"      ,"h_R_tg_ph_simc"      ,80,   -0.05,  0.05);
  tree_simc->Project("h_R_tg_ph_simc","h_yptar","");

  TH1D* h_L_y_data       = new TH1D("h_L_y_data"      ,"h_L_y_data [cm]"      ,80,   -6.,  6.);
  TH1D* h_L_y_simc       = new TH1D("h_L_y_simc"      ,"h_L_y_simc [cm]"      ,80,   -6.,  6.);
  tree_simc->Project("h_L_y_simc","e_yfp+0.807","");
  TH1D* h_L_x_data       = new TH1D("h_L_x_data"      ,"h_L_x_data [cm]"      ,80,   -80.,  80.);
  TH1D* h_L_x_simc       = new TH1D("h_L_x_simc"      ,"h_L_x_simc [cm]"      ,80,   -80.,  80.);
  tree_simc->Project("h_L_x_simc","e_xfp","");
  TH1D* h_L_th_data       = new TH1D("h_L_th_data"      ,"h_L_th_data"      ,80,   -0.2,  0.2);
  TH1D* h_L_th_simc       = new TH1D("h_L_th_simc"      ,"h_L_th_simc"      ,80,   -0.2,  0.2);
  tree_simc->Project("h_L_th_simc","e_xpfp","");
  TH1D* h_L_ph_data       = new TH1D("h_L_ph_data"      ,"h_L_ph_data"      ,80,   -0.05,  0.05);
  TH1D* h_L_ph_simc       = new TH1D("h_L_ph_simc"      ,"h_L_ph_simc"      ,80,   -0.05,  0.05);
  tree_simc->Project("h_L_ph_simc","e_ypfp","");
  TH1D* h_L_tg_th_data       = new TH1D("h_L_tg_th_data"      ,"h_L_tg_th_data"      ,80,   -0.1,  0.1);
  TH1D* h_L_tg_th_simc       = new TH1D("h_L_tg_th_simc"      ,"h_L_tg_th_simc"      ,80,   -0.1,  0.1);
  tree_simc->Project("h_L_tg_th_simc","e_xptar","");
  TH1D* h_L_tg_ph_data       = new TH1D("h_L_tg_ph_data"      ,"h_L_tg_ph_data"      ,80,   -0.05,  0.05);
  TH1D* h_L_tg_ph_simc       = new TH1D("h_L_tg_ph_simc"      ,"h_L_tg_ph_simc"      ,80,   -0.05,  0.05);
  tree_simc->Project("h_L_tg_ph_simc","e_yptar","");

  TH2D* h_L_tg_ph_y_data       = new TH2D("h_L_tg_ph_y_data"      ,"h_L_tg_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  TH2D* h_L_tg_ph_y_simc       = new TH2D("h_L_tg_ph_y_simc"      ,"h_L_tg_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  tree_simc->Project("h_L_tg_ph_y_simc","e_ypfp:-1.*e_yfp-0.807","");
  TH2D* h_R_tg_ph_y_data       = new TH2D("h_R_tg_ph_y_data"      ,"h_R_tg_ph_y_data"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  TH2D* h_R_tg_ph_y_simc       = new TH2D("h_R_tg_ph_y_simc"      ,"h_R_tg_ph_y_simc"   ,80,  -6., 6.   ,80,   -0.05,  0.05);
  tree_simc->Project("h_R_tg_ph_y_simc","h_ypfp:-1.*h_yfp-0.516","");
  TH2D* h_L_tg_th_x_data       = new TH2D("h_L_tg_th_x_data"      ,"h_L_tg_th_x_data"   ,80,  -80., 80.   ,80,   -0.1,  0.1);
  TH2D* h_L_tg_th_x_simc       = new TH2D("h_L_tg_th_x_simc"      ,"h_L_tg_th_x_simc"   ,80,  -80., 80.   ,80,   -0.1,  0.1);
  //tree_simc->Project("h_L_tg_th_x_simc","e_xpfp:e_xfp","");
  tree_simc->Project("h_L_tg_th_x_simc","e_xpfp:e_xfp","e_xpfp<0.17*e_xfp/100.+0.025&&e_xpfp>0.17*e_xfp/100.-0.035&&e_xpfp<0.40*e_xfp/100.+0.130");
  TH2D* h_R_tg_th_x_data       = new TH2D("h_R_tg_th_x_data"      ,"h_R_tg_th_x_data"   ,80,  -80., 80.   ,80,   -0.1,  0.1);
  TH2D* h_R_tg_th_x_simc       = new TH2D("h_R_tg_th_x_simc"      ,"h_R_tg_th_x_simc"   ,80,  -80., 80.   ,80,   -0.1,  0.1);
  //tree_simc->Project("h_R_tg_th_x_simc","h_xpfp:h_xfp","");
  tree_simc->Project("h_R_tg_th_x_simc","h_xpfp:h_xfp","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130");
  TH2D* h_R_y_vz_data       = new TH2D("h_R_y_vz_data"      ,"h_R_y_vz_data"   ,80,   -6.,  6.,80,  -15., 15.   );
  TH2D* h_R_y_vz_simc       = new TH2D("h_R_y_vz_simc"      ,"h_R_y_vz_simc"   ,80,   -6.,  6.,80,  -15., 15.   );
  tree_simc->Project("h_R_y_vz_simc","zposi:h_yfp","");                                                         
  TH2D* h_L_y_vz_data       = new TH2D("h_L_y_vz_data"      ,"h_L_y_vz_data"   ,80,   -6.,  6.,80,  -15., 15.   );
  TH2D* h_L_y_vz_simc       = new TH2D("h_L_y_vz_simc"      ,"h_L_y_vz_simc"   ,80,   -6.,  6.,80,  -15., 15.   );
  tree_simc->Project("h_L_y_vz_simc","zposi:e_yfp","");
  TH2D* h_R_ph_vz_data       = new TH2D("h_R_ph_vz_data"      ,"h_R_ph_vz_data"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  TH2D* h_R_ph_vz_simc       = new TH2D("h_R_ph_vz_simc"      ,"h_R_ph_vz_simc"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  tree_simc->Project("h_R_ph_vz_simc","zposi:h_ypfp","");                                                         
  TH2D* h_L_ph_vz_data       = new TH2D("h_L_ph_vz_data"      ,"h_L_ph_vz_data"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  TH2D* h_L_ph_vz_simc       = new TH2D("h_L_ph_vz_simc"      ,"h_L_ph_vz_simc"   ,80,   -0.05,  0.05,80,  -15., 15.   );
  tree_simc->Project("h_L_ph_vz_simc","zposi:e_ypfp","");
  TH2D* h_R_p_x_data       = new TH2D("h_R_p_x_data"      ,"h_R_p_x_data"   ,80,   1700.,  1950.,80,  -80., 80.  );
  TH2D* h_R_p_x_simc       = new TH2D("h_R_p_x_simc"      ,"h_R_p_x_simc"   ,80,   1700.,  1950.,80,  -80., 80.  );
  tree_simc->Project("h_R_p_x_simc","h_xfp:Rp_rec","");
  TH2D* h_L_p_x_data       = new TH2D("h_L_p_x_data"      ,"h_L_p_x_data"   ,80,   1950.,  2250.,80,  -80., 80.  );
  TH2D* h_L_p_x_simc       = new TH2D("h_L_p_x_simc"      ,"h_L_p_x_simc"   ,80,   1950.,  2250.,80,  -80., 80.  );
  tree_simc->Project("h_L_p_x_simc","e_xfp:Lp_rec","");

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

  tree_simc->Project("h_pepk_simc","Lp_rec:Rp_rec","");
  TH1D* h_pe_simc = new TH1D("h_pe_simc", "p_{e'}" ,200,1980.,2220.);
  tree_simc->Project("h_pe_simc","Lp_rec","");
  TH1D* h_pk_simc = new TH1D("h_pk_simc", "p_{K}" ,200,1720.,1940.);
  tree_simc->Project("h_pk_simc","Rp_rec","");
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",100,-0.15,0.15,100,-0.15,0.15);

  TH1D* h_pe_simc_FPcut = new TH1D("h_pe_simc_FPcut", "p_{e'} (FP cut)" ,200,1980.,2220.);
  tree_simc->Project("h_pe_simc_FPcut","Lp_rec","e_xpfp<0.17*e_xfp/100.+0.025&&e_xpfp>0.17*e_xfp/100.-0.035&&e_xpfp<0.40*e_xfp/100.+0.130");
  TH1D* h_pk_simc_FPcut = new TH1D("h_pk_simc_FPcut", "p_{K} (FP cut)" ,200,1720.,1940.);
  tree_simc->Project("h_pk_simc_FPcut","Rp_rec","h_xpfp<0.17*h_xfp/100.+0.025&&h_xpfp>0.17*h_xfp/100.-0.035&&h_xpfp<0.40*h_xfp/100.+0.130");

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
		if(event_selection&&ct_cut&&abs(L_tr_tg_ph)>0.03)hmm_L_strict->Fill(mm*1000.);
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


		//int ebin = (int)((L_mom-1.8)/0.004);
		//if(ebin>=0 &&ebin<150){
		//Ng = Ng_table[ebin]*Ne;//
		//if(Ng!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
		//else cs=0.;
		//}else{cs=0.;} 
		int kbin = (int)((R_mom-1.5)/0.004);
		if(kbin>=0 &&kbin<150){
		RHRS = RHRS_table[kbin];//
		effDAQ = daq_table[nrun-111000];
		if(effDAQ==0.2)cout<<"Strange!!! DAQ Eff. of run"<<nrun<<" does not exist."<<endl;
		efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		//if(RHRS!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
		if(RHRS!=0.&&effDAQ!=0.)cs = 1./effDAQ/RHRS/100.;//[nb/sr]
		else cs=0.;
		}else{cs=0.;}
		double cs_temp = cs*labtocm;
		//cs_ave += cs;
		//cout<<"Ng="<<Ng<<endl;
		//cout<<"cs="<<cs<<endl;
		if(event_selection&&ct_cut&&theta_gk_cm*180./PI>8.)hcs_L_fom_best->Fill(mm,cs);
//		//if(event_selection_nocut&&ct_cut)hcs_L_fom_nocut->SetBinContent(hmm_L_fom_best->FindBin(mm),cs);

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
			h_L_y_vz_data->Fill(L_tr_y*100.,L_tr_vz*100.);
			h_R_y_vz_data->Fill(R_tr_y*100.,R_tr_vz*100.);
			h_L_ph_vz_data->Fill(L_tr_ph,L_tr_vz*100.);
			h_R_ph_vz_data->Fill(R_tr_ph,R_tr_vz*100.);
			h_L_p_x_data->Fill(L_mom*1000.,L_tr_x*100.);
			h_R_p_x_data->Fill(R_mom*1000.,R_tr_x*100.);
			h_L_th_data->Fill(L_tr_th);
			h_L_ph_data->Fill(L_tr_ph);
			h_L_tg_th_data->Fill(L_tr_tg_th);
			h_L_tg_ph_data->Fill(L_tr_ph);
			h_L_tg_ph_y_data->Fill(-1.*L_tr_y*100.,L_tr_ph);
			h_R_tg_ph_y_data->Fill(-1.*R_tr_y*100.,R_tr_ph);
			h_L_tg_th_x_data->Fill(L_tr_x*100.,L_tr_th);
			h_R_tg_th_x_data->Fill(R_tr_x*100.,R_tr_th);
			//if(R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)h_pe->Fill(L_mom*1000.);
			//if(R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)h_pk->Fill(R_mom*1000.);
		}

		//if(abs(R_tr_vz-L_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		if(abs(R_tr_vz-L_tr_vz)<0.025)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		h_nltrack->Fill(NLtr);
		h_nrtrack->Fill(NRtr);

}//ENum
//	cout<<"nbunch="<<nbunch<<endl;
	//TCanvas* c1 = new TCanvas("c1","c1");
	//hmm_L_fom_best->Draw("");
	//TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
	//hmm_bg_fom_best->Sumw2();
	//hmm_bg_fom_best->Scale(1./nbunch);
	//hmm_bg_fom_best->SetLineColor(kGreen);
	//hmm_bg_fom_best->Draw("same");

	//TCanvas* c2 = new TCanvas("c2","c2");
	//hmm_wobg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.,-1.);
	//hmm_wobg_fom_best->Draw("");

	//TCanvas* c3 = new TCanvas("c3","c3");
	//hmm_L_strict->Draw("");
	//TH1F* hmm_bg_strict=(TH1F*)file_mea_mthesis->Get("hmm_mixed");
	//hmm_bg_strict->Sumw2();
	//hmm_bg_strict->Scale(1./nbunch);
	//hmm_bg_strict->SetLineColor(kGreen);
	//hmm_bg_strict->SetFillColor(kGreen);
	//hmm_bg_strict->SetMarkerColor(kGreen);
	//hmm_bg_strict->Draw("same");

	//TCanvas* c4 = new TCanvas("c4","c4");
	//hmm_wobg_strict->Add(hmm_L_strict,hmm_bg_strict,1.,-1.);
	//hmm_wobg_strict->Draw("");

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
	TCanvas* c10 = new TCanvas("c10","c10",900.,900.);
	c10->Divide(2,2);
	c10->cd(1);
	h_pe_simc->Draw("");
	c10->cd(2);
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
	h_pe_simc->Scale(2081./923474.);//rec
	h_pe_simc->Draw("same");
	c30->cd(2);
	h_pk->SetNdivisions(505);
	h_pk->Draw("");
	h_pk_simc->Scale(2226./1041430.);//rec
	//h_pk_simc->Scale(2286./1040300.);//orig
	//h_pk_simc->Scale(1863./1040300.);
	//h_pk_simc->Scale(2421./1.2e+6);
	h_pk_simc->Draw("same");
	c30->cd(3);
	h_pe->Draw("e");
	h_pe_simc->Draw("same");
	h_pe_simc_FPcut->Scale(2081./923474.);//rec
	h_pe_simc_FPcut->SetLineColor(kGreen);
	h_pe_simc_FPcut->Draw("same");
	c30->cd(4);
	h_pk->Draw("e");
	h_pk_simc->Draw("same");
	h_pk_simc_FPcut->Scale(2226./1041430.);//rec
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
	double dtemp, stemp;
	c50->Divide(3,2);
	c50->cd(1);
	stemp=h_R_y_simc->Integral();
	dtemp=h_R_y_data->Integral();
	h_R_y_simc->Scale(dtemp/stemp);
	h_R_y_data->Draw("e");
	h_R_y_data->Fit("gausn");
	h_R_y_simc->Draw("histsame");
	h_R_y_simc->Fit("gausn");
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
	h_L_y_data->Fit("gausn");
	h_L_y_simc->Draw("histsame");
	h_L_y_simc->Fit("gausn");
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

//cout << "h_R_y_simc = "<<h_R_y_simc->Integral()<<endl;
//cout << "h_R_y_data = "<<h_R_y_data->Integral()<<endl;
//cout << "h_R_x_simc = "<<h_R_x_simc->Integral()<<endl;
//cout << "h_R_x_data = "<<h_R_x_data->Integral()<<endl;
//cout << "h_R_th_simc = "<<h_R_th_simc->Integral()<<endl;
//cout << "h_R_ph_data = "<<h_R_ph_data->Integral()<<endl;
//cout << "h_R_tg_th_simc = "<<h_R_tg_th_simc->Integral()<<endl;
//cout << "h_R_tg_ph_data = "<<h_R_tg_ph_data->Integral()<<endl;
//cout << "h_L_y_simc = "<<h_L_y_simc->Integral()<<endl;
//cout << "h_L_y_data = "<<h_L_y_data->Integral()<<endl;
//cout << "h_L_x_simc = "<<h_L_x_simc->Integral()<<endl;
//cout << "h_L_x_data = "<<h_L_x_data->Integral()<<endl;
//cout << "h_L_th_simc = "<<h_L_th_simc->Integral()<<endl;
//cout << "h_L_ph_data = "<<h_L_ph_data->Integral()<<endl;
//cout << "h_L_tg_th_simc = "<<h_L_tg_th_simc->Integral()<<endl;
//cout << "h_L_tg_ph_data = "<<h_L_tg_ph_data->Integral()<<endl;
//
	TCanvas* c60 = new TCanvas("c60","c60",900.,900.);
	c60->Divide(2,2);
	c60->cd(1);
	h_L_tg_ph_y_simc->Draw("colz");
	c60->cd(2);
	h_L_tg_ph_y_data->Draw("colz");
	c60->cd(3);
	h_L_tg_th_x_simc->Draw("colz");
	c60->cd(4);
	h_L_tg_th_x_data->Draw("colz");
	TCanvas* c61 = new TCanvas("c61","c61",900.,900.);
	c61->Divide(2,2);
	c61->cd(1);
	h_R_tg_ph_y_simc->Draw("colz");
	c61->cd(2);
	h_R_tg_ph_y_data->Draw("colz");
	c61->cd(3);
	h_R_tg_th_x_simc->Draw("colz");
	c61->cd(4);
	h_R_tg_th_x_data->Draw("colz");
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

cout << "Well done!" << endl;
}//fit
