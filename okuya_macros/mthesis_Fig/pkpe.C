//--  Missing Mass   --//
//	vs. Q2
//	vs. theta_gk_cm
//
//K. Okuyama (Nov. 20, 2020)
//K. Okuyama (Dec. 28, 2020)//GeV-->MeV
//K. Okuyama (Jan. 07, 2021)
//
//This is taken over from MM_divide.C
//No array branch mode 
//
double mom_acpt(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double MT = par[0];//Target Mass
	double MM = par[1];//Missing Mass
	double theta_ee = par[2];
	double theta_epk = par[3];
	double MK = 494.;//MeV
	double Ee = 4312.;//MeV
	//double cosine = TMath::Cos(0.05);
	//double cosine = TMath::Cos(13.2*PI/180.);
	//double cosine2 = TMath::Cos(2.*13.2*PI/180.);
	//double cosine = TMath::Cos(0.23);
	//double cosine2 = TMath::Cos(0.46);
	double cosine = TMath::Cos(theta_ee);
	double cosine2 = TMath::Cos(theta_epk);
	double Ex=TMath::Sqrt(x0*x0+MK*MK);
	//double Qsq=2.*4318.*2200.*(1-cosine);
	double Qsq=0.476*1000000.;
	double Me=0.511;
	double pe=TMath::Sqrt(Ee*Ee-Me*Me);
	//double Eep=TMath::Sqrt(x0*x0+Me*Me);//Ee'
//	double Qsq=-(Ee*Ee-2*Ee*Eep+Eep*Eep)+pe*pe+x0*x0-2*pe*x0*TMath::Cos(13.2*PI/180.);
	
	//return 4300-(2*(MT+TMath::Sqrt(x0*x0+MK*MK))-TMath::Sqrt((MT+TMath::Sqrt(x0*x0+MK*MK))*((MT+TMath::Sqrt(x0*x0+MK*MK)))+2*x0*TMath::Cos(0.02)*(MM*MM-MT*MT-MK*MK+2*MT*TMath::Sqrt(x0*x0+MK*MK))))/2/x0/TMath::Cos(0.02);
//	return 4300.-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine);
	//return Ee-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine-Qsq);
	//cout<<"MM="<<MM<<endl;
	//cout<<"MT="<<MT<<endl;
	//cout<<"Ex="<<Ex<<endl;
	//cout<<"x0="<<x0<<endl;
	//cout<<"Ee="<<Ee<<endl;
	//cout<<"pe="<<pe<<endl;
	//cout<<"Qsq="<<Qsq<<endl;
	//return (-1.*MM*MM+MT*MT+Ex*Ex+2.*(MT-Ex)*Ee-2.*MT*Ex-x0*x0+x0*(pe-2.1)*cos(13.2*PI/180.))/(2.*(MT-Ex));
	//return Ee+(-1.*MM*MM+MT*MT+MK*MK-2*MT*Ex-Qsq)/(2*MT-2*Ex+2*x0*cos(3.*PI/180.));
	//double temp = (-1.*MM*MM+2.*Me*Me+MT*MT+Ex*Ex+2.*(MT-Ex)*Ee-2.*MT*Ex+2.*x0*pe*cos(13.2*PI/180.))/(2.*Ee+2.*(MT-Ex)-2.*pe*cos(13.2*PI/180.)+2.*x0*cos(26.4*PI/180.));
	//cout<<"pe'="<<temp<<endl;
	//return temp;
	 double a = 2.*(Ee+MT-Ex)-2.*pe*cosine+2.*x0*cosine2;
     double b = (Ee+MT-Ex)*(Ee+MT-Ex)-pe*pe-x0*x0+2.*x0*pe*cosine-MM*MM;
    // cout<<"a="<<a<<endl;
    // cout<<"b="<<b<<endl;
     return b/a;

	}

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

void pkpe(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
	//ACCBGの引き算はmea_hist.ccから
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_2020Nov.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  //TFile *file_mea_mthesis = new TFile("../MixedEventAnalysis/bgmea_mthesis.root","read");//MeV
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
  TH2F* h_mm_Qsq = new TH2F("h_mm_Qsq","",100,-100.,200.,50,0.2,0.8);
  h_mm_Qsq->GetYaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");
  h_mm_Qsq->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  TH2F* h_mm_W = new TH2F("h_mm_W","",100,-100.,200.,50,2.05,2.25);
  h_mm_W->GetYaxis()->SetTitle("W [GeV]");
  h_mm_W->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  TH2F* h_mm_theta_gk_cm = new TH2F("h_mm_theta_gk_cm","",100,-100.,200.,50,-5.,25.);
  h_mm_theta_gk_cm->GetYaxis()->SetTitle("#theta_{#gamma K}^{CM} [deg]");
  h_mm_theta_gk_cm->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  TH2F* h_mm_eps = new TH2F("h_mm_eps","",100,-100.,200.,50,0.72,0.8);
  h_mm_eps->GetYaxis()->SetTitle("#epsilon");
  h_mm_eps->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
  TH2F* h_pk_W = new TH2F("h_pk_W","",50,1750.,1950.,50,2.05,2.25);
  TH2F* h_pk_Qsq = new TH2F("h_pk_Qsq","",50,1750.,1950.,50,0.2,0.8);
  TH2F* h_pk_theta_gk_cm = new TH2F("h_pk_theta_gk_cm","",50,1750.,1950.,50,-5.,25.);
  TH2F* h_pk_eps = new TH2F("h_pk_eps","",50,1750.,1950.,50,0.72,0.8);
  TH1F* h_Qsq = new TH1F("h_Qsq","",500,0.2,0.8);
  TH1F* h_Qsq2 = new TH1F("h_Qsq2","",500,0.2,0.8);
  TH1F* h_theta_gk_cm = new TH1F("h_theta_gk_cm","",500,-5.,25.);
  TH1F* h_theta_gk_cm2 = new TH1F("h_theta_gk_cm2","",500,-5.,25.);
  TH1F* h_eps = new TH1F("h_eps","",500,0.72,0.8);
  TH1F* h_eps2 = new TH1F("h_eps2","",500,0.72,0.8);
  TH1F* h_W = new TH1F("h_W","",500,2.05,2.25);
  TH1F* h_W2 = new TH1F("h_W2","",500,2.05,2.25);
  TH2F* h_pkpe = new TH2F("h_pkpe","",60,1730.,1940.,60,1980.,2190.);


  TH2F* h_Qsq_eps = new TH2F("h_Qsq_eps","",50,0.2,0.8,50,0.72,0.8);
  h_Qsq_eps->GetYaxis()->SetTitle("#epsilon");
  h_Qsq_eps->GetXaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");
  TH2F* h_W_eps = new TH2F("h_W_eps","",50,2.05,2.25,50,0.72,0.8);
  h_W_eps->GetYaxis()->SetTitle("#epsilon");
  h_W_eps->GetXaxis()->SetTitle("W [GeV]");
  TH2F* h_theta_gk_cm_eps = new TH2F("h_theta_gk_cm_eps","",50,-5.,25.,50,0.72,0.8);
  h_theta_gk_cm_eps->GetYaxis()->SetTitle("#epsilon");
  h_theta_gk_cm_eps->GetXaxis()->SetTitle("#theta_{#gamma K}^{CM} [deg]");





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
 // TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "theta_gk_cm",1000,0.,0.3);
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
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",100,-0.15,0.15,100,-0.15,0.15);

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
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm = (gamma*pR_cm*pR_cm*(pR_cm*cos(theta_gk_cm)+beta*ER_cm))/(pow(sqrt(pR_cm*pR_cm*sin(theta_gk_cm)*sin(theta_gk_cm)+gamma*gamma*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)),3.));
//cout<<"labtocm="<<labtocm<<endl;
		double tan_lab1 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+beta*sqrt(MK*MK+pR_cm*pR_cm)/pR_cm));
		double tan_lab2 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+(omega*Mp-Qsq*Qsq)/(omega*Mp+Mp*Mp)));
		//if(tan_lab1!=tan_lab2)cout<<"tan1="<<atan(tan_lab1)<<", tan2="<<atan(tan_lab2)<<"theta_gk_lab="<<theta_gk_lab<<endl;
		
		double q2=Qsq+omega*omega;
		double eps=1/(1+2*(q2/Qsq)*tan(theta_ee/2)*tan(theta_ee/2));


// Qsq = -qsq = mom_g*mom_g - omega*omega
		//cout<<"Qsq="<<Qsq<<endl;
		//cout<<"q2="<<q2<<endl;
		//cout<<"mog_g^2="<<mom_g*mom_g<<endl;
		//cout<<"omega^2="<<omega*omega<<endl;

		if(event_selection&&ct_cut){
			hmm_L_fom_best->Fill(mm);
			hmm_L_strict->Fill(mm*1000.);
			h_mm_Qsq->Fill(mm*1000.,Qsq);
			h_mm_W->Fill(mm*1000.,W);
			h_mm_theta_gk_cm->Fill(mm*1000.,theta_gk_cm*180./PI);
			h_mm_eps->Fill(mm*1000.,eps);
			h_pk_Qsq->Fill(R_mom*1000.,Qsq);
			h_pk_W->Fill(R_mom*1000.,W);
			h_pk_theta_gk_cm->Fill(R_mom*1000.,theta_gk_cm*180./PI);
			h_pk_eps->Fill(R_mom*1000.,eps);
			h_Qsq_eps->Fill(Qsq,eps);
			h_W_eps->Fill(W,eps);
			h_theta_gk_cm_eps->Fill(theta_gk_cm*180./PI,eps);
			h_Qsq->Fill(Qsq);
			h_eps->Fill(eps);
			h_W->Fill(W);
			h_theta_gk_cm->Fill(theta_gk_cm*180./PI);
			if(L_mom>2.010&&L_mom<2.160&&R_mom>1.760&&R_mom<1.900){
				h_Qsq2->Fill(Qsq);
				h_eps2->Fill(eps);
				h_W2->Fill(W);
				h_theta_gk_cm2->Fill(theta_gk_cm*180./PI);
				h_pkpe->Fill(R_mom*1000.,L_mom*1000.);
			}
		}

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
		if(effDAQ==0.2)cout<<"Starange!!! DAQ Eff. of run"<<nrun<<" does not exist."<<endl;
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
			if(theta_gk_cm*180./PI>=8.)h_pepk->Fill(R_mom,L_mom);
		}

		//if(abs(R_tr_vz-L_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		if(abs(R_tr_vz-L_tr_vz)<0.025)h_zave->Fill((R_tr_vz+L_tr_vz)/2.);
		h_nltrack->Fill(NLtr);
		h_nrtrack->Fill(NRtr);

}//ENum

//	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
	hmm_bg_fom_best->Sumw2();
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_best->SetLineColor(kGreen);
	hmm_bg_fom_best->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2");
	hmm_wobg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.,-1.);
	hmm_wobg_fom_best->Draw("");

	TCanvas* c3 = new TCanvas("c3","c3");
h_pkpe->Draw("colz");
TF1* func_acpt = new TF1("func_acpt",mom_acpt,1700,2000,4);
func_acpt->SetNpx(600);
func_acpt->SetParameter(0,Mp*1000.);
func_acpt->SetParameter(1,ML*1000.);
func_acpt->SetParameter(2,0.23);
func_acpt->SetParameter(3,0.46);
func_acpt->SetLineColor(kAzure);
func_acpt->SetLineWidth(4);
func_acpt->Draw("same");
cout<<"Lambda(pK=1920)="<<func_acpt->Eval(1920.)<<endl;
cout<<"Lambda(pK=1910)="<<func_acpt->Eval(1910.)<<endl;
cout<<"Lambda(pK=1900)="<<func_acpt->Eval(1900.)<<endl;
cout<<"Lambda(pK=1890)="<<func_acpt->Eval(1890.)<<endl;
cout<<"Lambda(pK=1880)="<<func_acpt->Eval(1880.)<<endl;
cout<<"Lambda(pK=1870)="<<func_acpt->Eval(1870.)<<endl;
cout<<"Lambda(pK=1860)="<<func_acpt->Eval(1860.)<<endl;
cout<<"Lambda(pK=1850)="<<func_acpt->Eval(1850.)<<endl;
cout<<"Lambda(pK=1840)="<<func_acpt->Eval(1840.)<<endl;
TF1* func_acpt_a = new TF1("func_acpt_a",mom_acpt,1700,2000,4);
func_acpt_a->SetNpx(600);
func_acpt_a->SetParameter(0,Mp*1000.);
func_acpt_a->SetParameter(1,ML*1000.);
func_acpt_a->SetParameter(2,0.20);
func_acpt_a->SetParameter(3,0.40);
func_acpt_a->SetLineColor(kRed);
func_acpt_a->SetLineStyle(10);
func_acpt_a->SetLineWidth(4);
func_acpt_a->Draw("same");
cout<<"Lambda_a(pK=1900)="<<func_acpt_a->Eval(1900.)<<endl;
cout<<"Lambda_a(pK=1914.1)="<<func_acpt_a->Eval(1914.1)<<endl;
cout<<"Lambda_a(pK=1914.2)="<<func_acpt_a->Eval(1914.2)<<endl;
cout<<"Lambda_a(pK=1914.3)="<<func_acpt_a->Eval(1914.3)<<endl;
cout<<"Lambda_a(pK=1914.4)="<<func_acpt_a->Eval(1915.4)<<endl;
cout<<"Lambda_a(pK=1914.5)="<<func_acpt_a->Eval(1915.5)<<endl;
TF1* func_acpt_b = new TF1("func_acpt_b",mom_acpt,1700,2000,4);
func_acpt_b->SetNpx(600);
func_acpt_b->SetParameter(0,Mp*1000.);
func_acpt_b->SetParameter(1,ML*1000.);
func_acpt_b->SetParameter(2,0.26);
func_acpt_b->SetParameter(3,0.52);
func_acpt_b->SetLineColor(kRed);
func_acpt_b->SetLineStyle(10);
func_acpt_b->SetLineWidth(4);
func_acpt_b->Draw("same");
cout<<"Lambda_b(pK=1900)="<<func_acpt_b->Eval(1900.)<<endl;

TF1* func_acpt2 = new TF1("func_acpt2",mom_acpt,1700,2000,4);
func_acpt2->SetNpx(600);
func_acpt2->SetParameter(0,Mp*1000.);
func_acpt2->SetParameter(1,MS0*1000.);
func_acpt2->SetParameter(2,0.23);
func_acpt2->SetParameter(3,0.46);
func_acpt2->SetLineColor(kCyan);
func_acpt2->SetLineWidth(4);
func_acpt2->Draw("same");
TF1* func_acpt2_a = new TF1("func_acpt2_a",mom_acpt,1700,2000,4);
func_acpt2_a->SetNpx(600);
func_acpt2_a->SetParameter(0,Mp*1000.);
func_acpt2_a->SetParameter(1,MS0*1000.);
func_acpt2_a->SetParameter(2,0.20);
func_acpt2_a->SetParameter(3,0.40);
func_acpt2_a->SetLineColor(kRed);
func_acpt2_a->SetLineStyle(10);
func_acpt2_a->SetLineWidth(4);
func_acpt2_a->Draw("same");
TF1* func_acpt2_b = new TF1("func_acpt2_b",mom_acpt,1700,2000,4);
func_acpt2_b->SetNpx(600);
func_acpt2_b->SetParameter(0,Mp*1000.);
func_acpt2_b->SetParameter(1,MS0*1000.);
func_acpt2_b->SetParameter(2,0.26);
func_acpt2_b->SetParameter(3,0.52);
func_acpt2_b->SetLineColor(kRed);
func_acpt2_b->SetLineStyle(10);
func_acpt2_b->SetLineWidth(4);
func_acpt2_b->Draw("same");

TLine *tl1, *tl2, *tl3, *tl4;
tl1 = new TLine(1760.,2010.,1760.,2160.);
tl2 = new TLine(1900.,2010.,1900.,2160.);
tl3 = new TLine(1760.,2010.,1900.,2010.);
tl4 = new TLine(1760.,2160.,1900.,2160.);
	tl1->SetLineWidth(1);
	tl1->SetLineColor(kBlack);
	tl1->Draw("same");
	tl2->SetLineWidth(1);
	tl2->SetLineColor(kBlack);
	tl2->Draw("same");
	tl3->SetLineWidth(1);
	tl3->SetLineColor(kBlack);
	tl3->Draw("same");
	tl4->SetLineWidth(1);
	tl4->SetLineColor(kBlack);
	tl4->Draw("same");

//	c3->Print("./pdf/mm_tight.pdf");
//	c4->Print("./pdf/mm_wobg_tight.pdf");
//	c5->Print("./pdf/mm_Qsq.pdf");
//	c6->Print("./pdf/mm_W.pdf");
//	c7->Print("./pdf/mm_theta.pdf");
//	c8->Print("./pdf/mm_eps.pdf");

cout << "Well done!" << endl;
}//fit
