//--  Al Contamination   --//
//
//K. Okuyama (Oct. 25, 2020)
//K. Okuyama (Dec. 12, 2020)
//
//This is taken over from result.C
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
double F_gas( double *x, double *par)
{
//par[0]: Al Gauss sigma
//par[1]: Al second gauss strength
//par[2]: Al(front) second gauss pos. 
//par[3]: Al second gauss sigma 
//par[4]: pol2 coeff.1
//par[5]: pol2 coeff.2
//par[6]: total scale 

double gas = 0.;
int np = 2000;
for(int i=0;i<np;i++){
 double pol2 = par[4]*pow((x[0]-par[5]),2.)+1.;
 double step = -100.+(double)i/10.;
 if(step<-0.125*100.) pol2 = 0.;
 if(step> 0.125*100.) pol2 = 0.;

 double fland;
 fland = exp(-0.5*pow((x[0]-step)/par[0],2.))+par[1]*exp(-0.5*pow((x[0]+par[2]-step)/par[3],2.));
 gas += pol2 * fland;
}
gas = par[6]*gas;
return gas;
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

double g_f1, g_f2;//front Gauss
double g_r1, g_r2;//rear Gauss
g_f1=par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
g_f2=par[0]*par[3]*exp(-0.5*pow((x[0]-par[1]+par[4])/par[5],2.));
g_r1=par[6]*exp(-0.5*pow((x[0]-par[7])/par[2],2.));
g_r2=par[6]*par[3]*exp(-0.5*pow((x[0]-par[7]+par[4])/par[5],2.));
double g_front = g_f1 + g_f2;
double g_rear  = g_r1 + g_r2;

double gas = 0.;
int np = 2000;
for(int i=0;i<np;i++){
 double pol2 = par[8]*pow((x[0]-par[9]),2.)+1.;
 double step = -100.+(double)i/10.;
 if(step<-0.125*100.) pol2 = 0.;
 if(step> 0.125*100.) pol2 = 0.;

 double fland;
 fland = exp(-0.5*pow((x[0]-step)/par[2],2.))+par[3]*exp(-0.5*pow((x[0]+par[4]-step)/par[5],2.));
 gas += pol2 * fland;
}
gas = par[10]*gas;
double val =  g_front + g_rear + gas;
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

void vertex_fit3(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("../h2all_Lsingle.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file_dummy = new TFile("../dummy_Tkin/tritium_111325.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file_true = new TFile("../dummy_Tkin/tritium_111157.root","read");//input file of all H2 run(default: h2all4.root)
	//ACCBGの引き算はmea_hist.ccから
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  //TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_llccrr_Lsingle.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_csbase.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  double nbunch = 6000.;//effetive bunches (6 bunches x 1000 mixtures)
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
  TChain *chain_dummy = new TChain("T");
	//chain_dummy = (TChain*)file_dummy->Get("T");
  TChain *chain_true = new TChain("T");// = (TChain*)file_true->Get("T");
  //TChain *chain_dummy = new TChain("tree_dummy","");
  //TChain *chain_true = new TChain("tree_true","");
	chain_dummy->Add("../dummy_Tkin/tritium_111325.root");
	chain_dummy->Add("../dummy_Tkin/tritium_111325_1.root");
	chain_dummy->Add("../dummy_Tkin/tritium_111323.root");
	chain_dummy->Add("../dummy_Tkin/tritium_111324.root");
	chain_dummy->Add("../dummy_Tkin/tritium_111324_1.root");
	chain_true->Add("~/data_okuyama/rootfiles/nnL_temp/tritium_111157_1.root");
	chain_true->Add("~/data_okuyama/rootfiles/nnL_temp/tritium_111157_2.root");
	chain_true->Add("~/data_okuyama/rootfiles/nnL_temp/tritium_111157_3.root");
	chain_true->Add("~/data_okuyama/rootfiles/nnL_temp/tritium_111157_4.root");
	chain_true->Add("~/data_okuyama/rootfiles/nnL_temp/tritium_111157_5.root");
//	tree->Write();


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
  TH2F* h_zz_dummy  = new TH2F("h_zz_dummy","h_zz_dummy",100,-0.15,0.15,100,-0.15,0.15);

  TH1F* h_zave  = new TH1F("h_zave","Z-vertex (Ave.)",1000,-25.,25.);
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



  int ENum_dummy = chain_dummy->GetEntries();
cout<<"Entries(dummy): "<<ENum_dummy<<endl;
 // for(int i=0;i<ENum_dummy;i++){
//	chain_dummy->GetEntry(i);
	//if(abs(R_tr_vz2-L_tr_vz2)<0.025)h_zave_dummy->Fill((R_tr_vz2+L_tr_vz2)/2.);
	chain_dummy->Project("h_zave_dummy","(R.tr.vz+L.tr.vz)/2.","abs(R.tr.vz-L.tr.vz)<0.025","");
	chain_dummy->Project("h_zz_dummy","R.tr.vz:L.tr.vz","","");
//	}
  int ENum_true = chain_true->GetEntries();
cout<<"Entries(true): "<<ENum_true<<endl;
 // for(int i=0;i<ENum_true;i++){
//	chain_true->GetEntry(i);
	//if(abs(R_tr_vz3-L_tr_vz3)<0.025)h_zave_true->Fill((R_tr_vz3+L_tr_vz3)/2.);
	chain_true->Project("h_zave_true","(R.tr.vz+L.tr.vz)/2.","abs(R.tr.vz-L.tr.vz)<0.025","");
//	}

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
		if(abs(R_tr_vz-L_tr_vz)<0.025)h_zave->Fill((R_tr_vz+L_tr_vz)*100./2.);
		h_nltrack->Fill(NLtr);
		h_nrtrack->Fill(NRtr);

}//ENum
//	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
//	hcs_L_fom_best->Draw("");
//	TCanvas* c2 = new TCanvas("c2","c2");
//	TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
//	TH1F* hmm_Al_fom_nocut=(TH1F*)file->Get("hmm_Al_fom_best");
//	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
//	//TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
//	//TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
//	hcs_L_fom_best->Rebin(2);
//	TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut");
//	TH1F* hmm_Albg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_forAl");
//	//TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_nocut_result");
//    int fitmin = hmm_L_fom_best->FindBin(0.10);
//    int fitmax = hmm_L_fom_best->FindBin(0.15);
//    double num1 = hmm_L_fom_best->Integral(fitmin,fitmax);
//    double num2 = hmm_bg_fom_best->Integral(fitmin,fitmax);
//    double mixscale = num1/num2;
//	cout<<"hmm_L integral ="<<num1<<endl;
//	cout<<"hmm_bg integral ="<<num2<<endl;
//	cout<<"mixscale(mixed/original)="<<1/mixscale<<endl;
	hmm_bg_fom_best->Sumw2();
//	hmm_bg_fom_nocut->Sumw2();
//	hmm_Albg_fom_nocut->Sumw2();
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_best->SetLineColor(kGreen);
	hmm_bg_fom_best->Draw("same");
//	//cs	 = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
//	cs	 = 1./(effDAQ*RHRS*100.);//[nb/sr]
//	hmm_bg_fom_best->Scale(cs);
//	hmm_bg_fom_best->Rebin(2);
//	hmm_bg_fom_nocut->Scale(1./nbunch);
//	hmm_Albg_fom_nocut->Scale(1./nbunch);
//	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
//	hmm_wo_bg_fom_best->Add(hcs_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
//	hmm_wo_bg_fom_nocut->Add(hmm_L_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
//	//hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
//	//hmm_pi_wobg_fom_nocut->Add(hmm_pi_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
//	hmm_pi_wobg_fom_nocut->Add(hmm_Al_fom_nocut,hmm_Albg_fom_nocut,1.0,-1.0);
//
//
//	 double constL=0.;
//	 double meanL=0.;
//	 double sigL=0.;
//	 double constS=0.;
//	 double meanS=0.;
//	 double sigS=0.;
//	 double par1 = 0.;
//	 double par2 = 0.;
//	 double par3 = 0.;
//	 double par4 = 0.;
//	 double par5 = 0.;
//	 double par6 = 0.;
//	 double integralL_nocut = 0.;
//	 double integralS_nocut = 0.;
//	 n_L_nocut=0.;
//	 n_S_nocut=0.;
//	 double fmin = -0.05;
//	 double fmax =  0.15;
///****************************************/
///************BEST CUT********************/
///****************************************/
//cout<<"BEST CUT START"<<endl;
//
//	fmmbg_best=new TF1("fmmbg_best","pol4",min_mm,max_mm);
//	fmmbg_best->SetNpx(2000);
//	 fL_best=new TF1("fL_best","gaus(0)",min_mm,max_mm);
//	 fL_best->SetNpx(2000);
//	 fL_best->SetParameters(def_n_L,def_mean_L,def_sig_L);
//	 fL_best->SetParLimits(0,0.,100000);
//	 fL_best->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 fL_best->SetParLimits(2,0.,0.01);
//	 fS_best=new TF1("fS_best","gaus(0)",min_mm,max_mm);
//	 fS_best->SetNpx(2000);
//	 fS_best->SetParameters(def_n_S,def_mean_S,def_sig_S);
//	 fS_best->SetParLimits(0,0.,100000);
//	 fS_best->SetParLimits(1,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fS_best->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
//	 
//	 fmm_best=new TF1("fmm_best","gaus(0)+gaus(3)+pol4(6)",min_mm,max_mm);
//	 fmm_best->SetNpx(2000);
//	 fmm_best->SetTitle("Missing Mass (best cut)");
//	 fmm_best->SetParLimits(0,0.,1000000.);//positive
//	 fmm_best->SetParLimits(3,0.,1000000.);//positive
//		
//	 hmm_wo_bg_fom_best->Fit("fL_best","N","",def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 const_L_best=fL_best->GetParameter(0);
//	 mean_L_best=fL_best->GetParameter(1);
//	 sig_L_best=fL_best->GetParameter(2);
//	
//	 hmm_wo_bg_fom_best->Fit("fS_best","N","",def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 const_S_best=fS_best->GetParameter(0);
//	 mean_S_best=fS_best->GetParameter(1);
//	 sig_S_best=fS_best->GetParameter(2);
//
//
//	 constL=0.;
//	 meanL=0.;
//	 sigL=0.;
//	 constS=0.;
//	 meanS=0.;
//	 sigS=0.;
//	 par1 = 0.;
//	 par2 = 0.;
//	 par3 = 0.;
//	 par4 = 0.;
//	 par5 = 0.;
//	 par6 = 0.;
//	 double integralL_best = 0.;
//	 double integralS_best = 0.;
//	 n_L_best=0.;
//	 n_S_best=0.;
//
///*%%%%%%%%%%%%%%%%%%%%%%%%*/
///*%%    4th Polynomial	%%*/
///*%%%%%%%%%%%%%%%%%%%%%%%%*/
//	//--- w/ 4th Polynomial func.
//	 cout<<"4Poly MODE START"<<endl;
//	 fmm_best_4Poly=new TF1("fmm_best_4Poly",FMM_Res,fit_min_mm,fit_max_mm,14);
//   //par[0]=Width (scale) parameter of Landau density
//   //par[1]=Most Probable (MP, location) parameter of Landau density
//   //par[2]=Total area (integral -inf to inf, normalization constant)
//   //par[3]=Width (sigma) of convoluted Gaussian function
//   //par[4]=tau of exp function
//   //par[5]=Shift of Function Peak
//   //par[6]=Relative Strength
//   //
//   //
//   //
//   //par[0]=Width (scale) parameter of Landau density
//   //par[1]=Most Probable (MP, location) parameter of Landau density
//   //par[2]=Total area (integral -inf to inf, normalization constant)
//   //par[3]=Width (sigma) of convoluted Gaussian function
//	 fmmbg_best_4Poly=new TF1("fmmbg_best_4Poly","pol4",fit_min_mm,fit_max_mm);
//	 fmm_best_4Poly->SetNpx(20000);
//	 fmm_best_4Poly->SetTitle("Missing Mass (best)");
//	 fmm_best_4Poly->SetParLimits(2,0.,1000.);//positive
//	 fmm_best_4Poly->SetParLimits(9,0.,300.);//positive
//	 fmm_best_4Poly->SetParameter(0,0.0007);//Landau width
//	 fmm_best_4Poly->SetParameter(1,mean_L_best);
//	 fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 fmm_best_4Poly->SetParameter(2,1.5);//total scale
//	 fmm_best_4Poly->SetParameter(3,0.001);//sigma
//	 fmm_best_4Poly->SetParLimits(3,0.,0.01);
//	 fmm_best_4Poly->SetParameter(4,0.05);//att.
//	 fmm_best_4Poly->SetParLimits(4,0.005,0.08);
//	 fmm_best_4Poly->SetParameter(5,-0.004);//peak pos.
//	 fmm_best_4Poly->SetParLimits(5,-0.05,0.05);
//	 fmm_best_4Poly->SetParameter(6,0.6);//relative strength
//	 fmm_best_4Poly->SetParLimits(6,0.,1.5);//relative strength
//
//	 fmm_best_4Poly->SetParameter(7,0.0003);//Landau width
//	 fmm_best_4Poly->SetParameter(8,mean_S_best);//MPV
//	 fmm_best_4Poly->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fmm_best_4Poly->SetParameter(9,0.4);//total scale
//	 fmm_best_4Poly->SetParameter(10,0.0015);//sigma
//	 fmm_best_4Poly->SetParLimits(10,0.,0.01);
//	 fmm_best_4Poly->SetParameter(11,0.05);//att
//	 fmm_best_4Poly->SetParLimits(11,0.03,0.12);
//	 fmm_best_4Poly->SetParameter(12,0.080);//peak pos.
//	 //fmm_best_4Poly->SetParLimits(15,-0.085,-0.055);
//	 fmm_best_4Poly->SetParameter(13,0.6);
//	 fmm_best_4Poly->SetParLimits(13,0.,1.5);//relative strength
//
//	 hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
//	 double chisq = fmm_best_4Poly->GetChisquare();
//	 double dof  = fmm_best_4Poly->GetNDF();
//	 cout<<"chisq="<<chisq<<endl;
//	 cout<<"dof="<<dof<<endl;
//	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;
//
//
//	 TF1* fmm_Lambda_only = new TF1("fmm_Lambda_only",FMM_Response,fit_min_mm,fit_max_mm,7);
//	 TF1* fmm_Sigma_only  = new TF1("fmm_Sigma_only" ,FMM_Response, fit_min_mm,fit_max_mm,7);
////Lambda_only
//	 fmm_Lambda_only->SetNpx(20000);
//	 fmm_Lambda_only->SetParameter(0,fmm_best_4Poly->GetParameter(0));
//	 fmm_Lambda_only->SetParameter(1,fmm_best_4Poly->GetParameter(1));
//	 fmm_Lambda_only->SetParameter(2,fmm_best_4Poly->GetParameter(2));
//	 fmm_Lambda_only->SetParameter(3,fmm_best_4Poly->GetParameter(3));
//	 fmm_Lambda_only->SetParameter(4,fmm_best_4Poly->GetParameter(4));
//	 fmm_Lambda_only->SetParameter(5,fmm_best_4Poly->GetParameter(5));
//	 fmm_Lambda_only->SetParameter(6,fmm_best_4Poly->GetParameter(6));
////Sigma_only
//	 fmm_Sigma_only->SetNpx(20000);
//	 fmm_Sigma_only->SetParameter(0,fmm_best_4Poly->GetParameter(7));
//	 fmm_Sigma_only->SetParameter(1,fmm_best_4Poly->GetParameter(8));
//	 fmm_Sigma_only->SetParameter(2,fmm_best_4Poly->GetParameter(9));
//	 fmm_Sigma_only->SetParameter(3,fmm_best_4Poly->GetParameter(10));
//	 fmm_Sigma_only->SetParameter(4,fmm_best_4Poly->GetParameter(11));
//	 fmm_Sigma_only->SetParameter(5,fmm_best_4Poly->GetParameter(12));
//	 fmm_Sigma_only->SetParameter(6,fmm_best_4Poly->GetParameter(13));
//
//	double nofL = fmm_Lambda_only->Integral(fit_min_mm,fit_max_mm);
//	double nofL_old = fmm_Lambda_only->Integral(-0.006,0.006);
//	nofL = nofL/fit_bin_width;
//	nofL_old = nofL_old/fit_bin_width;
//	cout<<"Number of Lambda (TF1 Integral) = "<<nofL<<endl;
//	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old<<endl;
//	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(-0.006),hmm_wo_bg_fom_best->FindBin(0.006))<<endl;
//
//	double nofS = fmm_Sigma_only->Integral(fit_min_mm,fit_max_mm);
//	nofS = nofS/fit_bin_width;
//	cout<<"Number of Sigma (TF1 Integral) = "<<nofS<<endl;
//
//	 hmm_wo_bg_fom_best->Draw();
//	 fmm_best_4Poly->SetLineColor(kRed);
//	 fmm_best_4Poly->Draw("same");
//	 fmm_Lambda_only->SetLineColor(kAzure);
//	 fmm_Sigma_only->SetLineColor(kCyan);
//	 fmm_Lambda_only->Draw("same");
//	 fmm_Sigma_only->Draw("same");
//	 //fL_best->SetLineColor(kGreen);
//	 //fS_best->SetLineColor(kGreen);
//	 //fL_best->Draw("same");
//	 //fS_best->Draw("same");

	TCanvas* c3 = new TCanvas("c3","c3");
	c3->Divide(2,2);
	c3->cd(1);
	gklab_gkcm->Draw("colz");
	c3->cd(2);
	gklab_eklab->Draw("colz");
	c3->cd(3);
	eelab_eklab->Draw("colz");
	c3->cd(4);
	//eklab_gkcm->Draw("colz");
	ekcm_gkcm->Draw("colz");
	TCanvas* c4 = new TCanvas("c4","c4");
	c4->Divide(2,2);
	c4->cd(1);
	cos_gklab_gkcm->Draw("colz");
	c4->cd(2);
	cos_gklab_eklab->Draw("colz");
	c4->cd(3);
	cos_eelab_eklab->Draw("colz");
	c4->cd(4);
	//cos_eklab_gkcm->Draw("colz");
	cos_ekcm_gkcm->Draw("colz");
	//TCanvas* c4 = new TCanvas("c4","c4");
	//h_nltrack->Draw("");
	//TCanvas* c5 = new TCanvas("c5","c5");
	//h_nrtrack->Draw("");
	TCanvas* c6 = new TCanvas("c6","c6");
	h_pepk->Draw("colz");
	TCanvas* c7 = new TCanvas("c7","c7");
	c7->SetLogy(1);
	 TF1* fz_true_front = new TF1("fz_true_front","F_VZ",-0.15,0.15,13);
	 TF1* fz_true_back  = new TF1("fz_true_back ","F_VZ",-0.15,0.15,13);
	 TF1* fz_dummy_front = new TF1("fz_dummy_front","F_VZ",-0.15,0.15,13);
	 TF1* fz_dummy_back  = new TF1("fz_dummy_back ","F_VZ",-0.15,0.15,13);
	 TF1* fz_corr_front = new TF1("fz_corr_front","F_VZ",-0.15,0.15,13);
	 TF1* fz_corr_back  = new TF1("fz_corr_back ","F_VZ",-0.15,0.15,13);
	 //TF1* fz_corr  = new TF1("fz_corr","F_VZ",-0.15,0.15,17);
	 TF1* fz_corr  = new TF1("fz_corr","F_VZ2",-16.,16.,11);//total
	 TF1* fz_cell  = new TF1("fz_cell","F_VZ_cell",-16.,16.,8);//cell
	 TF1* fz_gas  = new TF1("fz_gas","F_gas",-25.,25.,7);//gas
//-------------------------------//
//    NEW Z FUNCTION             //
//-------------------------------//
//par[0]: Al front scale
//par[1]: Al front pos.
//par[2]: Al Gauss sigma
//par[3]: Al second gauss strength
//par[4]: Al(front) second gauss pos. 
//par[5]: Al second gauss sigma 
//par[6]: Al rear scale
//par[7]: Al rear pos. 
//par[8]: pol2 coeff.1
//par[9]: pol2 coeff.2
//par[10]: total scale 
    fz_corr->SetParameter(0,380.);
    fz_corr->SetParameter(1,-0.125*100.);
    fz_corr->SetParameter(2,0.005*100.);
    fz_corr->SetParameter(3,0.4);
    fz_corr->SetParameter(4,-0.002*100.);
    fz_corr->SetParameter(5,0.006*100.);
    fz_corr->SetParameter(6,337.);
    fz_corr->SetParameter(7,0.125*100.);
    fz_corr->SetParameter(8,-0.016/100.);
    fz_corr->SetParameter(9,-1.6*100.);
    fz_corr->SetParameter(10,0.9);

	fz_corr->SetNpx(20000);

//c7->DrawFrame(-0.25,0.1,0.25,100.);
//fz_corr->Draw("same");
	fz_cell->SetNpx(20000);
//   //==pol2gaus
//   //par[0]=x^2
//   //par[1]=x^1
//   //par[2]=x^0
//   //par[3]=Norm
//   //par[4]=Width (sigma) of convoluted Gaussian function
//   //==F_Voigt
//    // par[0] : area
//    // par[1] : location
//    // par[2] : gaussian sigma
//    // par[3] : lorentz fwhm
//	fz_corr->SetParameter(0,-200.);
//	//fz_corr->SetParLimits(0,-0.1,0.1);
//	fz_corr->SetParameter(1,-12.);
//	fz_corr->SetParameter(2,15.);
//	fz_corr->SetParameter(3,2.);
//	fz_corr->SetParameter(4,0.004);
////front
//	fz_corr->SetParameter(5,2.2);
//	fz_corr->SetParLimits(5,0.,10000.);
//	fz_corr->SetParameter(6,-0.125);
//	fz_corr->SetParameter(7,0.004);//sigma
//	fz_corr->SetParLimits(7,0.,10000.);
//	fz_corr->SetParameter(8,1.35);
//	fz_corr->SetParLimits(8,0.,10000.);
//	fz_corr->SetParameter(9,-0.125);
//	fz_corr->SetParameter(10,0.008);//sigma
//	fz_corr->SetParLimits(10,0.,10000.);
////back
//	fz_corr->SetParameter(11,1.);
//	fz_corr->SetParLimits(11,0.,10000.);
//	fz_corr->SetParameter(12,0.125);
//	fz_corr->SetParameter(13,0.004);//sigma
//	fz_corr->SetParLimits(13,0.,10000.);
//	fz_corr->SetParameter(14,2.);
//	fz_corr->SetParLimits(14,0.,10000.);
//	fz_corr->SetParameter(15,0.125);
//	fz_corr->SetParameter(16,0.008);//sigma
//	fz_corr->SetParLimits(16,0.,10000.);
	//h_zave->Draw("");
	h_zave_true->SetLineColor(kAzure);
	h_zave_true->Sumw2();
	//h_zave_dummy->Sumw2();
	h_zave_true->Scale(0.28*0.997442/40.6/0.946037);
	//h_zave_dummy->Scale(100./0.28/0.997442);
	h_zave_dummy->SetLineColor(kRed);
	h_zave_dummy->Draw("");
	h_zave_true->Draw("same");
	//h_zave_true->Fit("fz_true_front","","",-0.145,-0.105);
	//h_zave_true->Fit("fz_true_back ","","",0.105,0.145);
//	double Al_contami_true = fz_true_front->Integral(-0.1,0.1);
//	Al_contami_true+=fz_true_back->Integral(-0.1,0.1);
//	Al_contami_true/=0.0003;
//	double Al_contami_dummy = fz_dummy_front->Integral(-0.1,0.1);
//	Al_contami_dummy+=fz_dummy_back->Integral(-0.1,0.1);
//	Al_contami_dummy/=0.0003;
	//h_zave_dummy->Fit("fz_dummy_front","","",-0.145,-0.105);
	//h_zave_dummy->Fit("fz_dummy_back ","","",0.105,0.145);
	TCanvas* c8 = new TCanvas("c8","c8");
	c8->SetLogy(1);
	//h_zave->Fit("fz_corr_front","","",-0.145,-0.105);
	//h_zave->Fit("fz_corr_back ","","",0.105,0.145);
	TCanvas* c9 = new TCanvas("c9","c9");
	//h_zave_dummy->Draw("");
	h_zave->Fit("fz_corr","","",-15.,15.);
//	fz_front->SetParameter(0,fz_corr->GetParameter(5));
//	fz_front->SetParameter(1,fz_corr->GetParameter(6));
//	fz_front->SetParameter(2,fz_corr->GetParameter(7));
//	fz_front->SetParameter(3,fz_corr->GetParameter(8));
//	fz_front->SetParameter(4,fz_corr->GetParameter(9));
//	fz_front->SetParameter(5,fz_corr->GetParameter(10));
//	fz_back->SetParameter(0,fz_corr->GetParameter(11));
//	fz_back->SetParameter(1,fz_corr->GetParameter(12));
//	fz_back->SetParameter(2,fz_corr->GetParameter(13));
//	fz_back->SetParameter(3,fz_corr->GetParameter(14));
//	fz_back->SetParameter(4,fz_corr->GetParameter(15));
//	fz_back->SetParameter(5,fz_corr->GetParameter(16));
//	fz_front->SetLineColor(kViolet);
//	fz_back->SetLineColor(kViolet);
//	fz_front->SetLineStyle(2);
//	fz_back->SetLineStyle(2);
//	fz_front->Draw("same");
//	fz_back->Draw("same");
	fz_cell->SetParameter(0,fz_corr->GetParameter(0));
	fz_cell->SetParameter(1,fz_corr->GetParameter(1));
	fz_cell->SetParameter(2,fz_corr->GetParameter(2));
	fz_cell->SetParameter(3,fz_corr->GetParameter(3));
	fz_cell->SetParameter(4,fz_corr->GetParameter(4));
	fz_cell->SetParameter(5,fz_corr->GetParameter(5));
	fz_cell->SetParameter(6,fz_corr->GetParameter(6));
	fz_cell->SetParameter(7,fz_corr->GetParameter(7));
	fz_cell->SetLineColor(kViolet);
	fz_cell->SetLineStyle(2);
	fz_gas->SetLineColor(kCyan);
	fz_gas->SetLineStyle(2);

	fz_gas->SetParameter(0,fz_corr->GetParameter(2));
	fz_gas->SetParameter(1,fz_corr->GetParameter(3));
	fz_gas->SetParameter(2,fz_corr->GetParameter(4));
	fz_gas->SetParameter(3,fz_corr->GetParameter(5));
	fz_gas->SetParameter(4,fz_corr->GetParameter(8));
	fz_gas->SetParameter(5,fz_corr->GetParameter(9));
	fz_gas->SetParameter(6,fz_corr->GetParameter(10));
	
	fz_cell->Draw("same");
	fz_gas->Draw("same");
	fz_corr->Draw("same");
	//fz_corr->Draw("");
	 double chisq = fz_corr->GetChisquare();
	 double dof  = fz_corr->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

	double Al_contami_corr = fz_corr->Integral(-10.,10.);
	Al_contami_corr+=fz_corr->Integral(-10.,10.);
	Al_contami_corr/=0.05;
	cout<<"Al (gas)  = "<<Al_contami_corr<<endl;
	double Al_contami_cell = fz_cell->Integral(-10.,10.);
	Al_contami_cell+=fz_cell->Integral(-10.,10.);
	Al_contami_cell/=0.05;
	cout<<"Al (cell) = "<<Al_contami_cell<<endl;
	cout<<"Al Contamination = "<<Al_contami_cell*100./(Al_contami_cell+Al_contami_corr)<<endl;
	double Al_res1 = fz_corr->Integral(-16.,-14.5);
	double Al_res2 = fz_corr->Integral(14.5,16.);
	double Al_res_hist1 = h_zave->Integral(h_zave->FindBin(-16.),h_zave->FindBin(-14.5));
	double Al_res_hist2 = h_zave->Integral(h_zave->FindBin(14.5),h_zave->FindBin(16.));
	cout<<"Al (res from func.) = "<<(Al_res1+Al_res2)/0.05<<endl;
	cout<<"Al (res from hist.) = "<<(Al_res_hist1+Al_res_hist2)<<endl;
	cout<<"Al Contamination Systematice Error [%] = "<<((Al_res_hist1+Al_res_hist2)-((Al_res1+Al_res2)/0.05))*100./(Al_contami_cell+Al_contami_corr)<<endl;

	TH1F* h_den  = new TH1F("h_den","denominator",25,0.,25.);
	TH1F* h_num  = new TH1F("h_num","numerator",25,0.,25.);
	double den, num;
	h_den->SetBinContent(1,0.);
	h_num->SetBinContent(1,0.);
	for(int i=1;i<25;i++){
	//den = h_zave->Integral(h_zave->FindBin((double)i*(-0.01)),h_zave->FindBin((double)i*(0.01)));
	//den = h_zave->Integral(h_zave->FindBin(-25.),h_zave->FindBin(0.25));
	den = fz_gas->Integral(-25.,25.)/0.05;
	num = fz_gas->Integral((double)i*(-1.),(double)i)/0.05;
	h_den->SetBinContent(i+1,den);
	h_num->SetBinContent(i+1,num);
	}

	TH1F* h_den2  = new TH1F("h_den2","denominator",25,0.,25.);
	TH1F* h_num2  = new TH1F("h_num2","numerator",25,0.,25.);
	double den2, num2;
	h_den2->SetBinContent(1,0.);
	h_num2->SetBinContent(1,0.);
	for(int i=1;i<25;i++){
	//den = h_zave->Integral(h_zave->FindBin((double)i*(-0.01)),h_zave->FindBin((double)i*(0.01)));
	//den = h_zave->Integral(h_zave->FindBin(-25.),h_zave->FindBin(0.25));
	den2 = fz_corr->Integral((double)i*(-1.),(double)i)/0.05;
	num2 = fz_cell->Integral((double)i*(-1.),(double)i)/0.05;
	h_den2->SetBinContent(i+1,den2);
	h_num2->SetBinContent(i+1,num2);
	}
	cout << "TEfficiency! (Gas SR)" << endl;
	TEfficiency *pEff1;
	if(TEfficiency::CheckConsistency(*h_num,*h_den,"w")){
	pEff1 = new TEfficiency(*h_num,*h_den);
	}
	cout << "TEfficiency2! (Al Contamination)" << endl;
	TEfficiency *pEff2;
	if(TEfficiency::CheckConsistency(*h_num2,*h_den2,"w")){
	pEff2 = new TEfficiency(*h_num2,*h_den2);
	}
	TCanvas* c10 = new TCanvas("c10","c10");
	pEff1->Draw("");
	TCanvas* c11 = new TCanvas("c11","c11");
	pEff1->Draw("");
	pEff2->SetLineColor(kRed);
	pEff2->Draw("same");
//	h_zz_dummy->Draw("colz");

cout << "Well done!" << endl;
}//fit
