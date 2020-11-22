//----------------------------------//
//--  Fitting w/ Response func.   --//
//--  Differential Cross Section  --//
//----------------------------------//
//
//K. Okuyama (Nov. 21, 2020)
//
//This is taken over from result_2D.C
//Use 2D Acceptance Map (Z, pK)
//Lab => CM (event by event)
//No array branch mode 


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

void result_2D_Nov(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("h2all_2020Nov.root","read");//2020Nov updated
  //TFile *file = new TFile("h2all_Lsingle.root","read");//before 2020Nov
	//ACCBGの引き算はmea_hist.ccから
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_llccrr_new_new.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_2020Nov.root","read");// 2020/11/19 rootfile
  TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_csbase.root","read");// 2020/11/22 rootfile
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_llccrr_Lsingle.root","read");// 2020/11/19 rootfile
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_mthesis.root","read");//from h2all_Lsingle.root, HRS-L: Single-tracking
  double nbunch = 6000.;//effetive bunches (6 bunches x 5 mixtures)
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
//	tree->Write();


//---  DAQ Efficiency ---//
//H2 run (run111157~111222 & run111480~542)
	string daq_file = "./information/daq.dat";//DAQ Efficiency from ELOG
	int runnum;
	double daq_eff;
	double daq_eff_total=0.;
	int daq_eff_bin=0;
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
		daq_eff_total+=daq_eff;
		daq_eff_bin++;
	}
	cout<<"DAQ Efficiency (average)="<<daq_eff_total/(double)daq_eff_bin<<endl;

//----------------HRS-R Acceptance-----------------//

	int RHRS_bin;
	double RHRS_SIMC;
	double RHRS_table[150][5];//1.5<pk[GeV/c]<2.1, 150 partition --> 1bin=4MeV/c
	double RHRS_total=0.;
	int RHRS_total_bin=0;
/*----- -10 < z < -6 -----*/
	string AcceptanceR_table_z1 = "./information/RHRS_SIMC_z1.dat";//Acceptance Table (SIMC)
	string buf_z1;

	ifstream ifp_z1(AcceptanceR_table_z1.c_str(),ios::in);
	if (ifp_z1.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z1.c_str() << endl;
	while(1){
		getline(ifp_z1,buf_z1);
		if(buf_z1[0]=='#'){continue;}
		if(ifp_z1.eof())break;
		stringstream sbuf_z1(buf_z1);
		sbuf_z1 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][0] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z1 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -6 < z < -2 -----*/
	string AcceptanceR_table_z2 = "./information/RHRS_SIMC_z2.dat";//Acceptance Table (SIMC)
	string buf_z2;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z2(AcceptanceR_table_z2.c_str(),ios::in);
	if (ifp_z2.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z2.c_str() << endl;
	while(1){
		getline(ifp_z2,buf_z2);
		if(buf_z2[0]=='#'){continue;}
		if(ifp_z2.eof())break;
		stringstream sbuf_z2(buf_z2);
		sbuf_z2 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][1] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z2 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -2 < z < 2 -----*/
	string AcceptanceR_table_z3 = "./information/RHRS_SIMC_z3.dat";//Acceptance Table (SIMC)
	string buf_z3;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z3(AcceptanceR_table_z3.c_str(),ios::in);
	if (ifp_z3.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z3.c_str() << endl;
	while(1){
		getline(ifp_z3,buf_z3);
		if(buf_z3[0]=='#'){continue;}
		if(ifp_z3.eof())break;
		stringstream sbuf_z3(buf_z3);
		sbuf_z3 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][2] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z3 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 2 < z < 6 -----*/
	string AcceptanceR_table_z4 = "./information/RHRS_SIMC_z4.dat";//Acceptance Table (SIMC)
	string buf_z4;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z4(AcceptanceR_table_z4.c_str(),ios::in);
	if (ifp_z4.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z4.c_str() << endl;
	while(1){
		getline(ifp_z4,buf_z4);
		if(buf_z4[0]=='#'){continue;}
		if(ifp_z4.eof())break;
		stringstream sbuf_z4(buf_z4);
		sbuf_z4 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][3] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z4 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 6 < z < 10 -----*/
	string AcceptanceR_table_z5 = "./information/RHRS_SIMC_z5.dat";//Acceptance Table (SIMC)
	string buf_z5;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z5(AcceptanceR_table_z5.c_str(),ios::in);
	if (ifp_z5.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z5.c_str() << endl;
	while(1){
		getline(ifp_z5,buf_z5);
		if(buf_z5[0]=='#'){continue;}
		if(ifp_z5.eof())break;
		stringstream sbuf_z5(buf_z5);
		sbuf_z5 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][4] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z5 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;


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
 const double fit_max_mm=0.085;
 const double fmin_mm=-0.01;
 const double fmax_mm=0.12;
 //const double fmax_mm=10.;
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
  TH1F* hmm_L_fom_strict  = new TH1F("hmm_L_fom_strict","hmm_L_fom_strict",xbin,xmin,xmax);
  hmm_L_fom_strict->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_L_fom_strict->GetYaxis()->SetTitle("Counts / MeV");
  hmm_L_fom_strict->SetLineColor(1);
  TH1F* hcs_L_fom_best  = new TH1F("hcs_L_fom_best","hcs_L_fom_best",xbin,xmin,xmax);
  hcs_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} [GeV/c^{2}]");
  hcs_L_fom_best->GetYaxis()->SetTitle("d#sigma/d#Omega (C.M.F.) [nb/sr]");
  hcs_L_fom_best->SetLineColor(1);
  TH1F* hcs_L_fom_strict  = new TH1F("hcs_L_fom_strict","hcs_L_fom_strict",xbin,xmin,xmax);
  hcs_L_fom_strict->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} [GeV/c^{2}]");
  hcs_L_fom_strict->GetYaxis()->SetTitle("d#sigma/d#Omega (C.M.F.) [nb/sr]");
  hcs_L_fom_strict->SetLineColor(1);
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","RESULT (N_{#Lambda}/#Delta#Omega_{K}/#epsilon_{DAQ})",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");
  TH2F* gklab_gkcm = new TH2F("gklab_gkcm","#theta_{gk}^{lab} vs #theta_{gk}^{CM}",100,0.,0.15,100,0.,0.4);
  TH2F* gklab_eklab = new TH2F("gklab_eklab","#theta_{gk}^{lab} vs #theta_{ek}^{lab}",100,0.,0.15,100,0.15,0.3);
  TH2F* eelab_eklab = new TH2F("eelab_eklab","#theta_{ee}^{lab} vs #theta_{ek}^{lab}",100,0.15,0.3,100,0.15,0.3);
  TH2F* eklab_gkcm = new TH2F("eklab_gkcm","#theta_{ek}^{lab} vs #theta_{gk}^{CM}",100,0.15,0.3,100,0.,0.4);

  TH1D* h_theta_ee = new TH1D("h_theta_ee", "theta_ee",1000,0.1,0.35);
  TH1D* h_phi_ee = new TH1D("h_phi_ee", "phi_ee",1000,0.,PI);
  TH1D* h_theta_ek = new TH1D("h_theta_ek", "theta_ek",1000,0.1,0.35);
  TH1D* h_phi_ek = new TH1D("h_phi_ek", "phi_ek",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_g = new TH1D("h_theta_g", "theta_g",1000,0.1,0.35);
  TH1D* h_phi_g = new TH1D("h_phi_g", "phi_g",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_gk_lab = new TH1D("h_theta_gk_lab", "theta_gk_lab",1000,0.,0.2);
  TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "theta_gk_cm",100,0.,20.);
  TH1D* h_theta_gk_cm2 = new TH1D("h_theta_gk_cm2", "theta_gk_cm2",100,0.,20.);
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

/*-------------------*/
/*--Angle partition--*/
/*-------------------*/
  TH1F* hcs_L_cm2_1  = new TH1F("hcs_L_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm2_2  = new TH1F("hcs_L_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm3_1  = new TH1F("hcs_L_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm3_2  = new TH1F("hcs_L_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm3_3  = new TH1F("hcs_L_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm4_1  = new TH1F("hcs_L_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm4_2  = new TH1F("hcs_L_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm4_3  = new TH1F("hcs_L_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hcs_L_cm4_4  = new TH1F("hcs_L_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm2_1 = new TH1F("hmm_wo_bg_cm2_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm2_2 = new TH1F("hmm_wo_bg_cm2_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm3_1 = new TH1F("hmm_wo_bg_cm3_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm3_2 = new TH1F("hmm_wo_bg_cm3_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm3_3 = new TH1F("hmm_wo_bg_cm3_3","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm4_1 = new TH1F("hmm_wo_bg_cm4_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm4_2 = new TH1F("hmm_wo_bg_cm4_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm4_3 = new TH1F("hmm_wo_bg_cm4_3","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_cm4_4 = new TH1F("hmm_wo_bg_cm4_4","",xbin,xmin,xmax);

  TH1F* hcs_L_new_cm2_1  = new TH1F("hcs_L_new_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm2_2  = new TH1F("hcs_L_new_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_1  = new TH1F("hcs_L_new_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_2  = new TH1F("hcs_L_new_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_3  = new TH1F("hcs_L_new_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_1  = new TH1F("hcs_L_new_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_2  = new TH1F("hcs_L_new_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_3  = new TH1F("hcs_L_new_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_4  = new TH1F("hcs_L_new_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm2_1 = new TH1F("hmm_wo_bg_new_cm2_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm2_2 = new TH1F("hmm_wo_bg_new_cm2_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm3_1 = new TH1F("hmm_wo_bg_new_cm3_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm3_2 = new TH1F("hmm_wo_bg_new_cm3_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm3_3 = new TH1F("hmm_wo_bg_new_cm3_3","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm4_1 = new TH1F("hmm_wo_bg_new_cm4_1","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm4_2 = new TH1F("hmm_wo_bg_new_cm4_2","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm4_3 = new TH1F("hmm_wo_bg_new_cm4_3","",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_new_cm4_4 = new TH1F("hmm_wo_bg_new_cm4_4","",xbin,xmin,xmax);
  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;//best=loose
  bool event_selection_new = false;//strict=tight
  bool event_selection_nocut = false;
  bool cm2_angle1_cut=false;
  bool cm2_angle2_cut=false;
  bool cm3_angle1_cut=false;
  bool cm3_angle2_cut=false;
  bool cm3_angle3_cut=false;
  bool cm4_angle1_cut=false;
  bool cm4_angle2_cut=false;
  bool cm4_angle3_cut=false;
  bool cm4_angle4_cut=false;
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
		double effDAQ= 0.950;
		double efftr = 0.810;
		double effK	 = 0.170;
		double efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		double RHRS  = 0.005;
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
	double csL[xbin];



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

		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_new=true;
		else event_selection_new=false;

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
		if(event_selection_new&&ct_cut)hmm_L_fom_strict->Fill(mm);
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



		h_theta_gk_cm->Fill(theta_gk_cm*180./PI);
		if(L_mom>2.1)h_theta_gk_cm2->Fill(theta_gk_cm*180./PI);

		cm2_angle1_cut=false;
   		cm2_angle2_cut=false;
   		cm3_angle1_cut=false;
   		cm3_angle2_cut=false;
   		cm3_angle3_cut=false;
   		cm4_angle1_cut=false;
   		cm4_angle2_cut=false;
   		cm4_angle3_cut=false;
   		cm4_angle4_cut=false;
//======= CM Angle(gamma-K) ========//
		if(theta_gk_cm*180./PI<8.)cm2_angle1_cut=true;
		if(theta_gk_cm*180./PI>=8.)cm2_angle2_cut=true;
		if(theta_gk_cm*180./PI<6.)cm3_angle1_cut=true;
		if(theta_gk_cm*180./PI>=6. && theta_gk_cm*180./PI<10.)cm3_angle2_cut=true;
		if(theta_gk_cm*180./PI>=10.)cm3_angle3_cut=true;
		if(theta_gk_cm*180./PI<5.)cm4_angle1_cut=true;
		if(theta_gk_cm*180./PI>=5. && theta_gk_cm*180./PI<8.)cm4_angle2_cut=true;
		if(theta_gk_cm*180./PI>=8. && theta_gk_cm*180./PI<11.)cm4_angle3_cut=true;
		if(theta_gk_cm*180./PI>=11.)cm4_angle4_cut=true;
//======= CM Angle(gamma-K) ========//

		//int ebin = (int)((L_mom-1.8)/0.004);
		//if(ebin>=0 &&ebin<150){
		//Ng = Ng_table[ebin]*Ne;//
		//if(Ng!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
		//else cs=0.;
		//}else{cs=0.;} 
		int kbin = (int)((R_mom-1.5)/0.004);
		int zbin = (int)((((L_tr_vz+R_tr_vz)/2.)-0.1)/0.04);
		if(event_selection&&ct_cut&&kbin>=0 &&kbin<150){
		if((L_tr_vz+R_tr_vz)/2.>-0.1&&(L_tr_vz+R_tr_vz)/2.<-0.06)RHRS = RHRS_table[kbin][0];
		else if((L_tr_vz+R_tr_vz)/2.>-0.06&&(L_tr_vz+R_tr_vz)/2.<-0.02)RHRS = RHRS_table[kbin][1];
		else if((L_tr_vz+R_tr_vz)/2.>-0.02&&(L_tr_vz+R_tr_vz)/2.<0.02)RHRS = RHRS_table[kbin][2];
		else if((L_tr_vz+R_tr_vz)/2.>0.02&&(L_tr_vz+R_tr_vz)/2.<0.06)RHRS = RHRS_table[kbin][3];
		else if((L_tr_vz+R_tr_vz)/2.>0.06&&(L_tr_vz+R_tr_vz)/2.<0.1)RHRS = RHRS_table[kbin][4];
		else cout<<"Z Error"<<(L_tr_vz+R_tr_vz)/2.<<endl;

		effDAQ = daq_table[nrun-111000];
		if(effDAQ==0.2)cout<<"Starange!!! DAQ Eff. of run"<<nrun<<" does not exist."<<endl;
		efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		//if(RHRS!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
		if(RHRS!=0.&&effDAQ!=0.)cs = labtocm/effDAQ/RHRS/10.;//[nb/sr]
		else cs=0.;
		}else{cs=0.;}
		double cs_temp = cs*labtocm;
		//cs_ave += cs;
		//cout<<"Ng="<<Ng<<endl;
		//cout<<"cs="<<cs<<endl;
		if(event_selection&&ct_cut)hcs_L_fom_best->Fill(mm,cs);
		if(event_selection_new&&ct_cut)hcs_L_fom_strict->Fill(mm,cs);
		//if(event_selection&&ct_cut)hcs_L_fom_best->Fill(mm,cs);
//		//if(event_selection_nocut&&ct_cut)hcs_L_fom_nocut->SetBinContent(hmm_L_fom_best->FindBin(mm),cs);

		if(event_selection&&ct_cut){
			gklab_gkcm->Fill(theta_gk_lab,theta_gk_cm);
			gklab_eklab->Fill(theta_gk_lab,theta_ek);
			eelab_eklab->Fill(theta_ee,theta_ek);
			eklab_gkcm->Fill(theta_ek,theta_gk_cm);
			if(cm2_angle1_cut)hcs_L_cm2_1->Fill(mm,cs);
			if(cm2_angle2_cut)hcs_L_cm2_2->Fill(mm,cs);
			if(cm3_angle1_cut)hcs_L_cm3_1->Fill(mm,cs);
			if(cm3_angle2_cut)hcs_L_cm3_2->Fill(mm,cs);
			if(cm3_angle3_cut)hcs_L_cm3_3->Fill(mm,cs);
			if(cm4_angle1_cut)hcs_L_cm4_1->Fill(mm,cs);
			if(cm4_angle2_cut)hcs_L_cm4_2->Fill(mm,cs);
			if(cm4_angle3_cut)hcs_L_cm4_3->Fill(mm,cs);
			if(cm4_angle4_cut)hcs_L_cm4_4->Fill(mm,cs);
		}
		if(event_selection_new&&ct_cut){
			if(cm2_angle1_cut)hcs_L_new_cm2_1->Fill(mm,cs);
			if(cm2_angle2_cut)hcs_L_new_cm2_2->Fill(mm,cs);
			if(cm3_angle1_cut)hcs_L_new_cm3_1->Fill(mm,cs);
			if(cm3_angle2_cut)hcs_L_new_cm3_2->Fill(mm,cs);
			if(cm3_angle3_cut)hcs_L_new_cm3_3->Fill(mm,cs);
			if(cm4_angle1_cut)hcs_L_new_cm4_1->Fill(mm,cs);
			if(cm4_angle2_cut)hcs_L_new_cm4_2->Fill(mm,cs);
			if(cm4_angle3_cut)hcs_L_new_cm4_3->Fill(mm,cs);
			if(cm4_angle4_cut)hcs_L_new_cm4_4->Fill(mm,cs);
		}


}//ENum
	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
	//hcs_L_fom_best->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2");
	TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	TH1F* hmm_Al_fom_nocut=(TH1F*)file->Get("hmm_Al_fom_best");
	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	TH1F* hmm_bg_cm2_1=(TH1F*)file_mea->Get("hmm_mixacc_result_cm2_1");
	TH1F* hmm_bg_cm2_2=(TH1F*)file_mea->Get("hmm_mixacc_result_cm2_2");
	TH1F* hmm_bg_cm3_1=(TH1F*)file_mea->Get("hmm_mixacc_result_cm3_1");
	TH1F* hmm_bg_cm3_2=(TH1F*)file_mea->Get("hmm_mixacc_result_cm3_2");
	TH1F* hmm_bg_cm3_3=(TH1F*)file_mea->Get("hmm_mixacc_result_cm3_3");
	TH1F* hmm_bg_cm4_1=(TH1F*)file_mea->Get("hmm_mixacc_result_cm4_1");
	TH1F* hmm_bg_cm4_2=(TH1F*)file_mea->Get("hmm_mixacc_result_cm4_2");
	TH1F* hmm_bg_cm4_3=(TH1F*)file_mea->Get("hmm_mixacc_result_cm4_3");
	TH1F* hmm_bg_cm4_4=(TH1F*)file_mea->Get("hmm_mixacc_result_cm4_4");
	TH1F* hmm_bg_new_cm2_1=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm2_1");
	TH1F* hmm_bg_new_cm2_2=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm2_2");
	TH1F* hmm_bg_new_cm3_1=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm3_1");
	TH1F* hmm_bg_new_cm3_2=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm3_2");
	TH1F* hmm_bg_new_cm3_3=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm3_3");
	TH1F* hmm_bg_new_cm4_1=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm4_1");
	TH1F* hmm_bg_new_cm4_2=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm4_2");
	TH1F* hmm_bg_new_cm4_3=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm4_3");
	TH1F* hmm_bg_new_cm4_4=(TH1F*)file_mea->Get("hmm_mixacc_result_new_cm4_4");
	//TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
	//TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
	TH1F* hmm_bg_fom_strict=(TH1F*)file_mea->Get("hmm_mixacc_result_new");
	TH1F* hcs_bg_fom_strict=(TH1F*)file_mea->Get("hcs_mixacc_result_new");
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
	hmm_bg_fom_strict->Sumw2();
	hmm_bg_fom_nocut->Sumw2();
	hcs_bg_fom_strict->Sumw2();
	hmm_Albg_fom_nocut->Sumw2();
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_strict->Scale(1./nbunch);
	hcs_bg_fom_strict->Scale(1./nbunch);
	//cs	 = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
	effDAQ= 0.950;
	RHRS  = 0.0055;
	double labtocm = 0.126;
	cs	 = labtocm/(effDAQ*RHRS*10.);
	hmm_bg_fom_best->Scale(cs);
	hmm_bg_fom_strict->Scale(cs);
	hmm_bg_cm2_1->Scale(1./nbunch/cs);
	hmm_bg_cm2_2->Scale(1./nbunch/cs);
	hmm_bg_cm3_1->Scale(1./nbunch/cs);
	hmm_bg_cm3_2->Scale(1./nbunch/cs);
	hmm_bg_cm3_3->Scale(1./nbunch/cs);
	hmm_bg_cm4_1->Scale(1./nbunch/cs);
	hmm_bg_cm4_2->Scale(1./nbunch/cs);
	hmm_bg_cm4_3->Scale(1./nbunch/cs);
	hmm_bg_cm4_4->Scale(1./nbunch/cs);
	hmm_bg_new_cm2_1->Scale(1./nbunch/cs);
	hmm_bg_new_cm2_2->Scale(1./nbunch/cs);
	hmm_bg_new_cm3_1->Scale(1./nbunch/cs);
	hmm_bg_new_cm3_2->Scale(1./nbunch/cs);
	hmm_bg_new_cm3_3->Scale(1./nbunch/cs);
	hmm_bg_new_cm4_1->Scale(1./nbunch/cs);
	hmm_bg_new_cm4_2->Scale(1./nbunch/cs);
	hmm_bg_new_cm4_3->Scale(1./nbunch/cs);
	hmm_bg_new_cm4_4->Scale(1./nbunch/cs);
	hmm_bg_fom_nocut->Scale(1./nbunch);
	hmm_Albg_fom_nocut->Scale(1./nbunch);
	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
	//hmm_wo_bg_fom_best->Add(hcs_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
	


//MM spctrum to be fitted (2020/10/18)
//Choose one from the list below
//===CHANGE===//
//Loose Cut
	//hmm_wo_bg_fom_best->Add(hcs_L_fom_best,hmm_bg_fom_best,1.0,-1.0);//All
	//hmm_wo_bg_fom_best->Add(hcs_L_cm2_1,hmm_bg_cm2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm2_2,hmm_bg_cm2_2,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm3_1,hmm_bg_cm3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm3_2,hmm_bg_cm3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm3_3,hmm_bg_cm3_3,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm4_1,hmm_bg_cm4_1,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm4_2,hmm_bg_cm4_2,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm4_3,hmm_bg_cm4_3,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_cm4_4,hmm_bg_cm4_4,1.0,-1.0);//4 div.
	
//Tight Cut	
	//hmm_wo_bg_fom_best->Add(hcs_L_fom_strict,hmm_bg_fom_strict,1.0,-1.0);//All
	hmm_wo_bg_fom_best->Add(hcs_L_fom_strict,hcs_bg_fom_strict,1.0,-1.0);//All by hcs
	//hmm_wo_bg_fom_best->Add(hcs_L_fom_strict,hmm_bg_fom_strict,1.0,-1.0);//All
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm2_1,hmm_bg_new_cm2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm2_2,hmm_bg_new_cm2_2,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm3_1,hmm_bg_new_cm3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm3_2,hmm_bg_new_cm3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm3_3,hmm_bg_new_cm3_3,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm4_1,hmm_bg_new_cm4_1,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm4_2,hmm_bg_new_cm4_2,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm4_3,hmm_bg_new_cm4_3,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_best->Add(hcs_L_new_cm4_4,hmm_bg_new_cm4_4,1.0,-1.0);//4 div.
//===CHANGE===//

	
	for(int i=90;i<200;i++){
	if(hcs_L_cm3_1->GetBinContent(i)==0)cout<<"Empty bin at "<<((double)i*0.001-0.1)<<endl;
	}

	hmm_wo_bg_fom_nocut->Add(hmm_L_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	//hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
	//hmm_pi_wobg_fom_nocut->Add(hmm_pi_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	hmm_pi_wobg_fom_nocut->Add(hmm_Al_fom_nocut,hmm_Albg_fom_nocut,1.0,-1.0);
	//hmm_wo_bg_cm2_1->Add(hcs_L_cm2_1,hmm_bg_cm2_1,1.0,-1.0);
	//hmm_wo_bg_cm2_2->Add(hcs_L_cm2_2,hmm_bg_cm2_2,1.0,-1.0);
	//hmm_wo_bg_cm3_1->Add(hcs_L_cm3_1,hmm_bg_cm3_1,1.0,-1.0);
	//hmm_wo_bg_cm3_2->Add(hcs_L_cm3_2,hmm_bg_cm3_2,1.0,-1.0);
	//hmm_wo_bg_cm3_3->Add(hcs_L_cm3_3,hmm_bg_cm3_3,1.0,-1.0);
	//hmm_wo_bg_cm4_1->Add(hcs_L_cm4_1,hmm_bg_cm4_1,1.0,-1.0);
	//hmm_wo_bg_cm4_2->Add(hcs_L_cm4_2,hmm_bg_cm4_2,1.0,-1.0);
	//hmm_wo_bg_cm4_3->Add(hcs_L_cm4_3,hmm_bg_cm4_3,1.0,-1.0);
	//hmm_wo_bg_cm4_4->Add(hcs_L_cm4_4,hmm_bg_cm4_4,1.0,-1.0);


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
	 double fmin = -0.01;
	 double fmax =  0.12;
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
	 fmm_best_4Poly=new TF1("fmm_best_4Poly",FMM_Res,fmin_mm,fmax_mm,14);
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
	 fmmbg_best_4Poly=new TF1("fmmbg_best_4Poly","pol4",fmin_mm,fmax_mm);
	 fmm_best_4Poly->SetNpx(20000);
	 fmm_best_4Poly->SetTitle("Missing Mass (best)");
	 fmm_best_4Poly->SetParLimits(2,0.,1000.);//positive
	 fmm_best_4Poly->SetParLimits(9,0.,300.);//positive
	 fmm_best_4Poly->SetParameter(0,0.0007);//Landau width
	 fmm_best_4Poly->SetParameter(1,mean_L_best);
	 fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_best_4Poly->SetParameter(2,1.5);//total scale
	 fmm_best_4Poly->SetParameter(3,0.001);//sigma
	 fmm_best_4Poly->SetParLimits(3,0.,0.01);
	 fmm_best_4Poly->SetParameter(4,0.05);//att.
	 fmm_best_4Poly->SetParLimits(4,0.005,0.08);
	 fmm_best_4Poly->SetParameter(5,-0.004);//peak pos.
	 fmm_best_4Poly->SetParLimits(5,-0.05,0.05);
	 fmm_best_4Poly->SetParameter(6,0.6);//relative strength
	 fmm_best_4Poly->SetParLimits(6,0.,1.5);//relative strength

	 fmm_best_4Poly->SetParameter(7,0.0003);//Landau width
	 fmm_best_4Poly->SetParameter(8,mean_S_best);//MPV
	 fmm_best_4Poly->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_best_4Poly->SetParameter(9,0.4);//total scale
	 fmm_best_4Poly->SetParameter(10,0.0015);//sigma
	 fmm_best_4Poly->SetParLimits(10,0.,0.01);
	 fmm_best_4Poly->SetParameter(11,0.05);//att
	 fmm_best_4Poly->SetParLimits(11,0.03,0.12);
	 fmm_best_4Poly->SetParameter(12,0.080);//peak pos.
	 //fmm_best_4Poly->SetParLimits(15,-0.085,-0.055);
	 fmm_best_4Poly->SetParameter(13,0.6);
	 fmm_best_4Poly->SetParLimits(13,0.,1.5);//relative strength

	 hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
	 double chisq = fmm_best_4Poly->GetChisquare();
	 double dof  = fmm_best_4Poly->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;


	 TF1* fmm_Lambda_only = new TF1("fmm_Lambda_only",FMM_Response,fmin_mm,fmax_mm,7);
	 TF1* fmm_Sigma_only  = new TF1("fmm_Sigma_only" ,FMM_Response, fmin_mm,fmax_mm,7);
//Lambda_only
	 fmm_Lambda_only->SetNpx(20000);
	 fmm_Lambda_only->SetParameter(0,fmm_best_4Poly->GetParameter(0));
	 fmm_Lambda_only->SetParameter(1,fmm_best_4Poly->GetParameter(1));
	 fmm_Lambda_only->SetParameter(2,fmm_best_4Poly->GetParameter(2));
	 fmm_Lambda_only->SetParameter(3,fmm_best_4Poly->GetParameter(3));
	 fmm_Lambda_only->SetParameter(4,fmm_best_4Poly->GetParameter(4));
	 fmm_Lambda_only->SetParameter(5,fmm_best_4Poly->GetParameter(5));
	 fmm_Lambda_only->SetParameter(6,fmm_best_4Poly->GetParameter(6));
//Sigma_only
	 fmm_Sigma_only->SetNpx(20000);
	 fmm_Sigma_only->SetParameter(0,fmm_best_4Poly->GetParameter(7));
	 fmm_Sigma_only->SetParameter(1,fmm_best_4Poly->GetParameter(8));
	 fmm_Sigma_only->SetParameter(2,fmm_best_4Poly->GetParameter(9));
	 fmm_Sigma_only->SetParameter(3,fmm_best_4Poly->GetParameter(10));
	 fmm_Sigma_only->SetParameter(4,fmm_best_4Poly->GetParameter(11));
	 fmm_Sigma_only->SetParameter(5,fmm_best_4Poly->GetParameter(12));
	 fmm_Sigma_only->SetParameter(6,fmm_best_4Poly->GetParameter(13));

	double nofL = fmm_Lambda_only->Integral(fmin_mm,fmax_mm);
	double nofL_old = fmm_Lambda_only->Integral(-0.006,0.006);
	nofL = nofL/fit_bin_width;
	nofL_old = nofL_old/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(-0.006),hmm_wo_bg_fom_best->FindBin(0.006))<<endl;

	double nofS = fmm_Sigma_only->Integral(fmin_mm,fmax_mm);
	nofS = nofS/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS<<endl;

	 hmm_wo_bg_fom_best->Draw();
	 fmm_best_4Poly->SetLineColor(kRed);
	 fmm_best_4Poly->Draw("same");
	 fmm_Lambda_only->SetLineColor(kAzure);
	 fmm_Sigma_only->SetLineColor(kCyan);
	 fmm_Lambda_only->Draw("same");
	 fmm_Sigma_only->Draw("same");
	 //fL_best->SetLineColor(kGreen);
	 //fS_best->SetLineColor(kGreen);
	 //fL_best->Draw("same");
	 //fS_best->Draw("same");

	TCanvas* c3 = new TCanvas("c3","c3");
	c3->Divide(2,2);
	c3->cd(1);
	gklab_gkcm->Draw("colz");
	c3->cd(2);
	gklab_eklab->Draw("colz");
	c3->cd(3);
	eelab_eklab->Draw("colz");
	c3->cd(4);
	eklab_gkcm->Draw("colz");

	TCanvas* c4 = new TCanvas("c4","c4");
	//c4->Divide(2,2);
	//c4->cd(1);
	h_theta_gk_cm->Draw("");
	h_theta_gk_cm2->SetLineColor(kRed);
	h_theta_gk_cm2->Draw("same");

	TCanvas* c5 = new TCanvas("c5","c5");
	hmm_L_fom_best->SetLineColor(kBlack);
	hmm_L_fom_best->Draw("");
	hmm_bg_fom_best->SetLineColor(kGreen);
	hmm_bg_fom_best->Draw("same");
	TCanvas* c6 = new TCanvas("c6","c6");
	hmm_L_fom_strict->SetLineColor(kBlack);
	hmm_L_fom_strict->Draw("");
	hmm_bg_fom_strict->SetLineColor(kGreen);
	hmm_bg_fom_strict->Draw("same");
cout << "Well done!" << endl;
}//fit
