//----------------------------------//
//--  Fitting w/ Response func.   --//
//--  Differential Cross Section  --//
//----------------------------------//
//
//This is "result_2D_2023"
//
//K. Okuyama (Nov. 21, 2020)
//K. Okuyama (Dec. 14, 2020)//Mom cut
//K. Okuyama (May. 11, 2021)//Fitting w/ pion, Al
//K. Okuyama (Jan. 28, 2022)//Kaon SR event by event (pick up from PathLength/mom_calib.C)
//K. Okuyama (Feb. 1, 2022)//change it readable
//K. Okuyama (Oct. 7, 2023)//SIMC:4318.5-->4313 GeV;after Acc./MEA also updated
//
//This is taken over from result_2D.C
//	but directly taken over from result_2D_Jan.C
//Use 2D Acceptance Map (Z, pK)
//Lab => CM (event by event)
//Kaon survival ratio (event by event)
//subtracting Al and Pion contamination
//No array branch mode 
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetDecimals(2);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetDecimals(2);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//Voigt Function
double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

//Get PathLength from Cointime.
double GetPathLenfromCt(double *x, double *par){//Path Length from Cointime.
	double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
	double Mp = 0.938272046;         // proton       mass (GeV/c2)
	double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
	double mom = x[0]/1000.;//GeV/c
	double b = par[0];//intercept
	double a = par[1];//slope
	double tdiff = a*mom*1000.+b;//pol1
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double pathlen = betapi*betap*LightVelocity*tdiff/(betapi-betap);
	return pathlen;//m
}

//Get Kaon Survival Ratio from Cointime.
double GetSRfromCt(double *x, double *par){//Survival Ratio with Path Length (from Cointime.)
	double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
	double Mp = 0.938272046;         // proton       mass (GeV/c2)
	double MK = 0.493677;            // charged Kaon mass (GeV/c2)
	double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
	double mom = x[0]/1000.;//GeV/c
	double b = par[0];//intercept
	double a = par[1];//slope
	double tdiff = a*mom*1000.+b;//pol1
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double pathlen = betapi*betap*LightVelocity*tdiff/(betapi-betap);
	double sr = exp(-1.*pathlen*MK/(mom*3.71));//exp(-LM/(pct))
	return sr;
}

//Fitting function for Lambda/Sigma0 including radiative tail
//(Landau+Exp)*Gaus
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

  return par[2]*(val1+par[6]*val2)/(1.+par[6]);//N x (Landau*Gauss) + N' x (Exp*Gauss)

}//FMM_Response

//Global fitting function: Lambda+Sigma0
double FMM_Res( double *x, double *par ){

	return FMM_Response(x,par)+FMM_Response(x,&par[7]);

}

//Global fitting function: Lambda+Sigma0+B.G.(Al+pion)
double FMM_Res_nocut( double *x, double *par ){

	return FMM_Response(x,par)+FMM_Response(x,&par[7])+par[21]*par[14] * TMath::Voigt(x[0]-par[15],par[16],par[17],4)+(1.-par[21])*par[14] * TMath::Voigt(x[0]-par[18],par[19],par[20],4);

}

//Fitting function: B.G.(Al+pion)
double FMM_2BG( double *x, double *par ){

	return par[7]*par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4)+(1.-par[7])*par[0] * TMath::Voigt(x[0]-par[4],par[5],par[6],4);

}

void result_2D_2023Qsq2(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("h2all_2020Nov.root","read");//2020Nov updated
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_momDec.root","read");// 2020/12/14 Mom cut 
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_2021Jan.root","read");// 2021/1/4 Mom cut 
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_llccrr_new_2022effK.root","read");// 2022/1/30 Mom cut & effK in cs
  TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_llccrr_new_2023.root","read");// 2023/10/6 Mom cut & effK in cs & 4.313 GeV & nmix750
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea_llccrr_new_tripleDCS2022.root","read");// 2022/6/5 Mom cut & effK in cs && LHRS in cs
  double nbunch = 4500.;//effetive bunches (6 bunches x 750 mixtures)
  TTree *tree = (TTree*)file->Get("tree_out");

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
	double RHRS_table[100][10];//1.6<pk[GeV/c]<2.0, 100 partition --> 1bin=4MeV/c
							   //-10<Z-vertex <10,   10 partition --> 1bin=2cm
	double RHRS_total=0.;
	int RHRS_total_bin=0;
/*----- -10 < z < -8 -----*/
	string AcceptanceR_table_z0 = "./information/RHRS_SIMC2023_10_z0.dat";//Acceptance Table (SIMC)
	string buf_z0;

	ifstream ifp_z0(AcceptanceR_table_z0.c_str(),ios::in);
	if (ifp_z0.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z0.c_str() << endl;
	while(1){
		getline(ifp_z0,buf_z0);
		if(buf_z0[0]=='#'){continue;}
		if(ifp_z0.eof())break;
		stringstream sbuf_z0(buf_z0);
		sbuf_z0 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][0] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z0 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -8 < z < -6 -----*/
	string AcceptanceR_table_z1 = "./information/RHRS_SIMC2023_10_z1.dat";//Acceptance Table (SIMC)
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

		RHRS_table[RHRS_bin-1][1] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z1 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -6 < z < -4 -----*/
	string AcceptanceR_table_z2 = "./information/RHRS_SIMC2023_10_z2.dat";//Acceptance Table (SIMC)
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

		RHRS_table[RHRS_bin-1][2] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z2 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -4 < z < -2 -----*/
	string AcceptanceR_table_z3 = "./information/RHRS_SIMC2023_10_z3.dat";//Acceptance Table (SIMC)
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

		RHRS_table[RHRS_bin-1][3] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z3 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- -2 < z < 0 -----*/
	string AcceptanceR_table_z4 = "./information/RHRS_SIMC2023_10_z4.dat";//Acceptance Table (SIMC)
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

		RHRS_table[RHRS_bin-1][4] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z4 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 0 < z < 2 -----*/
	string AcceptanceR_table_z5 = "./information/RHRS_SIMC2023_10_z5.dat";//Acceptance Table (SIMC)
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

		RHRS_table[RHRS_bin-1][5] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z5 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 2 < z < 4 -----*/
	string AcceptanceR_table_z6 = "./information/RHRS_SIMC2023_10_z6.dat";//Acceptance Table (SIMC)
	string buf_z6;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z6(AcceptanceR_table_z6.c_str(),ios::in);
	if (ifp_z6.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z6.c_str() << endl;
	while(1){
		getline(ifp_z6,buf_z6);
		if(buf_z6[0]=='#'){continue;}
		if(ifp_z6.eof())break;
		stringstream sbuf_z6(buf_z6);
		sbuf_z6 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][6] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z6 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 4 < z < 6 -----*/
	string AcceptanceR_table_z7 = "./information/RHRS_SIMC2023_10_z7.dat";//Acceptance Table (SIMC)
	string buf_z7;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z7(AcceptanceR_table_z7.c_str(),ios::in);
	if (ifp_z7.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z7.c_str() << endl;
	while(1){
		getline(ifp_z7,buf_z7);
		if(buf_z7[0]=='#'){continue;}
		if(ifp_z7.eof())break;
		stringstream sbuf_z7(buf_z7);
		sbuf_z7 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][7] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z7 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 6 < z < 8 -----*/
	string AcceptanceR_table_z8 = "./information/RHRS_SIMC2023_10_z8.dat";//Acceptance Table (SIMC)
	string buf_z8;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z8(AcceptanceR_table_z8.c_str(),ios::in);
	if (ifp_z8.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z8.c_str() << endl;
	while(1){
		getline(ifp_z8,buf_z8);
		if(buf_z8[0]=='#'){continue;}
		if(ifp_z8.eof())break;
		stringstream sbuf_z8(buf_z8);
		sbuf_z8 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][8] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z8 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;
/*----- 8 < z < 10 -----*/
	string AcceptanceR_table_z9 = "./information/RHRS_SIMC2023_10_z9.dat";//Acceptance Table (SIMC)
	string buf_z9;
	RHRS_total=0.;
	RHRS_total_bin=0;

	ifstream ifp_z9(AcceptanceR_table_z9.c_str(),ios::in);
	if (ifp_z9.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << AcceptanceR_table_z9.c_str() << endl;
	while(1){
		getline(ifp_z9,buf_z9);
		if(buf_z9[0]=='#'){continue;}
		if(ifp_z9.eof())break;
		stringstream sbuf_z9(buf_z9);
		sbuf_z9 >> RHRS_bin >> RHRS_SIMC;
		//cout << RHRS_bin << ", " << RHRS_SIMC <<endl;

		RHRS_table[RHRS_bin-1][9] = RHRS_SIMC*0.001;//sr
		RHRS_total+=RHRS_SIMC;
		if(RHRS_SIMC!=0)RHRS_total_bin++;
	}
	cout<<"HRS-R Acceptance (z9 average)="<<RHRS_total/(double)RHRS_total_bin<<endl;

  TH2F* Acceptance_map = new TH2F("Acceptance_map","#Delta#Omega_{K}^{lab}(p_{K},Z)",100,1.6,2.0,10,-10.,10.);
	for(int i=0;i<100;i++){
		for(int j=0;j<10;j++){
			Acceptance_map->SetBinContent(i+1,j+1,RHRS_table[i][j]*1000.);
		}
	}
    
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
//default
 //const double fit_min_mm=-0.01;
 //const double fit_max_mm=0.085;
 const double fit_min_mm=-0.005;

//change
 //const double fit_max_mm=0.085;//Angle 2-div.
 const double fit_max_mm=0.105;//Full (Mom cut)
 
 const double fmin_mm=-0.05;//out of range
 const double fmax_mm=2.20;//out of range
 const double fmin_mm_inside=-0.05;//Full (Mom cut)

//change
 const double fmax_mm_inside=0.15;//Full (Mom cut)
 const double fmax_mm_insideS=0.15;//Full (Mom cut)
 //const double fmax_mm_inside=0.100;//Angle 2-div.(Mom cut)
 //const double fmax_mm_insideS=0.110;//Angle 2-div.(Mom cut)
 
 //const double fmin_mm=-0.01;
 //const double fmax_mm=0.12;
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
 TF1* fmm_strict_Lexp;
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
 TF1* f1_sr;//Survival Ratio



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
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","RESULT (N_{Y}/#Delta#Omega_{K}/#epsilon_{DAQ})",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_Al_wobg_fom_noZ_new  = new TH1F("hmm_Al_wobg_fom_noZ_new","hmm_Al_wobg_fom_new",xbin,xmin,xmax);
  TH1F* hmm_Al_fom_noZ_new  = new TH1F("hmm_Al_fom_noZ_new","hmm_Al_fom_noZ_new",xbin,xmin,xmax);
  TH1F* hmm_L_fom_noZ_new  = new TH1F("hmm_L_fom_noZ_new","hmm_L_fom_noZ_new",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_strict  = new TH1F("hmm_wo_bg_fom_strict","hmm_wo_bg_fom_strict",xbin,xmin,xmax);

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
  TH2D* h_pepk = new TH2D("h_pepk", "p_{K} vs p_{e'}" ,50,1720.,1940.,50,1980.,2220.);
  TH1D* h_pe = new TH1D("h_pe", "p_{e'}" ,200,1980.,2220.);
  TH1D* h_pk = new TH1D("h_pk", "p_{K}" ,200,1720.,1940.);
//added 2023/7/19 for D-thesis Chap.5
  TH2F* h_mm_Qsq = new TH2F("h_mm_Qsq","",100,-100.,200.,50,0.2,0.8);
  SetTH2(h_mm_Qsq, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "Q^{2} [(GeV/c)^{2}]", 0.4);
  TH2F* h_mm_W = new TH2F("h_mm_W","",100,-100.,200.,50,2.05,2.25);
  SetTH2(h_mm_W, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "W [GeV/c]", 0.4);
  TH2F* h_mm_theta_gk_cm = new TH2F("h_mm_theta_gk_cm","",100,-100.,200.,50,-5.,25.);
  SetTH2(h_mm_theta_gk_cm, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "#theta_{#gamma K}^{c.m.} [deg]", 0.4);
  TH2F* h_mm_phi_gk = new TH2F("h_mm_phi_gk","",100,-100.,200.,50,0.,360.);
  SetTH2(h_mm_phi_gk, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "#phi_{#gamma K} [deg]", 0.4);
  TH2F* h_mm_eps = new TH2F("h_mm_eps","",100,-100.,200.,50,0.74,0.8);
  SetTH2(h_mm_eps, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "Transverse polarization #varepsilon", 0.4);
  
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

  TH1F* hmm_L_new_cm2_1  = new TH1F("hmm_L_new_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm2_2  = new TH1F("hmm_L_new_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm3_1  = new TH1F("hmm_L_new_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm3_2  = new TH1F("hmm_L_new_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm3_3  = new TH1F("hmm_L_new_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm4_1  = new TH1F("hmm_L_new_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm4_2  = new TH1F("hmm_L_new_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm4_3  = new TH1F("hmm_L_new_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_cm4_4  = new TH1F("hmm_L_new_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm2_1  = new TH1F("hcs_L_new_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm2_2  = new TH1F("hcs_L_new_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_1  = new TH1F("hcs_L_new_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_2  = new TH1F("hcs_L_new_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm3_3  = new TH1F("hcs_L_new_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_1  = new TH1F("hcs_L_new_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_2  = new TH1F("hcs_L_new_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_3  = new TH1F("hcs_L_new_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hcs_L_new_cm4_4  = new TH1F("hcs_L_new_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hmm_L_new_Qsq2_1  = new TH1F("hmm_L_new_Qsq2_1","Q^{2}<0.5",xbin,xmin,xmax);
  TH1F* hmm_L_new_Qsq2_2  = new TH1F("hmm_L_new_Qsq2_2","Q^{2}>0.5",xbin,xmin,xmax);
  TH1F* hmm_L_new_Qsq3_1  = new TH1F("hmm_L_new_Qsq3_1","Q^{2}<0.45",xbin,xmin,xmax);
  TH1F* hmm_L_new_Qsq3_2  = new TH1F("hmm_L_new_Qsq3_2","0.45<Q^{2}<0.55",xbin,xmin,xmax);
  TH1F* hmm_L_new_Qsq3_3  = new TH1F("hmm_L_new_Qsq3_3","Q^{2}>0.55",xbin,xmin,xmax);
  TH1F* hcs_L_new_Qsq2_1  = new TH1F("hcs_L_new_Qsq2_1","Q^{2}<0.5",xbin,xmin,xmax);
  TH1F* hcs_L_new_Qsq2_2  = new TH1F("hcs_L_new_Qsq2_2","Q^{2}>0.5",xbin,xmin,xmax);
  TH1F* hcs_L_new_Qsq3_1  = new TH1F("hcs_L_new_Qsq3_1","Q^{2}<0.45",xbin,xmin,xmax);
  TH1F* hcs_L_new_Qsq3_2  = new TH1F("hcs_L_new_Qsq3_2","0.45<Q^{2}<0.55",xbin,xmin,xmax);
  TH1F* hcs_L_new_Qsq3_3  = new TH1F("hcs_L_new_Qsq3_3","Q^{2}>0.55",xbin,xmin,xmax);
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
  bool event_selection_strict = false;
  bool event_selection_noZ_new = false;
  bool event_selection_momL1 = false;//VP Flux Syst.
  bool event_selection_momL2 = false;//VP Flux Syst.
  bool event_selection_momL3 = false;//VP Flux Syst.
  bool event_selection_momL4 = false;//VP Flux Syst.
  bool event_selection_momL5 = false;//VP Flux Syst.
  bool event_selection_momL6 = false;//VP Flux Syst.
  bool event_selection_momS1 = false;//VP Flux Syst.
  bool event_selection_momS2 = false;//VP Flux Syst.
  bool event_selection_momS3 = false;//VP Flux Syst.
  bool event_selection_momS4 = false;//VP Flux Syst.
  bool event_selection_momS5 = false;//VP Flux Syst.
  bool event_selection_momS6 = false;//VP Flux Syst.
  bool cm2_angle1_cut=false;
  bool cm2_angle2_cut=false;
  bool cm3_angle1_cut=false;
  bool cm3_angle2_cut=false;
  bool cm3_angle3_cut=false;
  bool cm4_angle1_cut=false;
  bool cm4_angle2_cut=false;
  bool cm4_angle3_cut=false;
  bool cm4_angle4_cut=false;
  bool Qsq2_1_cut=false;
  bool Qsq2_2_cut=false;
  bool Qsq3_1_cut=false;
  bool Qsq3_2_cut=false;
  bool Qsq3_3_cut=false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum

/***********************************/
/**	What parameter do you choose? **/
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
		double effK	 = 0.140;
		double efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		double RHRS  = 0.005;
		//double LHRS  = 0.006;
		double Charge= 4.6486;//[C]
		double ee	 = 1.602*pow(10.,-19);
		double Ne	 = Charge/ee;//Num of e
		//double Ng	 = 2.48*pow(10.,-6)*Ne;//Geant4
		double Ng	 = 1.94*pow(10.,-6)*Ne;//SIMC
		double cs	 = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
cout<<"typical efficiency="<<efficiency<<endl;
cout<<"Ntar(H2)="<<ntar_h2<<endl;
cout<<"Ne="<<Ne<<endl;
cout<<"Ng="<<Ng<<endl;
cout<<"cs="<<cs<<endl;
	double csL[xbin];



  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  tree->Draw(">>elist" , "fabs(ct_orig)<1.006");//ctsum (does NOT dintinguish #track)
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  int ENum = elist->GetN(); 
  int ENum_strict = 0;
  int ENum_strict_cm2_1 = 0;
  int ENum_strict_cm2_2 = 0;
  int ENum_strict_cm3_1 = 0;
  int ENum_strict_cm3_2 = 0;
  int ENum_strict_cm3_3 = 0;
  int ENum_strict_cm4_1 = 0;
  int ENum_strict_cm4_2 = 0;
  int ENum_strict_cm4_3 = 0;
  int ENum_strict_cm4_4 = 0;
  int ENum_strict_Qsq2_1 = 0;
  int ENum_strict_Qsq2_2 = 0;
  int ENum_strict_Qsq3_1 = 0;
  int ENum_strict_Qsq3_2 = 0;
  int ENum_strict_Qsq3_3 = 0;
  double ENum_strict_cs = 0;
  double ENum_strict_cs_cm2_1 = 0;
  double ENum_strict_cs_cm2_2 = 0;
  double ENum_strict_cs_cm3_1 = 0;
  double ENum_strict_cs_cm3_2 = 0;
  double ENum_strict_cs_cm3_3 = 0;
  double ENum_strict_cs_cm4_1 = 0;
  double ENum_strict_cs_cm4_2 = 0;
  double ENum_strict_cs_cm4_3 = 0;
  double ENum_strict_cs_cm4_4 = 0;
  double ENum_strict_cs_Qsq2_1 = 0;
  double ENum_strict_cs_Qsq2_2 = 0;
  double ENum_strict_cs_Qsq3_1 = 0;
  double ENum_strict_cs_Qsq3_2 = 0;
  double ENum_strict_cs_Qsq3_3 = 0;
cout<<"Entries in Cointime gate: "<<ENum<<endl;
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

      
//Cut Condition Definitions
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

		if(fabs(ct)<1.006&&L_mom>2.01&&L_mom<2.16&&R_mom>1.76&&R_mom<1.9)ct_cut=true;
		//if(fabs(ct)<0.2)ct_cut=true;
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
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_noZ_new=true;
		else event_selection_noZ_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_strict=true;
		else event_selection_strict=false;

//VP Flux Systematic Error
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.89)event_selection_momL1=true;
		else event_selection_momL1=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.88)event_selection_momL2=true;
		else event_selection_momL2=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.87)event_selection_momL3=true;
		else event_selection_momL3=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.86)event_selection_momL4=true;
		else event_selection_momL4=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.85)event_selection_momL5=true;
		else event_selection_momL5=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p<1.84)event_selection_momL6=true;
		else event_selection_momL6=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.77)event_selection_momS1=true;
		else event_selection_momS1=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.78)event_selection_momS2=true;
		else event_selection_momS2=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.79)event_selection_momS3=true;
		else event_selection_momS3=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.80)event_selection_momS4=true;
		else event_selection_momS4=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.81)event_selection_momS5=true;
		else event_selection_momS5=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_tr_p>1.82)event_selection_momS6=true;
		else event_selection_momS6=false;
//----------------------------Cut Condition Definition


//-----Kinematics Calculation
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
		//if(event_selection_new&&ct_cut)hmm_L_fom_strict->Fill(mm);
		if(event_selection_strict&&ct_cut){
			hmm_L_fom_strict->Fill(mm);
			ENum_strict++;
			//if(mm<-0.03||mm>0.12)cout<<"mm="<<mm<<", L_mom="<<L_mom<<", R_mom="<<R_mom<<", L_tr_tg_th="<<L_tr_tg_th<<", R_tr_tg_th="<<R_tr_tg_th<<", L_tr_tg_ph="<<L_tr_tg_ph<<", R_tr_tg_ph="<<R_tr_tg_ph<<endl;
		}
		//if(event_selection_momS6&&ct_cut)hmm_L_fom_strict->Fill(mm);


		if(event_selection_noZ_new&&ct_cut)hmm_L_fom_noZ_new->Fill(mm);

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
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		double q2=Qsq+omega*omega;
		double eps=1/(1+2*(q2/Qsq)*tan(theta_ee/2)*tan(theta_ee/2));
//Phi dependence
		TVector3 G_3vec = G_4vec.Vect();
		TVector3 L_3vec = L_4vec.Vect();
		TVector3 R_3vec = R_4vec.Vect();
		TVector3 l_3vec = G_3vec.Cross(L_3vec);
		TVector3 r_3vec = G_3vec.Cross(R_3vec);
		TVector3 s_3vec = l_3vec.Cross(r_3vec);//for sgn(sign(phi_k))
		TVector3 axis_3vec;
		axis_3vec.SetXYZ(0.,-1.,0.);//along VP flux
		double sgn = s_3vec*axis_3vec;
		double phi_k_cos = (l_3vec*r_3vec)/l_3vec.Mag()/r_3vec.Mag();
		double phi_k;
		if(sgn>0.){phi_k = acos(phi_k_cos);}
		else{phi_k = 2*PI-acos(phi_k_cos);}
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
		Qsq2_1_cut=false;
    	Qsq2_2_cut=false;
    	Qsq3_1_cut=false;
    	Qsq3_2_cut=false;
    	Qsq3_3_cut=false;
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
		if(Qsq<0.5)Qsq2_1_cut=true;
		if(Qsq>=0.5)Qsq2_2_cut=true;
		if(Qsq<0.45)Qsq3_1_cut=true;
		if(Qsq>=0.45&&Qsq<0.55)Qsq3_2_cut=true;
		if(Qsq>=0.55)Qsq3_3_cut=true;
//======= CM Angle(gamma-K) ========//

		//int ebin = (int)((L_mom-1.8)/0.004);
		//if(ebin>=0 &&ebin<150){
		//Ng = Ng_table[ebin]*Ne;//
		//if(Ng!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
		//else cs=0.;
		//}else{cs=0.;} 
		//int kbin = (int)((R_mom-1.5)/0.004);
		//int zbin = (int)((((L_tr_vz+R_tr_vz)/2.)+0.1)/0.02);
		int kbin = (int)((R_mom-1.6)*100./0.4);
		if(event_selection&&ct_cut&&kbin>=0 &&kbin<100){
		if((L_tr_vz+R_tr_vz)/2.>=-0.10&&(L_tr_vz+R_tr_vz)/2.<-0.08)RHRS = RHRS_table[kbin][0];
		else if((L_tr_vz+R_tr_vz)/2.>=-0.08&&(L_tr_vz+R_tr_vz)/2.<-0.06)RHRS = RHRS_table[kbin][1];
		else if((L_tr_vz+R_tr_vz)/2.>=-0.06&&(L_tr_vz+R_tr_vz)/2.<-0.04)RHRS = RHRS_table[kbin][2];
		else if((L_tr_vz+R_tr_vz)/2.>=-0.04&&(L_tr_vz+R_tr_vz)/2.<-0.02)RHRS = RHRS_table[kbin][3];
		else if((L_tr_vz+R_tr_vz)/2.>=-0.02&&(L_tr_vz+R_tr_vz)/2.<0.)RHRS = RHRS_table[kbin][4];
		else if((L_tr_vz+R_tr_vz)/2.>=0.&&(L_tr_vz+R_tr_vz)/2.<0.02)RHRS = RHRS_table[kbin][5];
		else if((L_tr_vz+R_tr_vz)/2.>=0.02&&(L_tr_vz+R_tr_vz)/2.<0.04)RHRS = RHRS_table[kbin][6];
		else if((L_tr_vz+R_tr_vz)/2.>=0.04&&(L_tr_vz+R_tr_vz)/2.<0.06)RHRS = RHRS_table[kbin][7];
		else if((L_tr_vz+R_tr_vz)/2.>=0.06&&(L_tr_vz+R_tr_vz)/2.<0.08)RHRS = RHRS_table[kbin][8];
		else if((L_tr_vz+R_tr_vz)/2.>=0.08&&(L_tr_vz+R_tr_vz)/2.<0.10)RHRS = RHRS_table[kbin][9];
		else cout<<"Z Error"<<(L_tr_vz+R_tr_vz)/2.<<endl;

		effDAQ = daq_table[nrun-111000];
		if(effDAQ==0.2)cout<<"Starange!!! DAQ Eff. of run"<<nrun<<" does not exist."<<endl;
		efficiency = effAC*effZ*effFP*effch2*effct*effDAQ*efftr*effK;
		//if(RHRS!=0.)cs = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]

		//PathLength calculation from Cointime(Jan. 30, 2022)
		f1_sr = new TF1("f1_sr", GetSRfromCt, 1730.,1930.,2);
		f1_sr->SetParameter(0,31.2307);
		f1_sr->SetParameter(1,-0.0110408);
		effK = f1_sr->Eval(R_mom*1000.);
		//cout<<"effK="<<effK<<endl;
		//cout<<"mom="<<R_mom<<endl;
		//-------------

		//if(RHRS!=0.&&effDAQ!=0.)cs = labtocm/effDAQ/RHRS;//[nb/sr]
		if(RHRS!=0.&&effDAQ!=0.)cs = labtocm/effDAQ/RHRS/effK;//[nb/sr]
		else cs=0.;
		}else{cs=0.;}
		double cs_temp = cs*labtocm;
		//cs_ave += cs;
		//cout<<"Ng="<<Ng<<endl;
		//cout<<"cs="<<cs<<endl;
		if(event_selection&&ct_cut)h_pepk->Fill(R_mom*1000.,L_mom*1000.);
		if(event_selection&&ct_cut)h_pe->Fill(L_mom*1000.);
		if(event_selection&&ct_cut)h_pk->Fill(R_mom*1000.);
		if(event_selection&&ct_cut)hcs_L_fom_best->Fill(mm,cs);


		if(event_selection_new&&ct_cut){
			hcs_L_fom_strict->Fill(mm,cs);
			if(cs>0.)ENum_strict_cs+=cs/150.;
			h_mm_Qsq->Fill(mm*1000.,Qsq);
			h_mm_W->Fill(mm*1000.,W);
			h_mm_theta_gk_cm->Fill(mm*1000.,theta_gk_cm*180./PI);
			h_mm_phi_gk->Fill(mm*1000.,phi_k*180./PI);
			h_mm_eps->Fill(mm*1000.,eps);
		}
		//if(event_selection_momS6&&ct_cut)hcs_L_fom_strict->Fill(mm,cs);
		//
		//
		//if(event_selection&&ct_cut)hcs_L_fom_best->Fill(mm,cs);
//		//if(event_selection_nocut&&ct_cut)hcs_L_fom_nocut->SetBinContent(hmm_L_fom_best->FindBin(mm),cs);

		if(event_selection&&ct_cut){
			gklab_gkcm->Fill(theta_gk_lab,theta_gk_cm);
			gklab_eklab->Fill(theta_gk_lab,theta_ek);
			eelab_eklab->Fill(theta_ee,theta_ek);
			eklab_gkcm->Fill(theta_ek,theta_gk_cm);
			if(cm2_angle1_cut){hcs_L_cm2_1->Fill(mm,cs);}
			if(cm2_angle2_cut){hcs_L_cm2_2->Fill(mm,cs);}
			if(cm3_angle1_cut){hcs_L_cm3_1->Fill(mm,cs);}
			if(cm3_angle2_cut){hcs_L_cm3_2->Fill(mm,cs);}
			if(cm3_angle3_cut){hcs_L_cm3_3->Fill(mm,cs);}
			if(cm4_angle1_cut){hcs_L_cm4_1->Fill(mm,cs);}
			if(cm4_angle2_cut){hcs_L_cm4_2->Fill(mm,cs);}
			if(cm4_angle3_cut){hcs_L_cm4_3->Fill(mm,cs);}
			if(cm4_angle4_cut){hcs_L_cm4_4->Fill(mm,cs);}
		}
		if(event_selection_new&&ct_cut){
			if(cm2_angle1_cut){hmm_L_new_cm2_1->Fill(mm);ENum_strict_cm2_1++;}
			if(cm2_angle2_cut){hmm_L_new_cm2_2->Fill(mm);ENum_strict_cm2_2++;}
			if(cm3_angle1_cut){hmm_L_new_cm3_1->Fill(mm);ENum_strict_cm3_1++;}
			if(cm3_angle2_cut){hmm_L_new_cm3_2->Fill(mm);ENum_strict_cm3_2++;}
			if(cm3_angle3_cut){hmm_L_new_cm3_3->Fill(mm);ENum_strict_cm3_3++;}
			if(cm4_angle1_cut){hmm_L_new_cm4_1->Fill(mm);ENum_strict_cm4_1++;}
			if(cm4_angle2_cut){hmm_L_new_cm4_2->Fill(mm);ENum_strict_cm4_2++;}
			if(cm4_angle3_cut){hmm_L_new_cm4_3->Fill(mm);ENum_strict_cm4_3++;}
			if(cm4_angle4_cut){hmm_L_new_cm4_4->Fill(mm);ENum_strict_cm4_4++;}
			if(cm2_angle1_cut){hcs_L_new_cm2_1->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm2_1+=cs*2./150.;}
			if(cm2_angle2_cut){hcs_L_new_cm2_2->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm2_2+=cs*2./150.;}
			if(cm3_angle1_cut){hcs_L_new_cm3_1->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm3_1+=cs*3./150.;}
			if(cm3_angle2_cut){hcs_L_new_cm3_2->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm3_2+=cs*3./150.;}
			if(cm3_angle3_cut){hcs_L_new_cm3_3->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm3_3+=cs*3./150.;}
			if(cm4_angle1_cut){hcs_L_new_cm4_1->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm4_1+=cs*4./150.;}
			if(cm4_angle2_cut){hcs_L_new_cm4_2->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm4_2+=cs*4./150.;}
			if(cm4_angle3_cut){hcs_L_new_cm4_3->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm4_3+=cs*4./150.;}
			if(cm4_angle4_cut){hcs_L_new_cm4_4->Fill(mm,cs);if(cs>0.)ENum_strict_cs_cm4_4+=cs*4./150.;}
			if(Qsq2_1_cut){hmm_L_new_Qsq2_1->Fill(mm);ENum_strict_Qsq2_1++;}
			if(Qsq2_2_cut){hmm_L_new_Qsq2_2->Fill(mm);ENum_strict_Qsq2_2++;}
			if(Qsq3_1_cut){hmm_L_new_Qsq3_1->Fill(mm);ENum_strict_Qsq3_1++;}
			if(Qsq3_2_cut){hmm_L_new_Qsq3_2->Fill(mm);ENum_strict_Qsq3_2++;}
			if(Qsq3_3_cut){hmm_L_new_Qsq3_3->Fill(mm);ENum_strict_Qsq3_3++;}
			if(Qsq2_1_cut){hcs_L_new_Qsq2_1->Fill(mm,cs);if(cs>0.)ENum_strict_cs_Qsq2_1+=cs*2./150.;}
			if(Qsq2_2_cut){hcs_L_new_Qsq2_2->Fill(mm,cs);if(cs>0.)ENum_strict_cs_Qsq2_2+=cs*2./150.;}
			if(Qsq3_1_cut){hcs_L_new_Qsq3_1->Fill(mm,cs);if(cs>0.)ENum_strict_cs_Qsq3_1+=cs*3./150.;}
			if(Qsq3_2_cut){hcs_L_new_Qsq3_2->Fill(mm,cs);if(cs>0.)ENum_strict_cs_Qsq3_2+=cs*3./150.;}
			if(Qsq3_3_cut){hcs_L_new_Qsq3_3->Fill(mm,cs);if(cs>0.)ENum_strict_cs_Qsq3_3+=cs*3./150.;}
		}

		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01)&&ct_cut)hmm_Al_fom_noZ_new->Fill(mm);//Al selection
		//if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01)&&ct_cut)hmm_Al_fom_noZ_new->Fill(mm);//Al selection

}//ENum



//Analysis Part: fitting and fitting...
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
	//TH1F* hmm_bg_fom_strict=(TH1F*)file_mea->Get("hmm_mixacc_result_momS6");
	//TH1F* hcs_bg_fom_strict=(TH1F*)file_mea->Get("hcs_mixacc_result_momS6");
	
	TH1F* hcs_bg_cm2_1=(TH1F*)file_mea->Get("hcs_mixacc_result_cm2_1");
	TH1F* hcs_bg_cm2_2=(TH1F*)file_mea->Get("hcs_mixacc_result_cm2_2");
	TH1F* hcs_bg_cm3_1=(TH1F*)file_mea->Get("hcs_mixacc_result_cm3_1");
	TH1F* hcs_bg_cm3_2=(TH1F*)file_mea->Get("hcs_mixacc_result_cm3_2");
	TH1F* hcs_bg_cm3_3=(TH1F*)file_mea->Get("hcs_mixacc_result_cm3_3");
	TH1F* hcs_bg_cm4_1=(TH1F*)file_mea->Get("hcs_mixacc_result_cm4_1");
	TH1F* hcs_bg_cm4_2=(TH1F*)file_mea->Get("hcs_mixacc_result_cm4_2");
	TH1F* hcs_bg_cm4_3=(TH1F*)file_mea->Get("hcs_mixacc_result_cm4_3");
	TH1F* hcs_bg_cm4_4=(TH1F*)file_mea->Get("hcs_mixacc_result_cm4_4");
	TH1F* hcs_bg_new_cm2_1=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm2_1");
	TH1F* hcs_bg_new_cm2_2=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm2_2");
	TH1F* hcs_bg_new_cm3_1=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm3_1");
	TH1F* hcs_bg_new_cm3_2=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm3_2");
	TH1F* hcs_bg_new_cm3_3=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm3_3");
	TH1F* hcs_bg_new_cm4_1=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm4_1");
	TH1F* hcs_bg_new_cm4_2=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm4_2");
	TH1F* hcs_bg_new_cm4_3=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm4_3");
	TH1F* hcs_bg_new_cm4_4=(TH1F*)file_mea->Get("hcs_mixacc_result_new_cm4_4");
	TH1F* hcs_bg_new_Qsq2_1=(TH1F*)file_mea->Get("hcs_mixacc_result_new_Qsq2_1");
	TH1F* hcs_bg_new_Qsq2_2=(TH1F*)file_mea->Get("hcs_mixacc_result_new_Qsq2_2");
	TH1F* hcs_bg_new_Qsq3_1=(TH1F*)file_mea->Get("hcs_mixacc_result_new_Qsq3_1");
	TH1F* hcs_bg_new_Qsq3_2=(TH1F*)file_mea->Get("hcs_mixacc_result_new_Qsq3_2");
	TH1F* hcs_bg_new_Qsq3_3=(TH1F*)file_mea->Get("hcs_mixacc_result_new_Qsq3_3");
	TH1F* hmm_bg_new_Qsq2_1=(TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq2_1");
	TH1F* hmm_bg_new_Qsq2_2=(TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq2_2");
	TH1F* hmm_bg_new_Qsq3_1=(TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq3_1");
	TH1F* hmm_bg_new_Qsq3_2=(TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq3_2");
	TH1F* hmm_bg_new_Qsq3_3=(TH1F*)file_mea->Get("hmm_mixacc_result_new_Qsq3_3");
	TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut");
	TH1F* hmm_Albg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_forAl");
	//TH1F* hmm_bg_fom_nocut=(TH1F*)file_mea->Get("hmm_mixacc_nocut_result");
	TH1F* hmm_Albg_fom_noZ_new=(TH1F*)file_mea->Get("hmm_mixacc_result_nocut_new_forAl");//strict AC
    int fitmin = hmm_L_fom_best->FindBin(0.10);
    int fitmax = hmm_L_fom_best->FindBin(0.15);
    double num1 = hmm_L_fom_best->Integral(fitmin,fitmax);
    double num2 = hmm_bg_fom_best->Integral(fitmin,fitmax);
    double mixscale = num1/num2;
	cout<<"hmm_L integral ="<<num1<<endl;
	cout<<"hmm_bg integral ="<<num2<<endl;
	cout<<"mixscale(mixed/original)="<<1/mixscale<<endl;
	//hmm_bg_fom_best->Sumw2();
	//hmm_bg_fom_strict->Sumw2();
	//hmm_bg_fom_nocut->Sumw2();
	//hcs_bg_fom_strict->Sumw2();
	//hcs_bg_new_cm2_1->Sumw2();
	//hcs_bg_new_cm2_2->Sumw2();
	//hcs_bg_new_cm3_1->Sumw2();
	//hcs_bg_new_cm3_2->Sumw2();
	//hcs_bg_new_cm3_3->Sumw2();
	//hcs_bg_new_cm4_1->Sumw2();
	//hcs_bg_new_cm4_2->Sumw2();
	//hcs_bg_new_cm4_3->Sumw2();
	//hcs_bg_new_cm4_4->Sumw2();
	//hcs_bg_new_Qsq2_1->Sumw2();
	//hcs_bg_new_Qsq2_2->Sumw2();
	//hcs_bg_new_Qsq3_1->Sumw2();
	//hcs_bg_new_Qsq3_2->Sumw2();
	//hcs_bg_new_Qsq3_3->Sumw2();
	hmm_Albg_fom_nocut->Sumw2();
	hmm_Albg_fom_noZ_new->Scale(1./nbunch);
	hmm_bg_fom_best->Scale(1./nbunch);
	hmm_bg_fom_strict->Scale(1./nbunch);
	hcs_bg_fom_strict->Scale(1./nbunch);
	hcs_bg_new_cm2_1->Scale(1./nbunch);
	hcs_bg_new_cm2_2->Scale(1./nbunch);
	hcs_bg_new_cm3_1->Scale(1./nbunch);
	hcs_bg_new_cm3_2->Scale(1./nbunch);
	hcs_bg_new_cm3_3->Scale(1./nbunch);
	hcs_bg_new_cm4_1->Scale(1./nbunch);
	hcs_bg_new_cm4_2->Scale(1./nbunch);
	hcs_bg_new_cm4_3->Scale(1./nbunch);
	hcs_bg_new_cm4_4->Scale(1./nbunch);
	hcs_bg_new_Qsq2_1->Scale(1./nbunch);
	hcs_bg_new_Qsq2_2->Scale(1./nbunch);
	hcs_bg_new_Qsq3_1->Scale(1./nbunch);
	hcs_bg_new_Qsq3_2->Scale(1./nbunch);
	hcs_bg_new_Qsq3_3->Scale(1./nbunch);
	//cs	 = pow(10.,33.)/(ntar_h2*efficiency*RHRS*Ng);//[nb/sr]
	effDAQ= 0.950;
	RHRS  = 0.0055;
	effK= 0.14;
	double labtocm = 0.126;
//comment out Feb. 1, why was this needed?
//	//cs	 = labtocm/(effDAQ*RHRS*10.*effK);
//	cs	 = labtocm/(effDAQ*RHRS*effK);
//	hmm_bg_fom_best->Scale(cs);
////	hmm_bg_fom_strict->Scale(cs);
//	hmm_bg_cm2_1->Scale(1./nbunch/cs);
//	hmm_bg_cm2_2->Scale(1./nbunch/cs);
//	hmm_bg_cm3_1->Scale(1./nbunch/cs);
//	hmm_bg_cm3_2->Scale(1./nbunch/cs);
//	hmm_bg_cm3_3->Scale(1./nbunch/cs);
//	hmm_bg_cm4_1->Scale(1./nbunch/cs);
//	hmm_bg_cm4_2->Scale(1./nbunch/cs);
//	hmm_bg_cm4_3->Scale(1./nbunch/cs);
//	hmm_bg_cm4_4->Scale(1./nbunch/cs);
	hmm_bg_new_cm2_1->Scale(1./nbunch);
	hmm_bg_new_cm2_2->Scale(1./nbunch);
	hmm_bg_new_cm3_1->Scale(1./nbunch);
	hmm_bg_new_cm3_2->Scale(1./nbunch);
	hmm_bg_new_cm3_3->Scale(1./nbunch);
	hmm_bg_new_cm4_1->Scale(1./nbunch);
	hmm_bg_new_cm4_2->Scale(1./nbunch);
	hmm_bg_new_cm4_3->Scale(1./nbunch);
	hmm_bg_new_cm4_4->Scale(1./nbunch);
	hmm_bg_new_Qsq2_1->Scale(1./nbunch);
	hmm_bg_new_Qsq2_2->Scale(1./nbunch);
	hmm_bg_new_Qsq3_1->Scale(1./nbunch);
	hmm_bg_new_Qsq3_2->Scale(1./nbunch);
	hmm_bg_new_Qsq3_3->Scale(1./nbunch);
	hmm_bg_fom_nocut->Scale(1./nbunch);
	hmm_Albg_fom_nocut->Scale(1./nbunch);
	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
	//hmm_wo_bg_fom_best->Add(hcs_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
	


//MM spctrum to be fitted (2020/10/18)
//Choose one from the list below
	//hmm_wo_bg_fom_strict->Add(hcs_L_fom_strict,hcs_bg_fom_strict,1.0,-1.0);//All by hcs
	//hmm_wo_bg_fom_strict->Scale(1./150.);
	//hmm_wo_bg_fom_strict->Add(hmm_L_fom_strict,hmm_bg_fom_strict,1.0,-1.0);//All by hmm
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm2_1,hcs_bg_new_cm2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm2_2,hcs_bg_new_cm2_2,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Scale(2./150.);
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm3_1,hcs_bg_new_cm3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm3_2,hcs_bg_new_cm3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm3_3,hcs_bg_new_cm3_3,1.0,-1.0);//3 div.
//===CHANGE===//
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_Qsq2_1,hcs_bg_new_Qsq2_1,1.0,-1.0);//2 div.
	hmm_wo_bg_fom_strict->Add(hcs_L_new_Qsq2_2,hcs_bg_new_Qsq2_2,1.0,-1.0);//2 div.
	hmm_wo_bg_fom_strict->Scale(2./150.);
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq2_1,hmm_bg_new_Qsq2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq2_2,hmm_bg_new_Qsq2_2,1.0,-1.0);//2 div.
//============//
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_Qsq3_1,hcs_bg_new_Qsq3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_Qsq3_2,hcs_bg_new_Qsq3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_Qsq3_3,hcs_bg_new_Qsq3_3,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm4_1,hcs_bg_new_cm4_1,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm4_2,hcs_bg_new_cm4_2,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm4_3,hcs_bg_new_cm4_3,1.0,-1.0);//4 div.
	//hmm_wo_bg_fom_strict->Add(hcs_L_new_cm4_4,hcs_bg_new_cm4_4,1.0,-1.0);//4 div.
//MM
	//hmm_wo_bg_fom_strict->Add(hmm_L_fom_strict,hmm_bg_fom_strict,1.0,-1.0);//All
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_cm2_1,hmm_bg_new_cm2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_cm2_2,hmm_bg_new_cm2_2,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_cm3_1,hmm_bg_new_cm3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_cm3_2,hmm_bg_new_cm3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_cm3_3,hmm_bg_new_cm3_3,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq2_1,hmm_bg_new_Qsq2_1,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq2_2,hmm_bg_new_Qsq2_2,1.0,-1.0);//2 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq3_1,hmm_bg_new_Qsq3_1,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq3_2,hmm_bg_new_Qsq3_2,1.0,-1.0);//3 div.
	//hmm_wo_bg_fom_strict->Add(hmm_L_new_Qsq3_3,hmm_bg_new_Qsq3_3,1.0,-1.0);//3 div.

	
	//for(int i=90;i<200;i++){
	//if(hcs_L_new_cm2_1->GetBinContent(i)==0){
	//cout<<"Empty bin at "<<((double)i*0.001-0.1)<<endl;
	//hmm_wo_bg_fom_best->SetBinContent(i,0);
	//hmm_wo_bg_fom_best->SetBinError(i,0);
	//}
	//}

	hmm_wo_bg_fom_nocut->Add(hmm_L_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);

//For pion B.G. function fitting
	hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
	//hmm_pi_wobg_fom_nocut->Add(hmm_pi_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
	
//For Al B.G. function fitting
	hmm_Al_wobg_fom_noZ_new->Add(hmm_Al_fom_noZ_new,hmm_Albg_fom_noZ_new,1.0,-1.0);
	//hmm_Al_wobg_fom_noZ_new->Add(hmm_Al_fom_noZ_new,hmm_Albg_fom_nocut,1.0,-1.0);

	//hmm_pi_wobg_fom_nocut->Add(hmm_Al_fom_nocut,hmm_Albg_fom_nocut,1.0,-1.0);
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

//Al B.G. fitting
	 TF1 *fAl_new=new TF1("fAl_new",F_Voigt,fit_min_mm,fit_max_mm,4);
	 fAl_new->SetNpx(2000);
	 fAl_new->SetTitle("Al selected (new)");
	 fAl_new->SetParameters(1.,0.05,0.04,0.0001);
	 fAl_new->SetParLimits(0,0.,10000.);
	 fAl_new->SetLineColor(kRed);
	 hmm_Al_wobg_fom_noZ_new->Fit("fAl_new","I","",0.02,0.10);
	 double Al_par0 = fAl_new->GetParameter(0);
	 double Al_par1 = fAl_new->GetParameter(1);
	 double Al_par2 = fAl_new->GetParameter(2);
	 double Al_par3 = fAl_new->GetParameter(3);

//pion B.G. fitting
	 TF1 *fpion=new TF1("fpion",F_Voigt,fit_min_mm,fit_max_mm,4);
	 fpion->SetNpx(2000);
	 fpion->SetTitle("Pion selected");
	 fpion->SetParameters(3.,0.05,0.04,0.001);
	 fpion->SetParLimits(0,0.,10000.);
	 fpion->SetLineColor(kRed);
	 hmm_pi_wobg_fom_best->Fit("fpion","I","",0.,0.1);
	 double pion_par0 = fpion->GetParameter(0);
	 double pion_par1 = fpion->GetParameter(1);
	 double pion_par2 = fpion->GetParameter(2);
	 double pion_par3 = fpion->GetParameter(3);

/*%%%%%%%%%%%%%%%%*/
/*%%    Strict	%%*/
/*%%%%%%%%%%%%%%%%*/
	 cout<<"strict cut START"<<endl;
	 cout<<"(Landau+Exp)*(Gauss) MODE START"<<endl;
	 fmm_strict_Lexp=new TF1("fmm_strict_Lexp",FMM_Res_nocut,fmin_mm,fmax_mm,22);
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

//change
int fit_flag = 3;
	//1: Fixed from SIMC (most reliable)
	//2: Free			 (chi-square is the best)
	//3: Fixed from data (out of acceptance)

//Free Fit
#if 0
//mean_L_best=-0.001;
//mean_S_best=0.076;
//	 fmm_strict_Lexp->SetParLimits(2,0.,1000.);//positive
//	 fmm_strict_Lexp->SetParLimits(9,0.,300.);//positive
//	 fmm_strict_Lexp->SetParameter(0,0.0007);//Landau width
//	 fmm_strict_Lexp->SetParameter(1,mean_L_best);
//	 fmm_strict_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
//	 fmm_strict_Lexp->SetParameter(2,1.5);//total scale
//	 fmm_strict_Lexp->SetParameter(3,0.001);//sigma
//	 fmm_strict_Lexp->SetParLimits(3,0.,0.01);
//	 fmm_strict_Lexp->SetParameter(4,0.05);//att.
//	 fmm_strict_Lexp->SetParLimits(4,0.005,0.08);
//	 fmm_strict_Lexp->SetParameter(5,-0.004);//peak pos.
//	 fmm_strict_Lexp->SetParLimits(5,-0.05,0.05);
//	 fmm_strict_Lexp->SetParameter(6,0.6);//relative strength
//	 fmm_strict_Lexp->SetParLimits(6,0.,1.5);//relative strength
//
//	 fmm_strict_Lexp->SetParameter(7,0.0003);//Landau width
//	 fmm_strict_Lexp->SetParameter(8,mean_S_best);//MPV
//	 fmm_strict_Lexp->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fmm_strict_Lexp->SetParameter(9,0.4);//total scale
//	 fmm_strict_Lexp->SetParameter(10,0.0015);//sigma
//	 fmm_strict_Lexp->SetParLimits(10,0.,0.01);
//	 fmm_strict_Lexp->SetParameter(11,0.05);//att
//	 fmm_strict_Lexp->SetParLimits(11,0.03,0.12);
//	 fmm_strict_Lexp->SetParameter(12,0.080);//peak pos.
//	 //fmm_strict_Lexp->SetParLimits(15,-0.085,-0.055);
//	 fmm_strict_Lexp->SetParameter(13,0.6);
//	 fmm_strict_Lexp->SetParLimits(13,0.,1.5);//relative strength
//	 //fmm_strict_Lexp->SetParameter(14,500.);//scale
//	 //fmm_strict_Lexp->SetParLimits(14,0.,1000000.);//scale
//	 //fmm_strict_Lexp->SetParameter(15,0.05);//mean
//	 //fmm_strict_Lexp->SetParameter(16,0.04);//Gsigma
//	 //fmm_strict_Lexp->SetParameter(17,0.01);//Lfwhm
	 fmm_strict_Lexp->SetParameter(0,0.000658);
	 fmm_strict_Lexp->SetParameter(1,-0.001258);
	 fmm_strict_Lexp->SetParameter(2,1.42394);
	 fmm_strict_Lexp->SetParameter(3,0.00108);
	 fmm_strict_Lexp->SetParameter(4,0.06493);
	 fmm_strict_Lexp->SetParameter(5,0.003);
	 fmm_strict_Lexp->SetParameter(6,0.692);
	 fmm_strict_Lexp->SetParameter(7,0.0001579);
	 fmm_strict_Lexp->SetParameter(8,0.0761140);
	 fmm_strict_Lexp->SetParameter(9,0.551341);
	 fmm_strict_Lexp->SetParameter(10,0.00168);
	 fmm_strict_Lexp->SetParameter(11,0.1107);
	 fmm_strict_Lexp->SetParameter(12,-0.08095);
	 fmm_strict_Lexp->SetParameter(13,1.499);
	 fmm_strict_Lexp->FixParameter(14,(double)ENum_strict*0.021*0.001);//scale(1.8%(pion)+0.3%(Al)) //B.G. ratio
	 fmm_strict_Lexp->FixParameter(15,Al_par1);//mean
	 fmm_strict_Lexp->FixParameter(16,Al_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(17,Al_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(18,pion_par1);//mean
	 fmm_strict_Lexp->FixParameter(19,pion_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(20,pion_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(21,0.3/1.8);//Al vs Pi
//	 fmm_strict_Lexp->SetParameter(14,6.);//scale
//e	 fmm_strict_Lexp->SetParLimits(14,0.,1000000.);//scale
#endif


//Fix Fit
#if 1
//Lambda//
	 fmm_strict_Lexp->SetParLimits(2,0.,1000.);//positive
	 fmm_strict_Lexp->SetParLimits(9,0.,300.);//positive
	 fmm_strict_Lexp->SetParLimits(5,-0.005,0.005);//peak pos.
	 fmm_strict_Lexp->SetParLimits(12,-0.085,-0.065);//peak pos.
	 //fmm_strict_Lexp->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_strict_Lexp->SetParLimits(3,0.,0.01);
	 //fmm_strict_Lexp->SetParLimits(4,0.005,0.08);
	 //fmm_strict_Lexp->SetParLimits(5,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_strict_Lexp->SetParLimits(6,0.,1.5);//relative strength

//default
	 //fmm_strict_Lexp->SetParameter(0,0.0007);//Landau width
	 //fmm_strict_Lexp->SetParameter(1,mean_L_strict);
	 //fmm_strict_Lexp->SetParameter(2,48.);//total scale
	 //fmm_strict_Lexp->SetParameter(3,0.001);//sigma
	 //fmm_strict_Lexp->SetParameter(4,0.05);//att.
	 //fmm_strict_Lexp->SetParameter(5,-1.*mean_L_strict);//peak pos.
	 //fmm_strict_Lexp->SetParameter(6,0.6);//relative strength
//default
	 //fmm_strict_Lexp->FixParameter(0,0.000658270);//Landau width
	 ////fmm_strict_Lexp->SetParLimits(0,0.00061,0.00062);//Landau width
	 //fmm_strict_Lexp->FixParameter(1,-0.00125805);
	 ////fmm_strict_Lexp->SetParLimits(1,-0.0014,-0.0013);
	 //fmm_strict_Lexp->SetParameter(2,49./2);//total scale
	 ////fmm_strict_Lexp->SetParLimits(2,20.,30.);//total scale
	 //fmm_strict_Lexp->FixParameter(3,0.00108656);//sigma
	 ////fmm_strict_Lexp->SetParLimits(3,0.0011,0.0012);//sigma
	 //fmm_strict_Lexp->FixParameter(4,0.0649346);//att.
	 ////fmm_strict_Lexp->SetParLimits(4,0.05775,0.05777);//att.
	 //fmm_strict_Lexp->FixParameter(5,0.003);//peak pos.
	 ////fmm_strict_Lexp->SetParLimits(5,-0.00187,-0.00186);//peak pos.
	 //fmm_strict_Lexp->FixParameter(6,0.692215);//relative strength
	 ////fmm_strict_Lexp->SetParLimits(6,0.577,0.578);//relative strength
switch(fit_flag){
case 1: //Fixed from SIMC
//simc_fitting2023
	 fmm_strict_Lexp->FixParameter(0,0.000337527);//Landau width
	 fmm_strict_Lexp->FixParameter(1,-0.000399571);
	 fmm_strict_Lexp->SetParameter(2,1.24759);//total scale
	 //fmm_strict_Lexp->SetParLimits(2,20.,30.);//total scale
	 fmm_strict_Lexp->FixParameter(3,0.00129240);//sigma
	 fmm_strict_Lexp->FixParameter(4,0.0499249);//att.
	 fmm_strict_Lexp->FixParameter(5,0.000548956);//peak pos.
	 fmm_strict_Lexp->FixParameter(6,0.716364);//relative strength

	 fmm_strict_Lexp->FixParameter(7,0.000263497);//Landau width
	 fmm_strict_Lexp->FixParameter(8,0.0755465);
	 fmm_strict_Lexp->SetParameter(9,0.304532);//total scale
	 //fmm_strict_Lexp->SetParLimits(9,15.,20.);//total scale
	 fmm_strict_Lexp->FixParameter(10,0.00138720);//sigma
	 fmm_strict_Lexp->FixParameter(11,0.0145474);//att.
	 fmm_strict_Lexp->FixParameter(12,-0.0712964);//peak pos.
	 fmm_strict_Lexp->FixParameter(13,0.289446);//relative strength
	 //fmm_strict_Lexp->FixParameter(13,0.8);//relative strength
//simc_fitting2022kai
	 //fmm_strict_Lexp->FixParameter(0,0.000392691);//Landau width
	 //fmm_strict_Lexp->FixParameter(1,-0.000294567);
	 //fmm_strict_Lexp->SetParameter(2,1.30797);//total scale
	 ////fmm_strict_Lexp->SetParLimits(2,20.,30.);//total scale
	 //fmm_strict_Lexp->FixParameter(3,0.00146711);//sigma
	 //fmm_strict_Lexp->FixParameter(4,0.0566266);//att.
	 //fmm_strict_Lexp->FixParameter(5,-0.000345746);//peak pos.
	 //fmm_strict_Lexp->FixParameter(6,0.584759);//relative strength
	 //fmm_strict_Lexp->FixParameter(7,0.000343969);//Landau width
	 //fmm_strict_Lexp->FixParameter(8,0.0756934);
	 //fmm_strict_Lexp->SetParameter(9,0.304364);//total scale
	 ////fmm_strict_Lexp->SetParLimits(9,15.,20.);//total scale
	 //fmm_strict_Lexp->FixParameter(10,0.00140370);//sigma
	 //fmm_strict_Lexp->FixParameter(11,0.0158351);//att.
	 //fmm_strict_Lexp->FixParameter(12,-0.0712093);//peak pos.
	 //fmm_strict_Lexp->FixParameter(13,0.276257);//relative strength
	 ////fmm_strict_Lexp->FixParameter(13,0.8);//relative strength
	 break;
	 
case 2: //simc_fitting2022kai (Free)
//simc_fitting2023
	 fmm_strict_Lexp->SetParameter(0,0.000337527);//Landau width
	 fmm_strict_Lexp->SetParameter(1,-0.000399571);
	 fmm_strict_Lexp->SetParameter(2,1.24759);//total scale
	 //fmm_strict_Lexp->SetParLimits(2,20.,30.);//total scale
	 fmm_strict_Lexp->FixParameter(3,0.00129240);//sigma
	 fmm_strict_Lexp->SetParameter(4,0.0499249);//att.
	 fmm_strict_Lexp->SetParameter(5,0.000548956);//peak pos.
	 fmm_strict_Lexp->SetParameter(6,0.716364);//relative strength

	 fmm_strict_Lexp->SetParameter(7,0.000263497);//Landau width
	 fmm_strict_Lexp->SetParameter(8,0.0755465);
	 fmm_strict_Lexp->SetParameter(9,0.304532);//total scale
	 fmm_strict_Lexp->FixParameter(10,0.00138720);//sigma
	 fmm_strict_Lexp->SetParameter(11,0.0145474);//att.
	 fmm_strict_Lexp->SetParameter(12,-0.0712964);//peak pos.
	 fmm_strict_Lexp->SetParameter(13,0.289446);//relative strength
	 break;

case 3: //Fixed from old fit
	 fmm_strict_Lexp->FixParameter(0,0.000658270);//Landau width
	 fmm_strict_Lexp->FixParameter(1,-0.00125805);
	 fmm_strict_Lexp->SetParameter(2,49./2);//total scale
	 fmm_strict_Lexp->FixParameter(3,0.00108656);//sigma
	 fmm_strict_Lexp->FixParameter(4,0.0649346);//att.
	 fmm_strict_Lexp->FixParameter(5,0.003);//peak pos.
	 fmm_strict_Lexp->FixParameter(6,0.692215);//relative strength
	 fmm_strict_Lexp->FixParameter(7,0.000157968);//Landau width
	 fmm_strict_Lexp->FixParameter(8,0.0761140);//MPV
	 fmm_strict_Lexp->SetParameter(9,17./2.);//total scale
	 fmm_strict_Lexp->FixParameter(10,0.00168394);//sigma
	 fmm_strict_Lexp->FixParameter(11,0.110712);//att
	 fmm_strict_Lexp->FixParameter(12,-0.0809590);//peak pos.
	 fmm_strict_Lexp->FixParameter(13,1.49999);
	 break;
}


//change ENum???
	 //fmm_strict_Lexp->FixParameter(14,(double)ENum_strict_cs*0.021*0.001);//scale(1.8%(pion)+0.3%(Al)) //B.G. ratio
	 //fmm_strict_Lexp->FixParameter(14,(double)ENum_strict*0.021*0.001);//scale(1.8%(pion)+0.3%(Al)) //B.G. ratio
	 //fmm_strict_Lexp->FixParameter(14,(double)ENum_strict_Qsq2_2*0.021*0.001);//scale(1.8%(pion)+0.3%(Al)) //B.G. ratio
	 fmm_strict_Lexp->FixParameter(14,(double)ENum_strict_cs_Qsq2_2*0.021*0.001);//scale(1.8%(pion)+0.3%(Al)) //B.G. ratio
	 fmm_strict_Lexp->FixParameter(15,Al_par1);//mean
	 fmm_strict_Lexp->FixParameter(16,Al_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(17,Al_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(18,pion_par1);//mean
	 fmm_strict_Lexp->FixParameter(19,pion_par2);//Gsigma
	 fmm_strict_Lexp->FixParameter(20,pion_par3);//Lfwhm
	 fmm_strict_Lexp->FixParameter(21,0.3/1.8);//Al vs Pi
//	 fmm_strict_Lexp->FixParameter(21,0.1);//Al vs Pi
//	 fmm_strict_Lexp->SetParLimits(21,0.0,0.2);
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
#endif

	 TH1F* h_Resulting  = new TH1F("h_Resulting","h_Resulting",xbin,xmin,xmax);
	 double val, eval;
	 for(int i=0;i<xbin;i++){
		val  = hmm_wo_bg_fom_strict->GetBinContent(i+1);
		eval = hmm_wo_bg_fom_strict->GetBinError(i+1);
		if(val>=0){
		h_Resulting->SetBinContent(i+1,val);
		h_Resulting->SetBinError(i+1,eval);
		}else{continue;}
	 }//Filling to the resulting histogram for chi2 fitting
	 //change
	 h_Resulting->Fit("fmm_strict_Lexp","LLI","",fit_min_mm,fit_max_mm);//Total fitting (full) w/ 4Poly BG
	 //hmm_wo_bg_fom_strict->Fit("fmm_strict_Lexp","I","",fit_min_mm,fit_max_mm);//Total fitting (full) w/ 4Poly BG
	 //hmm_wo_bg_fom_strict->Fit("fmm_strict_Lexp","LLI","",fit_min_mm,fit_max_mm);//Total fitting (div.) w/ 4Poly BG
	 double chisq_strict = fmm_strict_Lexp->GetChisquare();
	 double dof_strict  = fmm_strict_Lexp->GetNDF();
	 cout<<"chisq_strict="<<chisq_strict<<endl;
	 cout<<"dof="<<dof_strict<<endl;
	 cout<<"Reduced chi-square = "<<chisq_strict/dof_strict<<endl;


	 TF1* fmm_Lambda_only_strict = new TF1("fmm_Lambda_only_strict",FMM_Response,fmin_mm,fmax_mm,7);
	 TF1* fmm_Sigma_only_strict  = new TF1("fmm_Sigma_only_strict" ,FMM_Response, fmin_mm,fmax_mm,7);
	 TF1* fmm_bg_only_strict  = new TF1("fmm_bg_only_strict" ,FMM_2BG, fmin_mm,fmax_mm,8);
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

	double nofL_strict_inside = fmm_Lambda_only_strict->Integral(fmin_mm_inside,fmax_mm_inside);
	double nofL_strict = fmm_Lambda_only_strict->Integral(fmin_mm,fmax_mm);
	double nofL_old_strict = fmm_Lambda_only_strict->Integral(-0.006,0.006);
	double nofL_bg_strict = fmm_bg_only_strict->Integral(-0.006,0.006);
	nofL_strict_inside = nofL_strict_inside/fit_bin_width;
	nofL_strict = nofL_strict/fit_bin_width;
	nofL_old_strict = nofL_old_strict/fit_bin_width;
	nofL_bg_strict = nofL_bg_strict/fit_bin_width;
	double nofL_old_hist_strict=hmm_wo_bg_fom_strict->Integral(hmm_wo_bg_fom_strict->FindBin(-0.006),hmm_wo_bg_fom_strict->FindBin(0.006)-nofL_bg_strict);
	cout<<"Number of Lambda (TF1 Integral, inside) = "<<nofL_strict_inside<<endl;
	cout<<"Number of Lambda (TF1 Integral, out of acceptance) = "<<nofL_strict<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old_strict<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<nofL_old_hist_strict<<endl;

	double nofS_strict_inside = fmm_Sigma_only_strict->Integral(fmin_mm_inside,fmax_mm_insideS);
	double nofS_strict = fmm_Sigma_only_strict->Integral(fmin_mm,fmax_mm);
	double nofS_old_strict = fmm_Sigma_only_strict->Integral(def_mean_S-0.008,def_mean_S+0.008);
	double nofS_bg_strict = fmm_bg_only_strict->Integral(def_mean_S-0.008,def_mean_S+0.008);
	nofS_strict_inside = nofS_strict_inside/fit_bin_width;
	nofS_strict = nofS_strict/fit_bin_width;
	nofS_old_strict = nofS_old_strict/fit_bin_width;
	nofS_bg_strict = nofS_bg_strict/fit_bin_width;
	double nofS_old_hist_strict=hmm_wo_bg_fom_strict->Integral(hmm_wo_bg_fom_strict->FindBin(def_mean_S-0.008),hmm_wo_bg_fom_strict->FindBin(def_mean_S+0.008))-nofS_bg_strict;
	cout<<"Number of Sigma (TF1 Integral, inside) = "<<nofS_strict_inside<<endl;
	cout<<"Number of Sigma (TF1 Integral, out of acceptance) = "<<nofS_strict<<endl;
	cout<<"Number of Sigma w/o radiative tail (TF1 Integral) = "<<nofS_old_strict<<endl;
	cout<<"Number of Sigma w/o radiative tail (TH1F Integral) = "<<nofS_old_hist_strict<<endl;

/***************************************/
/************Display********************/
/***************************************/
	 bool displayON = false;//OFF
	 //bool displayON = true;//ON

	 TCanvas* c2plus = new TCanvas("c2plus","Resulting Missing Mass Spectrum (if hmm)");
	 TH1D* hf2 = (TH1D*)c2plus->DrawFrame(-0.1,-20.,0.2,180.);
	 hf2->GetXaxis()->SetTitle("Missing Mass - M_{#Lambda} [MeV/c^{2}]");
	 hf2->GetYaxis()->SetTitle("Counts/(MeV/c^{2})");
	 hf2->GetXaxis()->SetTitleOffset(1.10);
	 hf2->GetYaxis()->SetTitleOffset(1.10);
	 hf2->GetXaxis()->SetLabelOffset(100);
	 hf2->GetYaxis()->SetLabelOffset(100);
	 TGaxis* ax_strict = new TGaxis(-0.1,-20.,0.2,-20.,-100.,200.,510);
	 TGaxis* ay_strict = new TGaxis(-0.1,-20.,-0.1,180.,-20.,180.,510);
	 if(displayON){
	 hmm_wo_bg_fom_strict->Draw("same");
	 ax_strict->Draw("same");
	 ay_strict->Draw("same");
	 }else{hmm_wo_bg_fom_strict->Draw("");}

	 fmm_strict_Lexp->SetLineColor(kRed);
	 fmm_Lambda_only_strict->SetLineColor(kAzure);
	 fmm_Sigma_only_strict->SetLineColor(kCyan);
	 fmm_Lambda_only_strict->SetFillStyle(3004);
	 fmm_Lambda_only_strict->SetFillColor(kAzure);
	 fmm_Lambda_only_strict->SetLineWidth(2);
	 fmm_Sigma_only_strict->SetFillStyle(3005);
	 fmm_Sigma_only_strict->SetFillColor(kCyan);
	 fmm_Sigma_only_strict->SetLineWidth(2);
	 fmm_bg_only_strict->SetLineColor(kOrange);
	 fmm_bg_only_strict->SetFillStyle(3017);
	 fmm_bg_only_strict->SetFillColor(kOrange);
	 fmm_bg_only_strict->SetLineWidth(2);
	 fmm_bg_only_strict->Draw("same");
	 fmm_Lambda_only_strict->Draw("same");
	 fmm_Sigma_only_strict->Draw("same");
	 fmm_strict_Lexp->Draw("same");
	 TLegend *tl_strict = new TLegend(0.55,0.65,0.85,0.85,Form("#chi^{2}/ndf=%d/%d",(int)chisq_strict,(int)dof_strict));
	 tl_strict->SetBorderSize(0);
	 //tl_strict->Draw("same");
	 //fL_strict->SetLineColor(kGreen);
	 //fS_strict->SetLineColor(kGreen);
	 //fL_strict->Draw("same");
	 //fS_strict->Draw("same");
	 if(displayON)c2plus->Print("result_2D_2022_temp.pdf");


#if 0
/****************************************/
/************BEST CUT********************/
/****************************************/
cout<<"BEST CUT START"<<endl;
	TCanvas* c3 = new TCanvas("c3","c3");

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

//Lambda//
	 fmm_best_4Poly->SetParLimits(2,0.,1000.);//positive
	 fmm_best_4Poly->SetParLimits(9,0.,300.);//positive
	 //fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_best_4Poly->SetParLimits(3,0.,0.01);
	 //fmm_best_4Poly->SetParLimits(4,0.005,0.08);
	 //fmm_best_4Poly->SetParLimits(5,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 //fmm_best_4Poly->SetParLimits(6,0.,1.5);//relative strength

//default
	 //fmm_best_4Poly->SetParameter(0,0.0007);//Landau width
	 //fmm_best_4Poly->SetParameter(1,mean_L_best);
	 //fmm_best_4Poly->SetParameter(2,48.);//total scale
	 //fmm_best_4Poly->SetParameter(3,0.001);//sigma
	 //fmm_best_4Poly->SetParameter(4,0.05);//att.
	 //fmm_best_4Poly->SetParameter(5,-1.*mean_L_best);//peak pos.
	 //fmm_best_4Poly->SetParameter(6,0.6);//relative strength
//default
	 fmm_best_4Poly->FixParameter(0,0.000658270);//Landau width
	 //fmm_best_4Poly->SetParLimits(0,0.00061,0.00062);//Landau width
	 fmm_best_4Poly->FixParameter(1,-0.00125805);
	 //fmm_best_4Poly->SetParLimits(1,-0.0014,-0.0013);
	 fmm_best_4Poly->SetParameter(2,49./2);//total scale
	 //fmm_best_4Poly->SetParLimits(2,20.,30.);//total scale
	 fmm_best_4Poly->FixParameter(3,0.00108656);//sigma
	 //fmm_best_4Poly->SetParLimits(3,0.0011,0.0012);//sigma
	 fmm_best_4Poly->FixParameter(4,0.0649346);//att.
	 //fmm_best_4Poly->SetParLimits(4,0.05775,0.05777);//att.
	 fmm_best_4Poly->FixParameter(5,0.003);//peak pos.
	 //fmm_best_4Poly->SetParLimits(5,-0.00187,-0.00186);//peak pos.
	 fmm_best_4Poly->FixParameter(6,0.692215);//relative strength
	 //fmm_best_4Poly->SetParLimits(6,0.577,0.578);//relative strength


//SigmaZ//
	 //fmm_best_4Poly->SetParLimits(8,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 //fmm_best_4Poly->SetParLimits(10,0.,0.01);
	 //fmm_best_4Poly->SetParLimits(11,0.03,0.12);
	 //fmm_best_4Poly->SetParLimits(12,-1.*def_mean_S-def_sig_S,-1.*def_mean_S+def_sig_S);
	 //fmm_best_4Poly->SetParLimits(13,0.,1.5);//relative strength

//default
	 //fmm_best_4Poly->SetParameter(7,0.0003);//Landau width
	 //fmm_best_4Poly->SetParameter(8,mean_S_best);//MPV
	 //fmm_best_4Poly->SetParameter(9,14.);//total scale
	 //fmm_best_4Poly->SetParameter(10,0.0015);//sigma
	 //fmm_best_4Poly->SetParameter(11,0.05);//att
	 //fmm_best_4Poly->SetParameter(12,-1.*def_mean_S);//peak pos.
	 //fmm_best_4Poly->SetParameter(13,0.6);
//default
	 fmm_best_4Poly->FixParameter(7,0.000157968);//Landau width
	 //fmm_best_4Poly->SetParLimits(7,0.000300,0.003100);//Landau width
	 fmm_best_4Poly->FixParameter(8,0.0761140);//MPV
	 //fmm_best_4Poly->SetParLimits(8,0.075,0.076);//MPV
	 fmm_best_4Poly->SetParameter(9,17./2.);//total scale
	 //fmm_best_4Poly->SetParLimits(9,3.,10.);//total scale
	 fmm_best_4Poly->FixParameter(10,0.00168394);//sigma
	 //fmm_best_4Poly->SetParLimits(10,0.0012,0.0013);//sigma
	 fmm_best_4Poly->FixParameter(11,0.110712);//att
	 //fmm_best_4Poly->SetParLimits(11,0.025,0.035);//att
	 fmm_best_4Poly->FixParameter(12,-0.0809590);//peak pos.
	 //fmm_best_4Poly->SetParLimits(12,-0.705,-0.704);//peak pos.
	 fmm_best_4Poly->FixParameter(13,1.49999);
	 //fmm_best_4Poly->SetParLimits(13,0.645,0.655);

	 hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","LL","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
	 //hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","","",fit_min_mm,fit_max_mm);//Total fitting w/ 4Poly BG
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
	 fmm_Lambda_only->SetParError(0,fmm_best_4Poly->GetParError(0));
	 fmm_Lambda_only->SetParError(1,fmm_best_4Poly->GetParError(1));
	 fmm_Lambda_only->SetParError(2,fmm_best_4Poly->GetParError(2));
	 fmm_Lambda_only->SetParError(3,fmm_best_4Poly->GetParError(3));
	 fmm_Lambda_only->SetParError(4,fmm_best_4Poly->GetParError(4));
	 fmm_Lambda_only->SetParError(5,fmm_best_4Poly->GetParError(5));
	 fmm_Lambda_only->SetParError(6,fmm_best_4Poly->GetParError(6));
//Sigma_only
	 fmm_Sigma_only->SetNpx(20000);
	 fmm_Sigma_only->SetParameter(0,fmm_best_4Poly->GetParameter(7));
	 fmm_Sigma_only->SetParameter(1,fmm_best_4Poly->GetParameter(8));
	 fmm_Sigma_only->SetParameter(2,fmm_best_4Poly->GetParameter(9));
	 fmm_Sigma_only->SetParameter(3,fmm_best_4Poly->GetParameter(10));
	 fmm_Sigma_only->SetParameter(4,fmm_best_4Poly->GetParameter(11));
	 fmm_Sigma_only->SetParameter(5,fmm_best_4Poly->GetParameter(12));
	 fmm_Sigma_only->SetParameter(6,fmm_best_4Poly->GetParameter(13));
	 fmm_Sigma_only->SetParError(0,fmm_best_4Poly->GetParError(7));
	 fmm_Sigma_only->SetParError(1,fmm_best_4Poly->GetParError(8));
	 fmm_Sigma_only->SetParError(2,fmm_best_4Poly->GetParError(9));
	 fmm_Sigma_only->SetParError(3,fmm_best_4Poly->GetParError(10));
	 fmm_Sigma_only->SetParError(4,fmm_best_4Poly->GetParError(11));
	 fmm_Sigma_only->SetParError(5,fmm_best_4Poly->GetParError(12));
	 fmm_Sigma_only->SetParError(6,fmm_best_4Poly->GetParError(13));

	double nofL = fmm_Lambda_only->Integral(fmin_mm,fmax_mm);
	double nofL2 = fmm_Lambda_only->Integral(-0.05,0.15);
	double nofL_old = fmm_Lambda_only->Integral(-0.006,0.006);
	nofL = nofL/fit_bin_width;
	nofL_old = nofL_old/fit_bin_width;
	cout<<"Number of Lambda (TF1 Integral) = "<<nofL<<endl;
	cout<<"Number of Lambda within (TF1 Integral) = "<<nofL2/fit_bin_width<<endl;
	cout<<"Number of Lambda w/o radiative tail (TF1 Integral) = "<<nofL_old<<endl;
	cout<<"Number of Lambda w/o radiative tail (TH1F Integral) = "<<hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(-0.006),hmm_wo_bg_fom_best->FindBin(0.006))<<endl;

	//double nofS = fmm_Sigma_only->Integral(-0.05,0.085);
	double nofS = fmm_Sigma_only->Integral(fmin_mm,fmax_mm);
	double nofS2 = fmm_Sigma_only->Integral(-0.05,0.15);
	nofS = nofS/fit_bin_width;
	cout<<"Number of Sigma (TF1 Integral) = "<<nofS<<endl;
	cout<<"Number of Sigma within (TF1 Integral) = "<<nofS2/fit_bin_width<<endl;

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


//	TCanvas* c3 = new TCanvas("c3","c3");
//	c3->Divide(2,2);
//	c3->cd(1);
//	gklab_gkcm->Draw("colz");
//	c3->cd(2);
//	gklab_eklab->Draw("colz");
//	c3->cd(3);
//	eelab_eklab->Draw("colz");
//	c3->cd(4);
//	eklab_gkcm->Draw("colz");
//
//	TCanvas* c4 = new TCanvas("c4","c4");
//	//c4->Divide(2,2);
//	//c4->cd(1);
//	h_theta_gk_cm->Draw("");
//	h_theta_gk_cm2->SetLineColor(kRed);
//	h_theta_gk_cm2->Draw("same");
//
//	TCanvas* c5 = new TCanvas("c5","c5");
//	hmm_L_fom_best->SetLineColor(kBlack);
//	hmm_L_fom_best->Draw("");
//	hmm_bg_fom_best->SetLineColor(kGreen);
//	hmm_bg_fom_best->Draw("same");
//
//	TCanvas* c6 = new TCanvas("c6","c6");
//	hmm_L_fom_strict->SetLineColor(kBlack);
//	hmm_L_fom_strict->Draw("");
//	hmm_bg_fom_strict->SetLineColor(kGreen);
//	hmm_bg_fom_strict->Draw("same");

	TCanvas* c8 = new TCanvas("c8","c8");
	hmm_L_new_Qsq2_1->SetLineColor(kBlack);
	hmm_L_new_Qsq2_1->Draw("");
	hmm_bg_new_Qsq2_1->SetLineColor(kGreen);
	hmm_bg_new_Qsq2_1->Draw("same");
	TCanvas* c9 = new TCanvas("c9","c9");
	hmm_L_new_Qsq2_2->SetLineColor(kBlack);
	hmm_L_new_Qsq2_2->Draw("");
	hmm_bg_new_Qsq2_2->SetLineColor(kGreen);
	hmm_bg_new_Qsq2_2->Draw("same");
	TCanvas* c10 = new TCanvas("c10","c10");
	hmm_L_new_Qsq3_1->SetLineColor(kBlack);
	hmm_L_new_Qsq3_1->Draw("");
	hmm_bg_new_Qsq3_1->SetLineColor(kGreen);
	hmm_bg_new_Qsq3_1->Draw("same");
	TCanvas* c11 = new TCanvas("c11","c11");
	hmm_L_new_Qsq3_2->SetLineColor(kBlack);
	hmm_L_new_Qsq3_2->Draw("");
	hmm_bg_new_Qsq3_2->SetLineColor(kGreen);
	hmm_bg_new_Qsq3_2->Draw("same");
	TCanvas* c12 = new TCanvas("c12","c12");
	hmm_L_new_Qsq3_3->SetLineColor(kBlack);
	hmm_L_new_Qsq3_3->Draw("");
	hmm_bg_new_Qsq3_3->SetLineColor(kGreen);
	hmm_bg_new_Qsq3_3->Draw("same");

	TCanvas* c7 = new TCanvas("c7","c7");
	//hmm_wo_bg_fom_best->Draw("");
	fmm_Lambda_only->SetFillStyle(3004);
	fmm_Lambda_only->SetFillColor(kAzure);
	fmm_Lambda_only->SetLineWidth(2);
	fmm_Sigma_only->SetFillStyle(3005);
	fmm_Sigma_only->SetFillColor(kCyan);
	fmm_Sigma_only->SetLineWidth(2);
	fmm_Lambda_only->Draw("");
	fmm_Sigma_only ->Draw("same");
	hmm_wo_bg_fom_best->Draw("same");
	//Acceptance_map->Draw("lego2z");
	
	TCanvas* c20 = new TCanvas("c20","c20");
	 hmm_Al_wobg_fom_noZ_new->Draw("");
	//c20->Divide(2,2);
	//c20->cd(1);
	//h_pe->Draw("");
	//c20->cd(2);
	//h_pk->Draw("");
	//c20->cd(3);
	//h_pepk->Draw("colz");
#endif
#if 0
	TCanvas* cQsq = new TCanvas("cQsq","cQsq");
	h_mm_Qsq->Draw("colz");
	cQsq->SetLeftMargin(0.14);
	cQsq->SetRightMargin(0.14);
	cQsq->SetTopMargin(0.14);
	cQsq->SetBottomMargin(0.14);
	cQsq->Modified();
	cQsq->Update();
	gPad->Modified();
	gPad->Update();
	
	TCanvas* cW = new TCanvas("cW","cW");
	h_mm_W->Draw("colz");
	cW->SetLeftMargin(0.14);
	cW->SetRightMargin(0.14);
	cW->SetTopMargin(0.14);
	cW->SetBottomMargin(0.14);
	cW->Modified();
	cW->Update();
	gPad->Modified();
	gPad->Update();

	TCanvas* cth = new TCanvas("cth","cth");
	h_mm_theta_gk_cm->Draw("colz");
	cth->SetLeftMargin(0.14);
	cth->SetRightMargin(0.14);
	cth->SetTopMargin(0.14);
	cth->SetBottomMargin(0.14);
	cth->Modified();
	cth->Update();
	gPad->Modified();
	gPad->Update();

	TCanvas* cph = new TCanvas("cph","cph");
	h_mm_phi_gk->Draw("colz");
	cph->SetLeftMargin(0.14);
	cph->SetRightMargin(0.14);
	cph->SetTopMargin(0.14);
	cph->SetBottomMargin(0.14);
	cph->Modified();
	cph->Update();
	gPad->Modified();
	gPad->Update();

	TCanvas* ceps = new TCanvas("ceps","ceps");
	h_mm_eps->Draw("colz");
	ceps->SetLeftMargin(0.14);
	ceps->SetRightMargin(0.14);
	ceps->SetTopMargin(0.14);
	ceps->SetBottomMargin(0.14);
	ceps->Modified();
	ceps->Update();
	gPad->Modified();
	gPad->Update();

	cQsq->Print("dthesis_Fig/pdf/h_mm_Qsq.pdf");
	cW->Print("dthesis_Fig/pdf/h_mm_W.pdf");
	cth->Print("dthesis_Fig/pdf/h_mm_th.pdf");
	cph->Print("dthesis_Fig/pdf/h_mm_ph.pdf");
	ceps->Print("dthesis_Fig/pdf/h_mm_eps.pdf");
#endif
	TCanvas* c30 = new TCanvas("c30","c30");
	c30->Divide(2,1);
	c30->cd(1);
	 hmm_Al_wobg_fom_noZ_new->Draw("");
	 fAl_new->Draw("same");
	c30->cd(2);
	 hmm_pi_wobg_fom_best->Draw("");
	 fpion->Draw("same");
cout << "Well done!" << endl;

}//fit
