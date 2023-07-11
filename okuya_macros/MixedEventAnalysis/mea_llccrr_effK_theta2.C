//----------------------------------//
//--  MIXED EVENT ANALYSIS        --//
//----------------------------------//
//
//K. Okuyama (Sep. 23, 2020)
//K. Okuyama (Nov. 19, 2020) CS factor
//K. Okuyama (Dec. 14, 2020) Mom cut 
//K. Okuyama (Jan.  3, 2021) VP Flux Syst. 
//K. Okuyama (Jan. 30, 2022) Kaon Survival Ratio 
//K. Okuyama (Jun. 11, 2022) elist modified, acc-map: 150bin-->100bin
//K. Okuyama (May   2, 2023) phi dependence
//		taken over from mea_llccrr_momDec.C
//
//taken over from mea_lcr.C
//using {left, right, center} bunch
double F1_path(double *x, double *par){//Path Length from Cointime.
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
double F1_sr(double *x, double *par){//Survival Ratio with Path Length (from Cointime.)
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

void mea_llccrr_effK_theta2(){
	string pdfname = "temp.pdf";
	string rootname= "bgmea_temp_theta2.root";
cout << "Output pdf file name is " << pdfname << endl;
cout << "Output root file name is " << rootname << endl;
  
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file (default: h2all2.root)
 // TFile *file = new TFile("../h2all5.root","read");//input file (default: h2all2.root)
  TFile *file_new = new TFile(rootname.c_str(),"recreate");//new root
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
//	tree->Write();
    
//---  DAQ Efficiency ---//
//H2 run (run111157~111222 & run111480~542)
	string daq_file = "../information/daq.dat";//DAQ Efficiency from ELOG
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
	string AcceptanceR_table_z0 = "../information/RHRS_SIMC100bin_10_z0_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z1 = "../information/RHRS_SIMC100bin_10_z1_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z2 = "../information/RHRS_SIMC100bin_10_z2_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z3 = "../information/RHRS_SIMC100bin_10_z3_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z4 = "../information/RHRS_SIMC100bin_10_z4_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z5 = "../information/RHRS_SIMC100bin_10_z5_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z6 = "../information/RHRS_SIMC100bin_10_z6_theta2.dat";//Acceptance Table (SIMC)
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
	string AcceptanceR_table_z7 = "../information/RHRS_SIMC100bin_10_z7_theta2.dat";//Acceptance Table (SIMC)
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
/*----- 2 < z < 4 -----*/
	string AcceptanceR_table_z8 = "../information/RHRS_SIMC100bin_10_z8_theta2.dat";//Acceptance Table (SIMC)
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
/*----- 2 < z < 4 -----*/
	string AcceptanceR_table_z9 = "../information/RHRS_SIMC100bin_10_z9_theta2.dat";//Acceptance Table (SIMC)
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

  TH2F* Acceptance_map = new TH2F("Acceptance_map","#Delta#Omega_{K}^{lab}(p_{K},Z)",100,1.9,2.3,10,-10.,10.);
	for(int i=0;i<100;i++){
		for(int j=0;j<10;j++){
			Acceptance_map->SetBinContent(i+1,j+1,RHRS_table[i][j]*1000.);
		}
	}

 double RHRS = 0.005;
 double effDAQ = 0.95;
 double effK = 0.14;
 double cs = 0.;
 double cs_lab = 0.;



//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
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
 int bin_mm=(max_mm-min_mm)/0.001; //Counts/1 MeV
 bin_mm=(int)bin_mm;

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;


//---------------------------------------//
//               Branch                  //
//---------------------------------------//

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
	int nrun;//RUN NUMBER
	
//	string branchname[]={"tr.ntrack_l","tr.ntrack_r","tr.Ls2_pad[100]"};
//	Int_t nbranch = sizeof(branchname)/sizeof(branchname[0]);
//	tree->SetBranchStatus("*",0);
//	for(Int_t ibranch;ibranch<nbranch;ibranch++){
//	tree->SetBranchStatus(branchname[ibranch].c_str(),1);
//	}

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

//ADD 2020/8/12
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
  TH1F* h2  = new TH1F("h2","",400,-20.,20.0);
  h2->GetXaxis()->SetTitle("coin time (ns)");
  h2->GetYaxis()->SetTitle("Counts / 100 ps");
  h2->GetXaxis()->SetRangeUser(-14.0,17.);
  TH1F* h3  = new TH1F("h3","",400,-20.,20.0);
  h3->GetXaxis()->SetTitle("coin time (ns)");
  h3->GetYaxis()->SetTitle("Counts / 100 ps");
  h3->GetXaxis()->SetRangeUser(-14.0,17.);
  TH1F* h4  = new TH1F("h4","",400,-20.,20.0);
  h4->GetXaxis()->SetTitle("coin time (ns)");
  h4->GetYaxis()->SetTitle("Counts / 100 ps");
  h4->GetXaxis()->SetRangeUser(-14.0,17.);
  double xmin = -0.1, xmax = 0.2; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_best  = new TH1F("hmm_best","hmm_best",xbin,xmin,xmax);
  hmm_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_best->GetYaxis()->SetTitle("Counts / MeV");
  hmm_best->SetLineColor(1);
  TH1F* hmm_nocut  = new TH1F("hmm_nocut","hmm_nocut",xbin,xmin,xmax);
  hmm_nocut->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_nocut->GetYaxis()->SetTitle("Counts / MeV");
  hmm_nocut->SetLineColor(1);
  TH1F* hm2  = new TH1F("hm2","Mix region 1",xbin,xmin,xmax);

  TH1F* hmm_acc  = new TH1F("hmm_acc","ACC (original)",xbin,xmin,xmax);
  TH1F* hmm_mixacc  = new TH1F("hmm_mixacc","ACC (mixed)",xbin,xmin,xmax);
  TH1F* hmm_mixacc2  = new TH1F("hmm_mixacc2","ACC (mixed)",xbin,xmin,xmax);
  TH1F* hmm_mixacc3  = new TH1F("hmm_mixacc3","ACC (mixed)",xbin,xmin,xmax);
  TF1* f1_sr;//Survival Ratio
/*-----------------*/
/*--RESULT OUTPUT--*/
/*-----------------*/
  TH1F* hmm_mixacc_result_best  = new TH1F("hmm_mixacc_result_best","MEA result (best)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_bestCT  = new TH1F("hmm_mixacc_result_bestCT","MEA result (best CT)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut  = new TH1F("hmm_mixacc_result_nocut","MEA result (no Z cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_new  = new TH1F("hmm_mixacc_result_nocut_new","MEA result (no Z, new AC cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new  = new TH1F("hmm_mixacc_result_new","MEA result (strict)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_newCT  = new TH1F("hmm_mixacc_result_newCT","MEA result (strict CT)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_woac  = new TH1F("hmm_mixacc_result_woac","MEA result (no AC cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_zdiff  = new TH1F("hmm_mixacc_result_zdiff","MEA result (Z diff cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_forAl  = new TH1F("hmm_mixacc_result_nocut_forAl","ACC (Al (best))",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_new_forAl  = new TH1F("hmm_mixacc_result_nocut_new_forAl","ACC (Al (strict))",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momcut  = new TH1F("hmm_mixacc_result_momcut","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL1  = new TH1F("hmm_mixacc_result_momL1","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL2  = new TH1F("hmm_mixacc_result_momL2","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL3  = new TH1F("hmm_mixacc_result_momL3","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL4  = new TH1F("hmm_mixacc_result_momL4","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL5  = new TH1F("hmm_mixacc_result_momL5","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL6  = new TH1F("hmm_mixacc_result_momL6","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL8  = new TH1F("hmm_mixacc_result_momL8","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momL9  = new TH1F("hmm_mixacc_result_momL9","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS1  = new TH1F("hmm_mixacc_result_momS1","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS2  = new TH1F("hmm_mixacc_result_momS2","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS3  = new TH1F("hmm_mixacc_result_momS3","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS4  = new TH1F("hmm_mixacc_result_momS4","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS5  = new TH1F("hmm_mixacc_result_momS5","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS6  = new TH1F("hmm_mixacc_result_momS6","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS8  = new TH1F("hmm_mixacc_result_momS8","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momS9  = new TH1F("hmm_mixacc_result_momS9","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL1  = new TH1F("hcs_mixacc_result_momL1","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL2  = new TH1F("hcs_mixacc_result_momL2","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL3  = new TH1F("hcs_mixacc_result_momL3","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL4  = new TH1F("hcs_mixacc_result_momL4","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL5  = new TH1F("hcs_mixacc_result_momL5","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL6  = new TH1F("hcs_mixacc_result_momL6","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL8  = new TH1F("hcs_mixacc_result_momL8","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momL9  = new TH1F("hcs_mixacc_result_momL9","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS1  = new TH1F("hcs_mixacc_result_momS1","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS2  = new TH1F("hcs_mixacc_result_momS2","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS3  = new TH1F("hcs_mixacc_result_momS3","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS4  = new TH1F("hcs_mixacc_result_momS4","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS5  = new TH1F("hcs_mixacc_result_momS5","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS6  = new TH1F("hcs_mixacc_result_momS6","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS8  = new TH1F("hcs_mixacc_result_momS8","MEA result (Mom. cut)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_momS9  = new TH1F("hcs_mixacc_result_momS9","MEA result (Mom. cut)",xbin,xmin,xmax);
// (NY/DAQ_Eff./RHRS)
  TH1F* hcs_mixacc_result_new  = new TH1F("hcs_mixacc_result_new","MEA result (new AC best)",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_lab  = new TH1F("hcs_mixacc_result_new_lab","MEA result (all lab)",xbin,xmin,xmax);
/*-------------------*/
/*--Angle partition--*/
/*-------------------*/
  TH1F* hmm_mixacc_result_cm2_1  = new TH1F("hmm_mixacc_result_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm2_2  = new TH1F("hmm_mixacc_result_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm3_1  = new TH1F("hmm_mixacc_result_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm3_2  = new TH1F("hmm_mixacc_result_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm3_3  = new TH1F("hmm_mixacc_result_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm4_1  = new TH1F("hmm_mixacc_result_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm4_2  = new TH1F("hmm_mixacc_result_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm4_3  = new TH1F("hmm_mixacc_result_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_cm4_4  = new TH1F("hmm_mixacc_result_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm2_1  = new TH1F("hmm_mixacc_result_new_cm2_1","(tight) #theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm2_2  = new TH1F("hmm_mixacc_result_new_cm2_2","(tight) #theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm3_1  = new TH1F("hmm_mixacc_result_new_cm3_1","(tight) #theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm3_2  = new TH1F("hmm_mixacc_result_new_cm3_2","(tight) 6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm3_3  = new TH1F("hmm_mixacc_result_new_cm3_3","(tight) #theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm4_1  = new TH1F("hmm_mixacc_result_new_cm4_1","(tight) #theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm4_2  = new TH1F("hmm_mixacc_result_new_cm4_2","(tight) 5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm4_3  = new TH1F("hmm_mixacc_result_new_cm4_3","(tight) 8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_cm4_4  = new TH1F("hmm_mixacc_result_new_cm4_4","(tight) #theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_Qsq2_1  = new TH1F("hmm_mixacc_result_new_Qsq2_1","(tight) Q^{2}<0.5 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_Qsq2_2  = new TH1F("hmm_mixacc_result_new_Qsq2_2","(tight) Q^{2}>0.5 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_Qsq3_1  = new TH1F("hmm_mixacc_result_new_Qsq3_1","(tight) Q^{2}<0.45 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_Qsq3_2  = new TH1F("hmm_mixacc_result_new_Qsq3_2","(tight) 0.45<Q^{2}<0.55 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_Qsq3_3  = new TH1F("hmm_mixacc_result_new_Qsq3_3","(tight) Q^{2}>0.55 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_phi1_1  = new TH1F("hmm_mixacc_result_new_phi1_1","(tight) |#phi'|<=#pi/2",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_phi2_1  = new TH1F("hmm_mixacc_result_new_phi2_1","(tight) |#phi'|<=#pi/4",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new_phi2_2  = new TH1F("hmm_mixacc_result_new_phi2_2","(tight) #pi/4<|#phi'|<=#pi/2",xbin,xmin,xmax);
//hmm => hcs(Lab) 2020/11/22
  TH1F* hcs_mixacc_result_cm2_1  = new TH1F("hcs_mixacc_result_cm2_1","#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm2_2  = new TH1F("hcs_mixacc_result_cm2_2","#theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm3_1  = new TH1F("hcs_mixacc_result_cm3_1","#theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm3_2  = new TH1F("hcs_mixacc_result_cm3_2","6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm3_3  = new TH1F("hcs_mixacc_result_cm3_3","#theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm4_1  = new TH1F("hcs_mixacc_result_cm4_1","#theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm4_2  = new TH1F("hcs_mixacc_result_cm4_2","5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm4_3  = new TH1F("hcs_mixacc_result_cm4_3","8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_cm4_4  = new TH1F("hcs_mixacc_result_cm4_4","#theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm2_1  = new TH1F("hcs_mixacc_result_new_cm2_1","(tight) #theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm2_2  = new TH1F("hcs_mixacc_result_new_cm2_2","(tight) #theta_{#gamma K}^{CM}>8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm3_1  = new TH1F("hcs_mixacc_result_new_cm3_1","(tight) #theta_{#gamma K}^{CM}<6 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm3_2  = new TH1F("hcs_mixacc_result_new_cm3_2","(tight) 6<#theta_{#gamma K}^{CM}<10 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm3_3  = new TH1F("hcs_mixacc_result_new_cm3_3","(tight) #theta_{#gamma K}^{CM}>10 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm4_1  = new TH1F("hcs_mixacc_result_new_cm4_1","(tight) #theta_{#gamma K}^{CM}<5 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm4_2  = new TH1F("hcs_mixacc_result_new_cm4_2","(tight) 5<#theta_{#gamma K}^{CM}<8 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm4_3  = new TH1F("hcs_mixacc_result_new_cm4_3","(tight) 8<#theta_{#gamma K}^{CM}<11 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_cm4_4  = new TH1F("hcs_mixacc_result_new_cm4_4","(tight) #theta_{#gamma K}^{CM}>11 deg",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_Qsq2_1  = new TH1F("hcs_mixacc_result_new_Qsq2_1","(tight) Q^{2}<0.5 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_Qsq2_2  = new TH1F("hcs_mixacc_result_new_Qsq2_2","(tight) Q^{2}>0.5 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_Qsq3_1  = new TH1F("hcs_mixacc_result_new_Qsq3_1","(tight) Q^{2}<0.45 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_Qsq3_2  = new TH1F("hcs_mixacc_result_new_Qsq3_2","(tight) 0.45<Q^{2}<0.55 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_Qsq3_3  = new TH1F("hcs_mixacc_result_new_Qsq3_3","(tight) Q^{2}>0.55 (GeV/c)^{2}",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_phi1_1  = new TH1F("hcs_mixacc_result_new_phi1_1","(tight) |#phi'|<=#pi/2",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_phi2_1  = new TH1F("hcs_mixacc_result_new_phi2_1","(tight) |#phi'|<=#pi/4",xbin,xmin,xmax);
  TH1F* hcs_mixacc_result_new_phi2_2  = new TH1F("hcs_mixacc_result_new_phi2_2","(tight) #pi/4<|#phi'|<=#pi/2",xbin,xmin,xmax);
/*-----------------*/
/*--RESULT OUTPUT--*/
/*-----------------*/
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);
  h2->SetLineColor(1);
  h3->SetLineColor(1);
  h3->SetFillColor(1);
  h3->SetFillStyle(3001);
  h4->SetLineColor(9);


  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool event_selection = false;
  bool event_selection_new = false;
  bool event_selection_woac = false;
  bool event_selection_zdiff = false;
  bool event_selection_nocut = false;
  bool event_selection_nocut_new = false;
  bool event_selection_momcut = false;
  bool event_selection_momL1 = false;//VP Flux Syst.
  bool event_selection_momL2 = false;//VP Flux Syst.
  bool event_selection_momL3 = false;//VP Flux Syst.
  bool event_selection_momL4 = false;//VP Flux Syst.
  bool event_selection_momL5 = false;//VP Flux Syst.
  bool event_selection_momL6 = false;//VP Flux Syst.
  bool event_selection_momL8 = false;//VP Flux Syst.
  bool event_selection_momL9 = false;//VP Flux Syst.
  bool event_selection_momS1 = false;//VP Flux Syst.
  bool event_selection_momS2 = false;//VP Flux Syst.
  bool event_selection_momS3 = false;//VP Flux Syst.
  bool event_selection_momS4 = false;//VP Flux Syst.
  bool event_selection_momS5 = false;//VP Flux Syst.
  bool event_selection_momS6 = false;//VP Flux Syst.
  bool event_selection_momS8 = false;//VP Flux Syst.
  bool event_selection_momS9 = false;//VP Flux Syst.
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
  bool phi1_1_cut=false;
  bool phi2_1_cut=false;
  bool phi2_2_cut=false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  bool mix_region1 = false;
  bool mix_region2 = false;
  bool mix_region3 = false;
  bool mix_region4 = false;
  bool mix_region5 = false;
  bool mix_region6 = false;
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum


  int ENum=0;
 // tree->Draw(">>elist", "abs(ct_orig+7*2.0)<1.0||abs(ct_orig+2*2.0)<1.0||abs(ct_orig+1.0*2.0)<1.0||abs(ct_orig-4.0*2.0)<1.0||abs(ct_orig-5.0*2.0)<1.0||abs(ct_orig-6.0*2.0)<1.0");
  //tree->Draw(">>elist", "abs(ct_orig+9.0*2.012)<1.006||abs(ct_orig+8.0*2.012)<1.006||abs(ct_orig+2.0*2.012)<1.006||abs(ct_orig+1.0*2.012)<1.006||abs(ct_orig-6.0*2.012)<1.006||abs(ct_orig-7.0*2.012)<1.006&&Rp_c>1.760&&Rp_c<1.900&&2.010<Lp_c&&Lp_c<2.160");
  tree->Draw(">>elist", "(abs(ct_orig+9.0*2.012)<1.006||abs(ct_orig+8.0*2.012)<1.006||abs(ct_orig+2.0*2.012)<1.006||abs(ct_orig+1.0*2.012)<1.006||abs(ct_orig-6.0*2.012)<1.006||abs(ct_orig-7.0*2.012)<1.006)&&Rp_c>1.760&&Rp_c<1.900&&2.010<Lp_c&&Lp_c<2.160");//2022/6/11
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  ENum = elist->GetN(); 
  //ENum = tree->GetEntries();
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/100;
  if(ENum<100000)time_div=10000;
  int nmix = 750;//Num of mix
  double mass,mm;
  TLorentzVector L_4vec;
  TLorentzVector R_4vec;
  TLorentzVector B_4vec;
  TLorentzVector T_4vec;
  TLorentzVector G_4vec;
  TVector3 boost;
  TLorentzVector GT_4vec;
  TLorentzVector Missing;
  TLorentzVector L_4vec_saved, B_4vec_saved, T_4vec_saved;

	time_t start, end;
	start = time(NULL);
	time(&start);

int best=0;
int ac_new=0;
int woz=0;
int zdiff=0;
int woac=0;

//***************************//
//  MIXED! EVENT! ANALYSIS!  //
//***************************//
cout<<"MIXED! EVENT! ANALYSIS!"<<endl;
  for(int i=0;i<ENum;i++){
	    //tree->GetEntry(i);
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


	


		if(abs(ct+7.0*rf_bunch)<1.0) mix_region1 = true;
    	else mix_region1 = false;
    	if(abs(ct+2.0*rf_bunch)<1.0) mix_region2 = true;
    	else mix_region2 = false;
    	if(abs(ct+1.0*rf_bunch)<1.0) mix_region3 = true;
    	else mix_region3 = false;
    	if(abs(ct-4.0*rf_bunch)<1.0) mix_region4 = true;
    	else mix_region4 = false;
    	if(abs(ct-5.0*rf_bunch)<1.0) mix_region5 = true;
    	else mix_region5 = false;
    	if(abs(ct-6.0*rf_bunch)<1.0) mix_region6 = true;
    	else mix_region6 = false;
		//TEventList check
//		if(mix_region1==false&&mix_region2==false&&mix_region3==false&&mix_region4==false&&mix_region5==false&&mix_region6==false)cout<<"Weird!"<<endl;
		//tree_out (in h2all.cc) check
		if(L_Tr==false||R_Tr==false||L_FP==false||R_FP==false)cout<<"Weird! (Tr, FP)"<<endl;


		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
		else event_selection_nocut=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut_new=true;
		else event_selection_nocut_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_new=true;
		else event_selection_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&L_mom>2.12&&L_mom<2.18&&R_mom>1.81&&R_mom<1.88)event_selection_momcut=true;
		else event_selection_momcut=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_woac=true;
		else event_selection_woac=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_zdiff=true;
		else event_selection_zdiff=false;



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


		TLorentzVector L_4vec;//HRS-L
		TLorentzVector R_4vec;//HRS-R
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target(H2)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);

	    R_4vec.RotateX(  13.2/180.*PI );
	    L_4vec.RotateX( -13.2/180.*PI );

		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
	    mm=mass - mh;//shift by ML
		//if(event_selection)hmm_acc->Fill(mm);
		double before_mm=mm;
//cout<<"before_mm="<<before_mm<<endl;


		/*--Electron information is saved.--*/
		L_4vec_saved = L_4vec;
		B_4vec_saved = B_4vec;
		T_4vec_saved = T_4vec;
		L_tr_vz_saved  = L_tr_vz;
		/*--Electron information is saved.--*/
	

		
	for(int j=0 ; j<(int)nmix ; j++){
	  
	  int ENum_mixed = i+j;

	  if(ENum_mixed<ENum){
	    //tree->GetEntry(ENum_mixed);
	    tree->GetEntry(elist->GetEntry(i+j));
	  }
	  else {
	    //tree->GetEntry(ENum_mixed-ENum);
	    tree->GetEntry(elist->GetEntry(i+j-ENum));
	  }


		//Electron info.
		L_tr_vz=L_tr_vz_saved;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2<0.01 ) R_Tr = true;
        if( R_tr_th<0.17*R_tr_x+0.025
         && R_tr_th>0.17*R_tr_x-0.035
         && R_tr_th<0.40*R_tr_x+0.130 ) R_FP = true;


	// Change Only Kaon information
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
	    double R_px = R_pz * (R_tr_tg_th );
	    double R_py = R_pz * ( R_tr_tg_ph );
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
		R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
	    R_4vec.RotateX(  13.2/180.*PI );
	

	
	// Electron information is not changed.
		L_4vec = L_4vec_saved;
		B_4vec = B_4vec_saved;
		T_4vec = T_4vec_saved;

	// Event Selection Again (R_tr_vz is changed)/2020/11/21
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
		else event_selection_nocut=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut_new=true;
		else event_selection_nocut_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_new=true;
		else event_selection_new=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.89)event_selection_momL1=true;
		else event_selection_momL1=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.88)event_selection_momL2=true;
		else event_selection_momL2=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.87)event_selection_momL3=true;
		else event_selection_momL3=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.86)event_selection_momL4=true;
		else event_selection_momL4=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.85)event_selection_momL5=true;
		else event_selection_momL5=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.84)event_selection_momL6=true;
		else event_selection_momL6=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.91)event_selection_momL8=true;
		else event_selection_momL8=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom<1.92)event_selection_momL9=true;
		else event_selection_momL9=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.77)event_selection_momS1=true;
		else event_selection_momS1=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.78)event_selection_momS2=true;
		else event_selection_momS2=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.79)event_selection_momS3=true;
		else event_selection_momS3=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.80)event_selection_momS4=true;
		else event_selection_momS4=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.81)event_selection_momS5=true;
		else event_selection_momS5=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.82)event_selection_momS6=true;
		else event_selection_momS6=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.75)event_selection_momS8=true;
		else event_selection_momS8=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.74)event_selection_momS9=true;
		else event_selection_momS9=false;

		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&L_mom>2.12&&L_mom<2.18&&R_mom>1.81&&R_mom<1.88)event_selection_momcut=true;
		else event_selection_momcut=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_woac=true;
		else event_selection_woac=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_zdiff=true;
		else event_selection_zdiff=false;

        double mass_mixed,mm_mixed;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass_mixed = Missing.M();
	    mm_mixed=mass_mixed - mh;//shift by ML
		double after_mm=mm_mixed;
		//cout<<"mm(Before):mm(after)="<<before_mm<<":"<<after_mm<<" (i,j)=("<<i<<","<<j<<")"<<endl;
		
		
//////////////////////////		
// CM frame information //
//////////////////////////		
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%DAQ Eff. & HRS-R Acceptance %%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

		int kbin = (int)((R_mom-1.6)*100./0.4);
		int zbin = (int)((((L_tr_vz+R_tr_vz)/2.)-0.1)/0.04);
		if(event_selection&&kbin>=0 &&kbin<100){
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
		//PathLength (Jan. 30, 2022)
		f1_sr = new TF1("f1_sr", F1_sr, 1730.,1930.,2);
		f1_sr->SetParameter(0,31.2307);
		f1_sr->SetParameter(1,-0.0110408);
		effK = f1_sr->Eval(R_mom*1000.);
		//cout<<"effK="<<effK<<endl;
		//cout<<"mom="<<R_mom<<endl;
		//-------------
		if(RHRS!=0.&&effDAQ!=0.){
			cs = labtocm/effDAQ/RHRS/effK;//[nb/sr]
			cs_lab = 1./effDAQ/RHRS/effK;//[nb/sr]
		}else{cs=0.;cs_lab=0.;}
		}else{cs=0.;cs_lab=0.;}

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
    phi1_1_cut=false;
    phi2_1_cut=false;
    phi2_2_cut=false;
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
	if(abs(phi_k-PI)<=PI/2.)phi1_1_cut=true;
	if(abs(phi_k-PI)<=PI/4.)phi2_1_cut=true;
	if(abs(phi_k-PI)>PI/4.&&abs(phi_k-PI)<=PI/2.)phi2_2_cut=true;

			if(mix_region1)hm2->Fill(mm_mixed);
			if(event_selection)hmm_mixacc->Fill(mm_mixed);
			if(event_selection)hmm_mixacc2->Fill(mm_mixed);
			if(event_selection)hmm_mixacc3->Fill(mm_mixed);
			if(event_selection)hmm_mixacc_result_best->Fill(mm_mixed);
			if(event_selection)hmm_mixacc_result_bestCT->Fill(mm_mixed);
			if(event_selection_nocut)hmm_mixacc_result_nocut->Fill(mm_mixed);
			if(event_selection_nocut_new)hmm_mixacc_result_nocut_new->Fill(mm_mixed);
			if(event_selection_new)hmm_mixacc_result_new->Fill(mm_mixed);
			if(event_selection_new)hmm_mixacc_result_newCT->Fill(mm_mixed);
			if(event_selection_new)hcs_mixacc_result_new->Fill(mm_mixed,cs);
			if(event_selection_new)hcs_mixacc_result_new_lab->Fill(mm_mixed,cs_lab);
			if(event_selection_momcut)hmm_mixacc_result_momcut->Fill(mm_mixed);
			if(event_selection_momL1)hmm_mixacc_result_momL1->Fill(mm_mixed);
			if(event_selection_momL2)hmm_mixacc_result_momL2->Fill(mm_mixed);
			if(event_selection_momL3)hmm_mixacc_result_momL3->Fill(mm_mixed);
			if(event_selection_momL4)hmm_mixacc_result_momL4->Fill(mm_mixed);
			if(event_selection_momL5)hmm_mixacc_result_momL5->Fill(mm_mixed);
			if(event_selection_momL6)hmm_mixacc_result_momL6->Fill(mm_mixed);
			if(event_selection_momL8)hmm_mixacc_result_momL8->Fill(mm_mixed);
			if(event_selection_momL9)hmm_mixacc_result_momL9->Fill(mm_mixed);
			if(event_selection_momS1)hmm_mixacc_result_momS1->Fill(mm_mixed);
			if(event_selection_momS2)hmm_mixacc_result_momS2->Fill(mm_mixed);
			if(event_selection_momS3)hmm_mixacc_result_momS3->Fill(mm_mixed);
			if(event_selection_momS4)hmm_mixacc_result_momS4->Fill(mm_mixed);
			if(event_selection_momS5)hmm_mixacc_result_momS5->Fill(mm_mixed);
			if(event_selection_momS6)hmm_mixacc_result_momS6->Fill(mm_mixed);
			if(event_selection_momS8)hmm_mixacc_result_momS8->Fill(mm_mixed);
			if(event_selection_momS9)hmm_mixacc_result_momS9->Fill(mm_mixed);
			if(event_selection_momL1)hcs_mixacc_result_momL1->Fill(mm_mixed,cs);
			if(event_selection_momL2)hcs_mixacc_result_momL2->Fill(mm_mixed,cs);
			if(event_selection_momL3)hcs_mixacc_result_momL3->Fill(mm_mixed,cs);
			if(event_selection_momL4)hcs_mixacc_result_momL4->Fill(mm_mixed,cs);
			if(event_selection_momL5)hcs_mixacc_result_momL5->Fill(mm_mixed,cs);
			if(event_selection_momL6)hcs_mixacc_result_momL6->Fill(mm_mixed,cs);
			if(event_selection_momL8)hcs_mixacc_result_momL8->Fill(mm_mixed,cs);
			if(event_selection_momL9)hcs_mixacc_result_momL9->Fill(mm_mixed,cs);
			if(event_selection_momS1)hcs_mixacc_result_momS1->Fill(mm_mixed,cs);
			if(event_selection_momS2)hcs_mixacc_result_momS2->Fill(mm_mixed,cs);
			if(event_selection_momS3)hcs_mixacc_result_momS3->Fill(mm_mixed,cs);
			if(event_selection_momS4)hcs_mixacc_result_momS4->Fill(mm_mixed,cs);
			if(event_selection_momS5)hcs_mixacc_result_momS5->Fill(mm_mixed,cs);
			if(event_selection_momS6)hcs_mixacc_result_momS6->Fill(mm_mixed,cs);
			if(event_selection_momS8)hcs_mixacc_result_momS8->Fill(mm_mixed,cs);
			if(event_selection_momS9)hcs_mixacc_result_momS9->Fill(mm_mixed,cs);
			if(event_selection_woac)hmm_mixacc_result_woac->Fill(mm_mixed);
			if(event_selection_zdiff)hmm_mixacc_result_zdiff->Fill(mm_mixed);
			if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01))hmm_mixacc_result_nocut_forAl->Fill(mm_mixed);//Al selection
			if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&(fabs(R_tr_vz-L_tr_vz)<0.025)&&(fabs(fabs(R_tr_vz+L_tr_vz)/2.-0.12)<0.01||fabs(fabs(R_tr_vz+L_tr_vz)/2.+0.12)<0.01))hmm_mixacc_result_nocut_new_forAl->Fill(mm_mixed);//Al selection
			if(event_selection&&cm2_angle1_cut)hmm_mixacc_result_cm2_1->Fill(mm_mixed);
			if(event_selection&&cm2_angle2_cut)hmm_mixacc_result_cm2_2->Fill(mm_mixed);
			if(event_selection&&cm3_angle1_cut)hmm_mixacc_result_cm3_1->Fill(mm_mixed);
			if(event_selection&&cm3_angle2_cut)hmm_mixacc_result_cm3_2->Fill(mm_mixed);
			if(event_selection&&cm3_angle3_cut)hmm_mixacc_result_cm3_3->Fill(mm_mixed);
			if(event_selection&&cm4_angle1_cut)hmm_mixacc_result_cm4_1->Fill(mm_mixed);
			if(event_selection&&cm4_angle2_cut)hmm_mixacc_result_cm4_2->Fill(mm_mixed);
			if(event_selection&&cm4_angle3_cut)hmm_mixacc_result_cm4_3->Fill(mm_mixed);
			if(event_selection&&cm4_angle4_cut)hmm_mixacc_result_cm4_4->Fill(mm_mixed);
			if(event_selection_new&&cm2_angle1_cut)hmm_mixacc_result_new_cm2_1->Fill(mm_mixed);
			if(event_selection_new&&cm2_angle2_cut)hmm_mixacc_result_new_cm2_2->Fill(mm_mixed);
			if(event_selection_new&&cm3_angle1_cut)hmm_mixacc_result_new_cm3_1->Fill(mm_mixed);
			if(event_selection_new&&cm3_angle2_cut)hmm_mixacc_result_new_cm3_2->Fill(mm_mixed);
			if(event_selection_new&&cm3_angle3_cut)hmm_mixacc_result_new_cm3_3->Fill(mm_mixed);
			if(event_selection_new&&cm4_angle1_cut)hmm_mixacc_result_new_cm4_1->Fill(mm_mixed);
			if(event_selection_new&&cm4_angle2_cut)hmm_mixacc_result_new_cm4_2->Fill(mm_mixed);
			if(event_selection_new&&cm4_angle3_cut)hmm_mixacc_result_new_cm4_3->Fill(mm_mixed);
			if(event_selection_new&&cm4_angle4_cut)hmm_mixacc_result_new_cm4_4->Fill(mm_mixed);
			if(event_selection_new&&Qsq2_1_cut)hmm_mixacc_result_new_Qsq2_1->Fill(mm_mixed);
			if(event_selection_new&&Qsq2_2_cut)hmm_mixacc_result_new_Qsq2_2->Fill(mm_mixed);
			if(event_selection_new&&Qsq3_1_cut)hmm_mixacc_result_new_Qsq3_1->Fill(mm_mixed);
			if(event_selection_new&&Qsq3_2_cut)hmm_mixacc_result_new_Qsq3_2->Fill(mm_mixed);
			if(event_selection_new&&Qsq3_3_cut)hmm_mixacc_result_new_Qsq3_3->Fill(mm_mixed);
			if(event_selection_new&&phi1_1_cut)hmm_mixacc_result_new_phi1_1->Fill(mm_mixed);
			if(event_selection_new&&phi2_1_cut)hmm_mixacc_result_new_phi2_1->Fill(mm_mixed);
			if(event_selection_new&&phi2_2_cut)hmm_mixacc_result_new_phi2_2->Fill(mm_mixed);

//hmm => hcs(Lab) 2020/11/22
			if(event_selection&&cm2_angle1_cut)hcs_mixacc_result_cm2_1->Fill(mm_mixed,cs);
			if(event_selection&&cm2_angle2_cut)hcs_mixacc_result_cm2_2->Fill(mm_mixed,cs);
			if(event_selection&&cm3_angle1_cut)hcs_mixacc_result_cm3_1->Fill(mm_mixed,cs);
			if(event_selection&&cm3_angle2_cut)hcs_mixacc_result_cm3_2->Fill(mm_mixed,cs);
			if(event_selection&&cm3_angle3_cut)hcs_mixacc_result_cm3_3->Fill(mm_mixed,cs);
			if(event_selection&&cm4_angle1_cut)hcs_mixacc_result_cm4_1->Fill(mm_mixed,cs);
			if(event_selection&&cm4_angle2_cut)hcs_mixacc_result_cm4_2->Fill(mm_mixed,cs);
			if(event_selection&&cm4_angle3_cut)hcs_mixacc_result_cm4_3->Fill(mm_mixed,cs);
			if(event_selection&&cm4_angle4_cut)hcs_mixacc_result_cm4_4->Fill(mm_mixed,cs);
			if(event_selection_new&&cm2_angle1_cut)hcs_mixacc_result_new_cm2_1->Fill(mm_mixed,cs);
			if(event_selection_new&&cm2_angle2_cut)hcs_mixacc_result_new_cm2_2->Fill(mm_mixed,cs);
			if(event_selection_new&&cm3_angle1_cut)hcs_mixacc_result_new_cm3_1->Fill(mm_mixed,cs);
			if(event_selection_new&&cm3_angle2_cut)hcs_mixacc_result_new_cm3_2->Fill(mm_mixed,cs);
			if(event_selection_new&&cm3_angle3_cut)hcs_mixacc_result_new_cm3_3->Fill(mm_mixed,cs);
			if(event_selection_new&&cm4_angle1_cut)hcs_mixacc_result_new_cm4_1->Fill(mm_mixed,cs);
			if(event_selection_new&&cm4_angle2_cut)hcs_mixacc_result_new_cm4_2->Fill(mm_mixed,cs);
			if(event_selection_new&&cm4_angle3_cut)hcs_mixacc_result_new_cm4_3->Fill(mm_mixed,cs);
			if(event_selection_new&&cm4_angle4_cut)hcs_mixacc_result_new_cm4_4->Fill(mm_mixed,cs);
			if(event_selection_new&&Qsq2_1_cut)hcs_mixacc_result_new_Qsq2_1->Fill(mm_mixed,cs);
			if(event_selection_new&&Qsq2_2_cut)hcs_mixacc_result_new_Qsq2_2->Fill(mm_mixed,cs);
			if(event_selection_new&&Qsq3_1_cut)hcs_mixacc_result_new_Qsq3_1->Fill(mm_mixed,cs);
			if(event_selection_new&&Qsq3_2_cut)hcs_mixacc_result_new_Qsq3_2->Fill(mm_mixed,cs);
			if(event_selection_new&&Qsq3_3_cut)hcs_mixacc_result_new_Qsq3_3->Fill(mm_mixed,cs);
			if(event_selection_new&&phi1_1_cut)hcs_mixacc_result_new_phi1_1->Fill(mm_mixed,cs);
			if(event_selection_new&&phi2_1_cut)hcs_mixacc_result_new_phi2_1->Fill(mm_mixed,cs);
			if(event_selection_new&&phi2_2_cut)hcs_mixacc_result_new_phi2_2->Fill(mm_mixed,cs);
//cout<<"mm="<<mm<<endl;
//cout<<"mm_mixed="<<mm_mixed<<endl;
	}//Mix (j loop)
	if(event_selection)best++;
	if(event_selection_new)ac_new++;
	if(event_selection_nocut)woz++;
	if(event_selection_zdiff)zdiff++;
	if(event_selection_woac)woac++;

}//ENum (i loop)

cout<<"best="<< best <<endl;;
cout<<"ac_new="<< ac_new <<endl;;
cout<<"woz="<< woz <<endl;;
cout<<"zdiff="<< zdiff <<endl;;
cout<<"woac="<< woac <<endl;;

/*************************************/
/*Mixed Event Analysis was completed!*/
/*************************************/
//Following code is just for some checks using true-kaon event
//e.g.) filling MM histogram, comparing MM & BG(MEA), etc.
  tree->Draw(">>elist_kaon", "abs(ct_orig)<1.006");
  TEventList *elist_kaon = (TEventList*)gROOT->FindObject("elist_kaon");
  int ENum_kaon = elist_kaon->GetN(); 
cout<<"Entries(kaon): "<<ENum_kaon<<endl;
  for(int i=0;i<ENum_kaon;i++){
  tree->GetEntry(elist_kaon->GetEntry(i));

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


		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
		else event_selection_nocut=false;
//----MM Calc. again----//
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
	    double R_px = R_pz * (R_tr_tg_th );
	    double R_py = R_pz * ( R_tr_tg_ph );

	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(( L_tr_tg_th ), 2.0) + pow(( L_tr_tg_ph),2.0));
	    double L_px = L_pz * ( L_tr_tg_th );
	    double L_py = L_pz * ( L_tr_tg_ph );


	    double B_E =sqrt(B_mom*B_mom + Me*Me);
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
	    double L_E =sqrt(L_mom*L_mom + Me*Me);


		TLorentzVector L_4vec;
		TLorentzVector R_4vec;
		TLorentzVector B_4vec;
		TLorentzVector T_4vec;
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);

	    R_4vec.RotateX(  13.2/180.*PI );
	    L_4vec.RotateX( -13.2/180.*PI );

		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		double mass_true, mm_true;
		mass_true = Missing.M();
	    mm_true=mass_true - mh;//shift by ML
	
		if(event_selection)hmm_best->Fill(mm_true);
		if(event_selection_nocut)hmm_nocut->Fill(mm_true);

  }//ENum_kaon



cout << "creating new rootfile ... "<<endl;
file_new->Write();
//file_new->Close();

cout << "Well done!" << endl;
cout << "bgmea_temp.root was successfully created."<<endl;
cout << "Please change filename."<<endl;
}//mea_hist
