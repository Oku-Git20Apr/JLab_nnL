//-- FP & Chi2 from small root  --//
//
//K. Okuyama (Feb. 17, 2023)

void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->SetTitleFont(42,"");
  h->SetTitleSize(0.04,"");

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNdivisions(505);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}
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
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void FP(){

  TFile *file = new TFile("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111500.root","read");//input file of all H2 run(default: h2all4.root)
  TChain *chain = new TChain("T");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111500.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111500_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111501.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111501_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111501_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111501_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111502.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111502_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111502_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111502_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111503.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111503_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111503_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111503_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111504.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111504_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111504_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111505.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111505_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111505_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111505_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111506.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111506_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111506_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111506_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111507.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111507_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111507_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111507_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111508.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111508_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111508_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111508_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111509.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111509_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111509_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111509_3.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111510.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111510_1.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111510_2.root");
  chain->Add("/data/41a/ELS/okuyama/rootfiles/nnL_small/tritium_111510_3.root");



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

  TH2F* h_xth_L  = new TH2F("h_xth_L","LHRS: FP X vs FP #theta",100,-0.2,0.2,100,-0.8,0.8);
  TH2F* h_xth_L_cut  = new TH2F("h_xth_L_cut","LHRS: FP X vs FP #theta (w/ cut)",100,-0.2,0.2,100,-0.8,0.8);
  TH2F* h_xth_R  = new TH2F("h_xth_R","RHRS: FP X vs FP #theta",100,-0.2,0.2,100,-0.8,0.8);
  TH2F* h_xth_R_cut  = new TH2F("h_xth_R_cut","RHRS: FP X vs FP #theta (w/ cut)",100,-0.2,0.2,100,-0.8,0.8);
  TH1F* h_chi2_L  = new TH1F("h_chi2_L","LHRS: #chi^{2}",1000,0.,0.01);
  TH1F* h_chi2_R  = new TH1F("h_chi2_R","RHRS: #chi^{2}",1000,0.,0.01);
  SetTH2(h_xth_L, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  SetTH2(h_xth_R, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  SetTH2(h_xth_L_cut, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  SetTH2(h_xth_R_cut, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  SetTH1(h_chi2_L, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);
  SetTH1(h_chi2_R, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);



	gStyle->SetPalette(1);
	int ENum = chain->GetEntries();
cout<<"Entries: "<<ENum<<endl;
	chain->Project("h_xth_L","L.tr.x:L.tr.th","","");
	chain->Project("h_xth_L_cut","L.tr.x:L.tr.th","L.tr.th<0.17*L.tr.x+0.025 && L.tr.th>0.17*L.tr.x-0.035 && L.tr.th<0.40*L.tr.x+0.130","");
	chain->Project("h_xth_R","R.tr.x:R.tr.th","","");
	chain->Project("h_xth_R_cut","R.tr.x:R.tr.th","R.tr.th<0.17*R.tr.x+0.025 && R.tr.th>0.17*R.tr.x-0.035 && R.tr.th<0.40*R.tr.x+0.130","");
	chain->Project("h_chi2_L","L.tr.chi2","","");
	chain->Project("h_chi2_R","R.tr.chi2","","");
	cout<<"ENum in TChain="<<chain->GetEntries()<<endl;
	cout<<"ENum FP_L (w/  cut)="<<h_xth_L_cut->Integral()<<endl;
	cout<<"ENum FP_L (w/o cut)="<<h_xth_L->GetEntries()<<endl;
	cout<<"Eff. (FP_L) = "<<h_xth_L_cut->Integral()/h_xth_L->GetEntries()<<endl;
	cout<<"ENum FP_R (w/  cut)="<<h_xth_R_cut->Integral()<<endl;
	cout<<"ENum FP_R (w/o cut)="<<h_xth_R->GetEntries()<<endl;
	cout<<"Eff. (FP_R) = "<<h_xth_R_cut->Integral()/h_xth_R->GetEntries()<<endl;
	cout<<"ENum Chi2_L (w/  cut)="<<h_chi2_L->Integral()<<endl;
	cout<<"ENum Chi2_L (w/o cut)="<<h_chi2_L->GetEntries()<<endl;
	cout<<"Eff. (Chi2_L) = "<<h_chi2_L->Integral()/h_chi2_L->GetEntries()<<endl;
	cout<<"ENum Chi2_R (w/  cut)="<<h_chi2_R->Integral()<<endl;
	cout<<"ENum Chi2_R (w/o cut)="<<h_chi2_R->GetEntries()<<endl;
	cout<<"Eff. (Chi2_R) = "<<h_chi2_R->Integral()/h_chi2_R->GetEntries()<<endl;
	TCanvas* c1 = new TCanvas("c1","c1", 800., 800.);
	h_xth_L->Draw("colz");
	c1->SetLeftMargin(0.14);
	c1->SetRightMargin(0.14);
	c1->SetTopMargin(0.14);
	c1->SetBottomMargin(0.14);
	c1->Modified();
	c1->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c2 = new TCanvas("c2","c2", 800., 800.);
	h_xth_L_cut->Draw("colz");
	c2->SetLeftMargin(0.14);
	c2->SetRightMargin(0.14);
	c2->SetTopMargin(0.14);
	c2->SetBottomMargin(0.14);
	c2->Modified();
	c2->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c3 = new TCanvas("c3","c3", 800., 800.);
	h_xth_R->Draw("colz");
	c3->SetLeftMargin(0.14);
	c3->SetRightMargin(0.14);
	c3->SetTopMargin(0.14);
	c3->SetBottomMargin(0.14);
	c3->Modified();
	c3->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c4 = new TCanvas("c4","c4", 800., 800.);
	h_xth_R_cut->Draw("colz");
	c4->SetLeftMargin(0.14);
	c4->SetRightMargin(0.14);
	c4->SetTopMargin(0.14);
	c4->SetBottomMargin(0.14);
	c4->Modified();
	c4->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c5 = new TCanvas("c5","c5", 800., 800.);
	h_chi2_L->Draw("");
	c5->SetLeftMargin(0.14);
	c5->SetRightMargin(0.14);
	c5->SetTopMargin(0.14);
	c5->SetBottomMargin(0.14);
	c5->Modified();
	c5->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c6 = new TCanvas("c6","c6", 800., 800.);
	h_chi2_R->Draw("");
	c6->SetLeftMargin(0.14);
	c6->SetRightMargin(0.14);
	c6->SetTopMargin(0.14);
	c6->SetBottomMargin(0.14);
	c6->Modified();
	c6->Update();
	gPad->Modified();
	gPad->Update();

//c1->Print("./pdf/LFP_xth.pdf");
//c2->Print("./pdf/LFP_xth_cut.pdf");
//c3->Print("./pdf/RFP_xth.pdf");
//c4->Print("./pdf/RFP_xth_cut.pdf");
//c5->Print("./pdf/Lchi2.pdf");
//c6->Print("./pdf/Rchi2.pdf");

cout << "Well done!" << endl;
}//fit
