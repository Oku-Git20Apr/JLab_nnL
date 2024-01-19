//-----------------------------//
//--  Full Kinematics        --//
//-----------------------------//
//
//K. Okuyama (Feb. 21, 2023)
//from kinematics.C
//array --> no array
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
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNdivisions(505);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.00);
  h->GetYaxis()->SetTitleSize(0.05);
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
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.00);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void kinematics(){
  
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file (default: h2all1.root)
  TTree *tree = (TTree*)file->Get("tree_out");

  TFile *file_ = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/okuyamacro/simcR.root","read");
  //TH1F* h_phi_k_simc=(TH1F*)file_->Get("h_phi_k");

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
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;


//---------------------------------------//
//               Branch                  //
//---------------------------------------//

	int nrun;
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

//---------------------------------------//
//             Histogram                 //
//---------------------------------------//

//scattered electrons, kaons
  TH1D* h_theta_ee = new TH1D("h_theta_ee", "",200,5.,20.);
  TH1D* h_phi_ee = new TH1D("h_phi_ee", "",200,0.,180.);
  TH1D* h_theta_ek = new TH1D("h_theta_ek", "",200,5.,20.);
  TH1D* h_phi_ek = new TH1D("h_phi_ek", "",200,210.,330.);
  SetTH1(h_theta_ee, "", "#theta_{ee'} [deg]", "Counts", kAzure, kRed, kRed);
  SetTH1(h_phi_ee, "", "#phi_{ee'} [deg]", "Counts", kAzure, kRed, kRed);
  SetTH1(h_theta_ek, "", "#theta_{eK} [deg]", "Counts", kAzure, kRed, kRed);
  SetTH1(h_phi_ek, "", "#phi_{eK} [deg]", "Counts", kAzure, kRed, kRed);

//virtual photon
  TH1D* h_theta_g = new TH1D("h_theta_g", "",200,5.,20.);
  SetTH1(h_theta_g, "", "#theta_{e#gamma} [deg]", "Counts", kAzure, kRed, kRed);
  TH1D* h_phi_g = new TH1D("h_phi_g", "",200,210.,330.);
  SetTH1(h_phi_g, "", "#phi_{e#gamma} [deg]", "Counts", kAzure, kRed, kRed);
  TH1D* h_theta_gk_lab = new TH1D("h_theta_gk_lab", "",70,0.,7.);
  SetTH1(h_theta_gk_lab, "", "#theta_{#gamma K}^{lab} [deg]", "Counts", kAzure, kRed, kRed);
  TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "",180,0.,18.);
  SetTH1(h_theta_gk_cm, "", "#theta_{#gamma K}^{c.m.} [deg]", "Counts", kAzure, kRed, kRed);
  TH1D* h_cos_gk_lab = new TH1D("h_cos_gk_lab", "",100,0.97,1.0);
  SetTH1(h_cos_gk_lab, "", "cos(#theta_{#gamma K}^{lab})", "Counts", kAzure, kBlack, kBlack);
  TH1D* h_cos_gk_cm = new TH1D("h_cos_gk_cm", "",100,0.8,1.0);
  SetTH1(h_cos_gk_cm, "", "cos(#theta_{#gamma K}^{c.m.})", "Counts", kAzure, kBlack, kBlack);
  TH1D* h_mom_g = new TH1D("h_mom_g", "",100,2.1,2.5);
  SetTH1(h_mom_g, "", "P_{#gamma} [GeV/c]", "Counts", kAzure, kRed, kRed);
  TH1D* h_phi_gk = new TH1D("h_phi_gk","",100,-PI,PI);
  SetTH1(h_phi_gk, "", "#phi_{#gamma K} [deg]", "Counts", kAzure, kRed, kRed);
  TH1D* h_phi_gk_cos = new TH1D("h_phi_gk_cos","",100,-1.1,1.1);
  SetTH1(h_phi_gk_cos, "", "cos(#phi_{#gamma K})", "Counts", kAzure, kRed, kRed);
//added for the through inqury to phi_gk
  TH2D* h_phi_pk = new TH2D("h_phi_pk", "" ,72,0.,360.,60,1.70,1.95);
  h_phi_pk->SetStats(0);
  SetTH2(h_phi_pk, "", "#phi_{#gamma K} [deg]", "pK [GeV/c]", 0.4);
  TH2D* h_phi_pe = new TH2D("h_phi_pe", "" ,72,0.,360.,60,2.00,2.18);
  h_phi_pe->SetStats(0);
  SetTH2(h_phi_pe, "", "#phi_{#gamma K} [deg]", "pe [GeV/c]", 0.4);
  TH2D* h_phi_z = new TH2D("h_phi_z", "" ,72,0.,360.,60,-15.,15.);
  h_phi_z->SetStats(0);
  SetTH2(h_phi_z, "", "#phi_{#gamma K} [deg]", "Z [cm]", 0.4);
  TH2D* h_phi_ct = new TH2D("h_phi_ct", "" ,72,0.,360.,60,-6.,6.);
  h_phi_ct->SetStats(0);
  SetTH2(h_phi_ct, "", "#phi_{#gamma K} [deg]", "Coin. time [ns]", 0.4);

//Q2, W
  TH1D* h_qsq = new TH1D("h_qsq", "",150,0.2,0.8);
  SetTH1(h_qsq, "", "Q^{2} [(GeV/c)^{2}]", "Counts", kAzure, kRed, kRed);
  TH1D* h_w = new TH1D("h_w", "",100,2.05,2.25);
  SetTH1(h_w, "", "W [GeV]", "Counts", kAzure, kRed, kRed);
  TH1D* h_labtocm = new TH1D("h_labtocm", "",100,0.0,0.25);
  SetTH1(h_labtocm, "", "(d#sigma/d#Omega)_{c.m.}/(d#sigma/d#Omega)_{lab}", "Counts", kAzure, kRed, kRed);
  TH2D* h_qw = new TH2D("h_qw", "" ,20,2.05,2.25,60,0.2,0.8);
  h_qw->SetStats(0);
  SetTH2(h_qw, "", "Q^{2} [(GeV/c)^{2}]", "W [GeV]", 0.4);
  TH2D* h_thph_ee = new TH2D("h_thph_ee", "" ,200,5.,20.,200,0.,180.);
  h_thph_ee->SetStats(0);
  SetTH2(h_thph_ee, "", "#theta_{ee'} [deg]", "#phi_{ee'} [deg]", 0.4);
  TH2D* h_thph_ek = new TH2D("h_thph_ek", "" ,200,5.,20.,200,210.,330.);
  h_thph_ek->SetStats(0);
  SetTH2(h_thph_ek, "", "#theta_{eK} [deg]", "#phi_{eK} [deg]", 0.4);
  TH2D* h_thph_g = new TH2D("h_thph_g", "" ,200,5.,20.,200,210.,330.);
  h_thph_g->SetStats(0);
  SetTH2(h_thph_g, "", "#theta_{e#gamma} [deg]", "#phi_{e#gamma} [deg]", 0.4);
  TH1D* h_pR_lab = new TH1D("h_pR_lab", "" ,100,1.7,1.95);
  SetTH1(h_pR_lab, "", "p_{K}^{lab} [GeV/c]", "Counts", kAzure, kRed, kRed);
  TH1D* h_pR_cm = new TH1D("h_pR_cm", "" ,50,0.5,0.75);
  SetTH1(h_pR_cm, "", "p_{K}^{c.m.} [GeV/c]", "Counts", kAzure, kRed, kRed);
  TH2D* h2_pR_lab_cm = new TH2D("h_pR_lab_cm", "" ,1000,1.7,1.95,1000,0.0,1.0);
  h2_pR_lab_cm->SetStats(0);
  SetTH2(h2_pR_lab_cm, "", "p_{K}^{lab} [GeV/c]", "p_{K}^{c.m.} [GeV/c]", 0.4);

  TH3D *h3_uni = new TH3D("h3_uni","(#theta_{#gamma K},#phi_{#gamma K})",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  h3_uni->SetFillColor(kBlack);
  TH3D *h3_uni_cm = new TH3D("h3_uni_cm","(#theta_{#gamma K}^{c.m.},#phi_{#gamma K})",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  h3_uni_cm->SetFillColor(kBlack);
  TH3D *h3_gk = new TH3D("h3_gk","(#theta_{#gamma K},#phi_{#gamma K})",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  h3_gk->SetFillColor(kRed);
  h3_gk->SetLineColor(kRed);
  h3_gk->SetMarkerColor(kRed);
  TH3D *h3_gk_cm = new TH3D("h3_gk_cm","(#theta_{#gamma K}^{c.m.},#phi_{#gamma K})",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  h3_gk_cm->SetFillColor(kRed);
  h3_gk_cm->SetLineColor(kRed);
  h3_gk_cm->SetMarkerColor(kRed);
  TH2D* h2_gk = new TH2D("h2_gk", "" ,50,0.,18.,50,0.,360.);
  h2_gk->SetStats(0);
  h2_gk->GetXaxis()->SetTitle("#theta_{#gamma K}^{c.m.} [deg]");
  h2_gk->GetYaxis()->SetTitle("#phi_{#gamma K} [deg]");

  TH2D* h_phi_Q2 = new TH2D("h_phi_Q2", "" ,72,0.,360.,60,0.2,0.8);
  h_phi_Q2->SetStats(0);
  h_phi_Q2->GetXaxis()->SetTitle("#phi_{#gamma K} [deg]");
  h_phi_Q2->GetYaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
  bool Lambda = false;
  bool Sigma = false;
  bool mix_region1 = false;
  bool mix_region2 = false;
  bool mix_region3 = false;
  bool mix_region4 = false;
  bool mix_region5 = false;
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum

  //int ENum=0;
  //ENum = tree->GetEntries();
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


	


	if(fabs(ct)<1.006)ct_cut=true;
	else ct_cut=false;
    if(abs(ct+5.0*rf_bunch)<1.0) mix_region1 = true;
    else mix_region1 = false;
    if(abs(ct-1.0*rf_bunch)<1.0) mix_region2 = true;
    else mix_region2 = false;
    if(abs(ct-2.0*rf_bunch)<1.0) mix_region3 = true;
    else mix_region3 = false;
    if(abs(ct-7.0*rf_bunch)<1.0) mix_region4 = true;
    else mix_region4 = false;
    if(abs(ct+4.0*rf_bunch)<1.0) mix_region5 = true;
    else mix_region5 = false;
	if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)event_selection=true;
	//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
	else event_selection=false;
	if(L_mom>2.092&&L_mom<2.160)Lambda=true; else Lambda=false;
	if(L_mom>2.010&&L_mom<2.108)Sigma =true; else Sigma =false;

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




	    double pL    = L_mom;//GeV
	    double pR    = R_mom;//GeV
		double theta = L_tr_tg_th;
		double theta_R = R_tr_tg_th;
		double phi = L_tr_tg_ph;
		double phi_R = R_tr_tg_ph;
		double phi0=13.2*PI/180;//rad
		double phi_L = L_4vec.Phi();//LHRS frame
		double phi_RHRS = R_4vec.Phi();//RHRS frame
		//test double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );
        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML
		//without Matrix & Energy Loss calibration
		double theta_ee = L_4vec.Theta();
		//test double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double phi_L = atan((phi*cos(phi0)+sin(phi0))/theta);//LHRS frame
		double phi_ee = L_4vec.Phi();//original frame
//cout<<"phi_ee="<<phi_ee<<endl;
		double phi_ek = R_4vec.Phi()+2*PI;//original frame
//cout<<"phi_ek="<<phi_ek<<endl;

		G_4vec = B_4vec - L_4vec;
		//double mom_g = G_4vec.Rho();
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
		TVector3 G_3vec = G_4vec.Vect();
		TVector3 L_3vec = L_4vec.Vect();
		TVector3 R_3vec = R_4vec.Vect();
		//TVector3 G_3vec;
		//TVector3 L_3vec;
		//TVector3 R_3vec;
		//G_3vec.SetXYZ(G_4vec.Px(),G_4vec.Py(),G_4vec.Pz());
		//L_3vec.SetXYZ(L_4vec.Px(),L_4vec.Py(),L_4vec.Pz());
		//R_3vec.SetXYZ(R_4vec.Px(),R_4vec.Py(),R_4vec.Pz());
		TVector3 l_3vec = G_3vec.Cross(L_3vec);
		TVector3 r_3vec = G_3vec.Cross(R_3vec);
		TVector3 s_3vec = l_3vec.Cross(r_3vec);//for sgn(sign(phi_gk))
		TVector3 axis_3vec;
		axis_3vec.SetXYZ(0.,-1.,0.);//along VP flux
		double sgn = s_3vec*axis_3vec;
		double phi_gk_cos = (l_3vec*r_3vec)/l_3vec.Mag()/r_3vec.Mag();
		double phi_gk;
		if(sgn>0.){phi_gk = acos(phi_gk_cos);}
		else{phi_gk = 2*PI-acos(phi_gk_cos);}
	
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

//cout<<"theta_gk_cm(TLorentz)="<<theta_gk_cm*180./PI<<endl;
		double n = MK/ML;
		double p_cm=sqrt(GT_4vec.Px()*GT_4vec.Px()+GT_4vec.Py()*GT_4vec.Py()+GT_4vec.Pz()*GT_4vec.Pz());
		double E_cm = GT_4vec.E();
//beta=2.3/(2.2+Mp);
//pR_cm=0.65;
//theta_gk_cm=0.12;
		//double gamma=1./sqrt(1-beta*beta);
		double ER_cm=sqrt(pR_cm*pR_cm+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm = (gamma*pR_cm*pR_cm*(pR_cm*cos(theta_gk_cm)+beta*ER_cm))/(pow(sqrt(pR_cm*pR_cm*sin(theta_gk_cm)*sin(theta_gk_cm)+gamma*gamma*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)),3.));
//cout<<"labtocm="<<labtocm<<endl;
		double tan_lab1 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+beta*sqrt(MK*MK+pR_cm*pR_cm)/pR_cm));
		double tan_lab2 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+(omega*Mp-Qsq*Qsq)/(omega*Mp+Mp*Mp)));
		//if(tan_lab1!=tan_lab2)cout<<"tan1="<<atan(tan_lab1)<<", tan2="<<atan(tan_lab2)<<"theta_gk_lab="<<theta_gk_lab<<endl;

		double xx = sin(theta_gk_lab)*cos(phi_gk);
		double yy = sin(theta_gk_lab)*sin(phi_gk);
		double zz = cos(theta_gk_lab);
		double xx_cm = sin(theta_gk_cm)*cos(phi_gk);
		double yy_cm = sin(theta_gk_cm)*sin(phi_gk);
		double zz_cm = cos(theta_gk_cm);


		//if(event_selection&&Lambda){//Lambda selection
		//if(event_selection&&Sigma){//Sigma0 selection
		//if(event_selection){//All data (w/ mom. cut)
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16){
		h_theta_ee ->Fill(theta_ee*180./PI);
		h_phi_ee ->Fill(phi_ee*180./PI);
		h_theta_ek ->Fill(theta_ek*180./PI);
		h_phi_ek ->Fill(phi_ek*180./PI);
		h_theta_g ->Fill(theta_g*180./PI);
		h_phi_g ->Fill(phi_g*180./PI);
		h_thph_ee ->Fill(theta_ee*180./PI,phi_ee*180./PI);
		h_thph_ek->Fill(theta_ek*180./PI,phi_ek*180./PI);
		h_thph_g->Fill(theta_g*180./PI,phi_g*180./PI);
		h_mom_g->Fill(mom_g);
		h_qsq->Fill(Qsq);
		h_w->Fill(W);
		h_labtocm->Fill(labtocm);
		h_qw->Fill(W,Qsq);
		h_theta_gk_lab->Fill(theta_gk_lab*180./PI);
		h_theta_gk_cm->Fill(theta_gk_cm*180./PI);
		h_cos_gk_lab->Fill(cos(theta_gk_lab));
		h_cos_gk_cm->Fill(cos(theta_gk_cm));
		h_pR_lab->Fill(pR);
		h_pR_cm->Fill(pR_cm);
		h2_pR_lab_cm->Fill(pR,pR_cm);
		h_phi_gk->Fill(phi_gk-PI);
		h_phi_gk_cos->Fill(phi_gk_cos);
		h_phi_pk->Fill(phi_gk*180./PI,pR);
		h_phi_pe->Fill(phi_gk*180./PI,pL);
		h_phi_z->Fill(phi_gk*180./PI,(R_tr_vz+L_tr_vz)*100./2.);
		h_phi_ct->Fill(phi_gk*180./PI,ct);
		h3_gk->Fill(xx,yy,zz);
		h3_gk_cm->Fill(xx_cm,yy_cm,zz_cm);
		h2_gk->Fill(theta_gk_cm*180./PI,phi_gk*180./PI);
		h_phi_Q2->Fill(phi_gk*180./PI,Qsq);
		}

}//ENum

//	TCanvas* c1 = new TCanvas("c1","c1");
//	c1->Divide(2,2);
//	c1->cd(1);
//	h_theta_ee->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c1->cd(2);
//	h_phi_ee->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c1->cd(3);
//	h_thph_ee->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c1->Modified();
//	c1->Update();
//	gPad->Modified();
//	gPad->Update();
//	TCanvas* c2 = new TCanvas("c2","c2");
//	c2->Divide(2,2);
//	c2->cd(1);
//	h_theta_ek->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c2->cd(2);
//	h_phi_ek->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c2->cd(3);
//	h_thph_ek->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c3 = new TCanvas("c3","c3");
//	c3->Divide(2,2);
//	c3->cd(1);
//	h_theta_g->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c3->cd(2);
//	h_phi_g->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c3->cd(3);
//	h_thph_g->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c4 = new TCanvas("c4","c4");
//	h_mom_g->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c5 = new TCanvas("c5","c5");
//	h_qsq->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c6 = new TCanvas("c6","c6");
//	c6->Divide(2,2);
//	c6->cd(1);
//	h_theta_gk_lab->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c6->cd(2);
//	h_theta_gk_cm->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c6->cd(3);
//	h_cos_gk_lab->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c6->cd(4);
//	h_cos_gk_cm->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c7 = new TCanvas("c7","c7");
//	c7->Divide(2,2);
//	c7->cd(1);
//	h_pR_lab->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c7->cd(2);
//	h_pR_cm->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c7->cd(3);
//	h_phi_gk->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	c7->cd(4);
//	h_phi_gk_cos->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c8 = new TCanvas("c8","c8");
//	h_w->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c9 = new TCanvas("c9","c9");
//	h_qw->Draw("colz");
//	c9->SetLeftMargin(0.14);
//	c9->SetRightMargin(0.14);
//	c9->SetTopMargin(0.14);
//	c9->SetBottomMargin(0.14);
//	c9->Modified();
//	c9->Update();
//	gPad->Modified();
//	gPad->Update();
//	TCanvas* c10 = new TCanvas("c10","c10");
//	h_labtocm->Draw("");
//	c10->SetLeftMargin(0.14);
//	c10->SetRightMargin(0.14);
//	c10->SetTopMargin(0.14);
//	c10->SetBottomMargin(0.14);
//	c10->Modified();
//	c10->Update();
//	gPad->Modified();
//	gPad->Update();
//	TCanvas* c11 = new TCanvas("c11","c11");
//	for(int i=0;i<10000;i++){
//		double theta_gen = acos(gRandom->Uniform(-1.,1.));
//		double phi_gen = 2.*PI*gRandom->Uniform(0.,1.);
//		double x_uni = sin(theta_gen)*cos(phi_gen);
//		double y_uni = sin(theta_gen)*sin(phi_gen);
//		double z_uni = cos(theta_gen);
//		h3_uni->Fill(x_uni,y_uni,z_uni);
//		h3_uni_cm->Fill(x_uni,y_uni,z_uni);
//	}
//	h3_uni->SetStats(0);
//	h3_gk->SetStats(0);
//	h3_uni->Draw("");
//	h3_gk->Draw("same");
//	TCanvas* c12 = new TCanvas("c12","c12");
//	h3_uni_cm->SetStats(0);
//	h3_gk_cm->SetStats(0);
//	h3_uni_cm->Draw("");
//	h3_gk_cm->Draw("same");
//	TCanvas* c13 = new TCanvas("c13","c13");
//	h2_gk->Draw("colz");
//	c13->SetLeftMargin(0.14);
//	c13->SetRightMargin(0.14);
//	c13->SetTopMargin(0.14);
//	c13->SetBottomMargin(0.14);
//	c13->Modified();
//	c13->Update();
//	gPad->Modified();
//	gPad->Update();
//
////theta,phi -- ee', eK, eg
//	TCanvas* c51 = new TCanvas("c51","c51");
//	h_theta_ee->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c52 = new TCanvas("c52","c52");
//	h_phi_ee->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c53 = new TCanvas("c53","c53");
//	h_thph_ee->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c54 = new TCanvas("c54","c54");
//	h_theta_ek->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c55 = new TCanvas("c55","c55");
//	h_phi_ek->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c56 = new TCanvas("c56","c56");
//	h_thph_ek->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c57 = new TCanvas("c57","c57");
//	h_theta_g->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c58 = new TCanvas("c58","c58");
//	h_phi_g->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c59 = new TCanvas("c59","c59");
//	h_thph_g->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//
////theta,phi -- gK -- lab,c.m.
//	TCanvas* c61 = new TCanvas("c61","c61");
//	h_theta_gk_lab->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c62 = new TCanvas("c62","c62");
//	h_theta_gk_cm->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c63 = new TCanvas("c63","c63");
//	h_phi_gk->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//
////pK -- lab,c.m.
//	TCanvas* c71 = new TCanvas("c71","c71");
//	h_pR_lab->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c72 = new TCanvas("c72","c72");
//	h_pR_cm->Draw("");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//	TCanvas* c73 = new TCanvas("c73","c73");
//	h2_pR_lab_cm->Draw("colz");
//	gPad->SetLeftMargin(0.14);
//	gPad->SetRightMargin(0.14);
//	gPad->SetTopMargin(0.14);
//	gPad->SetBottomMargin(0.14);
//
///*--- Print ---*/
//	string pdfname = "./pdf/kinematics_wMom.pdf";
//cout << "Output pdf file name is " << pdfname << endl;
//cout << "Print is starting" << endl;
//	c1->Print(Form("%s[",pdfname.c_str()));
//	c1->Print(Form("%s",pdfname.c_str()));
//	c2->Print(Form("%s",pdfname.c_str()));
//	c3->Print(Form("%s",pdfname.c_str()));
//	c4->Print(Form("%s",pdfname.c_str()));
//	c5->Print(Form("%s",pdfname.c_str()));
//	c6->Print(Form("%s",pdfname.c_str()));
//	c7->Print(Form("%s",pdfname.c_str()));
//	c8->Print(Form("%s",pdfname.c_str()));
//	c9->Print(Form("%s",pdfname.c_str()));
//	c10->Print(Form("%s",pdfname.c_str()));
//	c11->Print(Form("%s",pdfname.c_str()));
//	c12->Print(Form("%s",pdfname.c_str()));
//	c13->Print(Form("%s",pdfname.c_str()));
//	c13->Print(Form("%s]",pdfname.c_str()));
//
//	c5->Print("./pdf/Qsq_wMom.pdf");
//	c6->Print("./pdf/theta_gk_wMom.pdf");
//	c7->Print("./pdf/CM_wMom.pdf");
//	c8->Print("./pdf/W_wMom.pdf");
//	c9->Print("./pdf/QsqW_wMom.pdf");
//	c11->Print("./pdf/theta_lab_ph_3d_wMom.pdf");
//	c12->Print("./pdf/theta_cm_ph_3d_wMom.pdf");
//	c13->Print("./pdf/theta_cm_ph_2d_wMom.pdf");
//
////theta,phi -- ee', eK, eg
//	c51->Print("./pdf/theta_ee_wMom.pdf");
//	c52->Print("./pdf/phi_ee_wMom.pdf");
//	c53->Print("./pdf/theta_phi_ee_wMom.pdf");
//	c54->Print("./pdf/theta_ek_wMom.pdf");
//	c55->Print("./pdf/phi_ek_wMom.pdf");
//	c56->Print("./pdf/theta_phi_ek_wMom.pdf");
//	c57->Print("./pdf/theta_eg_wMom.pdf");
//	c58->Print("./pdf/phi_eg_wMom.pdf");
//	c59->Print("./pdf/theta_phi_eg_wMom.pdf");
//
////theta,phi -- gK -- lab,c.m.
//	c61->Print("./pdf/theta_gk_lab_wMom.pdf");
//	c62->Print("./pdf/theta_gk_cm_wMom.pdf");
//	c63->Print("./pdf/phi_gk_wMom.pdf");
////pK -- lab,c.m.
//	c71->Print("./pdf/pK_lab_wMom.pdf");
//	c72->Print("./pdf/pK_cm_wMom.pdf");
//	c73->Print("./pdf/pK_lab_cm_wMom.pdf");
//	c10->Print("./pdf/labtocm_wMom.pdf");
//cout << "Well done!" << endl;
	TCanvas* c63 = new TCanvas("c63","c63");
	h_phi_gk->Draw("");
	gPad->SetLeftMargin(0.14);
	gPad->SetRightMargin(0.14);
	gPad->SetTopMargin(0.14);
	gPad->SetBottomMargin(0.14);
	TCanvas* c64 = new TCanvas("c64","c64");
	c64->Divide(2,2);
	c64->cd(1);h_phi_pk->Draw("colz");
	gPad->SetLeftMargin(0.14);
	gPad->SetRightMargin(0.14);
	gPad->SetTopMargin(0.14);
	gPad->SetBottomMargin(0.14);
	c64->cd(2);h_phi_pe->Draw("colz");
	gPad->SetLeftMargin(0.14);
	gPad->SetRightMargin(0.14);
	gPad->SetTopMargin(0.14);
	gPad->SetBottomMargin(0.14);
	c64->cd(3);h_phi_z->Draw("colz");
	gPad->SetLeftMargin(0.14);
	gPad->SetRightMargin(0.14);
	gPad->SetTopMargin(0.14);
	gPad->SetBottomMargin(0.14);
	c64->cd(4);h_phi_ct->Draw("colz");
	gPad->SetLeftMargin(0.14);
	gPad->SetRightMargin(0.14);
	gPad->SetTopMargin(0.14);
	gPad->SetBottomMargin(0.14);

//	TCanvas* c65 = new TCanvas("c65","c65");
//	double val1,val2;
//	TH1D* h_phi_gk_cs = new TH1D("h_phi_gk_cs", "",100,0,2.*PI);
//	for(int i=0;i<100;i++){
//		val1 = h_phi_gk->GetBinContent(i+1);
//		val2 = h_phi_k_simc->GetBinContent(i+1);
//		val1 += PI;
//		val2 += PI;
//		if(val1!=0.&&val2!=0.){
//			h_phi_gk_cs->SetBinContent(i+1,1250.*val1/val2);
//			h_phi_gk_cs->SetBinError(i+1,1250.*sqrt(val1/val2/val2+val1*val1/val2/val2/val2));
//		}else{h_phi_gk_cs->SetBinContent(i+1,0.);}
//	}
//	h_phi_gk_cs->Draw("e");

	TCanvas* c66 = new TCanvas("c66","c66");
	h_phi_Q2->Draw("colz");
	c66->SetLeftMargin(0.14);
	c66->SetRightMargin(0.14);
	c66->SetTopMargin(0.14);
	c66->SetBottomMargin(0.14);
	c66->Modified();
	c66->Update();
	gPad->Modified();
	gPad->Update();
	c66->Print("./pdf/phi_Q2.pdf");
}//kinematics
