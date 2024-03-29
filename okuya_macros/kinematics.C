//-----------------------------//
//--  Full Kinematics        --//
//-----------------------------//
//
//K. Okuyama (Aug. 13, 2020)

#define MAX 100     // Maximum No. of Tracks
#define RS0 1      // No. of Segments of R-S0
#define RS2 16     // No. of Segments of R-S2
#define RA1 24     // No. of Segments of R-AC1
#define RA2 26     // No. of Segments of R-AC2
#define RCR 10     // No. of Segments of R-GC
#define LCL 10     // No. of Segments of L-GC
#define RPS 48     // No. of Segments of R-Pre-Shower
#define RSH 75     // No. of Segments of R-Shower
#define RF1TDC 64  // No. of ch of R-F1TDC
#define LS0 1      // No. of Segments of L-S0
#define LS2 16     // No. of Segments of L-S2
#define LF1TDC 64  // No. of ch of L-F1TDC


void kinematics(){
	string pdfname = "kinematics.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("h2all.root","read");//input file (default: h2all1.root)
  TFile *file_new = new TFile("kinematics.root","recreate");//new root
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
//	tree->Write();
    



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

	double L_tr_chi2[MAX];
	double L_tr_x[MAX], L_tr_y[MAX], L_tr_th[MAX], L_tr_ph[MAX];
	double L_tr_p[MAX];
	double L_tr_tg_th[MAX], L_tr_tg_ph[MAX];
	double L_tr_vz[MAX];
	double L_tr_vz_saved[MAX];
	double R_tr_chi2[MAX];
	double R_tr_x[MAX], R_tr_y[MAX], R_tr_th[MAX], R_tr_ph[MAX];
	double R_tr_p[MAX];
	double R_tr_tg_th[MAX], R_tr_tg_ph[MAX];
	double R_tr_vz[MAX];
	double L_mom[MAX], R_mom[MAX], B_mom; 
	double L_ene[MAX], R_ene[MAX], B_ene; 
	double ac1sum, ac2sum;//NPE SUM
	
//	string branchname[]={"tr.ntrack_l","tr.ntrack_r","tr.Ls2_pad[100]"};
//	Int_t nbranch = sizeof(branchname)/sizeof(branchname[0]);
//	tree->SetBranchStatus("*",0);
//	for(Int_t ibranch;ibranch<nbranch;ibranch++){
//	tree->SetBranchStatus(branchname[ibranch].c_str(),1);
//	}

	//tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("tr.ntrack_l",1);tree->SetBranchAddress("tr.ntrack_l",&NLtr);
	tree->SetBranchStatus("tr.ntrack_r",1);tree->SetBranchAddress("tr.ntrack_r",&NRtr);
	
  	tree->SetBranchStatus("ac1_npe_sum",1);  tree->SetBranchAddress("ac1_npe_sum", &ac1sum);
  	tree->SetBranchStatus("ac2_npe_sum",1);  tree->SetBranchAddress("ac2_npe_sum", &ac2sum);
  	tree->SetBranchStatus("Lp_c",1);  tree->SetBranchAddress("Lp_c", L_mom);
  	tree->SetBranchStatus("Rp_c",1);  tree->SetBranchAddress("Rp_c", R_mom);
  	tree->SetBranchStatus("Bp_c",1);  tree->SetBranchAddress("Bp_c", &B_mom);
  	tree->SetBranchStatus("ct_orig",1);  tree->SetBranchAddress("ct_orig", &ct);

//ADD 2020/8/12
  	tree->SetBranchStatus("L.tr.chi2",1);  tree->SetBranchAddress("L.tr.chi2", L_tr_chi2);
  	tree->SetBranchStatus("L.tr.x",1);  tree->SetBranchAddress("L.tr.x", L_tr_x);
  	tree->SetBranchStatus("L.tr.y",1);  tree->SetBranchAddress("L.tr.y", L_tr_y);
  	tree->SetBranchStatus("L.tr.th",1);  tree->SetBranchAddress("L.tr.th", L_tr_th);
  	tree->SetBranchStatus("L.tr.ph",1);  tree->SetBranchAddress("L.tr.ph", L_tr_ph);
  	tree->SetBranchStatus("L.tr.p",1);  tree->SetBranchAddress("L.tr.p", L_tr_p);
  	tree->SetBranchStatus("L.tr.tg_th",1);  tree->SetBranchAddress("L.tr.tg_th", L_tr_tg_th );
  	tree->SetBranchStatus("L.tr.tg_ph",1);  tree->SetBranchAddress("L.tr.tg_ph", L_tr_tg_ph );
  	tree->SetBranchStatus("L.tr.vz",1);  tree->SetBranchAddress("L.tr.vz", &L_tr_vz);

  	tree->SetBranchStatus("R.tr.chi2",1);  tree->SetBranchAddress("R.tr.chi2", R_tr_chi2);
	tree->SetBranchStatus("R.tr.x" ,1);  tree->SetBranchAddress("R.tr.x" , R_tr_x );
  	tree->SetBranchStatus("R.tr.y" ,1);  tree->SetBranchAddress("R.tr.y" , R_tr_y );
  	tree->SetBranchStatus("R.tr.th",1);  tree->SetBranchAddress("R.tr.th", R_tr_th);
  	tree->SetBranchStatus("R.tr.ph",1);  tree->SetBranchAddress("R.tr.ph", R_tr_ph);
  	tree->SetBranchStatus("R.tr.p",1);  tree->SetBranchAddress("R.tr.p", R_tr_p);
  	tree->SetBranchStatus("R.tr.tg_th",1);  tree->SetBranchAddress("R.tr.tg_th", R_tr_tg_th);
  	tree->SetBranchStatus("R.tr.tg_ph",1);  tree->SetBranchAddress("R.tr.tg_ph", R_tr_tg_ph);
  	tree->SetBranchStatus("R.tr.vz",1);  tree->SetBranchAddress("R.tr.vz", R_tr_vz);




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
  TH1F* hm1  = new TH1F("hm1","",xbin,xmin,xmax);
  hm1->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hm1->GetYaxis()->SetTitle("Counts / MeV");
  hm1->SetLineColor(1);
  TH1F* hm2   = (TH1F*)hm1->Clone("hm2");
  TH1F* hmm_acc  = new TH1F("hmm_acc","ACC (original)",xbin,xmin,xmax);
  TH1F* hmm_mixacc  = new TH1F("hmm_mixacc","ACC (mixed)",xbin,xmin,xmax);
  TH1F* hm4   = (TH1F*)hm1->Clone("hm4");

  TH2D* h_theta_eeek = new TH2D("h_theta_eeek", "theta_eeek",1000,0.1,0.35);
  TH1D* h_theta_ee = new TH1D("h_theta_ee", "theta_ee",1000,0.1,0.35);
  h_theta_ee->SetLineColor(kAzure);
  TH1D* h_phi_ee = new TH1D("h_phi_ee", "phi_ee",1000,0.,PI);
  TH1D* h_theta_ek = new TH1D("h_theta_ek", "theta_ek",1000,0.1,0.35);
  TH1D* h_phi_ek = new TH1D("h_phi_ek", "phi_ek",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_g = new TH1D("h_theta_g", "theta_g",1000,0.1,0.35);
  TH1D* h_phi_g = new TH1D("h_phi_g", "phi_g",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_gk_lab = new TH1D("h_theta_gk_lab", "#theta_{#gamma K}^{lab}",70,0.,7.);
  h_theta_gk_lab->GetXaxis()->SetTitle("#theta_{#gamma K}^{lab} [deg]");
  h_theta_gk_lab->GetYaxis()->SetTitle("Counts");
  h_theta_gk_lab->SetLineColor(kAzure);
  TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "#theta_{#gamma K}^{CM}",180,0.,18.);
  h_theta_gk_cm->GetXaxis()->SetTitle("#theta_{#gamma K}^{CM} [deg]");
  h_theta_gk_cm->GetYaxis()->SetTitle("Counts");
  h_theta_gk_cm->SetLineColor(kAzure);
  TH1D* h_cos_gk_lab = new TH1D("h_cos_gk_lab", "cos(#theta_{#gamma K}^{lab})",100,0.97,1.0);
  h_cos_gk_lab->GetXaxis()->SetTitle("cos(#theta_{#gamma K}^{lab})");
  h_cos_gk_lab->GetYaxis()->SetTitle("Counts");
  h_cos_gk_lab->SetLineColor(kAzure);
  TH1D* h_cos_gk_cm = new TH1D("h_cos_gk_cm", "cos(#theta_{#gamma K}^{CM})",100,0.8,1.0);
  h_cos_gk_cm->GetXaxis()->SetTitle("cos(#theta_{#gamma K}^{CM})");
  h_cos_gk_cm->GetYaxis()->SetTitle("Counts");
  h_cos_gk_cm->SetLineColor(kAzure);
  TH1D* h_mom_g = new TH1D("h_mom_g", "mom_g",1000,2.1,2.5);
  h_mom_g->SetLineColor(kAzure);
  TH1D* h_qsq = new TH1D("h_qsq", "Q^{2}",150,0.2,0.8);
  h_qsq->GetXaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");
  h_qsq->GetYaxis()->SetTitle("Counts");
  h_qsq->SetLineColor(kAzure);
  TH1D* h_w = new TH1D("h_w", "W",100,2.05,2.25);
  h_w->GetXaxis()->SetTitle("W [GeV]");
  h_w->GetYaxis()->SetTitle("Counts");
  h_w->SetLineColor(kAzure);
  TH1D* h_labtocm = new TH1D("h_labtocm", "labtocm",100,0.0,0.25);
  h_labtocm->GetXaxis()->SetTitle("(d#sigma/d#Omega)_{CM}/(d#sigma/d#Omega)_{lab}");
  h_labtocm->GetYaxis()->SetTitle("Counts");
  h_labtocm->SetLineColor(kAzure);
  TH2D* h_z_yfp = new TH2D("h_z_yfp", "z vs y(FP)" ,60,-0.06,0.06,60,-0.15,0.15);
  TH2D* h_z_phi = new TH2D("h_z_phi", "z vs phi" ,60,-0.06,0.06,60,-0.15,0.15);
  TH2D* h_yfp_phi = new TH2D("h_yfp_phi", "y(FP) vs phi" ,60,-0.06,0.06,60,-0.06,0.06);
  TH2D* h_qw = new TH2D("h_qw", "Q^{2}:W" ,20,2.05,2.25,60,0.2,0.8);
  h_qw->GetXaxis()->SetTitle("W [GeV]");
  h_qw->GetYaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]");
  TH2D* h_thph_ee = new TH2D("h_thph_ee", "theta_ee:phi_ee" ,1000,0.1,0.35,1000,PI/2-1.,PI/2+1.);
  TH2D* h_thph_ek = new TH2D("h_thph_ek", "theta_ek:phi_ek" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
  TH2D* h_thph_g = new TH2D("h_thph_g", "theta_g:phi_g" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_pR_lab = new TH1D("h_pR_lab", "p_{K}^{lab}" ,100,1.7,1.95);
  h_pR_lab->GetXaxis()->SetTitle("p_{K}^{lab} [GeV/c]");
  h_pR_lab->GetYaxis()->SetTitle("Counts");
  h_pR_lab->SetLineColor(kAzure);
  TH1D* h_pR_cm = new TH1D("h_pR_cm", "p_{K}^{CM}" ,50,0.5,0.75);
  h_pR_cm->GetXaxis()->SetTitle("p_{K}^{CM} [GeV/c]");
  h_pR_cm->GetYaxis()->SetTitle("Counts");
  h_pR_cm->SetLineColor(kAzure);
  TH1D* h_pL_cm = new TH1D("h_pL_cm", "h_pL_cm" ,1000,0.8,1.2);
  TH1D* h_pB_cm = new TH1D("h_pB_cm", "h_pB_cm" ,1000,1.6,2.0);
  TH2D* h2_pR_lab_cm = new TH2D("h_pR_lab_cm", "h_pR_lab_cm" ,1000,1.7,1.95,1000,0.0,1.0);
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);
  h2->SetLineColor(1);
  h3->SetLineColor(1);
  h3->SetFillColor(1);
  h3->SetFillStyle(3001);
  h4->SetLineColor(9);

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
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

	    //==============================//
	    //======  Initialization  ======//
	    //==============================//
		for(int j=0;j<MAX;j++){
    	}

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

      //int NLtr = (int)L_tr_n;  if(NLtr>MAX) NLtr = MAX;
      //int NRtr = (int)R_tr_n;  if(NRtr>MAX) NRtr = MAX;
      
	  NLtr = 1;
	  NRtr = 1;
      for(int lt=0;lt<NLtr;lt++){
        L_Tr = L_FP = false;
        if( L_tr_chi2[lt]<0.01 ) L_Tr = true;
        if( L_tr_th[lt]<0.17*L_tr_x[lt]+0.025
         && L_tr_th[lt]>0.17*L_tr_x[lt]-0.035
         && L_tr_th[lt]<0.40*L_tr_x[lt]+0.130 ) L_FP = true;
	
        for(int rt=0;rt<NRtr;rt++){
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2[rt]<0.01 ) R_Tr = true;
        if( R_tr_th[rt]<0.17*R_tr_x[rt]+0.025
         && R_tr_th[rt]>0.17*R_tr_x[rt]-0.035
         && R_tr_th[rt]<0.40*R_tr_x[rt]+0.130 ) R_FP = true;


	


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
	if(fabs(L_tr_vz[lt]-R_tr_vz[rt])<0.025&&fabs(R_tr_vz[rt]+L_tr_vz[lt])<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
	//if(fabs(L_tr_vz[lt]-R_tr_vz[rt])<0.025&&fabs(R_tr_vz[rt]+L_tr_vz[lt])<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
	else event_selection=false;

	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama

		
	    double R_pz = R_mom[lt]/sqrt(1.0*1.0 + pow((R_tr_tg_th[rt]), 2.0) + pow(( R_tr_tg_ph[rt]),2.0) );
	    double R_px = R_pz * (R_tr_tg_th[rt] );
	    double R_py = R_pz * ( R_tr_tg_ph[rt] );

	    double L_pz = L_mom[lt]/sqrt(1.0*1.0 + pow(( L_tr_tg_th[lt] ), 2.0) + pow(( L_tr_tg_ph[lt]),2.0));
	    double L_px = L_pz * ( L_tr_tg_th[lt] );
	    double L_py = L_pz * ( L_tr_tg_ph[lt] );


	    double B_E =sqrt(B_mom*B_mom + Me*Me);
	    double R_E =sqrt(R_mom[lt]*R_mom[lt] + MK*MK);
	    double L_E =sqrt(L_mom[lt]*L_mom[lt] + Me*Me);

		h_test->Fill(L_E);

		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);




	    double pL    = L_tr_p[lt];//GeV
	    double pR    = R_tr_p[rt];//GeV
		double theta = L_tr_tg_th[lt];
		double theta_R = R_tr_tg_th[rt];
		double phi = L_tr_tg_ph[lt];
		double phi_R = R_tr_tg_ph[rt];
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
		if(event_selection&&ct_cut)hm1->Fill(mm);
		if(mix_region1||mix_region2||mix_region3||mix_region4||mix_region5)hmm_acc->Fill(mm);
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



		if(event_selection){
		h_theta_ee ->Fill(theta_ee);
		h_phi_ee ->Fill(phi_ee);
		h_theta_ek ->Fill(theta_ek);
		h_phi_ek ->Fill(phi_ek);
		h_theta_g ->Fill(theta_g);
		h_phi_g ->Fill(phi_g);
		h_thph_ee ->Fill(theta_ee,phi_ee);
		h_thph_ek->Fill(theta_ek,phi_ek);
		h_thph_g->Fill(theta_g,phi_g);
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
		h_pL_cm->Fill(pL_cm);
		h_pB_cm->Fill(pB_cm);
		}
		if(fabs(L_tr_vz[lt]-R_tr_vz[rt])<0.025&&fabs(L_tr_vz[lt]+R_tr_vz[rt])<0.2){
		h_z_yfp->Fill(R_tr_y[rt],R_tr_vz[rt]); 
  		h_z_phi->Fill(R_tr_ph[rt],R_tr_vz[rt]); 
		h_yfp_phi->Fill(R_tr_y[rt],R_tr_ph[rt]); 
		}

		}//NRtr
	}//NLtr
}//ENum

	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
	h_theta_ee->Draw("");
	c1->cd(2);
	h_phi_ee->Draw("");
	c1->cd(3);
	h_thph_ee->Draw("colz");
	TCanvas* c2 = new TCanvas("c2","c2");
	c2->Divide(2,2);
	c2->cd(1);
	h_theta_ek->Draw("");
	c2->cd(2);
	h_phi_ek->Draw("");
	c2->cd(3);
	h_thph_ek->Draw("colz");
	TCanvas* c3 = new TCanvas("c3","c3");
	c3->Divide(2,2);
	c3->cd(1);
	h_theta_g->Draw("");
	c3->cd(2);
	h_phi_g->Draw("");
	c3->cd(3);
	h_thph_g->Draw("colz");
	TCanvas* c4 = new TCanvas("c4","c4");
	h_mom_g->Draw("");
	TCanvas* c5 = new TCanvas("c5","c5");
	h_qsq->Draw("");
	TCanvas* c6 = new TCanvas("c6","c6");
	c6->Divide(2,2);
	c6->cd(1);
	h_theta_gk_lab->Draw("");
	c6->cd(2);
	h_theta_gk_cm->Draw("");
	c6->cd(3);
	h_cos_gk_lab->Draw("");
	c6->cd(4);
	h_cos_gk_cm->Draw("");
	TCanvas* c7 = new TCanvas("c7","c7");
	c7->Divide(2,2);
	c7->cd(1);
	h_pR_lab->Draw("");
	c7->cd(2);
	h_pR_cm->Draw("");
	c7->cd(3);
	//h2_pR_lab_cm->Draw("colz");
	//h_pL_cm->Draw("");
	c7->cd(4);
	//h_pB_cm->Draw("");
	TCanvas* c8 = new TCanvas("c8","c8");
	h_w->Draw("");
	TCanvas* c9 = new TCanvas("c9","c9");
	h_qw->Draw("colz");
	TCanvas* c10 = new TCanvas("c10","c10");
	h_labtocm->Draw("");
	TCanvas* c11 = new TCanvas("c11","c11");
	c11->Divide(1,2);
	c11->cd(1);
	h_z_yfp->Draw("colz");
	c11->cd(2);
	h_z_phi->Draw("colz");
	TCanvas* c12 = new TCanvas("c12","c12");
	h_yfp_phi->Draw("colz");

/*--- Print ---*/
cout << "Print is starting" << endl;
	c1->Print(Form("%s[",pdfname.c_str()));
	c1->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s",pdfname.c_str()));
	c3->Print(Form("%s",pdfname.c_str()));
	c4->Print(Form("%s",pdfname.c_str()));
	c5->Print(Form("%s",pdfname.c_str()));
	c6->Print(Form("%s",pdfname.c_str()));
	c7->Print(Form("%s",pdfname.c_str()));
	c7->Print(Form("%s]",pdfname.c_str()));

	c5->Print("Qsq.pdf");
	c6->Print("theta_gk.pdf");
	c7->Print("CM.pdf");
	c8->Print("W.pdf");
	c9->Print("QsqW.pdf");

cout << "Well done!" << endl;
}//kinematics
