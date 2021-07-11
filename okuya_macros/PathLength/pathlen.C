//--  Path Length   --//
//Data vs. SIMC
//
//K. Okuyama (Jun. 14, 2021)
//
//This is taken over from pkpe.C
//No array branch mode 
//
void pathlen(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file = new TFile("../h2all_2021Mar.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_2020Nov.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TTree *tree = (TTree*)file->Get("tree_out");
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/radAl_June.root","read");
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");


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



//---Fitting Function---//
 TF1* fmmbg_nocut;
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
	double R_tr_pathl;
	double R_s2_trpath;

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
  	tree->SetBranchStatus("R.tr.pathl",1);  tree->SetBranchAddress("R.tr.pathl", &R_tr_pathl);
  	tree->SetBranchStatus("R.s2.trpath",1);  tree->SetBranchAddress("R.s2.trpath", &R_s2_trpath);




  TH1F* h_tr_pathl  = new TH1F("h_tr_pathl","h_tr_pathl",1000,25.,27.);
  h_tr_pathl->SetLineColor(kAzure);
  TH1F* h_s2_trpath  = new TH1F("h_s2_trpath","h_s2_trpath",1000,3.,3.6);
  h_s2_trpath->SetLineColor(kAzure);
  TH1F* h_path_ev  = new TH1F("h_path_ev","PathLength vs. #ev",5000,1.,5000.);
  h_path_ev->GetXaxis()->SetTitle("#ev");
  h_path_ev->GetYaxis()->SetTitle("L(Data)-L(SIMC) [cm]");
  h_path_ev->SetLineColor(kAzure);
  TH1F* h_path_dif  = new TH1F("h_path_dif","PathLength Diff.",1000,-50,50.);
  h_path_dif->GetXaxis()->SetTitle("L(Data)-L(SIMC) [cm]");
  h_path_dif->GetYaxis()->SetTitle("Counts");
  h_path_dif->SetLineColor(kAzure);
  TH2F* h_theta_gk_cm_eps = new TH2F("h_theta_gk_cm_eps","",50,-5.,25.,50,0.72,0.8);
  h_theta_gk_cm_eps->GetYaxis()->SetTitle("#epsilon");
  h_theta_gk_cm_eps->GetXaxis()->SetTitle("#theta_{#gamma K}^{CM} [deg]");

  TGraph* g_path = new TGraph();
	g_path->SetMarkerStyle(21);
	g_path->SetMarkerColor(kRed);
	g_path->SetMarkerSize(1.0);





  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool ct_cut_ctout = false;
  bool event_selection = false;
  bool event_selection_center = false;
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

		int evnum=1;


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
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.005&&abs(R_mom-1.82)<0.02&&abs(R_tr_tg_th)<0.03&&abs(R_tr_tg_ph)<0.01)event_selection_center=true;
		else event_selection_center=false;


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
		
		//if(event_selection&&event_selection_center){
		if(event_selection){
			h_tr_pathl->Fill(R_tr_pathl); 
  			h_s2_trpath->Fill(R_s2_trpath); 
			h_path_ev->SetBinContent(evnum,R_tr_pathl*100.-2740.+329.2-172.3);
			h_path_dif->Fill(R_tr_pathl*100.-2740.+329.2-172.3);
			g_path->SetPoint(evnum,evnum,R_tr_pathl*100.-2740.+329.2-172.3);
			evnum++;
		}
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


}//ENum

//	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
		h_path_ev->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2");
		h_path_dif->Draw("");
	TCanvas* c3 = new TCanvas("c3","c3");
		g_path->Draw("AP");
	TCanvas* c4 = new TCanvas("c4","c4");
		h_tr_pathl->Draw("");
	TCanvas* c5 = new TCanvas("c5","c5");
		h_s2_trpath->Draw("");


cout << "Well done!" << endl;
}
