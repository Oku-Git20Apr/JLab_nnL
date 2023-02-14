//--  Missing Mass   --//
//
//K. Okuyama (Nov. 17, 2020)
//K. Okuyama (Jan. 30, 2023)
//
//This is taken over from MM.C
//Comparison befor and after the MM tuning

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

void mm_tuning(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *file_w  = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file_wo = new TFile("../h2all_woMomCalib.root","read");//input file of all H2 run(default: h2all4.root)
  TTree *tree_w  = (TTree*)file_w ->Get("tree_out");
  TTree *tree_wo = (TTree*)file_wo->Get("tree_out");

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
//          Common Definition            //
//---------------------------------------//

  double xmin = -0.1, xmax = 0.2; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_L_strict_w  = new TH1F("hmm_L_strict_w","",xbin,xmin*1000.,xmax*1000.);
  SetTH1(hmm_L_strict_w, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "Counts/(MeV/c^{2})", kRed, kRed, kRed);
  TH1F* hmm_L_strict_wo  = new TH1F("hmm_L_strict_wo","",xbin,xmin*1000.,xmax*1000.);
  SetTH1(hmm_L_strict_wo, "", "Missing Mass - M_{#Lambda} [MeV/c^{2}]", "Counts/(MeV/c^{2})", kAzure, kRed, kRed);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass



//---------------------------------------//
//               Branch                  //
//---------------------------------------//

	int nrun_w, NLtr_w, NRtr_w;
	double ct_w;

	double L_tr_chi2_w;
	double L_tr_x_w, L_tr_y_w, L_tr_th_w, L_tr_ph_w;
	double L_tr_p_w;
	double L_tr_tg_th_w, L_tr_tg_ph_w;
	double L_tr_vz_w, L_tr_vz2_w, L_tr_vz3_w;
	double L_tr_vz_saved_w;
	double R_tr_chi2_w;
	double R_tr_x_w, R_tr_y_w, R_tr_th_w, R_tr_ph_w;
	double R_tr_p_w;
	double R_tr_tg_th_w, R_tr_tg_ph_w;
	double R_tr_vz_w, R_tr_vz2_w, R_tr_vz3_w;
	double L_mom_w, R_mom_w, B_mom_w; 
	double L_ene_w, R_ene_w, B_ene_w; 
	double ac1sum_w, ac2sum_w;//NPE SUM

	tree_w->SetBranchStatus("*",0);
	tree_w->SetBranchStatus("nrun",1);tree_w->SetBranchAddress("nrun",&nrun_w);
	tree_w->SetBranchStatus("tr.ntrack_l",1);tree_w->SetBranchAddress("tr.ntrack_l",&NLtr_w);
	tree_w->SetBranchStatus("tr.ntrack_r",1);tree_w->SetBranchAddress("tr.ntrack_r",&NRtr_w);
	
  	tree_w->SetBranchStatus("ac1_npe_sum",1);  tree_w->SetBranchAddress("ac1_npe_sum", &ac1sum_w);
  	tree_w->SetBranchStatus("ac2_npe_sum",1);  tree_w->SetBranchAddress("ac2_npe_sum", &ac2sum_w);
  	tree_w->SetBranchStatus("Lp_c",1);  tree_w->SetBranchAddress("Lp_c", &L_mom_w);
  	tree_w->SetBranchStatus("Rp_c",1);  tree_w->SetBranchAddress("Rp_c", &R_mom_w);
  	tree_w->SetBranchStatus("Bp_c",1);  tree_w->SetBranchAddress("Bp_c", &B_mom_w);
  	tree_w->SetBranchStatus("ct_orig",1);  tree_w->SetBranchAddress("ct_orig", &ct_w);

  	tree_w->SetBranchStatus("L.tr.chi2",1);  tree_w->SetBranchAddress("L.tr.chi2", &L_tr_chi2_w);
  	tree_w->SetBranchStatus("L.tr.x",1);  tree_w->SetBranchAddress("L.tr.x", &L_tr_x_w);
  	tree_w->SetBranchStatus("L.tr.y",1);  tree_w->SetBranchAddress("L.tr.y", &L_tr_y_w);
  	tree_w->SetBranchStatus("L.tr.th",1);  tree_w->SetBranchAddress("L.tr.th", &L_tr_th_w);
  	tree_w->SetBranchStatus("L.tr.ph",1);  tree_w->SetBranchAddress("L.tr.ph", &L_tr_ph_w);
  	tree_w->SetBranchStatus("L.tr.p",1);  tree_w->SetBranchAddress("L.tr.p", &L_tr_p_w);
  	tree_w->SetBranchStatus("L.tr.tg_th",1);  tree_w->SetBranchAddress("L.tr.tg_th", &L_tr_tg_th_w );
  	tree_w->SetBranchStatus("L.tr.tg_ph",1);  tree_w->SetBranchAddress("L.tr.tg_ph", &L_tr_tg_ph_w );
  	tree_w->SetBranchStatus("L.tr.vz",1);  tree_w->SetBranchAddress("L.tr.vz", &L_tr_vz_w);

  	tree_w->SetBranchStatus("R.tr.chi2",1);  tree_w->SetBranchAddress("R.tr.chi2", &R_tr_chi2_w);
	tree_w->SetBranchStatus("R.tr.x" ,1);  tree_w->SetBranchAddress("R.tr.x" , &R_tr_x_w );
  	tree_w->SetBranchStatus("R.tr.y" ,1);  tree_w->SetBranchAddress("R.tr.y" , &R_tr_y_w );
  	tree_w->SetBranchStatus("R.tr.th",1);  tree_w->SetBranchAddress("R.tr.th", &R_tr_th_w);
  	tree_w->SetBranchStatus("R.tr.ph",1);  tree_w->SetBranchAddress("R.tr.ph", &R_tr_ph_w);
  	tree_w->SetBranchStatus("R.tr.p",1);  tree_w->SetBranchAddress("R.tr.p", &R_tr_p_w);
  	tree_w->SetBranchStatus("R.tr.tg_th",1);  tree_w->SetBranchAddress("R.tr.tg_th", &R_tr_tg_th_w);
  	tree_w->SetBranchStatus("R.tr.tg_ph",1);  tree_w->SetBranchAddress("R.tr.tg_ph", &R_tr_tg_ph_w);
  	tree_w->SetBranchStatus("R.tr.vz",1);  tree_w->SetBranchAddress("R.tr.vz", &R_tr_vz_w);

  int ENum_w = tree_w->GetEntries(); 
cout<<"Entries: "<<ENum_w<<endl;
  int time_div_w=ENum_w/25;
  if(ENum_w<100000)time_div_w=10000;


	time_t start_w, end_w;
	start_w = time(NULL);
	time(&start_w);

  for(int i=0;i<ENum_w;i++){
	tree_w->GetEntry(i);

    if(i%time_div_w==0){
      end_w = time(NULL);
      time(&end_w);
      double diff_w = difftime(end_w,start_w);
      double esttime_w = diff_w * ENum_w / (i+1) - diff_w;
      cout<<i<<" / "<<ENum_w<<" ("<<i*100/ENum_w<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff_w,esttime_w)<<endl;
    }
      
        L_Tr = L_FP = false;
        if( L_tr_chi2_w<0.01 ) L_Tr = true;
        if( L_tr_th_w<0.17*L_tr_x_w+0.025
         && L_tr_th_w>0.17*L_tr_x_w-0.035
         && L_tr_th_w<0.40*L_tr_x_w+0.130 ) L_FP = true;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2_w<0.01 ) R_Tr = true;
        if( R_tr_th_w<0.17*R_tr_x_w+0.025
         && R_tr_th_w>0.17*R_tr_x_w-0.035
         && R_tr_th_w<0.40*R_tr_x_w+0.130 ) R_FP = true;

		if(fabs(ct_w)<1.006)ct_cut=true;
		else ct_cut=false;
		if(fabs(L_tr_vz_w-R_tr_vz_w)<0.025&&fabs(R_tr_vz_w+L_tr_vz_w)<0.2&&ac1sum_w<3.75&&ac2sum_w>3.&&ac2sum_w<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;

	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama
	    double R_pz_w = R_mom_w/sqrt(1.0*1.0 + pow((R_tr_tg_th_w), 2.0) + pow(( R_tr_tg_ph_w),2.0) );
	    double R_px_w = R_pz_w * ( R_tr_tg_th_w );
	    double R_py_w = R_pz_w * ( R_tr_tg_ph_w );

	    double L_pz_w = L_mom_w/sqrt(1.0*1.0 + pow(( L_tr_tg_th_w ), 2.0) + pow(( L_tr_tg_ph_w),2.0));
	    double L_px_w = L_pz_w * ( L_tr_tg_th_w );
	    double L_py_w = L_pz_w * ( L_tr_tg_ph_w );

	    double B_E_w =sqrt(B_mom_w*B_mom_w + Me*Me);
	    double R_E_w =sqrt(R_mom_w*R_mom_w + MK*MK);
	    double L_E_w =sqrt(L_mom_w*L_mom_w + Me*Me);

		TLorentzVector L_4vec_w;//Left
		TLorentzVector R_4vec_w;//Right
		TLorentzVector B_4vec_w;//Beam
		TLorentzVector T_4vec_w;//Target
		TLorentzVector G_4vec_w;//Gamma (Virtual Photon)
		L_4vec_w.SetPxPyPzE(L_px_w, L_py_w, L_pz_w, L_E_w);
        R_4vec_w.SetPxPyPzE(R_px_w, R_py_w, R_pz_w, R_E_w);
        B_4vec_w.SetPxPyPzE(0.0 ,  0.0,B_mom_w, B_E_w);
        T_4vec_w.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);

	    L_4vec_w.RotateX( -13.2/180.*PI );
	    R_4vec_w.RotateX(  13.2/180.*PI );

        double mass_w,mm_w;
		TLorentzVector Missing_w;
		Missing_w = B_4vec_w + T_4vec_w - L_4vec_w - R_4vec_w;
		mass_w = Missing_w.M();
	    mm_w=mass_w - mh;//shift by ML
		
		if(event_selection&&ct_cut)hmm_L_strict_w->Fill(mm_w*1000.);
		double theta_ee_w = L_4vec_w.Theta();
		double theta_ek_w = R_4vec_w.Theta();
		double phi_ee_w = L_4vec_w.Phi();//original frame
		double phi_ek_w = R_4vec_w.Phi()+2*PI;//original frame

		G_4vec_w = B_4vec_w - L_4vec_w;
		double mom_g_w=sqrt(G_4vec_w.Px()*G_4vec_w.Px()+G_4vec_w.Py()*G_4vec_w.Py()+G_4vec_w.Pz()*G_4vec_w.Pz());
		double Qsq_w = G_4vec_w.M()*G_4vec_w.M();
		double phi_g_w = G_4vec_w.Phi()+2*PI;
		double theta_g_w = G_4vec_w.Theta();
		double theta_gk_lab_w = G_4vec_w.Angle(R_4vec_w.Vect());
		double omega_w=G_4vec_w.E();
		double beta_w=mom_g_w/(omega_w+Mp);
	
		TVector3 boost_w;
		TLorentzVector GT_4vec_w;
		GT_4vec_w=G_4vec_w+T_4vec_w;
		boost_w=GT_4vec_w.BoostVector();
		R_4vec_w.Boost(-boost_w);
		L_4vec_w.Boost(-boost_w);
		B_4vec_w.Boost(-boost_w);
		double theta_gk_cm_w = G_4vec_w.Angle(R_4vec_w.Vect());
		double theta_ek_cm_w = R_4vec_w.Theta();
		double pR_cm_w=sqrt(R_4vec_w.Px()*R_4vec_w.Px()+R_4vec_w.Py()*R_4vec_w.Py()+R_4vec_w.Pz()*R_4vec_w.Pz());
		double pL_cm_w=sqrt(L_4vec_w.Px()*L_4vec_w.Px()+L_4vec_w.Py()*L_4vec_w.Py()+L_4vec_w.Pz()*L_4vec_w.Pz());
		double pB_cm_w=sqrt(B_4vec_w.Px()*B_4vec_w.Px()+B_4vec_w.Py()*B_4vec_w.Py()+B_4vec_w.Pz()*B_4vec_w.Pz());

		double n_w = MK/ML;
		double p_cm_w=sqrt(GT_4vec_w.Px()*GT_4vec_w.Px()+GT_4vec_w.Py()*GT_4vec_w.Py()+GT_4vec_w.Pz()*GT_4vec_w.Pz());
		double E_cm_w = GT_4vec_w.E();
		double gamma_w=1./sqrt(1-beta_w*beta_w);
		double ER_cm_w=sqrt(pR_cm_w*pR_cm_w+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm_w = (gamma_w*pR_cm_w*pR_cm_w*(pR_cm_w*cos(theta_gk_cm_w)+beta_w*ER_cm_w))/(pow(sqrt(pR_cm_w*pR_cm_w*sin(theta_gk_cm_w)*sin(theta_gk_cm_w)+gamma_w*gamma_w*(pR_cm_w*cos(theta_gk_cm_w)+beta_w*ER_cm_w)*(pR_cm_w*cos(theta_gk_cm_w)+beta_w*ER_cm_w)),3.));
		double tan_lab1_w = sin(theta_gk_cm_w)/(gamma_w*(cos(theta_gk_cm_w)+beta_w*sqrt(MK*MK+pR_cm_w*pR_cm_w)/pR_cm_w));
		double tan_lab2_w = sin(theta_gk_cm_w)/(gamma_w*(cos(theta_gk_cm_w)+(omega_w*Mp-Qsq_w*Qsq_w)/(omega_w*Mp+Mp*Mp)));
		//if(tan_lab1!=tan_lab2)cout<<"tan1="<<atan(tan_lab1)<<", tan2="<<atan(tan_lab2)<<"theta_gk_lab="<<theta_gk_lab<<endl;

}//ENum_w



	int nrun_wo, NLtr_wo, NRtr_wo;
	double ct_wo;

	double L_tr_chi2_wo;
	double L_tr_x_wo, L_tr_y_wo, L_tr_th_wo, L_tr_ph_wo;
	double L_tr_p_wo;
	double L_tr_tg_th_wo, L_tr_tg_ph_wo;
	double L_tr_vz_wo, L_tr_vz2_wo, L_tr_vz3_wo;
	double L_tr_vz_saved_wo;
	double R_tr_chi2_wo;
	double R_tr_x_wo, R_tr_y_wo, R_tr_th_wo, R_tr_ph_wo;
	double R_tr_p_wo;
	double R_tr_tg_th_wo, R_tr_tg_ph_wo;
	double R_tr_vz_wo, R_tr_vz2_wo, R_tr_vz3_wo;
	double L_mom_wo, R_mom_wo, B_mom_wo; 
	double L_ene_wo, R_ene_wo, B_ene_wo; 
	double ac1sum_wo, ac2sum_wo;//NPE SUM

	tree_wo->SetBranchStatus("*",0);
	tree_wo->SetBranchStatus("nrun",1);tree_wo->SetBranchAddress("nrun",&nrun_wo);
	tree_wo->SetBranchStatus("tr.ntrack_l",1);tree_wo->SetBranchAddress("tr.ntrack_l",&NLtr_wo);
	tree_wo->SetBranchStatus("tr.ntrack_r",1);tree_wo->SetBranchAddress("tr.ntrack_r",&NRtr_wo);
	
  	tree_wo->SetBranchStatus("ac1_npe_sum",1);  tree_wo->SetBranchAddress("ac1_npe_sum", &ac1sum_wo);
  	tree_wo->SetBranchStatus("ac2_npe_sum",1);  tree_wo->SetBranchAddress("ac2_npe_sum", &ac2sum_wo);
  	tree_wo->SetBranchStatus("Lp_c",1);  tree_wo->SetBranchAddress("Lp_c", &L_mom_wo);
  	tree_wo->SetBranchStatus("Rp_c",1);  tree_wo->SetBranchAddress("Rp_c", &R_mom_wo);
  	tree_wo->SetBranchStatus("Bp_c",1);  tree_wo->SetBranchAddress("Bp_c", &B_mom_wo);
  	tree_wo->SetBranchStatus("ct_orig",1);  tree_wo->SetBranchAddress("ct_orig", &ct_wo);

  	tree_wo->SetBranchStatus("L.tr.chi2",1);  tree_wo->SetBranchAddress("L.tr.chi2", &L_tr_chi2_wo);
  	tree_wo->SetBranchStatus("L.tr.x",1);  tree_wo->SetBranchAddress("L.tr.x", &L_tr_x_wo);
  	tree_wo->SetBranchStatus("L.tr.y",1);  tree_wo->SetBranchAddress("L.tr.y", &L_tr_y_wo);
  	tree_wo->SetBranchStatus("L.tr.th",1);  tree_wo->SetBranchAddress("L.tr.th", &L_tr_th_wo);
  	tree_wo->SetBranchStatus("L.tr.ph",1);  tree_wo->SetBranchAddress("L.tr.ph", &L_tr_ph_wo);
  	tree_wo->SetBranchStatus("L.tr.p",1);  tree_wo->SetBranchAddress("L.tr.p", &L_tr_p_wo);
  	tree_wo->SetBranchStatus("L.tr.tg_th",1);  tree_wo->SetBranchAddress("L.tr.tg_th", &L_tr_tg_th_wo );
  	tree_wo->SetBranchStatus("L.tr.tg_ph",1);  tree_wo->SetBranchAddress("L.tr.tg_ph", &L_tr_tg_ph_wo );
  	tree_wo->SetBranchStatus("L.tr.vz",1);  tree_wo->SetBranchAddress("L.tr.vz", &L_tr_vz_wo);

  	tree_wo->SetBranchStatus("R.tr.chi2",1);  tree_wo->SetBranchAddress("R.tr.chi2", &R_tr_chi2_wo);
	tree_wo->SetBranchStatus("R.tr.x" ,1);  tree_wo->SetBranchAddress("R.tr.x" , &R_tr_x_wo );
  	tree_wo->SetBranchStatus("R.tr.y" ,1);  tree_wo->SetBranchAddress("R.tr.y" , &R_tr_y_wo );
  	tree_wo->SetBranchStatus("R.tr.th",1);  tree_wo->SetBranchAddress("R.tr.th", &R_tr_th_wo);
  	tree_wo->SetBranchStatus("R.tr.ph",1);  tree_wo->SetBranchAddress("R.tr.ph", &R_tr_ph_wo);
  	tree_wo->SetBranchStatus("R.tr.p",1);  tree_wo->SetBranchAddress("R.tr.p", &R_tr_p_wo);
  	tree_wo->SetBranchStatus("R.tr.tg_th",1);  tree_wo->SetBranchAddress("R.tr.tg_th", &R_tr_tg_th_wo);
  	tree_wo->SetBranchStatus("R.tr.tg_ph",1);  tree_wo->SetBranchAddress("R.tr.tg_ph", &R_tr_tg_ph_wo);
  	tree_wo->SetBranchStatus("R.tr.vz",1);  tree_wo->SetBranchAddress("R.tr.vz", &R_tr_vz_wo);

  int ENum_wo = tree_wo->GetEntries(); 
cout<<"Entries: "<<ENum_wo<<endl;
  int time_div_wo=ENum_wo/25;
  if(ENum_wo<100000)time_div_wo=10000;


	time_t start_wo, end_wo;
	start_wo = time(NULL);
	time(&start_wo);

  for(int i=0;i<ENum_wo;i++){
	tree_wo->GetEntry(i);

    if(i%time_div_wo==0){
      end_wo = time(NULL);
      time(&end_wo);
      double diff_wo = difftime(end_wo,start_wo);
      double esttime_wo = diff_wo * ENum_wo / (i+1) - diff_wo;
      cout<<i<<" / "<<ENum_wo<<" ("<<i*100/ENum_wo<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff_wo,esttime_wo)<<endl;
    }
      
        L_Tr = L_FP = false;
        if( L_tr_chi2_wo<0.01 ) L_Tr = true;
        if( L_tr_th_wo<0.17*L_tr_x_wo+0.025
         && L_tr_th_wo>0.17*L_tr_x_wo-0.035
         && L_tr_th_wo<0.40*L_tr_x_wo+0.130 ) L_FP = true;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2_wo<0.01 ) R_Tr = true;
        if( R_tr_th_wo<0.17*R_tr_x_wo+0.025
         && R_tr_th_wo>0.17*R_tr_x_wo-0.035
         && R_tr_th_wo<0.40*R_tr_x_wo+0.130 ) R_FP = true;

		if(fabs(ct_wo)<1.006)ct_cut=true;
		else ct_cut=false;
		if(fabs(L_tr_vz_wo-R_tr_vz_wo)<0.025&&fabs(R_tr_vz_wo+L_tr_vz_wo)<0.2&&ac1sum_wo<3.75&&ac2sum_wo>3.&&ac2sum_wo<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;

	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama
	    double R_pz_wo = R_mom_wo/sqrt(1.0*1.0 + pow((R_tr_tg_th_wo), 2.0) + pow(( R_tr_tg_ph_wo),2.0) );
	    double R_px_wo = R_pz_wo * ( R_tr_tg_th_wo );
	    double R_py_wo = R_pz_wo * ( R_tr_tg_ph_wo );

	    double L_pz_wo = L_mom_wo/sqrt(1.0*1.0 + pow(( L_tr_tg_th_wo ), 2.0) + pow(( L_tr_tg_ph_wo),2.0));
	    double L_px_wo = L_pz_wo * ( L_tr_tg_th_wo );
	    double L_py_wo = L_pz_wo * ( L_tr_tg_ph_wo );

	    double B_E_wo =sqrt(B_mom_wo*B_mom_wo + Me*Me);
	    double R_E_wo =sqrt(R_mom_wo*R_mom_wo + MK*MK);
	    double L_E_wo =sqrt(L_mom_wo*L_mom_wo + Me*Me);

		TLorentzVector L_4vec_wo;//Left
		TLorentzVector R_4vec_wo;//Right
		TLorentzVector B_4vec_wo;//Beam
		TLorentzVector T_4vec_wo;//Target
		TLorentzVector G_4vec_wo;//Gamma (Virtual Photon)
		L_4vec_wo.SetPxPyPzE(L_px_wo, L_py_wo, L_pz_wo, L_E_wo);
        R_4vec_wo.SetPxPyPzE(R_px_wo, R_py_wo, R_pz_wo, R_E_wo);
        B_4vec_wo.SetPxPyPzE(0.0 ,  0.0,B_mom_wo, B_E_wo);
        T_4vec_wo.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);

	    L_4vec_wo.RotateX( -13.2/180.*PI );
	    R_4vec_wo.RotateX(  13.2/180.*PI );

        double mass_wo,mm_wo;
		TLorentzVector Missing_wo;
		Missing_wo = B_4vec_wo + T_4vec_wo - L_4vec_wo - R_4vec_wo;
		mass_wo = Missing_wo.M();
	    mm_wo=mass_wo - mh;//shift by ML
		
		if(event_selection&&ct_cut)hmm_L_strict_wo->Fill(mm_wo*1000.);
		double theta_ee_wo = L_4vec_wo.Theta();
		double theta_ek_wo = R_4vec_wo.Theta();
		double phi_ee_wo = L_4vec_wo.Phi();//original frame
		double phi_ek_wo = R_4vec_wo.Phi()+2*PI;//original frame

		G_4vec_wo = B_4vec_wo - L_4vec_wo;
		double mom_g_wo=sqrt(G_4vec_wo.Px()*G_4vec_wo.Px()+G_4vec_wo.Py()*G_4vec_wo.Py()+G_4vec_wo.Pz()*G_4vec_wo.Pz());
		double Qsq_wo = G_4vec_wo.M()*G_4vec_wo.M();
		double phi_g_wo = G_4vec_wo.Phi()+2*PI;
		double theta_g_wo = G_4vec_wo.Theta();
		double theta_gk_lab_wo = G_4vec_wo.Angle(R_4vec_wo.Vect());
		double omega_wo=G_4vec_wo.E();
		double beta_wo=mom_g_wo/(omega_wo+Mp);
	
		TVector3 boost_wo;
		TLorentzVector GT_4vec_wo;
		GT_4vec_wo=G_4vec_wo+T_4vec_wo;
		boost_wo=GT_4vec_wo.BoostVector();
		R_4vec_wo.Boost(-boost_wo);
		L_4vec_wo.Boost(-boost_wo);
		B_4vec_wo.Boost(-boost_wo);
		double theta_gk_cm_wo = G_4vec_wo.Angle(R_4vec_wo.Vect());
		double theta_ek_cm_wo = R_4vec_wo.Theta();
		double pR_cm_wo=sqrt(R_4vec_wo.Px()*R_4vec_wo.Px()+R_4vec_wo.Py()*R_4vec_wo.Py()+R_4vec_wo.Pz()*R_4vec_wo.Pz());
		double pL_cm_wo=sqrt(L_4vec_wo.Px()*L_4vec_wo.Px()+L_4vec_wo.Py()*L_4vec_wo.Py()+L_4vec_wo.Pz()*L_4vec_wo.Pz());
		double pB_cm_wo=sqrt(B_4vec_wo.Px()*B_4vec_wo.Px()+B_4vec_wo.Py()*B_4vec_wo.Py()+B_4vec_wo.Pz()*B_4vec_wo.Pz());

		double n_wo = MK/ML;
		double p_cm_wo=sqrt(GT_4vec_wo.Px()*GT_4vec_wo.Px()+GT_4vec_wo.Py()*GT_4vec_wo.Py()+GT_4vec_wo.Pz()*GT_4vec_wo.Pz());
		double E_cm_wo = GT_4vec_wo.E();
		double gamma_wo=1./sqrt(1-beta_wo*beta_wo);
		double ER_cm_wo=sqrt(pR_cm_wo*pR_cm_wo+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm_wo = (gamma_wo*pR_cm_wo*pR_cm_wo*(pR_cm_wo*cos(theta_gk_cm_wo)+beta_wo*ER_cm_wo))/(pow(sqrt(pR_cm_wo*pR_cm_wo*sin(theta_gk_cm_wo)*sin(theta_gk_cm_wo)+gamma_wo*gamma_wo*(pR_cm_wo*cos(theta_gk_cm_wo)+beta_wo*ER_cm_wo)*(pR_cm_wo*cos(theta_gk_cm_wo)+beta_wo*ER_cm_wo)),3.));
		double tan_lab1_wo = sin(theta_gk_cm_wo)/(gamma_wo*(cos(theta_gk_cm_wo)+beta_wo*sqrt(MK*MK+pR_cm_wo*pR_cm_wo)/pR_cm_wo));
		double tan_lab2_wo = sin(theta_gk_cm_wo)/(gamma_wo*(cos(theta_gk_cm_wo)+(omega_wo*Mp-Qsq_wo*Qsq_wo)/(omega_wo*Mp+Mp*Mp)));
		//if(tan_lab1!=tan_lab2)cout<<"tan1="<<atan(tan_lab1)<<", tan2="<<atan(tan_lab2)<<"theta_gk_lab="<<theta_gk_lab<<endl;

}//ENum_wo

	TCanvas* c3 = new TCanvas("c3","c3");
	hmm_L_strict_w->Draw("");
	hmm_L_strict_wo->Draw("same");
	
	c3->Print("./pdf/mm_tuning.pdf");

cout << "Well done!" << endl;
}//fit
