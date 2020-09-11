//-----------------------------//
//--  Cointime Histo.        --//
//-----------------------------//
//
//K. Okuyama (Sep. 5, 2020)

double fcoin_template( double *x, double *par , int shift, int num)
{
  return par[num] * TMath::Gaus(x[0],par[num+1]-2.0*shift,par[num+2]);//Lambda Gaussian
}

double fcoin_acc( double *x, double *par, int num)
{
	double val=0.;
	for(int i=0;i<=60;i++){
	val += par[num]*TMath::Gaus(x[0],par[num+1]+par[num+12]*i-20,par[num+2]);
	}
	return val;
}

double fcoin_total( double *x, double *par ){

	//return fcoin_template(x,par,-10,0)+expgaus2(x,par,6)+expgaus2(x,par,10);//+expgaus2(x,par,14);
	return fcoin_acc(x,par,0)+par[3]*TMath::Gaus(x[0],par[4],par[5])+par[6]*TMath::Gaus(x[0],par[7],par[8])+par[9]*TMath::Gaus(x[0],par[10],par[11]);

}

void hcoin_wide(){
	string pdfname = "temp.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("../h2all5.root","read");//input file of all H2 run(default: h2all4.root)
	//ACCBGの引き算はmea_hist.ccから
  TFile *file_mea = new TFile("../bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  double nbunch = 600.;//effetive bunches (6 bunches x 5 mixtures)
 // TTree *tree_old = (TTree*)file->Get("tree_out");
//cout<<"Please wait a moment. CloneTree() is working..."<<endl;
  //TTree *tree = tree_old->CloneTree();
  TTree *tree = (TTree*)file->Get("tree_out");
//	tree->Write();

    
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
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;


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

//---------------------------------------//
//               Branch                  //
//---------------------------------------//

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
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
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");

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

  TH1F* hcoin  = new TH1F("hcoin","",120000/56,-20.,100.);
  
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



  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  tree->Draw(">>elist" , "fabs(ct_orig)<1.0");//ctsum (does NOT dintinguish #track)
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
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


	

		bool ac_cut=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.)ac_cut=true;

		if(fabs(ct)<1)ct_cut=true;
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

		if(ac_cut)hcoin->Fill(ct);

//	    //===== Right Hand Coordinate ====//
//	    //th and phi are originally meant tan(theta) and tan(phi),
//	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama
//	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
//	    double R_px = R_pz * (R_tr_tg_th );
//	    double R_py = R_pz * ( R_tr_tg_ph );
//
//	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(( L_tr_tg_th ), 2.0) + pow(( L_tr_tg_ph),2.0));
//	    double L_px = L_pz * ( L_tr_tg_th );
//	    double L_py = L_pz * ( L_tr_tg_ph );
//
//	    double B_E =sqrt(B_mom*B_mom + Me*Me);
//	    double R_E =sqrt(R_mom*R_mom + MK*MK);
//	    double L_E =sqrt(L_mom*L_mom + Me*Me);
//
//		TLorentzVector L_4vec;//Left
//		TLorentzVector R_4vec;//Right
//		TLorentzVector B_4vec;//Beam
//		TLorentzVector T_4vec;//Target
//		TLorentzVector G_4vec;//Gamma (Virtual Photon)
//		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
//        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
//        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
//        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);
//
//	    L_4vec.RotateX( -13.2/180.*PI );
//	    R_4vec.RotateX(  13.2/180.*PI );
//
//        double mass,mm;
//		TLorentzVector Missing;
//		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
//		mass = Missing.M();
//	    mm=mass - mh;//shift by ML
//		if(event_selection&&ct_cut)hmm_L_fom_best->Fill(mm);
//		if(event_selection_nocut&&ct_cut)hmm_L_fom_nocut->Fill(mm);



}//ENum

	 TF1 *fcoin=new TF1("fcoin",fcoin_total,-20.,100.,13);
	 fcoin->SetNpx(2000);
	 fcoin->SetParameter(0,500.);
	 fcoin->SetParameter(1,0.);
	 fcoin->SetParameter(2,0.5);
	 fcoin->SetParameter(3,0.);//pi
	 fcoin->SetParameter(4,3.18);
	 fcoin->SetParLimits(4,2.5,3.5);
	 fcoin->SetParameter(5,0.3);
	 fcoin->SetParameter(6,0.);//k
	 fcoin->FixParameter(7,0.);
	 fcoin->SetParameter(8,0.3);
	 fcoin->SetParameter(9,0.);//p
	 fcoin->SetParameter(10,-8.);
	 fcoin->SetParLimits(10,-9.,-7.);
	 fcoin->SetParameter(11,0.7);
	 fcoin->FixParameter(12,2.0);
	 fcoin->SetLineColor(kAzure);
	 hcoin->Fit("fcoin","","",12.,60.);
	 double par0=fcoin->GetParameter(0);
	 double par1=fcoin->GetParameter(1);
	 double par2=fcoin->GetParameter(2);
	 double par12=fcoin->GetParameter(12);
	 fcoin->FixParameter(0,par0);
	 fcoin->FixParameter(1,par1);
	 fcoin->FixParameter(2,par2);
	 fcoin->FixParameter(12,par12);
	 fcoin->SetParameter(6,200.);//k
	 fcoin->SetParameter(7,0.);//k
	 fcoin->SetParLimits(7,-0.5,0.5);
	 fcoin->SetParameter(8,0.3);//k
	 hcoin->Fit("fcoin","","",-1.,1.);
	 par0=fcoin->GetParameter(0);
	 par1=fcoin->GetParameter(1);
	 par2=fcoin->GetParameter(2);
	 double par6=fcoin->GetParameter(6);
	 double par7=fcoin->GetParameter(7);
	 double par8=fcoin->GetParameter(8);
	 par12=fcoin->GetParameter(12);
	 fcoin->FixParameter(0,par0);
	 fcoin->FixParameter(1,par1);
	 fcoin->FixParameter(2,par2);
	 fcoin->FixParameter(12,par12);
	 fcoin->SetParameter(3,1200.);//pi
	 fcoin->SetParameter(4,3.18);
	 fcoin->SetParLimits(4,2.5,3.5);
	 fcoin->SetParameter(5,0.3);
	 fcoin->FixParameter(6,par6);//k
	 fcoin->FixParameter(7,par7);
	 fcoin->FixParameter(8,par8);
	 fcoin->SetParameter(9,800.);//p
	 fcoin->SetParameter(10,-8.);
	 fcoin->SetParLimits(10,-9.,-7.);
	 fcoin->SetParameter(11,0.7);
	 hcoin->Fit("fcoin","","",-20.,20.);
	 TF1 *fcoin_pi=new TF1("fcoin_pi","gaus",-20.,20.);
	 TF1 *fcoin_k=new TF1("fcoin_k","gaus",-20.,20.);
	 TF1 *fcoin_p=new TF1("fcoin_p","gaus",-20.,20.);
	 fcoin_pi->SetNpx(2000);
	 fcoin_k->SetNpx(2000);
	 fcoin_p->SetNpx(2000);
	 fcoin_pi->SetParameter(0,fcoin->GetParameter(3));
	 fcoin_pi->SetParameter(1,fcoin->GetParameter(4));
	 fcoin_pi->SetParameter(2,fcoin->GetParameter(5));
	 fcoin_k->SetParameter(0,fcoin->GetParameter(6));
	 fcoin_k->SetParameter(1,fcoin->GetParameter(7));
	 fcoin_k->SetParameter(2,fcoin->GetParameter(8));
	 fcoin_p->SetParameter(0,fcoin->GetParameter(9));
	 fcoin_p->SetParameter(1,fcoin->GetParameter(10));
	 fcoin_p->SetParameter(2,fcoin->GetParameter(11));

		TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
		hcoin->Draw("");
		TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
		hcoin->Draw("");
		fcoin_pi->SetLineColor(kOrange);
		fcoin_k->SetLineColor(kGreen);
		fcoin_p->SetLineColor(kRed);
		fcoin->Draw("same");
		fcoin_pi->Draw("same");
		fcoin_k->Draw("same");
		fcoin_p->Draw("same");
		
	cout<<"left bunch = "<<hcoin->Integral(hcoin->FindBin(-17.),hcoin->FindBin(-15.))<<endl;		
	cout<<"right bunch = "<<hcoin->Integral(hcoin->FindBin(75.),hcoin->FindBin(77.))<<endl;		
	cout<<"center bunch = "<<hcoin->Integral(hcoin->FindBin(-5.),hcoin->FindBin(-3.))<<endl;		

//	cout<<"nbunch="<<nbunch<<endl;
//	TCanvas* c1 = new TCanvas("c1","c1");
//	hmm_L_fom_best->Draw("");
//	TCanvas* c2 = new TCanvas("c2","c2");
//	TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
//	TH1F* hmm_Al_fom_nocut=(TH1F*)file->Get("hmm_Al_fom_best");
//	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
//	//TH1F* hmm_pi_fom_nocut=(TH1F*)file->Get("hmm_pi_fom_noZ");
//	//TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
//	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result_best");
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
//	hmm_bg_fom_best->Sumw2();
//	hmm_bg_fom_nocut->Sumw2();
//	hmm_Albg_fom_nocut->Sumw2();
//	hmm_bg_fom_best->Scale(1./nbunch);
//	hmm_bg_fom_nocut->Scale(1./nbunch);
//	hmm_Albg_fom_nocut->Scale(1./nbunch);
//	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
//	hmm_wo_bg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
//	hmm_wo_bg_fom_nocut->Add(hmm_L_fom_nocut,hmm_bg_fom_nocut,1.0,-1.0);
//	hmm_pi_wobg_fom_best->Add(hmm_pi_fom_best,hmm_bg_fom_best,1.0,-1.0);
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
//	 fmm_best_4Poly=new TF1("fmm_best_4Poly",FMM_Res,min_mm,max_mm,14);
//	 fmmbg_best_4Poly=new TF1("fmmbg_best_4Poly","pol4",min_mm,max_mm);
//	 fmm_best_4Poly->SetNpx(20000);
//	 fmm_best_4Poly->SetTitle("Missing Mass (best)");
//	 fmm_best_4Poly->SetParLimits(0,0.,1000.);//positive
//	 fmm_best_4Poly->SetParLimits(3,0.,300.);//positive
//	 fmm_best_4Poly->SetParameter(0,const_L_best*0.85);
//	 fmm_best_4Poly->SetParameter(1,mean_L_best);
//	 fmm_best_4Poly->SetParLimits(1,def_mean_L-def_sig_L*0.4,def_mean_L+def_sig_L*0.4);
//	 fmm_best_4Poly->SetParameter(2,sig_L_best);
//	 fmm_best_4Poly->SetParLimits(2,0.,0.01);
//	 fmm_best_4Poly->SetParameter(3,const_S_best*0.85);
//	 fmm_best_4Poly->SetParameter(4,mean_S_best);
//	 fmm_best_4Poly->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
//	 fmm_best_4Poly->SetParameter(5,sig_S_best);
//	 fmm_best_4Poly->SetParLimits(5,0.,0.003);
////subL
//	 fmm_best_4Poly->SetParameter(6,0.7);//scale
//	 fmm_best_4Poly->SetParLimits(6,0.,1.5);
//	 fmm_best_4Poly->SetParameter(7,0.02);//att.
//	 fmm_best_4Poly->SetParLimits(7,0.005,0.05);
//	 fmm_best_4Poly->SetParameter(8,0.02);//sigma
//	 fmm_best_4Poly->SetParLimits(8,0.001,0.1);
//	 fmm_best_4Poly->SetParameter(9,0.);//peak pos.
//	 fmm_best_4Poly->SetParLimits(9,-0.05,0.05);
////Sigma0
//	 fmm_best_4Poly->SetParameter(10,0.25);
//	 fmm_best_4Poly->SetParLimits(10,0.,1.5);
//	 fmm_best_4Poly->SetParameter(11,0.08);
//	 fmm_best_4Poly->SetParLimits(11,0.04,0.12);
//	 fmm_best_4Poly->SetParameter(12,0.01);
//	 fmm_best_4Poly->SetParLimits(12,0.001,0.01);
//	 fmm_best_4Poly->SetParameter(13,-0.067);
//	 fmm_best_4Poly->SetParLimits(13,-0.085,-0.055);
////mainL
////	 fmm_best_4Poly->SetParameter(14,0.7);
////	 fmm_best_4Poly->SetParLimits(14,0.3,1.5);
////	 fmm_best_4Poly->SetParameter(15,0.004);
////	 fmm_best_4Poly->SetParLimits(15,0.001,0.01);
////	 fmm_best_4Poly->SetParameter(16,0.002);
////	 fmm_best_4Poly->SetParLimits(16,0.,0.01);
////	 fmm_best_4Poly->SetParameter(17,0.0);
////	 fmm_best_4Poly->SetParLimits(17,-0.002,0.002);
//	 hmm_wo_bg_fom_best->Fit("fmm_best_4Poly","","",-0.05,0.12);//Total fitting w/ 4Poly BG
//	 double chisq = fmm_best_4Poly->GetChisquare();
//	 double dof  = fmm_best_4Poly->GetNDF();
//	 cout<<"chisq="<<chisq<<endl;
//	 cout<<"dof="<<dof<<endl;
//	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;
//	 hmm_wo_bg_fom_best->Draw();
//	 fmm_best_4Poly->SetLineColor(kRed);
//	 fL_best->SetLineColor(kGreen);
//	 fS_best->SetLineColor(kGreen);
//	 //fL_best->Draw("same");
//	 //fS_best->Draw("same");
//	 fmm_best_4Poly->Draw("same");
//
//cout << "Well done!" << endl;
}//fit
