//-----------------------------//
//--  Cointime Histo.        --//
//-----------------------------//
//
//K. Okuyama (Sep. 27, 2020)

double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }


double fcoin_template( double *x, double *par , int shift, int num)
{
  return par[num] * TMath::Gaus(x[0],par[num+1]-2.0*shift,par[num+2]);//Lambda Gaussian
}

double fcoin_acc( double *x, double *par)
{
	double val=0.;
	for(int i=0;i<=60;i++){
	//val += par[num]*TMath::Gaus(x[0],par[num+1]+par[num+12]*i-20,par[num+2]);
    val += par[0] * (1.-par[4]) * TMath::Gaus(x[0],par[1]+par[3]*i-20.,par[2]);
    val += par[0] * par[4] * TMath::Gaus(x[0],par[1]+par[3]*i-20.,par[5]);
	}
	return val;
}

double fcoin_total( double *x, double *par ){

	return fcoin_acc(x,par)+par[6]*TMath::Gaus(x[0],par[7],par[8])+par[9]*TMath::Gaus(x[0],par[10],par[11])+par[12]*TMath::Gaus(x[0],par[13],par[14])+par[15]*TMath::Gaus(x[0],par[16],par[17]);

}

void hcoin_pi_contami_G(){
//ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
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

  //TH1F* hcoin  = new TH1F("hcoin","",40000/56,-20.,20.);
  TH1F* hcoin  = new TH1F("hcoin","",120000/56,-20.,100.);
  TH1F* hcoin_bg  = new TH1F("hcoin_bg","",120000/56,-20.,100.);
  TH1F* hcoin_wobg  = new TH1F("hcoin_wobg","",120000/56,-20.,100.);
  TH1F* hcoin_new  = new TH1F("hcoin_new","",120000/56,-20.,100.);
  TH1F* hcoin_new_bg  = new TH1F("hcoin_new_bg","",120000/56,-20.,100.);
  TH1F* hcoin_newwobg  = new TH1F("hcoin_newwobg","",120000/56,-20.,100.);
  TH1F* hcoin_pi  = new TH1F("hcoin_pi","",120000/56,-20.,100.);
  TH1F* hcoin_pi_bg  = new TH1F("hcoin_pi_bg","",120000/56,-20.,100.);
  TH1F* hcoin_piwobg  = new TH1F("hcoin_piwobg","",120000/56,-20.,100.);
  TH1F* hcoin_p  = new TH1F("hcoin_p","",120000/56,-20.,100.);
  TH1F* hcoin_p_bg  = new TH1F("hcoin_p_bg","",120000/56,-20.,100.);
  TH1F* hcoin_pwobg  = new TH1F("hcoin_pwobg","",120000/56,-20.,100.);
  TH1F* hcoin_best  = new TH1F("hcoin_best","",120000/56,-20.,100.);
  TH1F* hcoin_best_bg  = new TH1F("hcoin_best_bg","",120000/56,-20.,100.);
  TH1F* hcoin_bestwobg  = new TH1F("hcoin_bestwobg","",120000/56,-20.,100.);
  TH1F* hcoin_strict  = new TH1F("hcoin_strict","",120000/56,-20.,100.);
  TH1F* hcoin_strict_bg  = new TH1F("hcoin_strict_bg","",120000/56,-20.,100.);
  TH1F* hcoin_strictwobg  = new TH1F("hcoin_strictwobg","",120000/56,-20.,100.);
  
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
  //tree->Draw(">>elist" , "fabs(ct_orig)<1.0");//ctsum (does NOT dintinguish #track)
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


	

		bool ac_cut=false;
		bool ac_cut_new=false;
		bool ac_cut_pi=false;
		bool ac_cut_p=false;
		bool best_cut=false;
		bool strict_cut=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.)ac_cut=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.)ac_cut_new=true;
		if(ac1sum>3.75&&ac2sum>10.)ac_cut_pi=true;
		if(ac1sum<3.75&&ac2sum<2.)ac_cut_p=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2)best_cut=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2)strict_cut=true;

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

		if(ac_cut)hcoin->Fill(ct);
		if(ac_cut_new)hcoin_new->Fill(ct);
		if(ac_cut_pi)hcoin_pi->Fill(ct);
		if(ac_cut_p)hcoin_p->Fill(ct);
		if(best_cut)hcoin_best->Fill(ct);
		if(strict_cut)hcoin_strict->Fill(ct);
		if(ac_cut&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_new&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_new_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_pi&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_pi_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_p&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_p_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(best_cut&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_best_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(strict_cut&&20.<ct && ct<100.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_strict_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime

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

//BG subtraction
     hcoin_bg->Scale(40./80.);
     hcoin_new_bg->Scale(40./80.);
     hcoin_pi_bg->Scale(40./80.);
     hcoin_p_bg->Scale(40./80.);
     hcoin_best_bg->Scale(40./80.);
     hcoin_strict_bg->Scale(40./80.);
     hcoin_wobg->Add(hcoin,hcoin_bg,1.0,-1.0);
     hcoin_newwobg->Add(hcoin_new,hcoin_new_bg,1.0,-1.0);
     hcoin_piwobg->Add(hcoin_pi,hcoin_pi_bg,1.0,-1.0);
     hcoin_pwobg->Add(hcoin_p,hcoin_p_bg,1.0,-1.0);
     hcoin_bestwobg->Add(hcoin_best,hcoin_best_bg,1.0,-1.0);
     hcoin_strictwobg->Add(hcoin_strict,hcoin_strict_bg,1.0,-1.0);

//
////Kaon Fitting
//	 double area=fcoin->GetParameter(0);
//	 double location=fcoin->GetParameter(1);
//	 double gsigma=fcoin->GetParameter(2);
//	 double lwidth=fcoin->GetParameter(3);
//	 double bwidth=fcoin->GetParameter(16);
//	 fcoin->FixParameter(0,area);
//	 fcoin->FixParameter(1,location);
//	 fcoin->FixParameter(2,gsigma);
//	 fcoin->FixParameter(3,lwidth);
//	 fcoin->FixParameter(16,bwidth);
//	 fcoin->FixParameter(4,fcoin_first->GetParameter(0));//pi
//	 fcoin->FixParameter(5,fcoin_first->GetParameter(1));
//	 fcoin->FixParameter(6,fcoin_first->GetParameter(2));
//	 fcoin->FixParameter(7,fcoin_first->GetParameter(3));
//	 fcoin->SetParameter(8,2000.);//k
//	 fcoin->SetParLimits(8,0.,100000.);
//	 fcoin->SetParameter(9,0.);//k
//	 fcoin->SetParLimits(9,-0.5,0.5);
//	 fcoin->SetParameter(10,fcoin_first->GetParameter(2));//k G Fix
//	 fcoin->SetParLimits(10,0.,0.6);
//	 fcoin->FixParameter(11,fcoin_first->GetParameter(3));//k L Fix
////L	 //fcoin->SetParLimits(11,0.,1.0);
//	 //hcoin->Fit("fcoin","","",-0.7,0.7);




//Pion Fitting
	 TF1 *fcoin_first= new TF1("fcoin_first","gaus",0.,6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_first->SetNpx(20000);
	 fcoin_first->SetParameter(0,15500.);
	 fcoin_first->SetParameter(1,3.18);
	 fcoin_first->SetParameter(2,0.3);
	 hcoin_piwobg->Fit("fcoin_first","","",2.,4.);
	 double pion_par0 = fcoin_first->GetParameter(0);
	 double pion_par1 = fcoin_first->GetParameter(1);
	 double pion_par2 = fcoin_first->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Pion Fitting after BG subtraction (fcoin_first) "<<endl;
	 double chisq = fcoin_first->GetChisquare();
	 double dof  = fcoin_first->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
	 TF1 *fcoin_firstp= new TF1("fcoin_firstp","gaus",-10,-6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_firstp->SetNpx(20000);
	 fcoin_firstp->SetParameter(0,1550.);
	 fcoin_firstp->SetParameter(1,-8.1);
	 fcoin_firstp->SetParameter(2,0.4);
	 hcoin_pwobg->Fit("fcoin_firstp","","",-10.,-6.);
	 double proton_par0 = fcoin_firstp->GetParameter(0);
	 double proton_par1 = fcoin_firstp->GetParameter(1);
	 double proton_par2 = fcoin_firstp->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Proton Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_firstp->GetChisquare();
	 dof  = fcoin_firstp->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Accidentals Fitting
	 TF1 *fcoin= new TF1("fcoin",fcoin_total,-20.,100.,18);
	 fcoin->SetNpx(20000);
	 fcoin->SetParameter(0,200.);//acc
	 fcoin->SetParameter(1,1.);//
	 fcoin->FixParameter(2,pion_par2);//
	 fcoin->SetParameter(3,2.012);//shift
	 fcoin->SetParameter(4,0.5);//relative strength
	 fcoin->SetParLimits(4,0.,1.);//
	 fcoin->FixParameter(5,proton_par2);//
	 fcoin->FixParameter(6,0.);//pi
	 fcoin->FixParameter(9,0.);//k
	 fcoin->FixParameter(12,0.);//p
	 fcoin->FixParameter(15,0.);//pi2
	 fcoin->SetLineColor(kCyan);
	 hcoin->Fit("fcoin","","",12.,60.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 double chisq2 = fcoin->GetChisquare();
	 double dof2  = fcoin->GetNDF();
	 cout<<"chisq="<<chisq2<<endl;
	 cout<<"dof="<<dof2<<endl;
	 cout<<"Reduced chi-square = "<<chisq2/dof2<<endl;
	 hcoin->Draw("");
	 fcoin->Draw("same");

//Global Fitting
	 fcoin->ReleaseParameter(0);
	 fcoin->ReleaseParameter(1);
	 fcoin->ReleaseParameter(2);
	 fcoin->ReleaseParameter(3);
	 fcoin->ReleaseParameter(4);
	 fcoin->ReleaseParameter(5);
	 fcoin->ReleaseParameter(6);
	 fcoin->ReleaseParameter(7);
	 fcoin->ReleaseParameter(8);
	 fcoin->ReleaseParameter(9);
	 fcoin->ReleaseParameter(10);
	 fcoin->ReleaseParameter(11);
	 fcoin->ReleaseParameter(12);
	 fcoin->ReleaseParameter(13);
	 fcoin->ReleaseParameter(14);
	 fcoin->ReleaseParameter(15);
	 fcoin->ReleaseParameter(16);
	 fcoin->ReleaseParameter(17);
	 fcoin->FixParameter(0,fcoin->GetParameter(0));
	 fcoin->FixParameter(1,fcoin->GetParameter(1));
	 fcoin->FixParameter(2,fcoin->GetParameter(2));
	 fcoin->FixParameter(3,fcoin->GetParameter(3));
	 fcoin->FixParameter(4,fcoin->GetParameter(4));
	 fcoin->FixParameter(5,fcoin->GetParameter(5));
	 fcoin->SetParameter(6,5000.);//pi scale
	 fcoin->SetParLimits(6,0.,100000.);//pi scale
	 fcoin->SetParameter(7,3.18);
	 fcoin->FixParameter(8,pion_par2);
	 fcoin->SetParameter(9,100.);//k scale
	 fcoin->SetParLimits(9,0.,50000.);//k scale
	 fcoin->SetParameter(10,0.);
	 fcoin->SetParLimits(10,-0.5,0.5);
	 fcoin->SetParameter(11,0.5);
	 fcoin->SetParLimits(11,0.,0.7);
	 fcoin->SetParameter(12,2000.);//p scale
	 fcoin->SetParLimits(12,0.,100000.);//p scale
	 fcoin->SetParameter(13,-7.9);
	 fcoin->FixParameter(14,proton_par2);
	 fcoin->SetParameter(15,200.);//pi2 scale
	 fcoin->SetParLimits(15,0.,100000.);//pi2 scale
	 fcoin->FixParameter(16,pion_par1);
	 fcoin->SetParameter(17,0.8);

	 hcoin->Fit("fcoin","","",-20.,20.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 double chisq3 = fcoin->GetChisquare();
	 double dof3  = fcoin->GetNDF();
	 cout<<"chisq="<<chisq3<<endl;
	 cout<<"dof="<<dof3<<endl;
	 cout<<"Reduced chi-square = "<<chisq3/dof3<<endl;
	 TF1 *fcoin_pi= new TF1("fcoin_pi","gaus(0)+gaus(3)",-20.,20.);
	 TF1 *fcoin_k= new TF1("fcoin_k","gaus",-20.,20.);
	 TF1 *fcoin_p= new TF1("fcoin_p","gaus",-20.,20.);
	 fcoin_pi->SetNpx(20000);
	 fcoin_k->SetNpx(20000);
	 fcoin_p->SetNpx(20000);
	 fcoin_pi->SetParameter(0,fcoin->GetParameter(6));
	 fcoin_pi->SetParameter(1,fcoin->GetParameter(7));
	 fcoin_pi->SetParameter(2,fcoin->GetParameter(8));
	 fcoin_pi->SetParameter(3,fcoin->GetParameter(15));
	 fcoin_pi->SetParameter(4,fcoin->GetParameter(16));
	 fcoin_pi->SetParameter(5,fcoin->GetParameter(17));
	 fcoin_k->SetParameter(0,fcoin->GetParameter(9));
	 fcoin_k->SetParameter(1,fcoin->GetParameter(10));
	 fcoin_k->SetParameter(2,fcoin->GetParameter(11));
	 fcoin_p->SetParameter(0,fcoin->GetParameter(12));
	 fcoin_p->SetParameter(1,fcoin->GetParameter(13));
	 fcoin_p->SetParameter(2,fcoin->GetParameter(14));
//
//cout<<"fcoin(4)="<<fcoin->GetParameter(4)<<endl;
//
		TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
		hcoin_piwobg->Draw("");
		fcoin_firstp->Draw("same");
		TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
		c3->SetLogy(0);
		TH1 *frame = c3->DrawFrame(-20.,1.,20.,2000.);
		frame->Draw("");
		hcoin->Draw("same");
		double ktegrated = fcoin_k->Integral(-1.006,1.006);
		cout<<"Kaon (-1.006ns<ct<1.006ns): "<<ktegrated/0.056<<endl;
		double pitegrated = fcoin_pi->Integral(-1.006,1.006);
		cout<<"Pion (-1.006ns<ct<1.006ns): "<<pitegrated/0.056<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated*100./(ktegrated+pitegrated)<<endl;
		fcoin_pi->SetLineColor(kOrange);
		fcoin_k->SetLineColor(kGreen);
		fcoin_p->SetLineColor(kRed);
		fcoin->Draw("same");
		fcoin_pi->Draw("same");
		fcoin_k->Draw("same");
		fcoin_p->Draw("same");
//		
//	cout<<"left bunch = "<<hcoin->Integral(hcoin->FindBin(-17.),hcoin->FindBin(-15.))<<endl;		
//	cout<<"right bunch = "<<hcoin->Integral(hcoin->FindBin(75.),hcoin->FindBin(77.))<<endl;		
//	cout<<"center bunch = "<<hcoin->Integral(hcoin->FindBin(-5.),hcoin->FindBin(-3.))<<endl;		


/*-------------*/
/*   Best Cut  */
/*-------------*/
cout<<"%%%%%%%%%%%%%%%%%%%%"<<endl;
cout<<"%%% Best Cut %%%%%%%"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%"<<endl;


		TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
//Pion Fitting
	 TF1 *fcoin_first_best= new TF1("fcoin_first_best","gaus",0.,6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_first_best->SetNpx(20000);
	 fcoin_first_best->SetParameter(0,15500.);
	 fcoin_first_best->SetParameter(1,3.18);
	 fcoin_first_best->SetParameter(2,0.3);
	 hcoin_piwobg->Fit("fcoin_first_best","","",2.,4.);
	 double pion_best_par0 = fcoin_first_best->GetParameter(0);
	 double pion_best_par1 = fcoin_first_best->GetParameter(1);
	 double pion_best_par2 = fcoin_first_best->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Pion Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_first_best->GetChisquare();
	 dof  = fcoin_first_best->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
	 TF1 *fcoin_firstp_best= new TF1("fcoin_firstp_best","gaus",-10,-6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_firstp_best->SetNpx(20000);
	 fcoin_firstp_best->SetParameter(0,1550.);
	 fcoin_firstp_best->SetParameter(1,-8.1);
	 fcoin_firstp_best->SetParameter(2,0.4);
	 hcoin_pwobg->Fit("fcoin_firstp_best","","",-10.,-6.);
	 double proton_best_par0 = fcoin_firstp_best->GetParameter(0);
	 double proton_best_par1 = fcoin_firstp_best->GetParameter(1);
	 double proton_best_par2 = fcoin_firstp_best->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Proton Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_firstp_best->GetChisquare();
	 dof  = fcoin_firstp_best->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Accidentals Fitting
	 TF1 *fcoin_best= new TF1("fcoin_best",fcoin_total,-20.,100.,18);
	 fcoin_best->SetNpx(20000);
	 fcoin_best->SetParameter(0,200.);//acc
	 fcoin_best->SetParameter(1,1.);//
	 fcoin_best->FixParameter(2,pion_best_par2);//
	 fcoin_best->SetParameter(3,2.012);//shift
	 fcoin_best->SetParameter(4,0.5);//relative strength
	 fcoin_best->SetParLimits(4,0.,1.);//
	 fcoin_best->FixParameter(5,proton_best_par2);//
	 fcoin_best->FixParameter(6,0.);//pi
	 fcoin_best->FixParameter(9,0.);//k
	 fcoin_best->FixParameter(12,0.);//p
	 fcoin_best->FixParameter(15,0.);//pi2
	 fcoin_best->SetLineColor(kCyan);
	 hcoin_best->Fit("fcoin_best","","",12.,60.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq2 = fcoin_best->GetChisquare();
	 dof2  = fcoin_best->GetNDF();
	 cout<<"chisq="<<chisq2<<endl;
	 cout<<"dof="<<dof2<<endl;
	 cout<<"Reduced chi-square = "<<chisq2/dof2<<endl;
	 hcoin_best->Draw("");
	 fcoin_best->Draw("same");

//Global Fitting
	 fcoin_best->ReleaseParameter(0);
	 fcoin_best->ReleaseParameter(1);
	 fcoin_best->ReleaseParameter(2);
	 fcoin_best->ReleaseParameter(3);
	 fcoin_best->ReleaseParameter(4);
	 fcoin_best->ReleaseParameter(5);
	 fcoin_best->ReleaseParameter(6);
	 fcoin_best->ReleaseParameter(7);
	 fcoin_best->ReleaseParameter(8);
	 fcoin_best->ReleaseParameter(9);
	 fcoin_best->ReleaseParameter(10);
	 fcoin_best->ReleaseParameter(11);
	 fcoin_best->ReleaseParameter(12);
	 fcoin_best->ReleaseParameter(13);
	 fcoin_best->ReleaseParameter(14);
	 fcoin_best->ReleaseParameter(15);
	 fcoin_best->ReleaseParameter(16);
	 fcoin_best->ReleaseParameter(17);
	 fcoin_best->FixParameter(0,fcoin_best->GetParameter(0));
	 fcoin_best->FixParameter(1,fcoin_best->GetParameter(1));
	 fcoin_best->FixParameter(2,fcoin_best->GetParameter(2));
	 fcoin_best->FixParameter(3,fcoin_best->GetParameter(3));
	 fcoin_best->FixParameter(4,fcoin_best->GetParameter(4));
	 fcoin_best->FixParameter(5,fcoin_best->GetParameter(5));
	 fcoin_best->SetParameter(6,5000.);//pi scale
	 fcoin_best->SetParLimits(6,0.,100000.);//pi scale
	 fcoin_best->SetParameter(7,3.18);
	 fcoin_best->FixParameter(8,pion_best_par2);
	 fcoin_best->SetParameter(9,100.);//k scale
	 fcoin_best->SetParLimits(9,0.,50000.);//k scale
	 fcoin_best->SetParameter(10,0.);
	 fcoin_best->SetParLimits(10,-0.5,0.5);
	 fcoin_best->SetParameter(11,0.5);
	 fcoin_best->SetParLimits(11,0.,0.7);
	 fcoin_best->SetParameter(12,2000.);//p scale
	 fcoin_best->SetParLimits(12,0.,100000.);//p scale
	 fcoin_best->SetParameter(13,-7.9);
	 fcoin_best->FixParameter(14,proton_best_par2);
	 fcoin_best->SetParameter(15,200.);//pi2 scale
	 fcoin_best->SetParLimits(15,0.,1000.);//pi2 scale
	 fcoin_best->FixParameter(16,pion_best_par1);
	 fcoin_best->SetParameter(17,0.8);

	 hcoin_best->Fit("fcoin_best","","",-20.,20.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq3 = fcoin_best->GetChisquare();
	 dof3  = fcoin_best->GetNDF();
	 cout<<"chisq="<<chisq3<<endl;
	 cout<<"dof="<<dof3<<endl;
	 cout<<"Reduced chi-square = "<<chisq3/dof3<<endl;
	 TF1 *fcoin_pi_best= new TF1("fcoin_pi_best","gaus(0)+gaus(3)",-20.,20.);
	 TF1 *fcoin_k_best= new TF1("fcoin_k_best","gaus",-20.,20.);
	 TF1 *fcoin_p_best= new TF1("fcoin_p_best","gaus",-20.,20.);
	 fcoin_pi_best->SetNpx(20000);
	 fcoin_k_best->SetNpx(20000);
	 fcoin_p_best->SetNpx(20000);
	 fcoin_pi_best->SetParameter(0,fcoin_best->GetParameter(6));
	 fcoin_pi_best->SetParameter(1,fcoin_best->GetParameter(7));
	 fcoin_pi_best->SetParameter(2,fcoin_best->GetParameter(8));
	 fcoin_pi_best->SetParameter(3,fcoin_best->GetParameter(15));
	 fcoin_pi_best->SetParameter(4,fcoin_best->GetParameter(16));
	 fcoin_pi_best->SetParameter(5,fcoin_best->GetParameter(17));
	 fcoin_k_best->SetParameter(0,fcoin_best->GetParameter(9));
	 fcoin_k_best->SetParameter(1,fcoin_best->GetParameter(10));
	 fcoin_k_best->SetParameter(2,fcoin_best->GetParameter(11));
	 fcoin_p_best->SetParameter(0,fcoin_best->GetParameter(12));
	 fcoin_p_best->SetParameter(1,fcoin_best->GetParameter(13));
	 fcoin_p_best->SetParameter(2,fcoin_best->GetParameter(14));
//
//cout<<"fcoin(4)="<<fcoin->GetParameter(4)<<endl;
//
		c4->SetLogy(0);
		TH1 *frame_best = c4->DrawFrame(-20.,1.,20.,2000.);
		frame_best->Draw("");
		hcoin_best->Draw("same");
		fcoin_best->Draw("same");
		double ktegrated_best = fcoin_k_best->Integral(-1.006,1.006);
		cout<<"Kaon (-1.006ns<ct<1.006ns): "<<ktegrated_best/0.056<<endl;
		double pitegrated_best = fcoin_pi_best->Integral(-1.006,1.006);
		cout<<"Pion (-1.006ns<ct<1.006ns): "<<pitegrated_best/0.056<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated_best*100./(ktegrated_best+pitegrated_best)<<endl;
		fcoin_pi_best->SetLineColor(kOrange);
		fcoin_k_best->SetLineColor(kGreen);
		fcoin_p_best->SetLineColor(kRed);
		fcoin_best->Draw("same");
		fcoin_pi_best->Draw("same");
		fcoin_k_best->Draw("same");
		fcoin_p_best->Draw("same");



/*------------------------*/
/*   New AC Cut (AC2<10)  */
/*------------------------*/
cout<<"%%%%%%%%%%%%%%%%%%%%%%"<<endl;
cout<<"%%% New AC Cut %%%%%%%"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%"<<endl;


		TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
//Pion Fitting
	 TF1 *fcoin_first_new= new TF1("fcoin_first_new","gaus",0.,6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_first_new->SetNpx(20000);
	 fcoin_first_new->SetParameter(0,15500.);
	 fcoin_first_new->SetParameter(1,3.18);
	 fcoin_first_new->SetParameter(2,0.3);
	 hcoin_piwobg->Fit("fcoin_first_new","","",2.,4.);
	 double pion_new_par0 = fcoin_first_new->GetParameter(0);
	 double pion_new_par1 = fcoin_first_new->GetParameter(1);
	 double pion_new_par2 = fcoin_first_new->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Pion Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_first_new->GetChisquare();
	 dof  = fcoin_first_new->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
	 TF1 *fcoin_firstp_new= new TF1("fcoin_firstp_new","gaus",-10,-6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_firstp_new->SetNpx(20000);
	 fcoin_firstp_new->SetParameter(0,1550.);
	 fcoin_firstp_new->SetParameter(1,-8.1);
	 fcoin_firstp_new->SetParameter(2,0.4);
	 hcoin_pwobg->Fit("fcoin_firstp_new","","",-10.,-6.);
	 double proton_new_par0 = fcoin_firstp_new->GetParameter(0);
	 double proton_new_par1 = fcoin_firstp_new->GetParameter(1);
	 double proton_new_par2 = fcoin_firstp_new->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Proton Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_firstp_new->GetChisquare();
	 dof  = fcoin_firstp_new->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Accidentals Fitting
	 TF1 *fcoin_new= new TF1("fcoin_new",fcoin_total,-20.,100.,18);
	 fcoin_new->SetNpx(20000);
	 fcoin_new->SetParameter(0,200.);//acc
	 fcoin_new->SetParameter(1,1.);//
	 fcoin_new->FixParameter(2,pion_new_par2);//
	 fcoin_new->SetParameter(3,2.012);//shift
	 fcoin_new->SetParameter(4,0.5);//relative strength
	 fcoin_new->SetParLimits(4,0.,1.);//
	 fcoin_new->FixParameter(5,proton_new_par2);//
	 fcoin_new->FixParameter(6,0.);//pi
	 fcoin_new->FixParameter(9,0.);//k
	 fcoin_new->FixParameter(12,0.);//p
	 fcoin_new->FixParameter(15,0.);//pi2
	 fcoin_new->SetLineColor(kCyan);
	 hcoin_new->Fit("fcoin_new","","",12.,60.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq2 = fcoin_new->GetChisquare();
	 dof2  = fcoin_new->GetNDF();
	 cout<<"chisq="<<chisq2<<endl;
	 cout<<"dof="<<dof2<<endl;
	 cout<<"Reduced chi-square = "<<chisq2/dof2<<endl;
	 hcoin_new->Draw("");
	 fcoin_new->Draw("same");

//Global Fitting
	 fcoin_new->ReleaseParameter(0);
	 fcoin_new->ReleaseParameter(1);
	 fcoin_new->ReleaseParameter(2);
	 fcoin_new->ReleaseParameter(3);
	 fcoin_new->ReleaseParameter(4);
	 fcoin_new->ReleaseParameter(5);
	 fcoin_new->ReleaseParameter(6);
	 fcoin_new->ReleaseParameter(7);
	 fcoin_new->ReleaseParameter(8);
	 fcoin_new->ReleaseParameter(9);
	 fcoin_new->ReleaseParameter(10);
	 fcoin_new->ReleaseParameter(11);
	 fcoin_new->ReleaseParameter(12);
	 fcoin_new->ReleaseParameter(13);
	 fcoin_new->ReleaseParameter(14);
	 fcoin_new->ReleaseParameter(15);
	 fcoin_new->ReleaseParameter(16);
	 fcoin_new->ReleaseParameter(17);
	 fcoin_new->FixParameter(0,fcoin_new->GetParameter(0));
	 fcoin_new->FixParameter(1,fcoin_new->GetParameter(1));
	 fcoin_new->FixParameter(2,fcoin_new->GetParameter(2));
	 fcoin_new->FixParameter(3,fcoin_new->GetParameter(3));
	 fcoin_new->FixParameter(4,fcoin_new->GetParameter(4));
	 fcoin_new->FixParameter(5,fcoin_new->GetParameter(5));
	 fcoin_new->SetParameter(6,5000.);//pi scale
	 fcoin_new->SetParLimits(6,0.,100000.);//pi scale
	 fcoin_new->SetParameter(7,3.18);
	 fcoin_new->FixParameter(8,pion_new_par2);
	 fcoin_new->SetParameter(9,100.);//k scale
	 fcoin_new->SetParLimits(9,0.,50000.);//k scale
	 fcoin_new->SetParameter(10,0.);
	 fcoin_new->SetParLimits(10,-0.5,0.5);
	 fcoin_new->SetParameter(11,0.5);
	 fcoin_new->SetParLimits(11,0.,0.7);
	 fcoin_new->SetParameter(12,2000.);//p scale
	 fcoin_new->SetParLimits(12,0.,100000.);//p scale
	 fcoin_new->SetParameter(13,-7.9);
	 fcoin_new->FixParameter(14,proton_new_par2);
	 fcoin_new->SetParameter(15,200.);//pi2 scale
	 fcoin_new->SetParLimits(15,0.,100000.);//pi2 scale
	 fcoin_new->FixParameter(16,pion_new_par1);
	 fcoin_new->SetParameter(17,0.8);

	 hcoin_new->Fit("fcoin_new","","",-20.,20.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq3 = fcoin_new->GetChisquare();
	 dof3  = fcoin_new->GetNDF();
	 cout<<"chisq="<<chisq3<<endl;
	 cout<<"dof="<<dof3<<endl;
	 cout<<"Reduced chi-square = "<<chisq3/dof3<<endl;
	 TF1 *fcoin_pi_new= new TF1("fcoin_pi_new","gaus(0)+gaus(3)",-20.,20.);
	 TF1 *fcoin_k_new= new TF1("fcoin_k_new","gaus",-20.,20.);
	 TF1 *fcoin_p_new= new TF1("fcoin_p_new","gaus",-20.,20.);
	 fcoin_pi_new->SetNpx(20000);
	 fcoin_k_new->SetNpx(20000);
	 fcoin_p_new->SetNpx(20000);
	 fcoin_pi_new->SetParameter(0,fcoin_new->GetParameter(6));
	 fcoin_pi_new->SetParameter(1,fcoin_new->GetParameter(7));
	 fcoin_pi_new->SetParameter(2,fcoin_new->GetParameter(8));
	 fcoin_pi_new->SetParameter(3,fcoin_new->GetParameter(15));
	 fcoin_pi_new->SetParameter(4,fcoin_new->GetParameter(16));
	 fcoin_pi_new->SetParameter(5,fcoin_new->GetParameter(17));
	 fcoin_k_new->SetParameter(0,fcoin_new->GetParameter(9));
	 fcoin_k_new->SetParameter(1,fcoin_new->GetParameter(10));
	 fcoin_k_new->SetParameter(2,fcoin_new->GetParameter(11));
	 fcoin_p_new->SetParameter(0,fcoin_new->GetParameter(12));
	 fcoin_p_new->SetParameter(1,fcoin_new->GetParameter(13));
	 fcoin_p_new->SetParameter(2,fcoin_new->GetParameter(14));
//
//cout<<"fcoin(4)="<<fcoin->GetParameter(4)<<endl;
//
		c5->SetLogy(0);
		TH1 *frame_new = c5->DrawFrame(-20.,1.,20.,2000.);
		frame_new->Draw("");
		hcoin_new->Draw("same");
		fcoin_new->Draw("same");
		double ktegrated_new = fcoin_k_new->Integral(-1.006,1.006);
		cout<<"Kaon (-1.006ns<ct<1.006ns): "<<ktegrated_new/0.056<<endl;
		double pitegrated_new = fcoin_pi_new->Integral(-1.006,1.006);
		cout<<"Pion (-1.006ns<ct<1.006ns): "<<pitegrated_new/0.056<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated_new*100./(ktegrated_new+pitegrated_new)<<endl;
		fcoin_pi_new->SetLineColor(kOrange);
		fcoin_k_new->SetLineColor(kGreen);
		fcoin_p_new->SetLineColor(kRed);
		fcoin_new->Draw("same");
		fcoin_pi_new->Draw("same");
		fcoin_k_new->Draw("same");
		fcoin_p_new->Draw("same");

/*----------------*/
/*   Strict Cut   */
/*----------------*/
cout<<"%%%%%%%%%%%%%%%%%%%%%%"<<endl;
cout<<"%%% Strict Cut %%%%%%%"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%"<<endl;

		TCanvas *c6 = new TCanvas("c6", "c6", 800, 800);
//Pion Fitting
	 TF1 *fcoin_first_strict= new TF1("fcoin_first_strict","gaus",0.,6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_first_strict->SetNpx(20000);
	 fcoin_first_strict->SetParameter(0,15500.);
	 fcoin_first_strict->SetParameter(1,3.18);
	 fcoin_first_strict->SetParameter(2,0.3);
	 hcoin_piwobg->Fit("fcoin_first_strict","","",2.,4.);
	 double pion_strict_par0 = fcoin_first_strict->GetParameter(0);
	 double pion_strict_par1 = fcoin_first_strict->GetParameter(1);
	 double pion_strict_par2 = fcoin_first_strict->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Pion Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_first_strict->GetChisquare();
	 dof  = fcoin_first_strict->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
	 TF1 *fcoin_firstp_strict= new TF1("fcoin_firstp_strict","gaus",-10,-6.);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
	 fcoin_firstp_strict->SetNpx(20000);
	 fcoin_firstp_strict->SetParameter(0,1550.);
	 fcoin_firstp_strict->SetParameter(1,-8.1);
	 fcoin_firstp_strict->SetParameter(2,0.4);
	 hcoin_pwobg->Fit("fcoin_firstp_strict","","",-10.,-6.);
	 double proton_strict_par0 = fcoin_firstp_strict->GetParameter(0);
	 double proton_strict_par1 = fcoin_firstp_strict->GetParameter(1);
	 double proton_strict_par2 = fcoin_firstp_strict->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Proton Fitting after BG subtraction (fcoin_first) "<<endl;
	 chisq = fcoin_firstp_strict->GetChisquare();
	 dof  = fcoin_firstp_strict->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Accidentals Fitting
	 TF1 *fcoin_strict= new TF1("fcoin_strict",fcoin_total,-20.,100.,18);
	 fcoin_strict->SetNpx(20000);
	 fcoin_strict->SetParameter(0,200.);//acc
	 fcoin_strict->SetParameter(1,1.);//
	 fcoin_strict->FixParameter(2,pion_strict_par2);//
	 fcoin_strict->SetParameter(3,2.012);//shift
	 fcoin_strict->SetParameter(4,0.5);//relative strength
	 fcoin_strict->SetParLimits(4,0.,1.);//
	 fcoin_strict->FixParameter(5,proton_strict_par2);//
	 fcoin_strict->FixParameter(6,0.);//pi
	 fcoin_strict->FixParameter(9,0.);//k
	 fcoin_strict->FixParameter(12,0.);//p
	 fcoin_strict->FixParameter(15,0.);//pi2
	 fcoin_strict->SetLineColor(kCyan);
	 hcoin_strict->Fit("fcoin_strict","","",12.,60.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq2 = fcoin_strict->GetChisquare();
	 dof2  = fcoin_strict->GetNDF();
	 cout<<"chisq="<<chisq2<<endl;
	 cout<<"dof="<<dof2<<endl;
	 cout<<"Reduced chi-square = "<<chisq2/dof2<<endl;
	 hcoin_strict->Draw("");
	 fcoin_strict->Draw("same");

//Global Fitting
	 fcoin_strict->ReleaseParameter(0);
	 fcoin_strict->ReleaseParameter(1);
	 fcoin_strict->ReleaseParameter(2);
	 fcoin_strict->ReleaseParameter(3);
	 fcoin_strict->ReleaseParameter(4);
	 fcoin_strict->ReleaseParameter(5);
	 fcoin_strict->ReleaseParameter(6);
	 fcoin_strict->ReleaseParameter(7);
	 fcoin_strict->ReleaseParameter(8);
	 fcoin_strict->ReleaseParameter(9);
	 fcoin_strict->ReleaseParameter(10);
	 fcoin_strict->ReleaseParameter(11);
	 fcoin_strict->ReleaseParameter(12);
	 fcoin_strict->ReleaseParameter(13);
	 fcoin_strict->ReleaseParameter(14);
	 fcoin_strict->ReleaseParameter(15);
	 fcoin_strict->ReleaseParameter(16);
	 fcoin_strict->ReleaseParameter(17);
	 fcoin_strict->FixParameter(0,fcoin_strict->GetParameter(0));
	 fcoin_strict->FixParameter(1,fcoin_strict->GetParameter(1));
	 fcoin_strict->FixParameter(2,fcoin_strict->GetParameter(2));
	 fcoin_strict->FixParameter(3,fcoin_strict->GetParameter(3));
	 fcoin_strict->FixParameter(4,fcoin_strict->GetParameter(4));
	 fcoin_strict->FixParameter(5,fcoin_strict->GetParameter(5));
	 fcoin_strict->SetParameter(6,5000.);//pi scale
	 fcoin_strict->SetParLimits(6,0.,100000.);//pi scale
	 fcoin_strict->SetParameter(7,3.18);
	 fcoin_strict->FixParameter(8,pion_strict_par2);
	 fcoin_strict->SetParameter(9,100.);//k scale
	 fcoin_strict->SetParLimits(9,0.,50000.);//k scale
	 fcoin_strict->SetParameter(10,0.);
	 fcoin_strict->SetParLimits(10,-0.5,0.5);
	 fcoin_strict->SetParameter(11,0.5);
	 fcoin_strict->SetParLimits(11,0.,0.7);
	 fcoin_strict->SetParameter(12,2000.);//p scale
	 fcoin_strict->SetParLimits(12,0.,100000.);//p scale
	 fcoin_strict->SetParameter(13,-7.9);
	 fcoin_strict->FixParameter(14,proton_strict_par2);
	 fcoin_strict->SetParameter(15,200.);//pi2 scale
	 fcoin_strict->SetParLimits(15,0.,100000.);//pi2 scale
	 fcoin_strict->FixParameter(16,pion_strict_par1);
	 fcoin_strict->SetParameter(17,0.8);

	 hcoin_strict->Fit("fcoin_strict","","",-20.,20.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [12,60] (fcoin) "<<endl;
	 chisq3 = fcoin_strict->GetChisquare();
	 dof3  = fcoin_strict->GetNDF();
	 cout<<"chisq="<<chisq3<<endl;
	 cout<<"dof="<<dof3<<endl;
	 cout<<"Reduced chi-square = "<<chisq3/dof3<<endl;
	 TF1 *fcoin_pi_strict= new TF1("fcoin_pi_strict","gaus(0)+gaus(3)",-20.,20.);
	 TF1 *fcoin_k_strict= new TF1("fcoin_k_strict","gaus",-20.,20.);
	 TF1 *fcoin_p_strict= new TF1("fcoin_p_strict","gaus",-20.,20.);
	 fcoin_pi_strict->SetNpx(20000);
	 fcoin_k_strict->SetNpx(20000);
	 fcoin_p_strict->SetNpx(20000);
	 fcoin_pi_strict->SetParameter(0,fcoin_strict->GetParameter(6));
	 fcoin_pi_strict->SetParameter(1,fcoin_strict->GetParameter(7));
	 fcoin_pi_strict->SetParameter(2,fcoin_strict->GetParameter(8));
	 fcoin_pi_strict->SetParameter(3,fcoin_strict->GetParameter(15));
	 fcoin_pi_strict->SetParameter(4,fcoin_strict->GetParameter(16));
	 fcoin_pi_strict->SetParameter(5,fcoin_strict->GetParameter(17));
	 fcoin_k_strict->SetParameter(0,fcoin_strict->GetParameter(9));
	 fcoin_k_strict->SetParameter(1,fcoin_strict->GetParameter(10));
	 fcoin_k_strict->SetParameter(2,fcoin_strict->GetParameter(11));
	 fcoin_p_strict->SetParameter(0,fcoin_strict->GetParameter(12));
	 fcoin_p_strict->SetParameter(1,fcoin_strict->GetParameter(13));
	 fcoin_p_strict->SetParameter(2,fcoin_strict->GetParameter(14));
//
//cout<<"fcoin(4)="<<fcoin->GetParameter(4)<<endl;
//
		c6->SetLogy(0);
		TH1 *frame_strict = c6->DrawFrame(-20.,1.,20.,2000.);
		frame_strict->Draw("");
		hcoin_strict->Draw("same");
		fcoin_strict->Draw("same");
		double ktegrated_strict = fcoin_k_strict->Integral(-1.006,1.006);
		cout<<"Kaon (-1.006ns<ct<1.006ns): "<<ktegrated_strict/0.056<<endl;
		double pitegrated_strict = fcoin_pi_strict->Integral(-1.006,1.006);
		cout<<"Pion (-1.006ns<ct<1.006ns): "<<pitegrated_strict/0.056<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated_strict*100./(ktegrated_strict+pitegrated_strict)<<endl;
		fcoin_pi_strict->SetLineColor(kOrange);
		fcoin_k_strict->SetLineColor(kGreen);
		fcoin_p_strict->SetLineColor(kRed);
		fcoin_strict->Draw("same");
		fcoin_pi_strict->Draw("same");
		fcoin_k_strict->Draw("same");
		fcoin_p_strict->Draw("same");

}//fit
