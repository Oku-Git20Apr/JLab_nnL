//---------------//
//--  Fitting  --//
//---------------//
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

double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

void fit(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("h2_10run.root","read");
  TFile *file = new TFile("h2all3.root","read");//input file of all H2 run(default: h2all3.root)
  TFile *file_mea = new TFile("z_MEA.root","read");//input file of BG(MEA) histo.(default: bgmea.root)
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
 const double def_n_L=20.; 
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_n_S=6.; 
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
 double n_pi_noZ, n_k_noZ, n_p_noZ, n_L_noZ, n_S_noZ;
 double const_pi_noZ, const_k_noZ, const_p_noZ, const_L_noZ, const_S_noZ;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_noZ, mean_k_noZ, mean_p_noZ, mean_L_noZ, mean_S_noZ;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_noZ, sig_k_noZ, sig_p_noZ, sig_L_noZ, sig_S_noZ;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];


//---Fitting Function---//
 TF1* fpi_noZ;
 TF1* fk_noZ;
 TF1* fp_noZ;
 TF1* fmmbg_noZ;
 TF1* fmm_pi_noZ;
 TF1* fL_noZ;
 TF1* fS_noZ;
 TF1* fmm_noZ;
 TF1* fcoin_noZ;
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

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct[100][100], ct_eff;
		for(int i=0;i<100;i++){
			for(int j=0;j<100;j++){
				ct[i][j]=-1000.;
				}
		}
  double ct_sum;

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

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("tr.ntrack_l",1);tree->SetBranchAddress("tr.ntrack_l",&NLtr);
	tree->SetBranchStatus("tr.ntrack_r",1);tree->SetBranchAddress("tr.ntrack_r",&NRtr);
	
  	tree->SetBranchStatus("ac1_npe_sum",1);  tree->SetBranchAddress("ac1_npe_sum", &ac1sum);
  	tree->SetBranchStatus("ac2_npe_sum",1);  tree->SetBranchAddress("ac2_npe_sum", &ac2sum);
  	tree->SetBranchStatus("Lp_c",1);  tree->SetBranchAddress("Lp_c", L_mom);
  	tree->SetBranchStatus("Rp_c",1);  tree->SetBranchAddress("Rp_c", R_mom);
  	tree->SetBranchStatus("Bp_c",1);  tree->SetBranchAddress("Bp_c", &B_mom);
  	tree->SetBranchStatus("ct_orig",1);  tree->SetBranchAddress("ct_orig", ct);
  	tree->SetBranchStatus("ct",1);  tree->SetBranchAddress("ct", &ct_sum);

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
  double xmin = -0.1, xmax = 0.2; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_L_fom_best  = new TH1F("hmm_L_fom_best","hmm_L_fom_best",xbin,xmin,xmax);
  hmm_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_L_fom_best->GetYaxis()->SetTitle("Counts / MeV");
  hmm_L_fom_best->SetLineColor(1);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin,xmin,xmax);
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
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum

  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  tree->Draw(">>elist" , "fabs(ct)<1.0||tr.ntrack_l>1||tr.ntrack_r>1");//ctsum (does NOT dintinguish #track)
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
	tree->GetEntry(elist->GetEntry(i));

    if(i%time_div==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (i+1) - diff;
      cout<<i<<" / "<<ENum<<" ("<<i*100/ENum<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }

      
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


	


		if(fabs(ct[lt][rt])<1)ct_cut=true;
		else ct_cut=false;
		//if(fabs(L_tr_vz[lt]-R_tr_vz[rt])<0.025&&fabs(R_tr_vz[rt]+L_tr_vz[lt])<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz[lt]-R_tr_vz[rt])<0.025&&fabs(R_tr_vz[rt]+L_tr_vz[lt])<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
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



		}//NRtr
	}//NLtr
}//ENum

	TCanvas* c1 = new TCanvas("c1","c1");
	hmm_L_fom_best->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2");
	TH1F* hmm_pi_fom_noZ=(TH1F*)file->Get("hmm_pi_fom_noZ");
	TH1F* hmm_pi_fom_best=(TH1F*)file->Get("hmm_pi_fom_best");
	//ACCBGの引き算はmea_hist.rootから
	TH1F* hmm_bg_fom_best=(TH1F*)file_mea->Get("hmm_mixacc_result");
    int fitmin = hmm_L_fom_best->FindBin(0.10);
    int fitmax = hmm_L_fom_best->FindBin(0.15);
    double num1 = hmm_L_fom_best->Integral(fitmin,fitmax);
    double num2 = hmm_bg_fom_best->Integral(fitmin,fitmax);
    double mixscale = num1/num2;
	cout<<"hmm_L integral ="<<num1<<endl;
	cout<<"hmm_bg integral ="<<num2<<endl;
	cout<<"mixscale(original/mixed)="<<mixscale<<endl;
	hmm_bg_fom_best->Scale(1./15.);
	//TH1F* hmm_wo_bg_fom_best = (TH1F*)hmm_L_fom_best->Clone("hmm_wo_bg_fom_best");
	hmm_wo_bg_fom_best->Add(hmm_L_fom_best,hmm_bg_fom_best,1.0,-1.0);
	fmmbg_noZ=new TF1("fmmbg_noZ","pol4",min_mm,max_mm);
	fmmbg_noZ->SetNpx(2000);
	fmm_pi_noZ=new TF1("fmm_pi_noZ",F_Voigt,min_mm,max_mm,4);
	fmm_pi_noZ->SetNpx(2000);
	 fL_noZ=new TF1("fL_noZ","gausn(0)",min_mm,max_mm);
	 fL_noZ->SetNpx(2000);
	 fL_noZ->SetParameters(def_n_L,def_mean_L,def_sig_L);
	 fL_noZ->SetParLimits(0,0.,100000);
	 fL_noZ->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fL_noZ->SetParLimits(2,0.,0.01);
	 fS_noZ=new TF1("fS_noZ","gausn(0)",min_mm,max_mm);
	 fS_noZ->SetNpx(2000);
	 fS_noZ->SetParameters(def_n_S,def_mean_S,def_sig_S);
	 fS_noZ->SetParLimits(0,0.,100000);
	 fS_noZ->SetParLimits(1,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fS_noZ->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
	 
	 fmm_noZ=new TF1("fmm_noZ","gausn(0)+gausn(3)+pol4(6)",min_mm,max_mm);
	 fmm_noZ->SetNpx(2000);
	 fmm_noZ->SetTitle("Missing Mass (best cut)");
	 fmm_noZ->SetParLimits(0,0.,1000000.);//positive
	 fmm_noZ->SetParLimits(3,0.,1000000.);//positive
		
	 hmm_wo_bg_fom_best->Fit("fL_noZ","","",def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 const_L_noZ=fL_noZ->GetParameter(0);
	 mean_L_noZ=fL_noZ->GetParameter(1);
	 sig_L_noZ=fL_noZ->GetParameter(2);
	 center_L=def_mean_L;
	 range_L=2*def_sig_L;
	
	 hmm_wo_bg_fom_best->Fit("fS_noZ","","",def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 const_S_noZ=fS_noZ->GetParameter(0);
	 mean_S_noZ=fS_noZ->GetParameter(1);
	 sig_S_noZ=fS_noZ->GetParameter(2);
	 center_S=def_mean_S;
	 range_S=2*def_sig_S;

cout<<"fmm_pi_noZ fit start"<<endl;
	fmm_pi_noZ->SetParameters(50,0.05,0.05,0.00006);
	hmm_pi_fom_noZ->Fit("fmm_pi_noZ","","",min_mm,max_mm);//1st Fit
	double fmmpi1=fmm_pi_noZ->GetParameter(1);
	double fmmpi2=fmm_pi_noZ->GetParameter(2);
	double fmmpi3=fmm_pi_noZ->GetParameter(3);
	fmm_pi_noZ->FixParameter(1,fmmpi1);
	fmm_pi_noZ->FixParameter(2,fmmpi2);
	fmm_pi_noZ->FixParameter(3,fmmpi3);
	fmm_pi_noZ->SetParLimits(0,0.,1000000);//positive
	hmm_wo_bg_fom_best->Fit("fmm_pi_noZ","","",-0.1,-0.02);//2nd Fit

	 fmm_noZ->SetParameter(0,const_L_noZ);
	 fmm_noZ->SetParameter(1,mean_L_noZ);
	 fmm_noZ->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
	 fmm_noZ->SetParameter(2,sig_L_noZ);
	 fmm_noZ->SetParLimits(2,0.,2*def_sig_L);
	// fmm_noZ->SetParameters(9,100);
	 fmm_noZ->SetParameter(3,const_S_noZ);
	 fmm_noZ->SetParameter(4,mean_S_noZ);
	 fmm_noZ->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
	 fmm_noZ->SetParameter(5,sig_S_noZ);
	 fmm_noZ->SetParLimits(5,0.,2*def_sig_S);
	 hmm_wo_bg_fom_best->Fit("fmm_noZ","Rq0","0",-0.05,0.15);
	double fmm_noZpar0 = fmm_noZ->GetParameter(0);cout<<"fmm_noZ[0]="<<fmm_noZpar0<<endl;//area(L)
	double fmm_noZpar1 = fmm_noZ->GetParameter(1);cout<<"fmm_noZ[1]="<<fmm_noZpar1<<endl;//mean(L)
	double fmm_noZpar2 = fmm_noZ->GetParameter(2);cout<<"fmm_noZ[2]="<<fmm_noZpar2<<endl;//sigma(L)
	double fmm_noZpar3 = fmm_noZ->GetParameter(3);cout<<"fmm_noZ[3]="<<fmm_noZpar3<<endl;//area(S)
	double fmm_noZpar4 = fmm_noZ->GetParameter(4);cout<<"fmm_noZ[4]="<<fmm_noZpar4<<endl;//mean(S)
	double fmm_noZpar5 = fmm_noZ->GetParameter(5);cout<<"fmm_noZ[5]="<<fmm_noZpar5<<endl;//sigma(S)
	double fmm_noZpar6 = fmm_noZ->GetParameter(6);cout<<"fmm_noZ[6]="<<fmm_noZpar6<<endl;//poly_const
	double fmm_noZpar7 = fmm_noZ->GetParameter(7);cout<<"fmm_noZ[7]="<<fmm_noZpar7<<endl;//poly_x
	double fmm_noZpar8 = fmm_noZ->GetParameter(8);cout<<"fmm_noZ[8]="<<fmm_noZpar8<<endl;//poly_x^2
	double fmm_noZpar9 = fmm_noZ->GetParameter(9);cout<<"fmm_noZ[9]="<<fmm_noZpar9<<endl;//poly_x^3
	double fmm_noZpar10 = fmm_noZ->GetParameter(10);cout<<"fmm_noZ[10]="<<fmm_noZpar10<<endl;//poly_x^4
	 fmmbg_noZ->SetParameters(fmm_noZpar6,fmm_noZpar7,fmm_noZpar8,fmm_noZpar9,fmm_noZpar10);
	
	 mean_L_noZ=def_mean_L;
	 mean_S_noZ=def_mean_S;
	 sig_L_noZ=def_sig_L;
	 sig_S_noZ=def_sig_S;
	 n_L_noZ=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_L-range_L),hmm_wo_bg_fom_best->FindBin(center_L+range_L));
	cout<<"before(L):: "<<n_L_noZ<<endl;
	 double integralL=fmmbg_noZ->Integral(center_L-range_L,center_L+range_L);
	 integralL=integralL/(2*range_L/(hmm_wo_bg_fom_best->FindBin(center_L+range_L)-hmm_wo_bg_fom_best->FindBin(center_L-range_L)));
	cout<<"integralL="<<integralL<<endl;
	 if(integralL>0)n_L_noZ=n_L_noZ-integralL;
	cout<<"after(L):: "<<n_L_noZ<<endl;
	 n_S_noZ=hmm_wo_bg_fom_best->Integral(hmm_wo_bg_fom_best->FindBin(center_S-range_S),hmm_wo_bg_fom_best->FindBin(center_S+range_S));
	cout<<"before(S):: "<<n_S_noZ<<endl;
	 double integralS=fmmbg_noZ->Integral(center_S-range_S,center_S+range_S);
	 integralS=integralS/(2*range_S/(hmm_wo_bg_fom_best->FindBin(center_S+range_S)-hmm_wo_bg_fom_best->FindBin(center_S-range_S)));
	cout<<"integralS="<<integralS<<endl;
	 if(integralS>0)n_S_noZ-=integralS;
	cout<<"after(S):: "<<n_S_noZ<<endl;
	 cout<<"n_L"<<n_L_noZ<<endl;
	cout<<"mean_L"<<mean_L_noZ<<endl;
	cout<<"sig_L"<<sig_L_noZ<<endl;
	 cout<<"n_S"<<n_S_noZ<<endl;
	cout<<"mean_S"<<mean_S_noZ<<endl;
	cout<<"sig_S"<<sig_S_noZ<<endl;
	
	hmm_L_fom_best->Draw();
	hmm_bg_fom_best->Draw("same");
	
		TCanvas* c3 = new TCanvas("c3","c3");
	fmmbg_noZ->SetLineColor(kGreen);
	fL_noZ->SetLineColor(kRed);
	fS_noZ->SetLineColor(kRed);
	fL_noZ->Draw("");
	fS_noZ->Draw("same");
	fmmbg_noZ->Draw("same");
		TCanvas* c4 = new TCanvas("c4","c4");
	hmm_wo_bg_fom_best->Draw("");
	fmmbg_noZ->Draw("same");
	fL_noZ->Draw("same");
	fS_noZ->Draw("same");
		TCanvas* c5 = new TCanvas("c5","c5");
	hmm_wo_bg_fom_best->Draw("");
	fmm_noZ->Draw("same");
		TCanvas* c6 = new TCanvas("c6","c6");
	//hmm_pi_fom_noZ->Draw("");
	hmm_wo_bg_fom_best->Draw("");
	fmm_pi_noZ->SetLineColor(kOrange);
	fmm_pi_noZ->Draw("same");
		TCanvas* c7 = new TCanvas("c7","c7");
	hmm_pi_fom_noZ->Draw("");
	fmm_pi_noZ->Draw("same");

/*--- Print ---*/
cout << "Print is starting" << endl;
	c1->Print(Form("%s[",pdfname.c_str()));
	c1->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s",pdfname.c_str()));
	c3->Print(Form("%s",pdfname.c_str()));
	c4->Print(Form("%s",pdfname.c_str()));
	c5->Print(Form("%s",pdfname.c_str()));
	c6->Print(Form("%s",pdfname.c_str()));
	c6->Print(Form("%s]",pdfname.c_str()));


cout << "Well done!" << endl;
}//fit
