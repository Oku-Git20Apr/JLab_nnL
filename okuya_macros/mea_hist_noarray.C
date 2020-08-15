//----------------------------------//
//--  MIXED EVENT ANALYSIS        --//
//----------------------------------//
//
//K. Okuyama (Aug. 11, 2020)
//
//add some necessary branch (Aug. 12, 2020)
//complete MEA algorithm (Aug. 12, 2020)
//delete matrix() and Calib() (Aug. 12, 2020)
//
//Converting to a PDF file becomes available. (Aug. 13, 2020)

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


void mea_hist_noarray(){
	string pdfname = "z_MEA.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("h2all4.root","read");//input file (default: h2all2.root)
  TFile *file_new = new TFile("z_MEA.root","recreate");//new root
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
  TH1F* hm2   = (TH1F*)hmm_best->Clone("hm2");
  TH1F* hmm_acc  = new TH1F("hmm_acc","ACC (original)",xbin,xmin,xmax);
  TH1F* hmm_mixacc  = new TH1F("hmm_mixacc","ACC (mixed)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_best  = new TH1F("hmm_mixacc_result_best","ACC (20 bunch mixed)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut  = new TH1F("hmm_mixacc_result_nocut","ACC (20 bunch mixed)",xbin,xmin,xmax);
  TH1F* hm4   = (TH1F*)hmm_best->Clone("hm4");
  
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
  bool event_selection = false;
  bool event_selection_nocut = false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  bool mix_region1 = false;
  bool mix_region2 = false;
  bool mix_region3 = false;
  bool mix_region4 = false;
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

  int ENum=0;
  tree->Draw(">>elist", "abs(ct_orig+2*2.0)<1.0||abs(ct_orig+1.0*2.0)<1.0||abs(ct_orig-3.0*2.0)<1.0||abs(ct_orig-4.0*2.0)<1.0");
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  ENum = elist->GetN(); 
  //ENum = tree->GetEntries();
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;
  int nmix = 5;//Num of mix
  double mass,mm;
  TLorentzVector L_4vec;
  TLorentzVector R_4vec;
  TLorentzVector B_4vec;
  TLorentzVector T_4vec;
  TLorentzVector Missing;
  TLorentzVector L_4vec_saved, B_4vec_saved, T_4vec_saved;

	    //==============================//
	    //======  Initialization  ======//
	    //==============================//
		for(int j=0;j<MAX;j++){
    	}

	time_t start, end;
	start = time(NULL);
	time(&start);

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


	


		if(abs(ct+2.0*rf_bunch)<1.0) mix_region1 = true;
    	else mix_region1 = false;
    	if(abs(ct+1.0*rf_bunch)<1.0) mix_region2 = true;
    	else mix_region2 = false;
    	if(abs(ct-3.0*rf_bunch)<1.0) mix_region3 = true;
    	else mix_region3 = false;
    	if(abs(ct-4.0*rf_bunch)<1.0) mix_region4 = true;
    	else mix_region4 = false;
		if(mix_region1==false&&mix_region2==false&&mix_region3==false&&mix_region4==false)cout<<"Weird!"<<endl;
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

		h_test->Fill(L_E);

		TLorentzVector L_4vec;
		TLorentzVector R_4vec;
		TLorentzVector B_4vec;
		TLorentzVector T_4vec;
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);


        //    TVector3 L_v, R_v, B_v;
	    //B_v.SetXYZ(0.0,0.0,B_p);
	    //L_v.SetXYZ(L_px, L_py, L_pz);
	    //R_v.SetXYZ(R_px, R_py, R_pz);	    
	    R_4vec.RotateX(  13.2/180.*PI );
	    L_4vec.RotateX( -13.2/180.*PI );

	    //======= W/ Matrix & Energy Loss calibraiton ============//
            TVector3 L_vc, R_vc, B_vc;
	    B_vc.SetXYZ(0.0,0.0,B_p);
	    L_vc.SetXYZ(L_px, L_py, L_pz);
	    R_vc.SetXYZ(R_px, R_py, R_pz);
	    R_vc.RotateX(  13.2/180.*PI );
	    L_vc.RotateX( -13.2/180.*PI );
	    double Eec =sqrt(B_p*B_p + Me*Me);
	    double R_Ec =sqrt(R_p*R_p + MK*MK);
	    double L_Ec =sqrt(L_p*L_p + Me*Me);
	    //========================================================//

		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML
		//without Matrix & Energy Loss calibration
		if((mix_region1||mix_region2||mix_region3||mix_region4)&&event_selection)hmm_acc->Fill(mm);




		/*--Electron information is saved.--*/
		L_4vec_saved = L_4vec;
		B_4vec_saved = B_4vec;
		T_4vec_saved = T_4vec;
		L_tr_vz_saved  = L_tr_vz;
		/*--Electron information is saved.--*/
		
		if(mix_region1||mix_region2||mix_region3||mix_region4){
	for(int j=0 ; j<(int)nmix ; j++){
	  
	  int ENum_mixed = i+j;

	  if(ENum_mixed<ENum){
	    //tree->GetEntry(ENum_mixed);
	    tree->GetEntry(elist->GetEntry(ENum_mixed));
	  }
	  else {
	    //tree->GetEntry(ENum_mixed-ENum);
	    tree->GetEntry(elist->GetEntry(ENum_mixed-ENum));
	  }

		//Re-evaluation: Dit this event truely exist in mix_region?
		if(abs(ct+2.0*rf_bunch)<1.0) mix_region1 = true;
    	else mix_region1 = false;
    	if(abs(ct+1.0*rf_bunch)<1.0) mix_region2 = true;
    	else mix_region2 = false;
    	if(abs(ct-3.0*rf_bunch)<1.0) mix_region3 = true;
    	else mix_region3 = false;
    	if(abs(ct-4.0*rf_bunch)<1.0) mix_region4 = true;
    	else mix_region4 = false;
		if(mix_region1==false&&mix_region2==false&&mix_region3==false&&mix_region4==false)cout<<"Weird!"<<endl;

		//Electron info.
		L_tr_vz=L_tr_vz_saved;
	
        R_Tr = R_FP = false;
        // FP and chi2 cuts
        if( R_tr_chi2<0.01 ) R_Tr = true;
        if( R_tr_th<0.17*R_tr_x+0.025
         && R_tr_th>0.17*R_tr_x-0.035
         && R_tr_th<0.40*R_tr_x+0.130 ) R_FP = true;

	if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
	else event_selection=false;
	if(event_selection_nocut){
	// Change Only Kaon information
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
	    double R_px = R_pz * (R_tr_tg_th );
	    double R_py = R_pz * ( R_tr_tg_ph );
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
		R_4vec.SetPxPyPzE(R_E, R_px, R_py, R_pz);
	    R_4vec.RotateX(  13.2/180.*PI );
	

	
	// Electron information is not changed.
		L_4vec = L_4vec_saved;
		B_4vec = B_4vec_saved;
		T_4vec = T_4vec_saved;

        double mass_mixed,mm_mixed;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass_mixed = sqrt( (B_ene + mt - L_ene - R_ene)*(B_ene + mt - L_ene - R_ene)-(B_mom - L_mom - R_mom)*(B_mom - L_mom - R_mom) );
	    mm_mixed=mass_mixed - mh;//shift by ML
		//without Matrix & Energy Loss calibration

		if(mix_region1||mix_region2||mix_region3||mix_region4){
			hmm_mixacc_result_nocut->Fill(mm);
		}
	  }//event_selection_nocut
//Cut condition
	if(event_selection){
	// Change Only Kaon information
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
	    double R_px = R_pz * (R_tr_tg_th );
	    double R_py = R_pz * ( R_tr_tg_ph );
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
		R_4vec.SetPxPyPzE(R_E, R_px, R_py, R_pz);
	    R_4vec.RotateX(  13.2/180.*PI );
	

	
	// Electron information is not changed.
		L_4vec = L_4vec_saved;
		B_4vec = B_4vec_saved;
		T_4vec = T_4vec_saved;

        double mass_mixed,mm_mixed;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass_mixed = sqrt( (B_ene + mt - L_ene - R_ene)*(B_ene + mt - L_ene - R_ene)-(B_mom - L_mom - R_mom)*(B_mom - L_mom - R_mom) );
	    mm_mixed=mass_mixed - mh;//shift by ML
		//without Matrix & Energy Loss calibration

		if(mix_region1)hm2->Fill(mm);
		if(mix_region1||mix_region2||mix_region3||mix_region4){
			hmm_mixacc->Fill(mm);
			hmm_mixacc_result_best->Fill(mm);
		}
	  }//event_selection
	}//Mix (j loop)
	}//mix_region_if

}//ENum (i loop)


/*************************************/
/*Mixed Event Analysis was completed!*/
/*************************************/
//Following code is just for some checks using true-kaon event
//e.g.) filling MM histogram, comparing MM & BG(MEA), etc.
  tree->Draw(">>elist_kaon", "abs(ct_orig)<1.0");
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
  event_selection=false;
  if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
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

		h_test->Fill(L_E);

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

	    //======= W/ Matrix & Energy Loss calibraiton ============//
            TVector3 L_vc, R_vc, B_vc;
	    B_vc.SetXYZ(0.0,0.0,B_p);
	    L_vc.SetXYZ(L_px, L_py, L_pz);
	    R_vc.SetXYZ(R_px, R_py, R_pz);
	    R_vc.RotateX(  13.2/180.*PI );
	    L_vc.RotateX( -13.2/180.*PI );
	    double Eec =sqrt(B_p*B_p + Me*Me);
	    double R_Ec =sqrt(R_p*R_p + MK*MK);
	    double L_Ec =sqrt(L_p*L_p + Me*Me);
	    //========================================================//

		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML
		if(event_selection)hmm_best->Fill(mm);
		if(event_selection_nocut)hmm_nocut->Fill(mm);
  }//ENum_kaon

  double dbin = (xmax-xmin)/(double)xbin;
  double minx=-0.15,maxx=-0.02;
  int fitmin = (minx-xmin)/dbin;
  int fitmax = (maxx-xmin)/dbin;
  double num1 = hmm_best->Integral(fitmin,fitmax);
  double num2 = hmm_mixacc->Integral(fitmin,fitmax);
  double mixscale = num1/num2;
  //cout << fitmin << "," << fitmax << ": "<< num1 << "/" << num2 << "= " << mixscale<< endl;
  cout<<"Information:"<<endl;
  cout<<"Scale adjustment: ["<<minx<<", "<<maxx<<"]"<<endl;
  cout<<"Mixscale="<<num1<<"/"<<num2<<"="<<mixscale<<endl;
  cout<<nmix<<" x 4 bunches"<<"= "<<nmix * 4.0 <<" (" << "effective scale="<<1/mixscale<<")"<<endl;
    
  cout << endl;
  TH1F* hmm_mixacc_clone   = (TH1F*)hmm_mixacc->Clone();
  //TH1F* h2_H_acc2 = (TH1F*)h2_H_acc->Clone();
  for(int i=0 ; i<xbin ; i++){
    hmm_mixacc_clone->SetBinContent( i+1, hmm_mixacc->GetBinContent(i+1)*mixscale);
    hmm_mixacc_clone->SetBinError( i+1, hmm_mixacc->GetBinError(i+1)*mixscale);
  }

  TCanvas* c1 = new TCanvas("c1","c1");
  hmm_best->Draw("");
  TCanvas* c2 = new TCanvas("c2","c2");
  hm2->Draw("");
  TCanvas* c3 = new TCanvas("c3","c3");
  hmm_mixacc->Draw("");
  TCanvas* c4 = new TCanvas("c4","c4");
  h_test->Draw("");
  TCanvas* c5 = new TCanvas("c5","BG + BG(MEA)");
  //h2->Draw();
  double ntemp1 = hmm_acc->GetEntries();
  double ntemp2 = hmm_mixacc->GetEntries();
  //cout << " SCALE: " << ntemp1/ntemp2/5. << endl;
//  hmm_acc->Scale(1./5.);
  hmm_mixacc->Draw("h");
  hmm_mixacc->SetMarkerStyle(1);
  hmm_mixacc->Scale(ntemp1/ntemp2);
cout<<"ntemp1/ntemp2="<<ntemp1<<"/"<<ntemp2<<"="<<ntemp1/ntemp2<<endl;
cout<<"hmm_mixacc was scaled."<<endl;
  hmm_mixacc->SetLineColor(kRed);
  hmm_acc->Draw("same");
  
  TCanvas* c6 = new TCanvas("c6","MM + BG(MEA)");
  hmm_best->Draw();
  hmm_mixacc_clone->SetLineColor(kAzure);
  hmm_mixacc_clone->Draw("same");


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


cout << "creating new rootfile ... "<<endl;
file_new->Write();
//file_new->Close();

cout << "Well done!" << endl;
}//mea_hist
