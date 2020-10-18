//----------------------------------//
//--  MIXED EVENT ANALYSIS        --//
//----------------------------------//
//
//K. Okuyama (Sep. 23, 2020)
//
//taken over from mea_lcr.C
//using {left, right, center} bunch

void mea_llccrr_new(){
	string pdfname = "temp.pdf";
	string rootname= "bgmea_llccrr_Lsingle.root";
cout << "Output pdf file name is " << pdfname << endl;
cout << "Output root file name is " << rootname << endl;
  
  TFile *file = new TFile("../h2all_Lsingle.root","read");//input file (default: h2all2.root)
 // TFile *file = new TFile("../h2all5.root","read");//input file (default: h2all2.root)
  TFile *file_new = new TFile(rootname.c_str(),"recreate");//new root
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
	
//	string branchname[]={"tr.ntrack_l","tr.ntrack_r","tr.Ls2_pad[100]"};
//	Int_t nbranch = sizeof(branchname)/sizeof(branchname[0]);
//	tree->SetBranchStatus("*",0);
//	for(Int_t ibranch;ibranch<nbranch;ibranch++){
//	tree->SetBranchStatus(branchname[ibranch].c_str(),1);
//	}

	tree->SetBranchStatus("*",0);
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
/*-----------------*/
/*--RESULT OUTPUT--*/
/*-----------------*/
  TH1F* hmm_mixacc_result_best  = new TH1F("hmm_mixacc_result_best","MEA result (best)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut  = new TH1F("hmm_mixacc_result_nocut","MEA result (no Z cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_new  = new TH1F("hmm_mixacc_result_nocut_new","MEA result (no Z, new AC cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_new  = new TH1F("hmm_mixacc_result_new","MEA result (new AC best)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_woac  = new TH1F("hmm_mixacc_result_woac","MEA result (no AC cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_zdiff  = new TH1F("hmm_mixacc_result_zdiff","MEA result (Z diff cut)",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_forAl  = new TH1F("hmm_mixacc_result_nocut_forAl","ACC (Al (best))",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_nocut_new_forAl  = new TH1F("hmm_mixacc_result_nocut_new_forAl","ACC (Al (strict))",xbin,xmin,xmax);
  TH1F* hmm_mixacc_result_momcut  = new TH1F("hmm_mixacc_result_momcut","MEA result (Mom. cut)",xbin,xmin,xmax);
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
  bool cm2_angle1_cut=false;
  bool cm2_angle2_cut=false;
  bool cm3_angle1_cut=false;
  bool cm3_angle2_cut=false;
  bool cm3_angle3_cut=false;
  bool cm4_angle1_cut=false;
  bool cm4_angle2_cut=false;
  bool cm4_angle3_cut=false;
  bool cm4_angle4_cut=false;
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
  tree->Draw(">>elist", "abs(ct_orig+9.0*2.012)<1.006||abs(ct_orig+8.0*2.012)<1.006||abs(ct_orig+2.0*2.012)<1.006||abs(ct_orig+1.0*2.012)<1.006||abs(ct_orig-6.0*2.012)<1.006||abs(ct_orig-7.0*2.012)<1.006");
  TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  ENum = elist->GetN(); 
  //ENum = tree->GetEntries();
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/100;
  if(ENum<100000)time_div=10000;
  int nmix = 1000;//Num of mix
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
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP&&L_tr_p>2.12&&L_tr_p<2.18&&R_tr_p>1.81&&R_tr_p<1.88)event_selection_momcut=true;
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
		if(event_selection)hmm_acc->Fill(mm);
		double before_mm=mm;


		/*--Electron information is saved.--*/
		L_4vec_saved = L_4vec;
		B_4vec_saved = B_4vec;
		T_4vec_saved = T_4vec;
		L_tr_vz_saved  = L_tr_vz;
		/*--Electron information is saved.--*/
	
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

	if(theta_gk_cm*180./PI<8.)cm2_angle1_cut=true;
	if(theta_gk_cm*180./PI>8.)cm2_angle2_cut=true;
	if(theta_gk_cm*180./PI<6.)cm3_angle1_cut=true;
	if(theta_gk_cm*180./PI>6. && theta_gk_cm*180./PI<10.)cm3_angle2_cut=true;
	if(theta_gk_cm*180./PI>10.)cm3_angle3_cut=true;
	if(theta_gk_cm*180./PI<5.)cm4_angle1_cut=true;
	if(theta_gk_cm*180./PI>5. && theta_gk_cm*180./PI<8.)cm4_angle2_cut=true;
	if(theta_gk_cm*180./PI>8. && theta_gk_cm*180./PI<11.)cm4_angle3_cut=true;
	if(theta_gk_cm*180./PI>11.)cm4_angle4_cut=true;

		
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

        double mass_mixed,mm_mixed;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass_mixed = Missing.M();
	    mm_mixed=mass_mixed - mh;//shift by ML
		double after_mm=mm_mixed;
		//cout<<"mm(Before):mm(after)="<<before_mm<<":"<<after_mm<<" (i,j)=("<<i<<","<<j<<")"<<endl;

			if(mix_region1)hm2->Fill(mm_mixed);
			if(event_selection)hmm_mixacc->Fill(mm_mixed);
			if(event_selection)hmm_mixacc2->Fill(mm_mixed);
			if(event_selection)hmm_mixacc3->Fill(mm_mixed);
			if(event_selection)hmm_mixacc_result_best->Fill(mm_mixed);
			if(event_selection_nocut)hmm_mixacc_result_nocut->Fill(mm_mixed);
			if(event_selection_nocut_new)hmm_mixacc_result_nocut_new->Fill(mm_mixed);
			if(event_selection_new)hmm_mixacc_result_new->Fill(mm_mixed);
			if(event_selection_momcut)hmm_mixacc_result_new->Fill(mm_mixed);
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
}//mea_hist
