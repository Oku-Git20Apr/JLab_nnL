//--  Kaon's survival ratio   --//
//
//K. Okuyama (Nov. 21, 2020)
//K. Okuyama (May. 14, 2020)//Kaon SR event by event re-try
//
//This is taken over from MM.C
//No array branch mode 
//


void kaon_SR(){

  //TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run
  //TTree *tree = (TTree*)file->Get("tree_out");
  TFile *file = new TFile("/data/41a/ELS/JLab/E12-17-003/root/tritium_111175_1.root","read");//TEST
  TTree *tree = (TTree*)file->Get("T");


	//ACCBGの引き算はmea_hist.ccから
  //TFile *file_mea = new TFile("./MixedEventAnalysis/bgmea6.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  TFile *file_mea = new TFile("../MixedEventAnalysis/bgmea_2020Nov.root","read");//input file of BG(MEA) histo.(default: bgmea3.root)
  double nbunch = 6000.;//effetive bunches (6 bunches x 5 mixtures)


    
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
	double R_s2_trpath, R_tr_pathl;
	double L_s2_trpath, L_tr_pathl;

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
//ADD for kaon_SR.C
    tree->SetBranchStatus("R.s2.trpath",1);  tree->SetBranchAddress("R.s2.trpath", &R_s2_trpath);
    tree->SetBranchStatus("R.tr.pathl",1);  tree->SetBranchAddress("R.tr.pathl", &R_tr_pathl);
    tree->SetBranchStatus("L.s2.trpath",1);  tree->SetBranchAddress("L.s2.trpath", &L_s2_trpath);
    tree->SetBranchStatus("L.tr.pathl",1);  tree->SetBranchAddress("L.tr.pathl", &L_tr_pathl);




  TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  h1->GetXaxis()->SetTitle("coin time (ns)");
  h1->GetYaxis()->SetTitle("Counts / 100 ps");
  h1->GetXaxis()->SetRangeUser(-14.0,17.);
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);
  TH1F* hR_T2S2  = new TH1F("hR_T2S2","hR_T2S2",1000,20.,30.);
  TH1F* hL_T2S2  = new TH1F("hL_T2S2","hL_T2S2",1000,20.,30.);
  TH1F* hR_T2FP  = new TH1F("hR_T2FP","hR_T2FP",1000,20.,30.);
  TH1F* hL_T2FP  = new TH1F("hL_T2FP","hL_T2FP",1000,20.,30.);
  TH1F* hR2_T2S2  = new TH1F("hR2_T2S2","hR2_T2S2",1000,20.,30.);
  TH1F* hL2_T2S2  = new TH1F("hL2_T2S2","hL2_T2S2",1000,20.,30.);
  TH1F* hR2_T2FP  = new TH1F("hR2_T2FP","hR2_T2FP",1000,20.,30.);
  TH1F* hL2_T2FP  = new TH1F("hL2_T2FP","hL2_T2FP",1000,20.,30.);
  TH1F* hR_mom  = new TH1F("hR_mom","hR_mom",1000,1.7,2.0);
  TH1F* hR2_mom  = new TH1F("hR2_mom","hR2_mom",1000,1.7,2.0);
  TH2F* hR_T2S2_mom  = new TH2F("hR_T2S2_mom","hR_T2S2_mom",100,25.,27.,100,1.7,2.);
  TH2F* hL_T2S2_mom  = new TH2F("hL_T2S2_mom","hL_T2S2_mom",100,25.,27.,100,1.9,2.2);
  TH2F* hR2_T2S2_mom  = new TH2F("hR2_T2S2_mom","hR2_T2S2_mom",100,25.,27.,100,1.76,1.9);
  TH2F* hL2_T2S2_mom  = new TH2F("hL2_T2S2_mom","hL2_T2S2_mom",100,25.,27.,100,1.9,2.2);

  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  TH2F* hR_mom_x  = new TH2F("hR_mom_x","hR_mom_x;Momentum [GeV/c];X(FP) [m]",100,1.76,1.9,100,-1.0,1.0);
  TH2F* hR_mom_th  = new TH2F("hR_mom_th","hR_mom_th;Momentum [GeV/c];X'(FP) [rad]",100,1.76,1.9,100,-0.15,0.15);
  TH2F* hR_len_x  = new TH2F("hR_len_x","hR_len_x;Path Length (T to S2) [m];X(FP) [m]",100,25.,27.,100,-1.0,1.0);
  TH2F* hR_len_th  = new TH2F("hR_len_th","hR_len_th;Path Length (T to S2) [m];X'(FP) [rad]",100,25.,27.,100,-0.15,0.15);
  TH2F* hR_lenfp_x  = new TH2F("hR_lenfp_x","hR_lenfp_x;Path Length (T to FP) [m];X(FP) [m]",100,21.,25.,100,-1.0,1.0);
  TH2F* hR_lenfp_th  = new TH2F("hR_lenfp_th","hR_lenfp_th;Path Length (T to FP) [m];X'(FP) [rad]",100,21.,25.,100,-0.15,0.15);

  TH2F* hR_mom_y  = new TH2F("hR_mom_y","hR_mom_y;Momentum [GeV/c];Y(FP) [m]",100,1.76,1.9,100,-0.1,0.1);
  TH2F* hR_mom_ph  = new TH2F("hR_mom_ph","hR_mom_ph;Momentum [GeV/c];Y'(FP) [rad]",100,1.76,1.9,100,-0.15,0.15);
  TH2F* hR_len_y  = new TH2F("hR_len_y","hR_len_y;Path Length (T to S2) [m];Y(FP) [m]",100,25.,27.,100,-0.1,0.1);
  TH2F* hR_len_ph  = new TH2F("hR_len_ph","hR_len_ph;Path Length (T to S2) [m];Y'(FP) [rad]",100,25.,27.,100,-0.15,0.15);

  TH2F* h2_kSR  = new TH2F("h2_kSR","h2_kSR",100,25.,27.,100,1.7,2.);
  int srbin = 80;
  TH1D* h_kSR  = new TH1D("h_kSR","h_kSR",srbin,1.65,2.);
  double kaonsr[srbin];
  double ave[srbin];//Ensemble average: \bar{X} 
  double dif[srbin];// (X_i - \bar{X})^2
  //double sigma[srbin];//sigma of ave. \sigma = \sqrt{\frac{\sum_i^n (X_i-\bar{X})^2}{n(n-1)}}
  int counts[srbin];
  for(int i=0;i<srbin;i++){//Initialization
	kaonsr[i]=0.;
	ave[i]=0.;
	counts[i]=0;
	dif[i]=0.;
  }

//---  Loading... ---//average of ksr
	string daq_file = "../information/ksr_tmp.dat";//input file
	int binnum;//running bin
	int nbin;//total # of bin
	double ave_tmp;
	double sig_tmp;
	for(int i=0;i<srbin;i++){
		ave[i]=0.;
	}

	string buf;
	ifstream ifp(daq_file.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << daq_file.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(buf=="Number of Bin"){
			getline(ifp,buf);
			stringstream sbinbuf(buf);
			sbinbuf >> nbin;
			//cout << nbin <<endl;
			if(nbin==srbin)continue;
			else{
				cout<<"Number of Bins is different!"<<endl;
				cout<<"srbin should be set "<<nbin<<endl;
				break;
			}
		}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> binnum >> ave_tmp;
		ave[binnum] = ave_tmp;
		cout << binnum<< ", " << ave[binnum] <<endl;
	}
	cout<<"Average of SR was successfully loaded from external file."<<endl;
	


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


	


		if(fabs(ct)<1.006)ct_cut=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		//if(ct_cut){
			hR_mom->Fill(R_tr_p);
			hR_T2S2->Fill(R_tr_pathl);
			hL_T2S2->Fill(L_tr_pathl);
			hR_T2FP->Fill(R_tr_pathl-R_s2_trpath);
			hL_T2FP->Fill(L_tr_pathl-L_s2_trpath);
			hL_T2S2_mom->Fill(L_tr_pathl,L_tr_p);
			hR_T2S2_mom->Fill(R_tr_pathl,R_tr_p);
			hR_mom_x->Fill(R_tr_p,R_tr_x);
			hR_mom_th->Fill(R_tr_p,R_tr_th);
			hR_len_x->Fill(R_tr_pathl,R_tr_x);
			hR_len_th->Fill(R_tr_pathl,R_tr_th);
			hR_lenfp_x->Fill(R_tr_pathl-R_s2_trpath,R_tr_x);
			hR_lenfp_th->Fill(R_tr_pathl-R_s2_trpath,R_tr_th);

			hR_mom_y->Fill(R_tr_p,R_tr_y);
			hR_mom_ph->Fill(R_tr_p,R_tr_ph);
			hR_len_y->Fill(R_tr_pathl,R_tr_y);
			hR_len_ph->Fill(R_tr_pathl,R_tr_ph);
			double ksr = exp(-1.*R_tr_pathl*MK/R_tr_p/3.713);
//cout<<"survival ratio = "<<ksr<<endl;
//cout<<"R_tr_pathl = "<<R_tr_pathl<<endl;
//cout<<"R_tr_p = "<<R_tr_p<<endl;
			if(R_tr_p>1.7&&R_tr_p<2.0&&R_tr_pathl>25.&&R_tr_pathl<27.){	
			kaonsr[(int)((R_tr_p-1.65)*srbin/0.35)] += ksr;
			counts[(int)((R_tr_p-1.65)*srbin/0.35)] ++;
			double ksrave = ave[(int)((R_tr_p-1.65)*srbin/0.35)];
			dif[(int)((R_tr_p-1.65)*srbin/0.35)] += (ksr - ksrave)*(ksr - ksrave);
			}
			h2_kSR->SetBinContent(h2_kSR->GetXaxis()->FindBin(R_tr_pathl),h2_kSR->GetYaxis()->FindBin(R_tr_p),ksr);
//cout<<"SR= "<<ksr<<endl;
	//	}
		if(ct_cut&&event_selection){
			hR2_mom->Fill(R_tr_p);
			hR2_T2S2->Fill(R_tr_pathl);
			hL2_T2S2->Fill(L_tr_pathl);
			hR2_T2FP->Fill(R_tr_pathl-R_s2_trpath);
			hL2_T2FP->Fill(L_tr_pathl-L_s2_trpath);
			hL2_T2S2_mom->Fill(L_tr_pathl,L_tr_p);
			hR2_T2S2_mom->Fill(R_tr_pathl,R_tr_p);
		}


}//ENum

double tmp,tmp2;
double ave_new[srbin];
  for(int i=0;i<srbin;i++){
	if(counts[i]!=0){
		tmp = kaonsr[i]/counts[i];//average
		ave_new[i] = tmp;
	}
	else tmp = 0.;
	if(counts[i]!=0&&counts[i]!=1) tmp2 = sqrt(dif[i]/counts[i]/(counts[i]-1));//sigma of ave.
	else tmp2 = 0.;
//cout<<"survival ratio is "<<tmp<<endl;
//cout<<"kaonsr["<<i<<"] = "<<kaonsr[i]<<endl;
//cout<<"counts["<<i<<"] = "<<counts[i]<<endl;
	h_kSR->SetBinContent(i+1,tmp);
	h_kSR->SetBinError(i+1,tmp2);
  }



//	cout<<"nbunch="<<nbunch<<endl;
	TCanvas* c1 = new TCanvas("c1","c1");
	c1->Divide(2,2);
	c1->cd(1);
			hR_T2S2->Draw("");
			cout<<"hR_T2S2="<<hR_T2S2->Integral()<<endl;
	c1->cd(2);
			hL_T2S2->Draw("");
			cout<<"hL_T2S2="<<hL_T2S2->Integral()<<endl;
	c1->cd(3);
			hR_T2FP->Draw("");
			cout<<"hR_T2FP="<<hR_T2FP->Integral()<<endl;
	c1->cd(4);
			hL_T2FP->Draw("");
			cout<<"hL_T2FP="<<hL_T2FP->Integral()<<endl;

//Momentum vs. Path Length
	TCanvas* c2 = new TCanvas("c2","c2");
	c2->Divide(2,2);
	c2->cd(1);
	hL_T2S2_mom->Draw("colz");
	c2->cd(2);
	hR_T2S2_mom->Draw("colz");
	c2->cd(3);
	hL2_T2S2_mom->Draw("colz");
	c2->cd(4);
	hR2_T2S2_mom->Draw("colz");

	TCanvas* c3 = new TCanvas("c3","c3");
	hR_mom->Draw("");
	TCanvas* c4 = new TCanvas("c4","c4");
	hR2_mom->Draw("");

	TCanvas* c5 = new TCanvas("c5","c5");
	h2_kSR->Draw("colz");

	TCanvas* c6 = new TCanvas("c6","c6");
	TH1F *frame = c6->DrawFrame(1.65,0.,2.,0.20);
	h_kSR->SetLineColor(kRed);
	h_kSR->SetLineWidth(2);
	h_kSR->Draw("esame");

	TCanvas* c7 = new TCanvas("c7","c7");
	c7->Divide(2,3);
	c7->cd(1);
	hR_mom_x->Draw("colz");
	c7->cd(2);
	hR_mom_th->Draw("colz");
	c7->cd(3);
	hR_len_x->Draw("colz");
	c7->cd(4);
	hR_len_th->Draw("colz");
	c7->cd(5);
	hR_lenfp_x->Draw("colz");
	c7->cd(6);
	hR_lenfp_th->Draw("colz");

	TCanvas* c8 = new TCanvas("c8","c8");
	c8->Divide(2,3);
	c8->cd(1);
	hR_mom_y->Draw("colz");
	c8->cd(2);
	hR_mom_ph->Draw("colz");
	c8->cd(3);
	hR_len_y->Draw("colz");
	c8->cd(4);
	hR_len_ph->Draw("colz");

			hR_mom_th->Fill(R_tr_p,R_tr_th);
			hR_len_x->Fill(R_tr_pathl,R_tr_x);
			hR_len_th->Fill(R_tr_pathl,R_tr_th);
			hR_len_x->Fill(R_tr_pathl-R_s2_trpath,R_tr_x);
			hR_len_th->Fill(R_tr_pathl-R_s2_trpath,R_tr_th);
	ofstream fout("../information/ksr_tmp.dat");
		fout<<"#from mthesis_Fig/kaon_SR.C"<<endl;
		fout<<"Number of Bin"<<endl;
		fout<<srbin<<endl;
		fout<<"#"<<endl;
		fout<<"#"<<endl;
		fout<<"#"<<endl;
		double vpflux_temp;
		double vpflux_total=0.;
	for(int i=0; i<srbin; i++){
		fout<<i<<" "<<ave_new[i]<<endl;
	}
	cout<<"Average of each bin was successfully written in ../information/ksr_tmp.dat."<<endl;

cout << "Well done!" << endl;
}//fit
