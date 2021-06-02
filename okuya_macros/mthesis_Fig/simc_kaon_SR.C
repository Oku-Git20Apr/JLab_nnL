//--  Kaon's survival ratio   --//
//
//K. Okuyama (May. 24, 2020)
//
//This is taken over from kaon_SR.C
//No array branch mode 
//


void simc_kaon_SR(){

  TFile *file = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/decay_L.root","read");//
  TTree *tree = (TTree*)file->Get("SNT");

    
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


	float R_mom; 
	float R_tr_pathl;
	float R_x;//h_xfp
	float R_y;//h_yfp
	float R_dx;//h_xpfp
	float R_dy;//h_ypfp
	float R_ksr;

	tree->SetBranchStatus("*",0);
  	tree->SetBranchStatus("Rp_rec",1);  tree->SetBranchAddress("Rp_rec", &R_mom);
    tree->SetBranchStatus("h_xfp",1);  tree->SetBranchAddress("h_xfp", &R_x);
    tree->SetBranchStatus("h_yfp",1);  tree->SetBranchAddress("h_yfp", &R_y);
    tree->SetBranchStatus("h_xpfp",1);  tree->SetBranchAddress("h_xpfp", &R_dx);
    tree->SetBranchStatus("h_ypfp",1);  tree->SetBranchAddress("h_ypfp", &R_dy);
    tree->SetBranchStatus("PathLength",1);  tree->SetBranchAddress("PathLength", &R_tr_pathl);
    tree->SetBranchStatus("SurvivalRatio",1);  tree->SetBranchAddress("SurvivalRatio", &R_ksr);




  //TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  //h1->GetXaxis()->SetTitle("coin time (ns)");
  //h1->GetYaxis()->SetTitle("Counts / 100 ps");
  //h1->GetXaxis()->SetRangeUser(-14.0,17.);
  //h1 ->SetLineColor(2);
  //h1->SetLineWidth(2);
  TH1F* hR_mom  = new TH1F("hR_mom","hR_mom",1000,1.7,2.0);
  TH1F* hR2_mom  = new TH1F("hR2_mom","hR2_mom",1000,1.7,2.0);
  TH1F* hR_T2S2  = new TH1F("hR_T2S2","hR_T2S2",1000,26.,29.);
  TH2F* hR_T2S2_mom  = new TH2F("hR_T2S2_mom","hR_T2S2_mom",100,25.,27.,100,1.7,2.);
  TH2F* hR2_T2S2_mom  = new TH2F("hR2_T2S2_mom","hR2_T2S2_mom",100,25.,27.,100,1.76,1.9);

  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  TH2F* hR_mom_x  = new TH2F("hR_mom_x","hR_mom_x;Momentum [GeV/c];X(FP) [m]",100,1.76,1.9,100,-1.0,1.0);
  TH2F* hR_mom_th  = new TH2F("hR_mom_th","hR_mom_th;Momentum [GeV/c];X'(FP) [rad]",100,1.76,1.9,100,-0.15,0.15);
  TH2F* hR_len_x  = new TH2F("hR_len_x","hR_len_x;Path Length (T to S2) [m];X(FP) [m]",100,25.,27.,100,-1.0,1.0);
  TH2F* hR_len_th  = new TH2F("hR_len_th","hR_len_th;Path Length (T to S2) [m];X'(FP) [rad]",100,25.,27.,100,-0.15,0.15);
  TH2F* hR_lenfp_x  = new TH2F("hR_lenfp_x","hR_lenfp_x;Path Length (T to FP) [m];X(FP) [m]",100,21.,25.,100,-1.0,1.0);
  TH2F* hR_lenfp_th  = new TH2F("hR_lenfp_th","hR_lenfp_th;Path Length (T to FP) [m];X'(FP) [rad]",100,21.,25.,100,-0.15,0.15);
  TH2F* hR_mom_y  = new TH2F("hR_mom_y","hR_mom_y;Momentum [GeV/c];Y(FP) [m]",100,1.76,1.9,100,-0.1,0.1);
  TH2F* hR_mom_ph  = new TH2F("hR_mom_ph","hR_mom_ph;Momentum [GeV/c];Y'(FP) [rad]",100,1.76,1.9,100,-0.25,0.25);
  TH2F* hR_len_y  = new TH2F("hR_len_y","hR_len_y;Path Length (T to S2) [m];Y(FP) [m]",100,25.,27.,100,-0.1,0.1);
  TH2F* hR_len_ph  = new TH2F("hR_len_ph","hR_len_ph;Path Length (T to S2) [m];Y'(FP) [rad]",100,25.,27.,100,-0.25,0.25);

  TH2F* h2_kSR  = new TH2F("h2_kSR","h2_kSR",100,25.,27.,100,1.7,2.);
  int srbin = 80;
  TH1D* h_kSR  = new TH1D("h_kSR","h_kSR",srbin,1.65,2.);
  TH1D* h_kSR2  = new TH1D("h_kSR2","h_kSR2 from SIMC",srbin,1.65,2.);
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
	string daq_file = "../information/simc_ksr_tmp.dat";//input file
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

	R_mom /= 1000.;// MeV/c-->GeV/c
	R_tr_pathl /= 100.;// cm-->m

    if(i%time_div==0){
      end = time(NULL);
      time(&end);
      double diff = difftime(end,start);
      double esttime = diff * ENum / (i+1) - diff;
      cout<<i<<" / "<<ENum<<" ("<<i*100/ENum<<"%) : "<<Form("%.0lf sec passed,  %.0lf sec left",diff,esttime)<<endl;
    }

			R_tr_pathl += -1.83*(1.+R_dx*R_dx+R_dy*R_dy);// 183cm (Z_aerogel)
			hR_mom->Fill(R_mom);
			hR_T2S2->Fill(R_tr_pathl);
			hR_T2S2_mom->Fill(R_tr_pathl,R_mom);

			hR_mom_x->Fill(R_mom,R_x*0.01);
			hR_mom_th->Fill(R_mom,R_dx);
			hR_len_x->Fill(R_tr_pathl,R_x*0.01);
			hR_len_th->Fill(R_tr_pathl,R_dx);
			hR_mom_y->Fill(R_mom,R_y*0.01);
			hR_mom_ph->Fill(R_mom,R_dy);
			hR_len_y->Fill(R_tr_pathl,R_y*0.01);
			hR_len_ph->Fill(R_tr_pathl,R_dy);

			double ksr = exp(-1.*R_tr_pathl*MK/R_mom/3.713);
//cout<<"survival ratio = "<<ksr<<endl;
//cout<<"R_tr_pathl = "<<R_tr_pathl<<endl;
//cout<<"R_mom = "<<R_mom<<endl;
			//ksr = R_ksr;
			if(R_mom>1.7&&R_mom<2.0&&R_tr_pathl>25.&&R_tr_pathl<29.){	
			kaonsr[(int)((R_mom-1.65)*srbin/0.35)] += ksr;
			counts[(int)((R_mom-1.65)*srbin/0.35)] ++;
			double ksrave = ave[(int)((R_mom-1.65)*srbin/0.35)];
			dif[(int)((R_mom-1.65)*srbin/0.35)] += (ksr - ksrave)*(ksr - ksrave);
			}
			h2_kSR->SetBinContent(h2_kSR->GetXaxis()->FindBin(R_tr_pathl),h2_kSR->GetYaxis()->FindBin(R_mom),ksr);
			h_kSR2->Fill(R_ksr);
//cout<<"SR= "<<ksr<<endl;
	


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
	hR_T2S2->Draw("");

//	Momentum vs. Path Length
	TCanvas* c2 = new TCanvas("c2","c2");
	hR_T2S2_mom->Draw("colz");

	TCanvas* c3 = new TCanvas("c3","c3");
	hR_mom->Draw("");

	TCanvas* c4 = new TCanvas("c4","c4");
	hR2_mom->Draw("");

	TCanvas* c5 = new TCanvas("c5","c5");
	h2_kSR->Draw("colz");

	TCanvas* c6 = new TCanvas("c6","c6");
	TH1F *frame = c6->DrawFrame(1.65,0.,2.,0.20);
	h_kSR->SetLineColor(kAzure);
	h_kSR->SetLineWidth(2);
	h_kSR->Draw("esame");

	TCanvas* c7 = new TCanvas("c7","c7");
	TH1F *frame2 = c7->DrawFrame(1.65,0.,2.,0.20);
	h_kSR2->SetLineColor(kRed);
	h_kSR2->SetLineWidth(2);
	h_kSR2->Draw("esame");

	TCanvas* c8 = new TCanvas("c8","c8");
	c8->Divide(2,3);
	c8->cd(1);
	hR_mom_x->Draw("colz");
	c8->cd(2);
	hR_mom_th->Draw("colz");
	c8->cd(3);
	hR_len_x->Draw("colz");
	c8->cd(4);
	hR_len_th->Draw("colz");
	c8->cd(5);
	hR_lenfp_x->Draw("colz");
	c8->cd(6);
	hR_lenfp_th->Draw("colz");

	TCanvas* c9 = new TCanvas("c9","c9");
	c9->Divide(2,3);
	c9->cd(1);
	hR_mom_y->Draw("colz");
	c9->cd(2);
	hR_mom_ph->Draw("colz");
	c9->cd(3);
	hR_len_y->Draw("colz");
	c9->cd(4);
	hR_len_ph->Draw("colz");

	ofstream fout("../information/simc_ksr_tmp.dat");
		fout<<"#from mthesis_Fig/simc_kaon_SR.C"<<endl;
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
	cout<<"Average of each bin was successfully written in ../information/simc_ksr_tmp.dat."<<endl;

cout << "Well done!" << endl;
}//fit
