//--  FP cut & Chi2 cut efficiency after Mom. cut   --//
//
//K. Okuyama (Feb. 17, 2023)
//
//This is taken over from mm_tuning.C
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
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void Chi2FP(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *file  = new TFile("../h2all_woFPChi2.root","read");
  TTree *tree  = (TTree*)file ->Get("tree_out");

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
  bool mom_cut = false;
  bool event_selection = false;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass



//---------------------------------------//
//               Branch                  //
//---------------------------------------//

	int nrun, NLtr, NRtr;
	double ct;

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

	TH2F* h_xth_L  = new TH2F("h_xth_L","LHRS: FP X vs FP #theta",100,-0.2,0.2,100,-0.8,0.8);
  	TH2F* h_xth_L_cut  = new TH2F("h_xth_L_cut","LHRS: FP X vs FP #theta (w/ cut)",100,-0.2,0.2,100,-0.8,0.8);
  	TH2F* h_xth_R  = new TH2F("h_xth_R","RHRS: FP X vs FP #theta",100,-0.2,0.2,100,-0.8,0.8);
  	TH2F* h_xth_R_cut  = new TH2F("h_xth_R_cut","RHRS: FP X vs FP #theta (w/ cut)",100,-0.2,0.2,100,-0.8,0.8);
  	TH1F* h_chi2_L  = new TH1F("h_chi2_L","LHRS: #chi^{2}",1000,0.,0.01);
  	TH1F* h_chi2_R  = new TH1F("h_chi2_R","RHRS: #chi^{2}",1000,0.,0.01);
  	TH1F* h_chi2_L_cut  = new TH1F("h_chi2_L_cut","LHRS: #chi^{2}",1000,0.,0.01);
  	TH1F* h_chi2_R_cut  = new TH1F("h_chi2_R_cut","RHRS: #chi^{2}",1000,0.,0.01);
  	SetTH2(h_xth_L, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  	SetTH2(h_xth_R, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  	SetTH2(h_xth_L_cut, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  	SetTH2(h_xth_R_cut, "", "X(FP) [m]", "#theta(FP) [rad]", 0.4);
  	SetTH1(h_chi2_L, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);
  	SetTH1(h_chi2_R, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);
  	SetTH1(h_chi2_L_cut, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);
  	SetTH1(h_chi2_R_cut, "", "#chi^{2}", "Counts", kAzure, kRed, kRed);

  int ENum = tree->GetEntries(); 
  long long int Nnocut = 0;
  long long int NFPLcut = 0;
  long long int NFPRcut = 0;
  long long int NTrLcut = 0;
  long long int NTrRcut = 0;
  long long int NFPcut = 0;
  long long int NTrcut = 0;
  long long int Nallcut = 0;
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(int i=0;i<ENum;i++){
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

		if(R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16)mom_cut=true;
		else mom_cut=false;

		if(fabs(ct)<1.006)ct_cut=true;
		else ct_cut=false;

	//w/o FP or Chi2 Cut
		if(mom_cut){
			h_xth_L->Fill(L_tr_x,L_tr_th);
			h_xth_R->Fill(R_tr_x,R_tr_th);
			h_chi2_L->Fill(L_tr_chi2);
			h_chi2_R->Fill(R_tr_chi2);
			Nnocut++;
			if(L_Tr){h_chi2_L_cut->Fill(L_tr_chi2);NTrLcut++;}
			if(R_Tr){h_chi2_R_cut->Fill(R_tr_chi2);NTrRcut++;}
			if(L_FP){h_xth_L_cut->Fill(L_tr_x,L_tr_th);NFPLcut++;}
			if(R_FP){h_xth_R_cut->Fill(R_tr_x,R_tr_th);NFPRcut++;}
			if(L_Tr&&R_Tr){NTrcut++;}
			if(L_FP&&R_FP){NFPcut++;}
			if(L_Tr&&R_Tr&&L_FP&&R_FP){Nallcut++;}
		}
			

}//ENum
	TH1F* he1_n = new TH1F("he1_n", "he1_n", 20, 0., 20.);//nominator
	TH1F* he1_d = new TH1F("he1_d", "he1_d", 20, 0., 20.);//denominator
	he1_n->SetLineColor(kRed);
	he1_d->SetLineColor(kAzure);

	he1_n->SetBinContent(1,NFPLcut);
	he1_d->SetBinContent(1,Nnocut);
	he1_n->SetBinContent(2,NFPRcut);
	he1_d->SetBinContent(2,Nnocut);
	he1_n->SetBinContent(3,NTrLcut);
	he1_d->SetBinContent(3,Nnocut);
	he1_n->SetBinContent(4,NTrRcut);
	he1_d->SetBinContent(4,Nnocut);
	he1_n->SetBinContent(5,NFPcut);
	he1_d->SetBinContent(5,Nnocut);
	he1_n->SetBinContent(6,NTrcut);
	he1_d->SetBinContent(6,Nnocut);
	he1_n->SetBinContent(7,Nallcut);
	he1_d->SetBinContent(7,Nnocut);

cout << "TEfficiency!" << endl;
cout<<"pEff1:"<<endl;
TEfficiency *pEff1;
if(TEfficiency::CheckConsistency(*he1_n,*he1_d,"w")){
pEff1 = new TEfficiency(*he1_n,*he1_d);
}
	cout<<"ENum FP_L (w/  cut)="<<NFPLcut<<endl;
	cout<<"ENum FP_L (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (FP_L) = "<<(double)NFPLcut/(double)Nnocut<<endl;
	cout<<"Eff. (FP_L) = "<<pEff1->GetEfficiency(1)<<",-"<<pEff1->GetEfficiencyErrorLow(1)<<",+"<<pEff1->GetEfficiencyErrorUp(1)<<endl;
	cout<<"ENum FP_R (w/  cut)="<<NFPRcut<<endl;
	cout<<"ENum FP_R (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (FP_R) = "<<(double)NFPRcut/(double)Nnocut<<endl;
	cout<<"Eff. (FP_R) = "<<pEff1->GetEfficiency(2)<<",-"<<pEff1->GetEfficiencyErrorLow(2)<<",+"<<pEff1->GetEfficiencyErrorUp(2)<<endl;
	cout<<"ENum Tr_L (w/  cut)="<<NTrLcut<<endl;
	cout<<"ENum Tr_L (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (Tr_L) = "<<(double)NTrLcut/(double)Nnocut<<endl;
	cout<<"Eff. (Tr_L) = "<<pEff1->GetEfficiency(4)<<",-"<<pEff1->GetEfficiencyErrorLow(3)<<",+"<<pEff1->GetEfficiencyErrorUp(3)<<endl;
	cout<<"ENum Tr_R (w/  cut)="<<NTrRcut<<endl;
	cout<<"ENum Tr_R (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (Tr_R) = "<<(double)NTrRcut/(double)Nnocut<<endl;
	cout<<"Eff. (Tr_R) = "<<pEff1->GetEfficiency(4)<<",-"<<pEff1->GetEfficiencyErrorLow(4)<<",+"<<pEff1->GetEfficiencyErrorUp(4)<<endl;
	cout<<"ENum FP (w/  cut)="<<NFPcut<<endl;
	cout<<"ENum FP (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (FP) = "<<(double)NFPcut/(double)Nnocut<<endl;
	cout<<"Eff. (FP) = "<<pEff1->GetEfficiency(5)<<",-"<<pEff1->GetEfficiencyErrorLow(5)<<",+"<<pEff1->GetEfficiencyErrorUp(5)<<endl;
	cout<<"ENum Tr (w/  cut)="<<NTrcut<<endl;
	cout<<"ENum Tr (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (Tr) = "<<(double)NTrcut/(double)Nnocut<<endl;
	cout<<"Eff. (Tr) = "<<pEff1->GetEfficiency(6)<<",-"<<pEff1->GetEfficiencyErrorLow(6)<<",+"<<pEff1->GetEfficiencyErrorUp(6)<<endl;
	cout<<"ENum All (w/  cut)="<<Nallcut<<endl;
	cout<<"ENum All (w/o cut)="<<Nnocut<<endl;
	cout<<"Eff. (All) = "<<(double)Nallcut/(double)Nnocut<<endl;
	cout<<"Eff. (All) = "<<pEff1->GetEfficiency(7)<<",-"<<pEff1->GetEfficiencyErrorLow(7)<<",+"<<pEff1->GetEfficiencyErrorUp(7)<<endl;

	TCanvas* c3 = new TCanvas("c3","c3");
	
	//c3->Print("./pdf/g.pdf");

cout << "Well done!" << endl;
}//fit
