//--  Beam Energy vs #event  --//
//
//K. Okuyama (Nov. 23, 2023)
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
  h->SetMarkerSize(0.005);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetDecimals(3);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void gc(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //TFile *file  = new TFile("../h2all_2020Nov.root","read");
  TFile *file  = new TFile("../h2all_2024test.root","read");
  TTree *tree  = (TTree*)file ->Get("tree_out");


//---------------------------------------//
//               Branch                  //
//---------------------------------------//

	int nevent;
	double ebeam; 
	double L_mom, R_mom, B_mom; 
	double gcsum, ac1sum, ac2sum;//NPE SUM
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

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("nev",1);tree->SetBranchAddress("nev",&nevent);
	tree->SetBranchStatus("Bp_before",1);tree->SetBranchAddress("Bp_before",&ebeam);

  	tree->SetBranchStatus("ac1_npe_sum",1);  tree->SetBranchAddress("ac1_npe_sum", &ac1sum);
  	tree->SetBranchStatus("ac2_npe_sum",1);  tree->SetBranchAddress("ac2_npe_sum", &ac2sum);
  	tree->SetBranchStatus("L.cer.asum_c",1);  tree->SetBranchAddress("L.cer.asum_c", &gcsum);
  	tree->SetBranchStatus("ct_orig",1);  tree->SetBranchAddress("ct_orig", &ct);
  	tree->SetBranchStatus("Lp_c",1);  tree->SetBranchAddress("Lp_c", &L_mom);
  	tree->SetBranchStatus("Rp_c",1);  tree->SetBranchAddress("Rp_c", &R_mom);
  	tree->SetBranchStatus("Bp_c",1);  tree->SetBranchAddress("Bp_c", &B_mom);

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
	

	TH2F* ac1_mom  = new TH2F("ac1_mom","AC1 vs Mom.",100,1700.,1900.,50,0.,15.);
  	SetTH2(ac1_mom, "", "Momentum [MeV]", "AC1 NPE [p.e.]", 0.4);
	TH2F* ac2_mom  = new TH2F("ac2_mom","AC2 vs Mom.",100,1700.,1900.,50,0.,25.);
  	SetTH2(ac2_mom, "", "Momentum [MeV]", "AC2 NPE [p.e.]", 0.4);
	TH2F* ac1_ct  = new TH2F("ac1_ct","AC1 vs Cointime",100,-20.,20.,50,0.,15.);
  	SetTH2(ac1_ct, "", "Coincidence Time [ns]", "AC1 NPE [p.e.]", 0.4);
	TH2F* ac2_ct  = new TH2F("ac2_ct","AC2 vs Mom.",100,-20.,20.,50,0.,25.);
  	SetTH2(ac2_ct, "", "Coincidence Time [ns]", "AC2 NPE [p.e.]", 0.4);
	TH2F* gc_ct  = new TH2F("gc_ct","GC vs Mom.",100,-20.,20.,100,0.,20000.);
  	SetTH2(gc_ct, "", "Coincidence Time [ns]", "GC [ch]", 0.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool event_selection = false;

  long int ENum = tree->GetEntries(); 
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(long int i=0;i<ENum;i++){
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
      
	if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP&&R_mom>1.76&&R_mom<1.90&&L_mom>2.01&&L_mom<2.16){
	ac1_mom->Fill(R_mom*1000.,ac1sum);
	ac2_mom->Fill(R_mom*1000.,ac2sum);
	if(ac2sum>3.&&ac2sum<10.)ac1_ct->Fill(ct,ac1sum);
	if(ac1sum<3.75)ac2_ct->Fill(ct,ac2sum);
	gc_ct->Fill(ct,gcsum);
	}
			
//if(nevent%100000==0)cout<<"nevent="<<nevent<<", "<<i<<endl;
		//hist->Fill((double)i, ebeam);

}//ENum

	TCanvas* c = new TCanvas("c","c",1800,1000);
	gc_ct->Draw("colz");
	//ac1_mom->Draw("colz");
	//TCanvas* c2 = new TCanvas("c2","c2",1800,1000);
	//ac2_mom->Draw("colz");
	//TCanvas* c3 = new TCanvas("c3","c3",1800,1000);
	//ac1_ct->Draw("colz");
	//TCanvas* c4 = new TCanvas("c4","c4",1800,1000);
	//ac2_ct->Draw("colz");
	//c->SetLeftMargin(0.14);
	//c->SetRightMargin(0.14);
	//c->SetTopMargin(0.14);
	//c->SetBottomMargin(0.14);
	//c->Modified();
	//c->Update();
	//gPad->Modified();
	//gPad->Update();
	//c->Print("./pdf/ebeam.pdf");

cout << "Well done!" << endl;
}//fit
