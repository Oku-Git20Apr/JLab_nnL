//-- Drawing Multifoil distribution  --//
//
//K. Okuyama (Feb. 14, 2023)
//
//This is taken over from Vertex_fig.C
//No array branch mode 
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
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetNdivisions(510);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void Multifoil_dist(){
  
  TFile *file = new TFile("../h2all_Multi_Hkin.root","read");
  TTree *tree = (TTree*)file->Get("tree_out");
    
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	double L_tr_vz, R_tr_vz;
	double L_tr_x, L_tr_th;
	double R_tr_x, R_tr_th;
	double L_tr_chi2, R_tr_chi2;
	tree->SetBranchStatus("*",0);
  	tree->SetBranchStatus("L.tr.vz",1);  tree->SetBranchAddress("L.tr.vz", &L_tr_vz);
  	tree->SetBranchStatus("R.tr.vz",1);  tree->SetBranchAddress("R.tr.vz", &R_tr_vz);
  	tree->SetBranchStatus("L.tr.x" ,1);  tree->SetBranchAddress("L.tr.x" , &L_tr_x);
  	tree->SetBranchStatus("L.tr.th",1);  tree->SetBranchAddress("L.tr.th", &L_tr_th);
	tree->SetBranchStatus("R.tr.x" ,1);  tree->SetBranchAddress("R.tr.x" , &R_tr_x );
  	tree->SetBranchStatus("R.tr.th",1);  tree->SetBranchAddress("R.tr.th", &R_tr_th);
  	tree->SetBranchStatus("L.tr.chi2",1);  tree->SetBranchAddress("L.tr.chi2", &L_tr_chi2);
  	tree->SetBranchStatus("R.tr.chi2",1);  tree->SetBranchAddress("R.tr.chi2", &R_tr_chi2);
    TH1F* h_zL  = new TH1F("h_zL","",1000.,-25.,25.);
    TH1F* h_zR  = new TH1F("h_zR","",1000.,-25.,25.);
    SetTH1(h_zL, "", "Z-vertex [cm]", "Counts", kAzure, kRed, kRed);
    SetTH1(h_zR, "", "Z-vertex [cm]", "Counts", kAzure, kRed, kRed);

    bool L_Tr = false;
    bool L_FP = false;
    bool R_Tr = false;
    bool R_FP = false;

    int ENum = tree->GetEntries(); 
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
		if(R_Tr && L_Tr && R_FP && L_FP){
			h_zL->Fill(L_tr_vz*100.);
			h_zR->Fill(R_tr_vz*100.);
		}

}//ENum
	TCanvas* c1 = new TCanvas("c1","c1",800,400);
	h_zL->Draw("");
	c1->SetLeftMargin(0.14);
	c1->SetRightMargin(0.14);
	c1->SetTopMargin(0.14);
	c1->SetBottomMargin(0.14);
	c1->Modified();
	c1->Update();
	gPad->Modified();
	gPad->Update();
	TCanvas* c2 = new TCanvas("c2","c2",800,400);
	h_zR->Draw("");
	c2->SetLeftMargin(0.14);
	c2->SetRightMargin(0.14);
	c2->SetTopMargin(0.14);
	c2->SetBottomMargin(0.14);
	c2->Modified();
	c2->Update();
	gPad->Modified();
	gPad->Update();

c1->Print("./pdf/Multifoil_HRSL.pdf");
c2->Print("./pdf/Multifoil_HRSR.pdf");

cout << "Well done!" << endl;
}//fit
