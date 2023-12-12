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

void ebeam(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *file  = new TFile("../h2all_2020Nov.root","read");
  TTree *tree  = (TTree*)file ->Get("tree_out");


//---------------------------------------//
//               Branch                  //
//---------------------------------------//

	int nevent;
	double ebeam; 

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("nev",1);tree->SetBranchAddress("nev",&nevent);
	tree->SetBranchStatus("Bp_before",1);tree->SetBranchAddress("Bp_before",&ebeam);
	

	TH2F* hist  = new TH2F("hist","Ebeam vs #event",13000,0.,13000000.,1000,4.312,4.314);
  	SetTH2(hist, "", "Num. of events", "Ee [GeV]", 0.4);

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
      
			
//if(nevent%100000==0)cout<<"nevent="<<nevent<<", "<<i<<endl;
		hist->Fill((double)i, ebeam);

}//ENum

	TCanvas* c = new TCanvas("c","c",1800,1000);
	hist->Draw("");
	c->SetLeftMargin(0.14);
	c->SetRightMargin(0.14);
	c->SetTopMargin(0.14);
	c->SetBottomMargin(0.14);
	c->Modified();
	c->Update();
	gPad->Modified();
	gPad->Update();
	c->Print("./pdf/ebeam.pdf");

cout << "Well done!" << endl;
}//fit
