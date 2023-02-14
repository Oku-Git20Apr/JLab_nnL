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
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.75);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

void hole_SS(){

  bool IsR = true;//true: HRS-R, false: HRS-L
  TFile *f1;
  if(IsR){f1=new TFile(Form("./Sieve/ang_RHRS_5th_0602.root"));
  }else{  f1=new TFile(Form("./Sieve/ang_LHRS_5th_0601.root"));}
  TTree* t1=(TTree*)f1->Get("T");
  bool rarm=false;
  string name;
  if(IsR){name = "RHRS";
  }else{  name = "LHRS";}
  int bin_ss =100;
  double min_ss = -10.0;
  double max_ss =  10.0;

  TH2F* hist_a=(TH2F*)f1->Get("h3_a");
  TH2F* hist_c=(TH2F*)f1->Get("h3_c");  
  hist_a->SetName("hist_a");
  hist_c->SetName("hist_c");

  TH2F* h_holeSS = new TH2F("hole_SS",Form("SS at %s ; Y at Sieve Slit [cm]; X at Sieve Slit [cm] ; Counts",name.c_str()),bin_ss,-5.0,5.0,bin_ss,-8.0,8.0);
  SetTH2(h_holeSS, "", "Y at Sieve Slit [cm]", "X at Sieve Slit [cm]", 0.4);
  h_holeSS->GetZaxis()->SetTitle("Counts");

  for(int i=1;i<=100;i++){
	for(int j=1;j<=100;j++){
		h_holeSS->SetBinContent(i,j,hist_c->GetBinContent(i,j));
	}
  }
  

  double a1,a2,gc;
  double ss_x,ss_y;
  double xp[100],yp[100];
  double trig5;
  
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("L.cer.asum_c",&gc);
  t1->SetBranchAddress("ss_x", &ss_x);
  t1->SetBranchAddress("ss_y", &ss_y);      
  t1->SetBranchAddress("DR.T5",&trig5);
  
  int ENum=t1->GetEntries();
  cout<<"Entries : "<<ENum<<endl;
  
  //TH2D* hss = new TH2D("hss",Form("SS at %s ; Y at Sieve Slit [cm]; X at Sieve Slit [cm] ; Counts",name.c_str()),bin_ss,min_ss,max_ss,bin_ss,min_ss,max_ss);
  TH2D* hss = new TH2D("hss",Form("SS at %s ; Y at Sieve Slit [cm]; X at Sieve Slit [cm] ; Counts",name.c_str()),bin_ss,-5.0,5.0,bin_ss,-8.0,8.0);
  SetTH2(hss, "", "Y at Sieve Slit [cm]", "X at Sieve Slit [cm]", 0.4);
  hss->GetZaxis()->SetTitle("Counts");

  int bin_ang = 100;
  double min_ang = -0.1;
  double max_ang =  0.1;
  
  TH2D* hang = new TH2D("hang",Form("Xp and Yp at target (%s) ; Y' [mrad] ; X' [mrad]",name.c_str()),bin_ang,min_ang,max_ang,bin_ang,min_ang,max_ang);

  // Fill Tree
  for(int i = 0; i<ENum;i++){
    
    t1->GetEntry(i);
    bool pid_cut = false;
    if(!IsR && gc > 5.) pid_cut = true;
    if(IsR && a1>1. && a2>2000.) pid_cut = true;

    if(pid_cut ){
      hss ->Fill(ss_y, ss_x);
      hang -> Fill(yp[0], xp[0]);
    }
    
  }

  

  
const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 8;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double refx[nsshole];
double refy[nsshole];
double selec_widthx = 0.60; // selection width in x (dispersive plane)
double selec_widthy = 0.45; // selection width in y   
 int nhole = 0;
   double ssy_cent_real[nrow];
   double ssx_cent_real[ncol];


  //===== Draw ======//
  TCanvas* c0 = new TCanvas("c0","c0");
  TCanvas* c1 = new TCanvas("c1","c1");  
  TCanvas* c2 = new TCanvas("c2","c2");
  TCanvas* c3 = new TCanvas("c3","c3");
  
  c0->cd();
  hist_a->Draw("colz");
  hist_a->SetStats(0);
  c1->cd();
  hist_c->Draw("colz");
  hist_c->SetStats(0);
  c2->cd();
  hss->Draw("colz");
  hss->SetStats(0);
  c3->cd();
  h_holeSS->Draw("colz");
  h_holeSS->SetStats(0);
  
  TMarker* mark[nsshole];
  

    for(int j=0; j<nrow; j++){
  for(int i=0; i<ncol ; i++){      
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      //=== correction error point ====//
      if(j==10 && i==2) refy[nhole]=-1.87452;
      if(j==9  && i==1) refy[nhole]=-2.49936;
      if(j==8  && i==0) refy[nhole]=-4.377388;      
      //===============================//
      //      cout<<"j : "<<j<<" ssx "<<ssx_cent_real[j]<<" i "<<i<<" ssy "<<ssy_cent_real[i]<<endl;
      
      //      mark[nhole] = new TMarker(refy[nhole],refx[nhole],28);
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],46);
      
      //      if(rarm) mark[10]->SetMarkerColor(2);
      //      else     mark[12]->SetMarkerColor(2);
      //      mark[43]->SetMarkerColor(2);
      

      //      mark[nhole]->SetMarkerColor(3);
      mark[nhole]->SetMarkerColor(1);
      c0->cd();
      mark[nhole]->Draw("same");
      c1->cd();
      mark[nhole]->Draw("same");
      c2->cd();
      mark[nhole]->Draw("same");
	  c2->SetLeftMargin(0.14);
	  c2->SetRightMargin(0.14);
	  c2->SetTopMargin(0.14);
	  c2->SetBottomMargin(0.14);
	  c2->Modified();
	  c2->Update();
	  gPad->Modified();
	  gPad->Update();
      c3->cd();
      mark[nhole]->Draw("same");
	  c3->SetLeftMargin(0.14);
	  c3->SetRightMargin(0.14);
	  c3->SetTopMargin(0.14);
	  c3->SetBottomMargin(0.14);
	  c3->Modified();
	  c3->Update();
	  gPad->Modified();
	  gPad->Update();

      
      nhole++;

    }
  }


    
    if(IsR){ mark[10]->SetMarkerColor(1); mark[10]->SetMarkerStyle(47);  mark[10]->SetMarkerSize(1.5); }
    else    { mark[12]->SetMarkerColor(1); mark[12]->SetMarkerStyle(47);  mark[12]->SetMarkerSize(1.5); }
    
      mark[43]->SetMarkerStyle(47);
      mark[43]->SetMarkerColor(1);
      mark[43]->SetMarkerSize(1.5);
      hist_a->Write();
      hist_c->Write();
      hss   ->Write();
      hang  ->Write();


    //==============================//
    //======= < PDF  > =============//
    //==============================//
    cout<<"Print is starting"<<endl;
    string ofpname;
    if(IsR){ ofpname = "./pdf/ss_HRSR.pdf";
    c3->Print(Form("%s",ofpname.c_str()));
    }else{   ofpname = "./pdf/ss_HRSL.pdf";
    c2->Print(Form("%s",ofpname.c_str()));}
}
