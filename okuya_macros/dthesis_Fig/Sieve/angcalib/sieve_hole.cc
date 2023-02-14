

void sieve_hole(){

  //  bool rarm=true;
  bool rarm=false;
  TFile*f1 ;
  //  if(rarm)f1=new TFile(Form("../rootfiles/angcalib/ang_RHRS_woRas.root"));
  //  else {f1=new TFile(Form("../rootfiles/angcalib/ang_LHRS_woRas.root"));}
  //  f1=new TFile(Form("../rootfiles/zcalib/zt_RHRS_sieve.root"));
  //f1=new TFile(Form("../rootfiles/angcalib/ang_LHRS_woRas.root"));
   f1=new TFile(Form("../ang_RHRS_5th_0602.root"));
  //if(rarm) f1=new TFile(Form("../../rootfiles/angcalib/ang_RHRS_woRas.root"));
  //else  f1=new TFile(Form("../../rootfiles/angcalib/test_LHRS_woRas_test.root"));
  //  f1=new TFile(Form("../angcalib/rootfile/ang_RHRS_5th_0830_0.root"));
  TTree* t1=(TTree*)f1->Get("T");

  TH2F* hist_a=(TH2F*)f1->Get("h3_a");
  TH2F* hist_c=(TH2F*)f1->Get("h3_c");  
  hist_a->SetName("hist_a");
  hist_c->SetName("hist_c");



  Double_t trig5;
  Double_t trig4;
  Double_t trig1;
  double ps_asum,sh_asum;
  double gs_asum;
  double Zt[100];
  double Rvz[100];
  double Lvz[100];
  double Zt_tuned[100];
  double ztR_opt[100];
  double a1,a2,gc;
  double ss_x,ss_y;
  double xp[100],yp[100];
  
  t1->SetBranchAddress("DR.T1", &trig1);
  t1->SetBranchAddress("DR.T4", &trig4);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("R.sh.asum_c", &sh_asum);  
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("ss_x", &ss_x);
  t1->SetBranchAddress("ss_y", &ss_y);      
  
  if(rarm==true){
    t1->SetBranchAddress("R.tr.vz_opt",ztR_opt);
    t1->SetBranchAddress("R.tr.vz_tuned",Zt);
    t1->SetBranchAddress("R.tr.tg_th_tuned", xp);
    t1->SetBranchAddress("R.tr.tg_ph_tuned", yp);    
  } else{
    t1->SetBranchAddress("L.cer.asum_c",&gc);
    t1->SetBranchAddress("L.tr.vz_opt",ztR_opt);
    t1->SetBranchAddress("L.tr.vz_tuned",Zt);
    t1->SetBranchAddress("L.tr.tg_th_tuned", xp);
    t1->SetBranchAddress("L.tr.tg_ph_tuned", yp);
  }


  
  int ENum=t1->GetEntries();
  cout<<"Entries : "<<ENum<<endl;
  bool tree=false;
  if(tree){

    string ofname;
    string dir = "hist/";
    if(rarm)    ofname="sieveRHRS_drow.root";
    else     ofname="sieveLHRS_drow.root";
    TFile *ofp = new TFile(Form("%s%s",dir.c_str(),ofname.c_str()),"recreate");
    TTree *newtree = t1->CloneTree();
    

  }
  
  string name;
  if(rarm ) name = "RHRS";
  else      name = "LHRS";
  int bin_ss =100;
  double min_ss = -10.0;
  double max_ss =  10.0;
  
  TH2D* hss = new TH2D("hss",Form("SS at %s ; SSY [cm]; SSX [cm] ; Counts",name.c_str()),bin_ss,min_ss,max_ss,bin_ss,min_ss,max_ss);

  int bin_ang = 100;
  double min_ang = -0.1;
  double max_ang =  0.1;
  
  TH2D* hang = new TH2D("hang",Form("Xp and Yp at target (%s) ; Y' [mrad] ; X' [mrad]",name.c_str()),bin_ang,min_ang,max_ang,bin_ang,min_ang,max_ang);

  // Fill Tree
  for(int i = 0; i<ENum;i++){
    
    t1->GetEntry(i);
    bool pid_cut = false;
    if(!rarm && gc > 2000.) pid_cut = true;
    if(rarm  && a1>2. && a2>200.) pid_cut = true;
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
  hang->Draw();
  hang->SetStats(0);
  
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
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],34);
      
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
      c3->cd();
      mark[nhole]->Draw("same");

      
      nhole++;

    }
  }


    
    if(rarm){ mark[10]->SetMarkerColor(1); mark[10]->SetMarkerStyle(22);  mark[10]->SetMarkerSize(1.5); }
    else    { mark[12]->SetMarkerColor(1); mark[12]->SetMarkerStyle(22);  mark[12]->SetMarkerSize(1.5); }
    
      mark[43]->SetMarkerStyle(22);
      mark[43]->SetMarkerColor(1);
      mark[43]->SetMarkerSize(1.5);
    if(tree){
      hist_a->Write();
      hist_c->Write();
      hss   ->Write();
      hang  ->Write();
    }


    //==============================//
    //======= < PDF  > =============//
    //==============================//
    string ofpname = "./pdf/";//
    if(rarm) ofpname + "sieveRHRS_draw.pdf";
    else     ofpname + "sieveLHRS_draw.pdf";
    
    cout<<"Print is starting"<<endl;
    cout<<"pdf name: "<<ofpname<<endl;
    c0->Print(Form("%s[",ofpname.c_str()));
    c0->Print(Form("%s",ofpname.c_str()));
    c1->Print(Form("%s",ofpname.c_str()));
    c2->Print(Form("%s",ofpname.c_str()));
    c3->Print(Form("%s",ofpname.c_str()));
    c0->Print(Form("%s]",ofpname.c_str()));
}
