
void test_coin(){
TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1223.root");
  double ct,Ra1sum,Ra2sum,mm,Lz,Rz,Rpz,Lpz,ct_acc;

  TTree *T=(TTree*)fin->Get("tree");
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&ct);
 // T->SetBranchAddress("ct",&coin_t);
 T->SetBranchStatus("ac1_sum",1);
 T->SetBranchAddress("ac1_sum",&Ra1sum);
 T->SetBranchStatus("ac2_sum",1);
 T->SetBranchAddress("ac2_sum",&Ra2sum);
 T->SetBranchStatus("mm",1);
 T->SetBranchAddress("mm",&mm);
 T->SetBranchStatus("Lz",1);    
 T->SetBranchAddress("Lz",&Lz);
 T->SetBranchStatus("Rz",1);    
 T->SetBranchAddress("Rz",&Rz);
 T->SetBranchStatus("Rp",1);
 T->SetBranchAddress("Rp",&Rpz);
 T->SetBranchStatus("Lp",1);
 T->SetBranchAddress("Lp",&Lpz);
 T->SetBranchStatus("ct_acc",1);
 T->SetBranchAddress("ct_acc",&ct_acc);


 bool vz_cut;
 bool ac_cut;
 bool ct_cut;
 bool acc_cut;
 double ac1_th=1.4;
 double ac2_th_b=3.6;
 double ac2_th_t=10.;

  double min_coin=-20.0;
  double max_coin=10.0;
  double bin_coin=(max_coin-min_coin)/0.056;


 bin_coin=(int)bin_coin;
 double min_mm=0.5;
 double max_mm=1.5;
 //double bin_mm=(max_mm-min_mm)/0.002;
 //bin_mm=(int)bin_mm;
 int bin_mm=500;
 TH1F* hcoin_ac=new TH1F("hcoin_ac","Coin AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hcoin_acc=new TH1F("hcoin_acc","Coin ACC AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hcoin_p=new TH1F("hcoin_p","Coin ACC AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hmm_ac=new TH1F("hmm_ac","MM AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 TH1F* hmm_ac_acc=new TH1F("hmm_ac_acc","MM ACC AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 TH1F* hmm_ac_p=new TH1F("hmm_ac_p","MM AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 double ac1_npe,ac2_npe;
 double rz,lz;
 int evnt=T->GetEntries();
 cout<<"Events is "<<evnt<<endl;
 //  evnt=10000;
 int i=0;

 for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   if(k==1000000*i){
   cout<<"k "<<k<<endl;
   i=i+1;}
   vz_cut=false;
   ac_cut=false;
   ct_cut=false;
   acc_cut=false;
   ac1_npe=Ra1sum;
   ac2_npe=Ra2sum;
   rz=Rz;
   lz=Lz;

      if((rz<0.1 && -0.1<rz) && (-0.1<lz && lz<0.1))vz_cut=true;
      if(ac1_npe<ac1_th && ac2_th_b<ac2_npe && ac2_npe<ac2_th_t)ac_cut=true;
      if(-1.0<ct && ct<1.0)ct_cut=true;
      if((-40<ct && ct<-20)||(20<ct && ct<40))acc_cut=true;
   //======== FIll Hist ==================//
      if(vz_cut && ac_cut){
	hcoin_ac->Fill(ct);
	hcoin_acc->Fill(ct_acc);}
   if(vz_cut&&ac_cut&&ct_cut)hmm_ac->Fill(mm);
   if(vz_cut&&ac_cut&&acc_cut)hmm_ac_acc->Fill(mm);


}



}
