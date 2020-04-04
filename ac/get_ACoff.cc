//const   int nAC=24;
const   int nAC=26;
void get_ACoff(){

  TH1D* hac1[nAC];
  TF1* fac1[nAC];
  double min_adc=1000.;
  double max_adc=10000.;
  int bin_adc=9000;
  
  for(int i=0;i<nAC;i++){
    
    
    hac1[i]=new TH1D(Form("hac1_%d",i),"",bin_adc,min_adc,max_adc);
    fac1[i]=new TF1(Form("fac1_%d",i),"gausn(0)",min_adc,max_adc);
    fac1[i]->SetLineColor(2);
    fac1[i]->SetNpx(2000);
  }


  TChain* T=new TChain("T");
  string root_file="/data2/AC/tritium_111114_all.root";
  T->Add(root_file.c_str());
  int ENum=  T->GetEntries();
  cout<<"Get Entries : "<<ENum<<endl;
  double adc[nAC];
  T->SetBranchStatus("*",0);
  //  T->SetBranchStatus("R.a1.a",1);
  //  T->SetBranchAddress("R.a1.a",adc);
  T->SetBranchStatus("R.a2.a",1);
  T->SetBranchAddress("R.a2.a",adc);  
  for(int k=0;k<ENum;k++){
    T->GetEntry(k);

    for(int i=0;i<nAC;i++)
      hac1[i]->Fill(adc[i]);

  }//end Fill
  

  double th[nAC];
  for(int i=0;i<nAC;i++){
    th[i]=hac1[i]->GetBinCenter(hac1[i]->GetMaximumBin());
    fac1[i]->SetParameter(1,th[i]);
    hac1[i]->Fit(Form("fac1_%d",i),"QR","QR",th[i]-100.,th[i]+100.);
    th[i]=fac1[i]->GetParameter(1);
    //    cout<<"i "<<i<<" threshold "<<th[i]<<endl;
    cout<<th[i]<<" ";
  }


  TCanvas* c[5];
  for(int j=0;j<3;j++){
    c[j]=new TCanvas(Form("c%d",j),Form("c%d",j));
    c[j]->Divide(3,2);
    for(int l= 0;l<8;l++){
      c[j]->cd(l+1);
      hac1[l+8*j]->Draw();
      fac1[l+8*j]->Draw("same");
    }
  }
}
