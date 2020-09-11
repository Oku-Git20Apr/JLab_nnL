void z_fitting_test(){


	string ofname = "z_fitting.dat";
	string pdfname = "z_fitting.pdf";
cout << " output file name is " << ofname << endl;
cout << " output pdf file name is " << pdfname << endl;
  
  TFile *file = new TFile("fout.root","read");
  TFile *file_new = new TFile("fout_new.root","recreate");
  TTree *tree = ((TTree*)file->Get("tree_out"))->CloneTree();
  file->Write();
 // TTree *tree = tree_out->CloneTree(); 
//  TTree *tree = (TTree*)file->Get("tree");

  //TChain *tree_new = new TChain("tree_new");
  //TChain *tree = new TChain("tree");
  //tree->Add("zeff.root");
 
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
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
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;
 
 double center_pi, center_k, center_p, center_L, center_S;
 double range_pi, range_k, range_p, range_L, range_S;
 double n_pi_noZ, n_k_noZ, n_p_noZ, n_L_noZ, n_S_noZ;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_noZ, mean_k_noZ, mean_p_noZ, mean_L_noZ, mean_S_noZ;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_noZ, sig_k_noZ, sig_p_noZ, sig_L_noZ, sig_S_noZ;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];


 TF1* fpi_noZ;
 TF1* fk_noZ;
 TF1* fp_noZ;
 TF1* fmmbg_noZ;
 TF1* fL_noZ;
 TF1* fS_noZ;
 TF1* fmm_noZ;
 TF1* fcoin_noZ;
 TF1* facc[100];
 TF1* fpi[100];
 TF1* fk[100];
 TF1* fp[100];
 TF1* fcoin[100];
 TF1* fL[100];
 TF1* fS[100];
 TF1* fmm[100];
 TF1* fmmbg[100];

 TH1F *hcoin_k_fom_noZ;
 TH1F *hcoin_bg_fom_noZ;
 TH1F *hcoin_wo_bg_fom_noZ;
 TH1F *hcoin_pi_noZ;
 TH1F *hmm_L_fom_noZ;
 TH1F *hmm_bg_fom_noZ;
 TH1F *hmm_wo_bg_fom_noZ;
 TH1F *hmm_pi_fom_noZ;
 TH1F *hmm_pibg_fom_noZ;
 TH1F *hmm_pi_wobg_fom_noZ;
 TH1F *hcoin_k_fom[100];
 TH1F *hcoin_bg_fom[100];
 TH1F *hcoin_wo_bg_fom[100];
 TH1F *hcoin_pi[100];
 TH1F *hmm_L_fom[100];
 TH1F *hmm_bg_fom[100];
 TH1F *hmm_wo_bg_fom[100];
 TH1F *hmm_pi_fom[100];
 TH1F *hmm_pibg_fom[100];
 TH1F *hmm_pi_wobg_fom[100];

 int nth=100;
 double zver[100];
 zver[0]=0.;
 for(int i=1; i<nth; i++) zver[i]=zver[i-1]+0.005;//SUM
TH1F* h_pisr1 = new TH1F("h_pisr1","Pion Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_ksr1 = new TH1F("h_ksr1","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_psr1 = new TH1F("h_psr1","Proton Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_Lsr1 = new TH1F("h_Lsr1","Lambda Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_Ssr1 = new TH1F("h_Ssr1","Signa Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_pitot1 = new TH1F("h_pitot1","Pion Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_ktot1 = new TH1F("h_ktot1","Kaon Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_ptot1 = new TH1F("h_ptot1","Proton Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_Ltot1 = new TH1F("h_Ltot1","Lambda Survival Ratio",nth-1,zver[0],zver[nth-1]);
TH1F* h_Stot1 = new TH1F("h_Stot1","Sigma Survival Ratio",nth-1,zver[0],zver[nth-1]);
//--------------------------------------------//
 h_pisr1->Sumw2();
 h_pisr1->SetLineColor(kRed);
 h_ksr1->Sumw2();
 h_ksr1->SetLineColor(kRed);
 h_psr1->Sumw2();
 h_psr1->SetLineColor(kRed);
 h_Lsr1->Sumw2();
 h_Lsr1->SetLineColor(kAzure);
 h_Ssr1->Sumw2();
 h_Ssr1->SetLineColor(kCyan);
 h_pitot1->Sumw2();
 h_pitot1->SetLineColor(kBlack);
 h_ktot1->Sumw2();
 h_ktot1->SetLineColor(kBlack);
 h_ptot1->Sumw2();
 h_ptot1->SetLineColor(kBlack);
 h_Ltot1->Sumw2();
 h_Ltot1->SetLineColor(kAzure);
 h_Stot1->Sumw2();
 h_Stot1->SetLineColor(kCyan);
//-----Main FOM histgram------//
 int h3_fom_bin_ac1 = 10;
 int h3_fom_bin_ac2l = 10;
 int h3_fom_bin_ac2u = 10;
 

	hcoin_k_fom_noZ=new TH1F("hcoin_k_fom_noZ","",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_bg_fom_noZ=new TH1F("hcoin_bg_fom_noZ","",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_wo_bg_fom_noZ=new TH1F("hcoin_wo_bg_fom_noZ","",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_pi_noZ=new TH1F("hcoin_pi_noZ","",bin_coin_c,min_coin_c,max_coin_c);

	hmm_L_fom_noZ=new TH1F("hmm_L_fom_noZ","",bin_mm,min_mm,max_mm);
	hmm_bg_fom_noZ=new TH1F("hmm_bg_fom_noZ","",bin_mm,min_mm,max_mm);
	hmm_wo_bg_fom_noZ=new TH1F("hmm_wo_bg_fom_noZ","",bin_mm,min_mm,max_mm);
	hmm_pi_fom_noZ=new TH1F("hmm_pi_fom_noZ","",bin_mm,min_mm,max_mm);
	hmm_pibg_fom_noZ=new TH1F("hmm_pibg_fom_noZ","",bin_mm,min_mm,max_mm);
	hmm_pi_wobg_fom_noZ=new TH1F("hmm_pi_wobg_fom_noZ","",bin_mm,min_mm,max_mm);
 for (int i=0;i<nth;i++){
	hcoin_k_fom[i]=new TH1F(Form("hcoin_k_fom[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_bg_fom[i]=new TH1F(Form("hcoin_bg_fom[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_wo_bg_fom[i]=new TH1F(Form("hcoin_wo_bg_fom[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_pi[i]=new TH1F(Form("hcoin_pi[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c);

	hmm_L_fom[i]=new TH1F(Form("hmm_L_fom[%d]",i),"",bin_mm,min_mm,max_mm);
	hmm_bg_fom[i]=new TH1F(Form("hmm_bg_fom[%d]",i),"",bin_mm,min_mm,max_mm);
	hmm_wo_bg_fom[i]=new TH1F(Form("hmm_wo_bg_fom[%d]",i),"",bin_mm,min_mm,max_mm);
	hmm_pi_fom[i]=new TH1F(Form("hmm_pi_fom[%d]",i),"",bin_mm,min_mm,max_mm);
	hmm_pibg_fom[i]=new TH1F(Form("hmm_pibg_fom[%d]",i),"",bin_mm,min_mm,max_mm);
	hmm_pi_wobg_fom[i]=new TH1F(Form("hmm_pi_wobg_fom[%d]",i),"",bin_mm,min_mm,max_mm);
  }

//---Cointime BG---//
  int ev = 0;
  int ENum = 0;
  double ct,mm,ct_den,mm_den,mmbg_den;
  double ct_eff[100], mm_eff[100], mmbg_eff[100];
	//tree_new->SetBranchStatus("*",0);
	//tree_new->SetBranchStatus("tr.ct_gb",1);tree_new->SetBranchAddress("tr.ct_gb",&ct);
	//tree_new->SetBranchStatus("tr.ct_eff[100]",1);tree_new->SetBranchAddress("tr.ct_eff[100]",&ct_eff);
	//tree_new->SetBranchStatus("tr.ct_den",1);tree_new->SetBranchAddress("tr.ct_den",&ct_den);
	//tree_new->SetBranchStatus("tr.mm_den",1);tree_new->SetBranchAddress("tr.mm_den",&mm_den);
	//tree_new->SetBranchStatus("tr.mm_eff[100]",1);tree_new->SetBranchAddress("tr.mm_eff[100]",&mm_eff);
	//tree_new->SetBranchStatus("tr.mmbg_den",1);tree_new->SetBranchAddress("tr.mmbg_den",&mmbg_den);
	//tree_new->SetBranchStatus("tr.mmbg_eff[100]",1);tree_new->SetBranchAddress("tr.mmbg_eff[100]",&mmbg_eff);
  ENum=tree->GetEntries();
  cout<<"Events: "<<ENum<<endl; 
	for(int k=0;k<ENum;k++){
		tree->GetEntry(k);

	if(k==ev*100000){
cout << "Event (Fill) : " << k << "/" << ENum << endl;
	ev += 1;
	}

//////////////////////////////////////No Z
				hcoin_k_fom_noZ->Fill(ct_den);
				if(20.<ct_den && ct_den<100.){
					double ct_den_ = ct_den;
				        while(1){
					  if(-20.<ct_den && ct_den<20.){
						 hcoin_bg_fom_noZ->Fill(ct_den);
						 hmm_bg_fom_noZ->Fill(mm_den);
						 hmm_pibg_fom_noZ->Fill(mm_den);
						 break;}
					       else if(ct_den<-20.){ct_den=ct_den+40.;}
					       else if(20.<ct_den){ct_den=ct_den-40.;}
					 }
					ct_den = ct_den_;
					}//cointime
					if(fabs(ct_den)<1.){
						hmm_L_fom_noZ->Fill(mm_den);
				}
					if(fabs(ct_den-3.05)<0.7){
						hmm_pi_fom_noZ->Fill(mm_den);//MM if pion
					}
//////////////////////////////////////No Z

//////////////////////////////////////Z Dependence
	for(int i=0;i<nth;i++){
				hcoin_k_fom[i]->Fill(ct_eff[i]);
				if(20.<ct_eff[i]&&ct_eff[i]<100.){
					double ct_eff_ = ct_eff[i];
//kazuki
				        while(1){
					  if(-20.<ct_eff[i] && ct_eff[i]<20.){
						 hcoin_bg_fom[i]->Fill(ct_eff[i]);
						 hmm_bg_fom[i]->Fill(mm_eff[i]);
						 hmm_pibg_fom[i]->Fill(mm_eff[i]);//MM if pion
						 break;}
					       else if(ct_eff[i]<-20.){ct_eff[i]=ct_eff[i]+40.;}
					       else if(20.<ct_eff[i]){ct_eff[i]=ct_eff[i]-40.;}
					 }
						  ct_eff[i] = ct_eff_;
				}
		//-------------------------------------------//



					if(fabs(ct_eff[i])<1.){
						hmm_L_fom[i]->Fill(mm_eff[i]);
					}//cointime
					if(fabs(ct_eff[i]-3.05)<0.7){
						// def_sig_pi=0.443; def_mean_pi=3.0;
						hmm_pi_fom[i]->Fill(mm_eff[i]);//MM if pion
					}
	}//for i	

//////////////////////////////////////Z Dependence

	}//ENum


//-----No Z Cut-----//
 hcoin_bg_fom_noZ->Scale(40./80.);
 hcoin_wo_bg_fom_noZ->Add(hcoin_k_fom_noZ,hcoin_bg_fom_noZ,1.0,-1.0);
 fp_noZ=new TF1("fp_noZ","gausn(0)",min_coin_c,max_coin_c);
 fp_noZ->SetNpx(2000);
 fpi_noZ =new TF1("fpi_noZ","gausn(0)+gausn(3)",def_mean_pi-3*def_sig_pi,def_mean_pi+6*def_sig_pi);
 fpi_noZ->SetNpx(2000);
 fk_noZ=new TF1("fk_noZ","gausn(0)",min_coin_c,max_coin_c);
 fk_noZ->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom_noZ->Fit("fp_noZ","Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
//n_p_noZ=fp_noZ->GetParameter(0);
 mean_p_noZ=fp_noZ->GetParameter(1);
 sig_p_noZ=fp_noZ->GetParameter(2);
 //center_p=mean_p_noZ;
 //range_p=2*sig_p_noZ;
 //range_p=1.435;
 center_p=def_mean_p;
 range_p=2*def_sig_p;
 n_p_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_p-range_p),hcoin_wo_bg_fom_noZ->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi_noZ->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 //fpi_noZ->SetParameters(1000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom_noZ->Fit("fpi_noZ","Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=fpi_noZ->Integral(def_mean_pi-5*def_sig_pi,def_mean_pi+5*def_sig_pi);
 n_pi_noZ=n_pi_noZ/((max_coin_c-min_coin_c)/bin_coin_c);
 mean_pi_noZ=fpi_noZ->GetParameter(1);
 sig_pi_noZ=fpi_noZ->GetParameter(2);
 hcoin_pi_noZ->FillRandom("fpi_noZ",n_pi_noZ*1000.);
 hcoin_pi_noZ->Scale(1./1000.);
 hcoin_wo_bg_fom_noZ->Add(hcoin_wo_bg_fom_noZ,hcoin_pi_noZ,1.,-1.);
 //center_pi=mean_pi_noZ;
 //range_pi=2*sig_pi_noZ;
 //range_pi=0.7;
 center_pi=def_mean_pi;
 range_pi=2*def_sig_pi;
 n_pi_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_pi-range_pi),hcoin_wo_bg_fom_noZ->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
 //hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom_noZ->Fit("fk_noZ","Rq0","0",-1,1);
 //n_k_noZ=fk_noZ->GetParameter(0);
 mean_k_noZ=fk_noZ->GetParameter(1);
// mean_k_noZ=0;
 sig_k_noZ=fk_noZ->GetParameter(2);
 //center_k=mean_k_noZ;
 center_k=0.;
 //range_k=2*sig_k_noZ;
 range_k=1.;

 n_k_noZ=hcoin_wo_bg_fom_noZ->Integral(hcoin_wo_bg_fom_noZ->FindBin(center_k-range_k),hcoin_wo_bg_fom_noZ->FindBin(center_k+range_k));

 //-----Fitting as a whole function---------//

 fcoin_noZ =new TF1("fcoin_noZ","gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin_noZ->SetNpx(2000);
 fcoin_noZ->SetTitle("Cointime w/o AC cut;Cointime [ns];Counts [1/56 ns]");
 fcoin_noZ->SetParameters(n_pi_noZ,mean_pi_noZ,sig_pi_noZ,n_k_noZ,mean_k_noZ,sig_k_noZ,n_p_noZ,mean_p_noZ,sig_p_noZ);
 //hcoin_wo_bg_fom_noZ->Fit("fcoin_woAC","Rq0","0",min_coin_c,max_coin_c);
 //n_pi_noZ=fcoin_noZ->GetParameter(0);//Npi_nocut
 //mean_pi_noZ=fcoin_noZ->GetParameter(1);
 //sig_pi_noZ=fcoin_noZ->GetParameter(2);
 //n_k_noZ=fcoin_noZ->GetParameter(3);//Nk_nocut
 //mean_k_noZ=fcoin_noZ->GetParameter(4);
 //sig_k_noZ=fcoin_noZ->GetParameter(5);
 //n_p_noZ=fcoin_noZ->GetParameter(6);//Np_nocut
 //mean_p_noZ=fcoin_noZ->GetParameter(7);
 //sig_p_noZ=fcoin_noZ->GetParameter(8);
cout<<"n_pi_noZ="<<n_pi_noZ<<"n_k_noZ="<<n_k_noZ<<"n_p_noZ="<<n_p_noZ
<<"mean_pi_noZ="<<mean_pi_noZ<<"sig_pi_noZ="<<sig_pi_noZ<<"mean_k_noZ="<<mean_k_noZ<<"sig_k_noZ="<<sig_k_noZ<<"mean_p_noZ="<<mean_p_noZ<<"sig_p_noZ="<<sig_p_noZ<<endl;

 

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom_noZ->Scale(2./80.);
 hmm_pibg_fom_noZ->Scale(0.7/80.);
 hmm_wo_bg_fom_noZ->Add(hmm_L_fom_noZ,hmm_bg_fom_noZ,1.0,-1.0);
 hmm_pi_wobg_fom_noZ->Add(hmm_pi_fom_noZ,hmm_pibg_fom_noZ,1.0,-1.0);


// fmmbg_noZ=new TF1("fmmbg_noZ","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg_noZ=new TF1("fmmbg_noZ",F_Voigt,min_mm,max_mm,4);
 fmmbg_noZ=new TF1("fmmbg_noZ","pol4",-0.05,0.15);
// fmmbg_noZ->SetParameters(5,0.05,0.05,0.01);
// fmmbg_noZ->SetParLimits(0,0.,100000.);//positive
// fmmbg_noZ->SetParLimits(3,0.,100.);//positive
 fmmbg_noZ->SetNpx(2000);
 //fmmbg_noZ->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg_noZ->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg_noZ->SetParameter(1,0.05);
// fmmbg_noZ->SetParameter(2,0.03);
 fL_noZ=new TF1("fL_noZ","gausn(0)",min_mm,max_mm);
 fL_noZ->SetNpx(2000);
 fL_noZ->SetParLimits(2,0.,0.01);
 fS_noZ=new TF1("fS_noZ","gausn(0)",min_mm,max_mm);
 fS_noZ->SetNpx(2000);
 fS_noZ->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm_noZ=new TF1("fmm_noZ","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
 fmm_noZ->SetNpx(2000);
 fmm_noZ->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm_noZ->SetParLimits(0,0.,1000000.);//positive
 fmm_noZ->SetParLimits(3,0.,1000000.);//positive
// fmm_noZ->SetParLimits(3,0.,100.);//positive//Voigt
// fmm_noZ->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm_noZ->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm_noZ->SetParameter(1,def_mean_L);
// fmm_noZ->SetParameter(4,def_mean_S);

 hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","0",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
 mean_L_noZ=fL_noZ->GetParameter(1);
 sig_L_noZ=fL_noZ->GetParameter(2);
 center_L=def_mean_L;
 range_L=2*def_sig_L;

 hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","0",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
 mean_S_noZ=fS_noZ->GetParameter(1);
 sig_S_noZ=fS_noZ->GetParameter(2);
 center_S=def_mean_S;
 range_S=2*def_sig_S;


 //------- Fitting ----------//

cout<<"fmmbg fit start"<<endl;
// hmm_bg_fom_noZ->Fit("fmmbg_noZ","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom_noZ->Fit("fmmbg_noZ","Rq0","",-0.05,0.15);
//double fmmbga = fmmbg_noZ->GetParameter(0);
//double fmmbgb = fmmbg_noZ->GetParameter(1);
//double fmmbgc = fmmbg_noZ->GetParameter(2);
//double fmmbgd = fmmbg_noZ->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm_noZ->FixParameter(1,b);//mean
// fmm_noZ->FixParameter(2,c);//sigma
// fmm_noZ->FixParameter(3,d);//lg
//double e = fmmbg_noZ->GetParameter(4);
//double f = fmmbg_noZ->GetParameter(5);
// fmm_noZ->SetParameter(0,a);
// fmm_noZ->SetParameter(1,b);
// fmm_noZ->SetParameter(2,c);
// fmm_noZ->SetParameter(3,d);
// fmm_noZ->SetParameter(4,e);
// fmm_noZ->SetParameter(5,f);
// fmm_noZ->SetParameter(6,500);
 fmm_noZ->SetParameter(1,def_mean_L);
 fmm_noZ->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm_noZ->SetParameter(2,def_sig_L);
 fmm_noZ->SetParLimits(2,0.,2*def_sig_L);
// fmm_noZ->SetParameters(9,100);
 fmm_noZ->SetParameter(4,def_mean_S);
 fmm_noZ->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm_noZ->SetParameter(5,def_sig_S);
 fmm_noZ->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom_noZ->Fit("fmm_noZ","Rq0","0",-0.05,0.15);
double fmm_noZpar0 = fmm_noZ->GetParameter(0);cout<<"fmm_noZ[0]="<<fmm_noZpar0<<endl;//area(L)
double fmm_noZpar1 = fmm_noZ->GetParameter(1);cout<<"fmm_noZ[1]="<<fmm_noZpar1<<endl;//mean(L)
double fmm_noZpar2 = fmm_noZ->GetParameter(2);cout<<"fmm_noZ[2]="<<fmm_noZpar2<<endl;//sigma(L)
double fmm_noZpar3 = fmm_noZ->GetParameter(3);cout<<"fmm_noZ[3]="<<fmm_noZpar3<<endl;//area(S)
double fmm_noZpar4 = fmm_noZ->GetParameter(4);cout<<"fmm_noZ[4]="<<fmm_noZpar4<<endl;//mean(S)
double fmm_noZpar5 = fmm_noZ->GetParameter(5);cout<<"fmm_noZ[5]="<<fmm_noZpar5<<endl;//sigma(S)
double fmm_noZpar6 = fmm_noZ->GetParameter(6);cout<<"fmm_noZ[6]="<<fmm_noZpar6<<endl;//poly_const
double fmm_noZpar7 = fmm_noZ->GetParameter(7);cout<<"fmm_noZ[7]="<<fmm_noZpar7<<endl;//poly_x
double fmm_noZpar8 = fmm_noZ->GetParameter(8);cout<<"fmm_noZ[8]="<<fmm_noZpar8<<endl;//poly_x^2
double fmm_noZpar9 = fmm_noZ->GetParameter(9);cout<<"fmm_noZ[9]="<<fmm_noZpar9<<endl;//poly_x^3
double fmm_noZpar10 = fmm_noZ->GetParameter(10);cout<<"fmm_noZ[10]="<<fmm_noZpar10<<endl;//poly_x^4
 fmmbg_noZ->SetParameters(fmm_noZpar6,fmm_noZpar7,fmm_noZpar8,fmm_noZpar9,fmm_noZpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fL_noZ","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L_noZ=fL_noZ->GetParameter(0);
// mean_L_noZ=fL_noZ->GetParameter(1);
// sig_L_noZ=fL_noZ->GetParameter(2);
 mean_L_noZ=def_mean_L;
 mean_S_noZ=def_mean_S;
 sig_L_noZ=def_sig_L;
 sig_S_noZ=def_sig_S;
 n_L_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_L-range_L),hmm_wo_bg_fom_noZ->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L_noZ<<endl;
 double integralL=fmmbg_noZ->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom_noZ->FindBin(center_L+range_L)-hmm_wo_bg_fom_noZ->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L_noZ=n_L_noZ-integralL;
cout<<"after(L):: "<<n_L_noZ<<endl;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,5)-pow(mean_L_noZ-2*sig_L_noZ,5))*a/5;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,4)-pow(mean_L_noZ-2*sig_L_noZ,4))*b/4;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,3)-pow(mean_L_noZ-2*sig_L_noZ,3))*c/3;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,2)-pow(mean_L_noZ-2*sig_L_noZ,2))*d/2;
//n_L_noZ-=(pow(mean_L_noZ+2*sig_L_noZ,1)-pow(mean_L_noZ-2*sig_L_noZ,1))*e;
//
 n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(center_S-range_S),hmm_wo_bg_fom_noZ->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S_noZ<<endl;
 double integralS=fmmbg_noZ->Integral(center_S-range_S,center_S+range_S);
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom_noZ->FindBin(center_S+range_S)-hmm_wo_bg_fom_noZ->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
 if(integralS>0)n_S_noZ-=integralS;
cout<<"after(S):: "<<n_S_noZ<<endl;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,5)-pow(mean_S_noZ-2*sig_S_noZ,5))*a/5;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,4)-pow(mean_S_noZ-2*sig_S_noZ,4))*b/4;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,3)-pow(mean_S_noZ-2*sig_S_noZ,3))*c/3;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,2)-pow(mean_S_noZ-2*sig_S_noZ,2))*d/2;
//n_S_noZ-=(pow(mean_S_noZ+2*sig_S_noZ,1)-pow(mean_S_noZ-2*sig_S_noZ,1))*e;

 cout<<"n_L"<<n_L_noZ<<endl;
cout<<"mean_L"<<mean_L_noZ<<endl;
cout<<"sig_L"<<sig_L_noZ<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom_noZ->Fit("fS_noZ","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S_noZ=fS_noZ->GetParameter(0);
// mean_S_noZ=fS_noZ->GetParameter(1);
// mean_S_noZ=def_mean_S;
// sig_S_noZ=def_sig_S;
// sig_S_noZ=fS_noZ->GetParameter(2);
// n_S_noZ=hmm_wo_bg_fom_noZ->Integral(hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ-2*sig_S_noZ),hmm_wo_bg_fom_noZ->FindBin(mean_S_noZ+2*sig_S_noZ));
 cout<<"n_S"<<n_S_noZ<<endl;
cout<<"mean_S"<<mean_S_noZ<<endl;
cout<<"sig_S"<<sig_S_noZ<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//
 //----------------------------------------------------------------//

	for(int i=0;i<nth;i++){

			
//-----Background subtraction-----//
//------------cointime------------//
//cout<<"BG subtraction cointime"<<endl;
 hcoin_bg_fom[i]->Scale(40./80.);
 hcoin_wo_bg_fom[i]->Add(hcoin_k_fom[i],hcoin_bg_fom[i],1.0,-1.0);

 fp[i]=new TF1(Form("fp[%d]",i),"gausn(0)",min_coin_c,max_coin_c);
 fp[i]->SetNpx(2000);
 fpi[i] =new TF1(Form("fpi[%d]",i),"gausn(0)+gausn(3)",min_coin_c,max_coin_c);
 fpi[i]->SetNpx(2000);
 fk[i]=new TF1(Form("fk[%d]",i),"gausn(0)",min_coin_c,max_coin_c);
 fk[i]->SetNpx(2000);
//cout<<"fp fit start"<<endl;
 hcoin_wo_bg_fom[i]->Fit(Form("fp[%d]",i),"Rq0","0",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 //n_p[i]=fp[i]->GetParameter(0);
 mean_p[i]=fp[i]->GetParameter(1);
 sig_p[i]=fp[i]->GetParameter(2);
 n_p[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_p-range_p),hcoin_wo_bg_fom[i]->FindBin(center_p+range_p));
//cout<<"fpi fit start"<<endl;
 fpi[i]->SetParameters(10000.,def_mean_pi,def_sig_pi,2000.,def_mean_pi,def_sig_pi);
 hcoin_wo_bg_fom[i]->Fit(Form("fpi[%d]",i),"Rq0","0",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i]=fpi[i]->Integral(min_coin_c,max_coin_c);
 n_pi[i]=n_pi[i]/((max_coin_c-min_coin_c)/bin_coin_c);
 //n_pi[i]=n_pi[i]/(2*range_pi/(hcoin_wo_bg_fom[i]->FindBin(center_pi+range_pi)-hcoin_wo_bg_fom[i]->FindBin(center_pi-range_pi)));
 mean_pi[i]=fpi[i]->GetParameter(1);
cout << mean_pi[i] << endl;
 sig_pi[i]=fpi[i]->GetParameter(2);
 hcoin_pi[i]->FillRandom(Form("fpi[%d]",i),n_pi[i]);
 //hcoin_wo_bg_fom[i]->Add(hcoin_wo_bg_fom[i],hcoin_pi[i],1.,-1.);
n_pi[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_pi-range_pi),hcoin_wo_bg_fom[i]->FindBin(center_pi+range_pi));
//cout<<"fk fit start"<<endl;
// hcoin_wo_bg_fom[i]->Fit(Form("fk[%d][%d][%d]",i,j,l),"Rq0","0",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_wo_bg_fom[i]->Fit(Form("fk[%d]",i),"Rq0","0",-1,1);
// n_k[i]=fk[i]->GetParameter(0);
 mean_k[i]=fk[i]->GetParameter(1);
 sig_k[i]=fk[i]->GetParameter(2);
 n_k[i]=hcoin_wo_bg_fom[i]->Integral(hcoin_wo_bg_fom[i]->FindBin(center_k-range_k),hcoin_wo_bg_fom[i]->FindBin(center_k+range_k));

 // n_k[i][th1][1]=hcoin_k_fom[i][th1]->Integral(hcoin_k_fom[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_k_fom[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 fcoin[i] =new TF1(Form("fcoin[%d]",i),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i]->SetNpx(2000);
 fcoin[i]->SetParameters(n_pi[i],mean_pi[i],sig_pi[i],n_k[i],mean_k[i],sig_k[i],n_p[i],mean_p[i],sig_p[i]);
 //hcoin_wo_bg_fom[i]->Fit(Form("fcoin[%d][%d][%d]",i,j,l),"Rq0","0",min_coin_c,max_coin_c);
 //n_pi[i]=fcoin[i]->GetParameter(0);//Npi_nocut
 //mean_pi[i]=fcoin[i]->GetParameter(1);
 //sig_pi[i]=fcoin[i]->GetParameter(2);
 //n_k[i]=fcoin[i]->GetParameter(3);//Nk_nocut
 //mean_k[i]=fcoin[i]->GetParameter(4);
 //sig_k[i]=fcoin[i]->GetParameter(5);
 //n_p[i]=fcoin[i]->GetParameter(6);//Np_nocut
 //mean_p[i]=fcoin[i]->GetParameter(7);
 //sig_p[i]=fcoin[i]->GetParameter(8);

cout<<"n_pi["<<i<<"]="<<n_pi[i]<<endl;
cout<<"n_k["<<i<<"]="<<n_k[i]<<endl;
cout<<"n_p["<<i<<"]="<<n_p[i]<<endl;

//----------------------------------------------//
//--	Missing Mass  Start     ----------------//
//----------------------------------------------//
 hmm_bg_fom[i]->Scale(2./80.);
 hmm_pibg_fom[i]->Scale(0.7/80.);
 hmm_wo_bg_fom[i]->Add(hmm_L_fom[i],hmm_bg_fom[i],1.0,-1.0);
 hmm_pi_wobg_fom[i]->Add(hmm_pi_fom[i],hmm_pibg_fom[i],1.0,-1.0);


// fmmbg[i]=new TF1("fmmbg[i]","gausn(0)+gausn(3)",min_mm,max_mm);
 //fmmbg[i]=new TF1("fmmbg[i]",F_Voigt,min_mm,max_mm,4);
 fmmbg[i]=new TF1("fmmbg[i]","pol4",-0.05,0.15);
// fmmbg[i]->SetParameters(5,0.05,0.05,0.01);
// fmmbg[i]->SetParLimits(0,0.,100000.);//positive
// fmmbg[i]->SetParLimits(3,0.,100.);//positive
 fmmbg[i]->SetNpx(2000);
 //fmmbg[i]->SetParameters(100,0.05,0.03,10,0.05,0.03);//test.list
// fmmbg[i]->SetParameters(300,0.05,1.2,30,0.1,0.02);//small.list
// fmmbg[i]->SetParameter(1,0.05);
// fmmbg[i]->SetParameter(2,0.03);
 fL[i]=new TF1("fL[i]","gausn(0)",min_mm,max_mm);
 fL[i]->SetNpx(2000);
 fL[i]->SetParLimits(2,0.,0.01);
 fS[i]=new TF1("fS[i]","gausn(0)",min_mm,max_mm);
 fS[i]->SetNpx(2000);
 fS[i]->SetParLimits(2,0.,0.01);//sigma limit in order not to mix with bg
 
 fmm[i]=new TF1("fmm[i]","gausn(0)+gausn(3)+pol4(6)",-0.05,0.15);
 fmm[i]->SetNpx(2000);
 fmm[i]->SetTitle("Missing Mass w/o AC cut;Coin time [ns];Counts [1/56 ns]");
 fmm[i]->SetParLimits(0,0.,1000000.);//positive
 fmm[i]->SetParLimits(3,0.,1000000.);//positive
// fmm[i]->SetParLimits(3,0.,100.);//positive//Voigt
// fmm[i]->SetParLimits(4,0.,100000.);//positive//Voigt
// fmm[i]->SetParLimits(7,0.,100000.);//positive//Voigt
// fmm[i]->SetParameter(1,def_mean_L);
// fmm[i]->SetParameter(4,def_mean_S);



 //------- Fitting ----------//

//cout<<"fmmbg fit start"<<endl;
// mmbg_eff[i]->Fit("fmmbg[i]","Rq0","",min_mm,max_mm);
// hmm_pi_wobg_fom[i]->Fit("fmmbg[i]","Rq0","",-0.05,0.15);
//double fmmbga = fmmbg[i]->GetParameter(0);
//double fmmbgb = fmmbg[i]->GetParameter(1);
//double fmmbgc = fmmbg[i]->GetParameter(2);
//double fmmbgd = fmmbg[i]->GetParameter(3);
//cout<<"0:1:2:3:4:5="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<endl;//"::"<<e<<"::"<<f<<endl;
// fmm[i]->FixParameter(1,b);//mean
// fmm[i]->FixParameter(2,c);//sigma
// fmm[i]->FixParameter(3,d);//lg
//double e = fmmbg[i]->GetParameter(4);
//double f = fmmbg[i]->GetParameter(5);
// fmm[i]->SetParameter(0,a);
// fmm[i]->SetParameter(1,b);
// fmm[i]->SetParameter(2,c);
// fmm[i]->SetParameter(3,d);
// fmm[i]->SetParameter(4,e);
// fmm[i]->SetParameter(5,f);
// fmm[i]->SetParameter(6,500);
 fmm[i]->SetParameter(1,def_mean_L);
 fmm[i]->SetParLimits(1,def_mean_L-def_sig_L,def_mean_L+def_sig_L);
 fmm[i]->SetParameter(2,def_sig_L);
 fmm[i]->SetParLimits(2,0.,2*def_sig_L);
// fmm[i]->SetParameters(9,100);
 fmm[i]->SetParameter(4,def_mean_S);
 fmm[i]->SetParLimits(4,def_mean_S-def_sig_S,def_mean_S+def_sig_S);
 fmm[i]->SetParameter(5,def_sig_S);
 fmm[i]->SetParLimits(5,0.,2*def_sig_S);
 hmm_wo_bg_fom[i]->Fit("fmm[i]","Rq0","0",-0.05,0.15);
//double fmmpar0 = fmm[i]->GetParameter(0);//cout<<"fmm[0]="<<fmmpar0<<endl;//area(L)
//double fmmpar1 = fmm[i]->GetParameter(1);//cout<<"fmm[1]="<<fmmpar1<<endl;//mean(L)
//double fmmpar2 = fmm[i]->GetParameter(2);//cout<<"fmm[2]="<<fmmpar2<<endl;//sigma(L)
//double fmmpar3 = fmm[i]->GetParameter(3);//cout<<"fmm[3]="<<fmmpar3<<endl;//area(S)
//double fmmpar4 = fmm[i]->GetParameter(4);//cout<<"fmm[4]="<<fmmpar4<<endl;//mean(S)
//double fmmpar5 = fmm[i]->GetParameter(5);//cout<<"fmm[5]="<<fmmpar5<<endl;//sigma(S)
mean_L[i]=fmm[i]->GetParameter(1);
sig_L[i]=fmm[i]->GetParameter(2);
mean_S[i]=fmm[i]->GetParameter(4);
sig_S[i]=fmm[i]->GetParameter(5);
double fmmpar6 = fmm[i]->GetParameter(6);//cout<<"fmm[6]="<<fmmpar6<<endl;//poly_const
double fmmpar7 = fmm[i]->GetParameter(7);//cout<<"fmm[7]="<<fmmpar7<<endl;//poly_x
double fmmpar8 = fmm[i]->GetParameter(8);//cout<<"fmm[8]="<<fmmpar8<<endl;//poly_x^2
double fmmpar9 = fmm[i]->GetParameter(9);//cout<<"fmm[9]="<<fmmpar9<<endl;//poly_x^3
double fmmpar10 = fmm[i]->GetParameter(10);//cout<<"fmm[10]="<<fmmpar10<<endl;//poly_x^4
 fmmbg[i]->SetParameters(fmmpar6,fmmpar7,fmmpar8,fmmpar9,fmmpar10);
//cout<<"0:1:2:3:4:5(as a total func)="<<a<<"::"<<b<<"::"<<c<<"::"<<d<<"::"<<e<<"::"<<f<<endl;
//
//cout<<"fL fit start"<<endl;
// hmm_wo_bg_fom[i]->Fit("fL[i]","Rq0","",def_mean_L-3*def_sig_L,def_mean_L+3*def_sig_L);
// n_L[i]=fL[i]->GetParameter(0);
// mean_L[i]=fL[i]->GetParameter(1);
// sig_L[i]=fL[i]->GetParameter(2);
 n_L[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(center_L-range_L),hmm_wo_bg_fom[i]->FindBin(center_L+range_L));
cout<<"before(L):: "<<n_L[i]<<endl;
 double integralL=fmmbg[i]->Integral(center_L-range_L,center_L+range_L);
 integralL=integralL/(2*range_L/(hmm_wo_bg_fom[i]->FindBin(center_L+range_L)-hmm_wo_bg_fom[i]->FindBin(center_L-range_L)));
cout<<"integralL="<<integralL<<endl;
 if(integralL>0)n_L[i]=n_L[i]-integralL;
else cout<<"negative BG: ("<<i<<")"<<endl;
cout<<"after(L):: "<<n_L[i]<<endl;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],5)-pow(mean_L[i]-2*sig_L[i],5))*a/5;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],4)-pow(mean_L[i]-2*sig_L[i],4))*b/4;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],3)-pow(mean_L[i]-2*sig_L[i],3))*c/3;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],2)-pow(mean_L[i]-2*sig_L[i],2))*d/2;
//n_L[i]-=(pow(mean_L[i]+2*sig_L[i],1)-pow(mean_L[i]-2*sig_L[i],1))*e;
//
 n_S[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(center_S-range_S),hmm_wo_bg_fom[i]->FindBin(center_S+range_S));
cout<<"before(S):: "<<n_S[i]<<endl;
 double integralS=fmmbg[i]->Integral(center_S-range_S,center_S+range_S);
 integralS=integralS/(2*range_S/(hmm_wo_bg_fom[i]->FindBin(center_S+range_S)-hmm_wo_bg_fom[i]->FindBin(center_S-range_S)));
cout<<"integralS="<<integralS<<endl;
 if(integralS>0)n_S[i]-=integralS;
else cout<<"negative BG: ("<<i<<")"<<endl;
cout<<"after(S):: "<<n_S[i]<<endl;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],5)-pow(mean_S[i]-2*sig_S[i],5))*a/5;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],4)-pow(mean_S[i]-2*sig_S[i],4))*b/4;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],3)-pow(mean_S[i]-2*sig_S[i],3))*c/3;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],2)-pow(mean_S[i]-2*sig_S[i],2))*d/2;
//n_S[i]-=(pow(mean_S[i]+2*sig_S[i],1)-pow(mean_S[i]-2*sig_S[i],1))*e;

 cout<<"n_L"<<n_L[i]<<endl;
cout<<"mean_L"<<mean_L[i]<<endl;
cout<<"sig_L"<<sig_L[i]<<endl;

////cout<<"fS fit start"<<endl;
// hmm_wo_bg_fom[i]->Fit("fS[i]","Rq0","",def_mean_S-3*def_sig_S,def_mean_S+3*def_sig_S);
//// n_S[i]=fS[i]->GetParameter(0);
// mean_S[i]=fS[i]->GetParameter(1);
// mean_S[i]=def_mean_S;
// sig_S[i]=def_sig_S;
// sig_S[i]=fS[i]->GetParameter(2);
// n_S[i]=hmm_wo_bg_fom[i]->Integral(hmm_wo_bg_fom[i]->FindBin(mean_S[i]-2*sig_S[i]),hmm_wo_bg_fom[i]->FindBin(mean_S[i]+2*sig_S[i]));
 cout<<"n_S"<<n_S[i]<<endl;
cout<<"mean_S"<<mean_S[i]<<endl;
cout<<"sig_S"<<sig_S[i]<<endl;
//----------------------------------------------//
//--	Missing Mass  End     ------------------//
//----------------------------------------------//

	
	
	if(n_pi[i]>0.){}else{n_pi[i]=1.;}
	if(n_k[i]>0.){}else{n_k[i]=1.;}
	if(n_p[i]>0.){}else{n_p[i]=1.;}
	if(n_L[i]>0.){}else{n_L[i]=1.;}
	if(n_S[i]>0.){}else{n_S[i]=1.;}




	if(n_pi[i]>n_pi_noZ){
		for(int fill=0;fill<n_pi_noZ;fill++){h_pisr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_pi[i];fill++){h_pisr1->Fill(zver[i]+0.0025);}}

	if(n_k[i]>n_k_noZ){
		for(int fill=0;fill<n_k_noZ;fill++){h_ksr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_k[i];fill++){h_ksr1->Fill(zver[i]+0.0025);}}
	if(n_p[i]>n_p_noZ){
		for(int fill=0;fill<n_p_noZ;fill++){h_psr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_p[i];fill++){h_psr1->Fill(zver[i]+0.0025);}}

	if(n_L[i]>n_L_noZ){
		for(int fill=0;fill<n_L_noZ;fill++){h_Lsr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_L[i];fill++){h_Lsr1->Fill(zver[i]+0.0025);}}
	if(n_S[i]>n_S_noZ){
		for(int fill=0;fill<n_S_noZ;fill++){h_Ssr1->Fill(zver[i]+0.0025);}
	}
	else {for(int fill=0;fill<n_S[i];fill++){h_Ssr1->Fill(zver[i]+0.0025);}}

	for(int fill=0;fill<n_pi_noZ;fill++) h_pitot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_k_noZ;fill++) h_ktot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_p_noZ;fill++) h_ptot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_L_noZ;fill++) h_Ltot1->Fill(zver[i]+0.0025);
	for(int fill=0;fill<n_S_noZ;fill++) h_Stot1->Fill(zver[i]+0.0025);





	

	}//for i
	

 return 0;

}


