bool RHRS=true;
//bool RHRS=false;
const  int nfoil =10;
const double selection_width = 0.01; // event selection width for z


void makehist_ang(){



double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125};
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};


  

  string rname ="./ang_LHRS_5th_0601.root";
  if(RHRS)rname ="./ang_RHRS_5th_0602.root";
  TFile* f = new TFile(rname.c_str());
  TTree* T =(TTree*)f->Get("T");

  //==== Set Branch ======//
  double ss_x,ss_y;
  int zfoil;
  double theta[100],phi[100],theta_c[100],phi_c[100],vz[100],vz_c[100];
  double Lcer_asum_c;
  double Ra1_asum_c,Ra2_asum_c;
  double DRT1,DRT5,DRT4;
  T->SetBranchAddress("DR.T1",&DRT1);
  T->SetBranchAddress("DR.T4",&DRT4);
  T->SetBranchAddress("DR.T5",&DRT5);
  T->SetBranchAddress("ss_x",&ss_x);
  T->SetBranchAddress("ss_y",&ss_y);
  T->SetBranchAddress("zfoil",&zfoil);
  T->SetBranchAddress("L.cer.asum_c",&Lcer_asum_c);
  T->SetBranchAddress("R.a1.asum_c",&Ra1_asum_c);
  T->SetBranchAddress("R.a2.asum_c",&Ra2_asum_c);
  if(RHRS){
  T->SetBranchAddress("R.tr.vz",vz);
  T->SetBranchAddress("R.tr.vz_tuned",vz_c);
  T->SetBranchAddress("R.tr.tg_th",theta);
  T->SetBranchAddress("R.tr.tg_th_tuned",theta_c);
  T->SetBranchAddress("R.tr.tg_ph",phi);
  T->SetBranchAddress("R.tr.tg_ph_tuned",phi_c);
  }else{
  T->SetBranchAddress("L.tr.vz",vz);
  T->SetBranchAddress("L.tr.vz_tuned",vz_c);
  T->SetBranchAddress("L.tr.tg_th",theta);
  T->SetBranchAddress("L.tr.tg_th_tuned",theta_c);
  T->SetBranchAddress("L.tr.tg_ph",phi);
  T->SetBranchAddress("L.tr.tg_ph_tuned",phi_c);
  }


  int ENum = T->GetEntries();
  

  
  

  //=== out root file =======//
  string rname1 ="LHRS_angclaib.root";
  if(RHRS) rname1 ="RHRS_angclaib.root";
  rname1 ="test.root";
  TFile* f1= new TFile(rname1.c_str(),"recreate");
  TTree* Tnew =new TTree("T","reacreate");
  bool LPID, RPID;
  Tnew->Branch("zfoil",&zfoil,"zfoil/I");
  Tnew->Branch("ss_x",&ss_x,"ss_x/D");
  Tnew->Branch("ss_y",&ss_y,"ss_y/D");
  Tnew->Branch("vz",vz,"vz[100]/D");
  Tnew->Branch("vz_c",vz_c,"vz_c[100]/D");
  Tnew->Branch("theta",theta,"theta[100]/D");
  Tnew->Branch("theta_c",theta_c,"theta_c[100]/D");
  Tnew->Branch("phi",phi,"phi[100]/D");
  Tnew->Branch("phi_c",phi_c,"phi_c[100]/D");
  Tnew->Branch("LPID",&LPID,"LPID/B");
  Tnew->Branch("RPID",&RPID,"RPID/B");
  Tnew->Branch("DRT1",&DRT1,"DRT1/D");
  Tnew->Branch("DRT4",&DRT4,"DRT4/D");
  Tnew->Branch("DRT5",&DRT5,"DRT5/D");

  //  TH2F* hss =(TH2F*)f->Get("h3_c");
  //  hss->SetName("hss");
  //  hss->Write();

  TH2F* hss_z[10];
  TH2F* hang_z[10];
  TH2F* hss;
  TH1F* hz[10];
  TH1F* hz_select;


  double min_z = -0.2;
  double max_z =  0.2;
  int    bin_z =  1000;
  double min_ss = -10.;
  double max_ss =  10.;
  int    bin_ss = 100;
  double min_ang = -0.1;
  double max_ang =  0.1;
  int    bin_ang =  200.;

  hz_select = new TH1F("hz_select","",bin_z,min_z,max_z);
  hss       = new TH2F("hss","Sieve Slit Pattern ; ss_y [cm] ; ss_x [cm]",bin_ss,min_ss,max_ss,bin_ss,min_ss,max_ss);
  
  for(int i=0;i<nfoil;i++){
    hss_z[i] = new TH2F(Form("hss_z%d",i),Form("SS pattern with foil %d select; ss_y [cm] ; ss_x [cm]",i),bin_ss,min_ss,max_ss,bin_ss,min_ss,max_ss);
    hang_z[i] = new TH2F(Form("hang_z%d",i),Form("Theta vs Phi with foil %d select; phi [rad] ; theta [cm]",i),bin_ang,min_ang,max_ang,bin_ang,min_ang,max_ang);
    
  }

  
  /*
  TH1F* hz_select =(TH1F*)f->Get("hz_all");
  hz_select->SetName("hz_select");
  hz_select->Write();
  for(int i=0;i<10;i++){
    hss_z[i] =(TH2F*)f->Get(Form("hss_%d",i));
    hss_z[i]->SetName(Form("hss_%d",i));
    hss_z[i]->Write();
    hang_z[i] =(TH2F*)f->Get(Form("hang_%d",i));
    hang_z[i]->SetName(Form("hang_%d",i));
    hang_z[i]->Write();    
    hz[i] =(TH1F*)f->Get(Form("hz_%d",i));
    hz[i]->SetName(Form("hz_%d",i));
    hz[i]->Write();      }
  */
  
  for(int i=0;i<ENum;i++){

    T->GetEntry(i);
    LPID = false;
    RPID = false;
    if(Lcer_asum_c>2000.) LPID =true;
    if(Ra1_asum_c> 1. && Ra2_asum_c > 2000.) RPID = true;

    if((RHRS && RPID) || (!RHRS && LPID)){
     
      for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<vz_c[0] && vz_c[0]<fcent[j]+selection_width){
	  zfoil=j;
	  hss_z[j]  -> Fill(ss_y,ss_x);
	  hang_z[j] -> Fill(phi_c[0],theta_c[0]);
	  
	}
      }
      hss->Fill(ss_y,ss_x);
      Tnew->Fill();
    }
    
  }// ent loop

  Tnew->Write();

  for(int i=0;i<nfoil;i++){

    hss_z[i]->Write();
    hang_z[i] ->Write();
  }
  hss->Write();
  
    

}
