//-----------------------------//
//--  Path Length Calc.      --//
//--	from Cointime        --//
//-----------------------------//

//K. Okuyama (Jun. 22, 2021)

double PI=4.*atan(1.);

double F_Voigt( double *x, double *par )
  {
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
    double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    return val;
  }

double Dgauss( double *x, double *par )
  {
    // par[0] : scale1 
    // par[1] : location1
    // par[2] : sigma1
    // par[3] : scale2 
    // par[4] : location2
    // par[5] : sigma2
    double val = par[0] * TMath::Gaus(x[0],par[1],par[2]);
    val += par[3] * TMath::Gaus(x[0],par[4],par[5]);
    return val;
  }


double fcoin_template( double *x, double *par , int shift, int num)
{
  return par[num] * TMath::Gaus(x[0],par[num+1]-2.0*shift,par[num+2]);//Lambda Gaussian
}

double fcoin_acc( double *x, double *par)
{
	double val=0.;
	for(int i=0;i<=80;i++){
	//val += par[num]*TMath::Gaus(x[0],par[num+1]+par[num+12]*i-20,par[num+2]);
    val += par[0] * (1.-par[5]) * TMath::Voigt(x[0]-3.179-par[1]-par[4]*(double)i+30.,par[2],par[3],4);
    //val += par[0] * par[5] * TMath::Voigt(x[0]+7.925-par[1]-par[4]*i+30.,par[6],par[7],4);
    val += par[0] * par[5] * 0.6 * TMath::Gaus(x[0],-8.303+par[1]+par[4]*(double)i-30.,par[6])/(sqrt(2.*PI)*par[6]);
    val += par[0] * par[5] * 0.4 * TMath::Gaus(x[0],-7.463+par[1]+par[4]*(double)i-30.,par[7])/(sqrt(2.*PI)*par[7]);
	}
	return val;
}

//double fcoin_total( double *x, double *par ){
//
//	if(reject_flag&&x[0]>-12.&&x[0]<10.){TF1::RejectPoint();return 0;}
//	//return fcoin_template(x,par,-10,0)+expgaus2(x,par,6)+expgaus2(x,par,10);//+expgaus2(x,par,14);
//	else return fcoin_acc(x,par)+par[8]*TMath::Voigt(x[0]-par[9],par[10],par[11],4)+par[12]*TMath::Voigt(x[0]-par[13],par[14],par[15],4)+par[16]*TMath::Gaus(x[0],par[17],par[18])+par[19]*TMath::Gaus(x[0],par[20],par[21]);
//
//}

void path_from_ct(){
//ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
	string pdfname = "temp.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("../h2all5.root","read");//input file of all H2 run(default: h2all4.root)
  TFile *file = new TFile("../h2all_2020Nov.root","read");//input file of all H2 run(default: h2all4.root)
  TTree *tree = (TTree*)file->Get("tree_out");

    
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);



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
 int bin_coin_c=(int)((max_coin_c-min_coin_c)/tdc_time);
 const double min_mm=-0.1;//GeV/c^2
 const double max_mm=0.2;//GeV/c^2
 int bin_mm=(max_mm-min_mm)/0.002; //Counts/2 MeV
 bin_mm=(int)bin_mm;


//---For Efficiency Result---//
 double center_pi, center_k, center_p, center_L, center_S;
 double range_pi, range_k, range_p, range_L, range_S;
 double n_pi_nocut, n_k_nocut, n_p_nocut, n_L_nocut, n_S_nocut;
 double const_pi_nocut, const_k_nocut, const_p_nocut, const_L_nocut, const_S_nocut;
 double n_pi_best, n_k_best, n_p_best, n_L_best, n_S_best;
 double const_pi_best, const_k_best, const_p_best, const_L_best, const_S_best;
 double n_pi[100], n_k[100],n_p[100], n_L[100], n_S[100];
 double mean_pi_nocut, mean_k_nocut, mean_p_nocut, mean_L_nocut, mean_S_nocut;
 double mean_pi_best, mean_k_best, mean_p_best, mean_L_best, mean_S_best;
 double mean_pi[100], mean_k[100],mean_p[100], mean_L[100], mean_S[100];
 double sig_pi_nocut, sig_k_nocut, sig_p_nocut, sig_L_nocut, sig_S_nocut;
 double sig_pi_best, sig_k_best, sig_p_best, sig_L_best, sig_S_best;
 double sig_pi[100], sig_k[100],sig_p[100], sig_L[100], sig_S[100];

//---------------------------------------//
//               Branch                  //
//---------------------------------------//

 int NLtr, NRtr, Ls2_pad[100], Rs2_pad[100];
 double ct, ct_eff;

	double L_tr_chi2;
	double L_tr_x, L_tr_y, L_tr_th, L_tr_ph;
	double L_tr_p;
	double L_tr_tg_th, L_tr_tg_ph;
	double L_tr_vz;
	double L_tr_vz_saved;
	double R_tr_chi2;
	double R_tr_x, R_tr_y, R_tr_th, R_tr_ph;
	double R_tr_p;
	double R_tr_tg_th, R_tr_tg_ph;
	double R_tr_vz;
	double L_mom, R_mom, B_mom; 
	double L_ene, R_ene, B_ene; 
	double ac1sum, ac2sum;//NPE SUM

	tree->SetBranchStatus("*",0);
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


  //TH1F* hcoin  = new TH1F("hcoin","",40000/56,-20.,20.);
  TH1F* hcoin  = new TH1F("hcoin","",130000/56,-30.,100.);
  

  TH1F* h_test  = new TH1F("h_test","",1000,1.8,2.4);

  bool L_Tr = false;
  bool L_FP = false;
  bool R_Tr = false;
  bool R_FP = false;
  bool ct_cut = false;
  bool event_selection = false;
  bool event_selection_nocut = false;
  double z_par[100], ac_par[100], ct_par[100];
  double z2_par[100][100], ac2_par[100][100];
  double rf_bunch=2.0;//ns (RF bunch structure)
  const double kcenter = 0.0;
  double mh = ML;//hypernuclei
  double mt = Mp;//target mass
  double B_p, L_p, R_p;//Momentum


  //tree->Draw(">>elist" , "fabs(ct_orig[0][0])<1.0");
  //tree->Draw(">>elist" , "fabs(ct_orig)<1.0");//ctsum (does NOT dintinguish #track)
  //TEventList *elist = (TEventList*)gROOT->FindObject("elist");
  //int ENum = elist->GetN(); 
  int ENum = tree->GetEntries(); 
cout<<"Entries: "<<ENum<<endl;
  int time_div=ENum/25;
  if(ENum<100000)time_div=10000;


	time_t start, end;
	start = time(NULL);
	time(&start);

  for(int i=0;i<ENum;i++){
	//tree->GetEntry(elist->GetEntry(i));
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


	

		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;
		
		if(abs(R_mom-1.90)<0.02)hcoin->Fill(ct);

}//ENum

	double chisq = 0.;
	double dof = 0.;


		TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
//Pion Fitting
cout<<"Pion Fitting is performed as a first step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	 TF1 *fcoin_pion= new TF1("fcoin_pion","gausn",1.,5.);
	 fcoin_pion->SetNpx(20000);
	 fcoin_pion->SetParameter(0,30000.);
	 fcoin_pion->SetParameter(1,3.18);
	 fcoin_pion->SetParameter(2,0.4);
	 hcoin->Fit("fcoin_pion","","",2.,4.);
	 double pion_par0 = fcoin_pion->GetParameter(0);
	 double pion_par1 = fcoin_pion->GetParameter(1);
	 double pion_par2 = fcoin_pion->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Pion Fitting after BG subtraction (fcoin_pion) "<<endl;
	 cout<<"Pion Fitting w/o BG subtraction (fcoin_pion) "<<endl;
	 chisq = fcoin_pion->GetChisquare();
	 dof  = fcoin_pion->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
cout<<"Proton Fitting is performed as a next step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	 TF1 *fcoin_proton= new TF1("fcoin_proton","gausn",-9,-7.);
	 fcoin_proton->SetNpx(20000);
	 fcoin_proton->SetParameter(0,5000.);
	 fcoin_proton->SetParameter(1,-8.0);
	 fcoin_proton->SetParameter(2,0.4);
	 hcoin->Fit("fcoin_proton","","",-8.,-6.);
	 double proton_par0 = fcoin_proton->GetParameter(0);
	 double proton_par1 = fcoin_proton->GetParameter(1);
	 double proton_par2 = fcoin_proton->GetParameter(2);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Proton Fitting after BG subtraction (fcoin_proton) "<<endl;
	 cout<<"Proton Fitting w/o BG subtraction (fcoin_proton) "<<endl;
	 chisq = fcoin_proton->GetChisquare();
	 dof  = fcoin_proton->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;


	TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
	c5->DrawFrame(-20.,0.,20.,40000.);
	hcoin->SetLineColor(kAzure);
	hcoin->Draw("same");
	fcoin_pion->SetLineColor(kOrange);
	fcoin_pion->Draw("same");
	fcoin_proton->SetLineColor(kRed);
	fcoin_proton->Draw("same");

	double mom = 1.90;//GeV/c
	double betapi = mom/sqrt(Mpi*Mpi+mom*mom);
	double betap  = mom/sqrt(Mp*Mp+mom*mom);
	double t_diff = proton_par1-pion_par1;
	double pathlen = betapi*betap*LightVelocity*t_diff/(betap-betapi);

	double t_diff_sig = sqrt(proton_par2*proton_par2+pion_par2*pion_par2);
	double t_diff_sig2 = sqrt(proton_par2*proton_par2/proton_par0*0.056+pion_par2*pion_par2/pion_par0*0.056);
	double pathlen_sig = betapi*betap*LightVelocity*t_diff_sig/(betap-betapi);
	double pathlen_sig2 = betapi*betap*LightVelocity*t_diff_sig2/(betap-betapi);
	cout<<"Npi="<<pion_par0/0.056<<endl;
	cout<<"Np="<<proton_par0/0.056<<endl;
	cout<<"betapi="<<betapi<<endl;
	cout<<"betap="<<betap<<endl;
	cout<<"t_diff="<<abs(t_diff)<<" [ns]"<<endl;
	cout<<"pathlen="<<pathlen<<"+/-"<<abs(pathlen_sig)<<" [m]"<<endl;
	cout<<"pathlen2="<<pathlen<<"+/-"<<abs(pathlen_sig2)<<" [m]"<<endl;

}//path_from_ct
