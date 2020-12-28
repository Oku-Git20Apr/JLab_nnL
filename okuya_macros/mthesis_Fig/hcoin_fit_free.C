//-----------------------------//
//--  Cointime Histo.        --//
//-----------------------------//
//
//K. Okuyama (Dec. 5, 2020)
//K. Okuyama (Dec. 7, 2020)

#include "TF1.h"
#include "TH1.h"
bool reject_flag = false;
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

double fcoin_total( double *x, double *par ){

	if(reject_flag&&x[0]>-12.&&x[0]<10.){TF1::RejectPoint();return 0;}
	//return fcoin_template(x,par,-10,0)+expgaus2(x,par,6)+expgaus2(x,par,10);//+expgaus2(x,par,14);
	else return fcoin_acc(x,par)+par[8]*TMath::Voigt(x[0]-par[9],par[10],par[11],4)+par[12]*TMath::Voigt(x[0]-par[13],par[14],par[15],4)+par[16]*TMath::Gaus(x[0],par[17],par[18])+par[19]*TMath::Gaus(x[0],par[20],par[21]);

}

void hcoin_fit_free(){
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




  TH1F* h1  = new TH1F("h1","",400,-20.,20.0);
  h1->GetXaxis()->SetTitle("coin time (ns)");
  h1->GetYaxis()->SetTitle("Counts / 100 ps");
  h1->GetXaxis()->SetRangeUser(-14.0,17.);
  double xmin = -0.1, xmax = 0.2; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_L_fom_best  = new TH1F("hmm_L_fom_best","hmm_L_fom_best",xbin,xmin,xmax);
  hmm_L_fom_best->GetXaxis()->SetTitle("M_{x} - M_{#Lambda} (GeV/c^{2})");
  hmm_L_fom_best->GetYaxis()->SetTitle("Counts / MeV");
  hmm_L_fom_best->SetLineColor(1);
  TH1F* hmm_L_fom_nocut  = new TH1F("hmm_L_fom_nocut","hmm_L_fom_nocut",xbin,xmin,xmax);
//  TH1F* hmm_bg_fom_best  = new TH1F("hmm_bg_fom_best","hmm_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_best  = new TH1F("hmm_wo_bg_fom_best","hmm_wo_bg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_wo_bg_fom_nocut  = new TH1F("hmm_wo_bg_fom_nocut","hmm_wo_bg_fom_nocut",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_best  = new TH1F("hmm_pi_wobg_fom_best","hmm_pi_wobg_fom_best",xbin,xmin,xmax);
  TH1F* hmm_pi_wobg_fom_nocut  = new TH1F("hmm_pi_wobg_fom_nocut","hmm_pi_wobg_fom_nocut",xbin,xmin,xmax);
  TH1F* hm2   = (TH1F*)hmm_L_fom_best->Clone("hm2");
  TH1F* hm4   = (TH1F*)hmm_L_fom_best->Clone("hm4");

  TH1D* h_theta_ee = new TH1D("h_theta_ee", "theta_ee",1000,0.1,0.35);
  TH1D* h_phi_ee = new TH1D("h_phi_ee", "phi_ee",1000,0.,PI);
  TH1D* h_theta_ek = new TH1D("h_theta_ek", "theta_ek",1000,0.1,0.35);
  TH1D* h_phi_ek = new TH1D("h_phi_ek", "phi_ek",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_g = new TH1D("h_theta_g", "theta_g",1000,0.1,0.35);
  TH1D* h_phi_g = new TH1D("h_phi_g", "phi_g",1000,3*PI/2-1.,3*PI/2+1.);
  TH1D* h_theta_gk_lab = new TH1D("h_theta_gk_lab", "theta_gk_lab",1000,0.,0.2);
  TH1D* h_theta_gk_cm = new TH1D("h_theta_gk_cm", "theta_gk_cm",1000,0.,0.3);
  TH1D* h_cos_gk_lab = new TH1D("h_cos_gk_lab", "cos_gk_lab",1000,0.97,1.0);
  TH1D* h_cos_gk_cm = new TH1D("h_cos_gk_cm", "cos_gk_cm",1000,0.8,1.0);
  TH1D* h_mom_g = new TH1D("h_mom_g", "mom_g",1000,1.8,2.5);
  TH1D* h_qsq = new TH1D("h_qsq", "Q^2",1000,0.,0.8);
  TH1D* h_w = new TH1D("h_w", "W",1000,0.,0.8);
  TH2D* h_thph_ee = new TH2D("h_thph_ee", "theta_ee:phi_ee" ,1000,0.1,0.35,1000,PI/2-1.,PI/2+1.);
  TH2D* h_thph_ek = new TH2D("h_thph_ek", "theta_ek:phi_ek" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);
  TH2D* h_thph_g = new TH2D("h_thph_g", "theta_g:phi_g" ,1000,0.1,0.35,1000,3*PI/2-1.,3*PI/2+1.);

  //TH1F* hcoin  = new TH1F("hcoin","",40000/56,-20.,20.);
  TH1F* hcoin  = new TH1F("hcoin","",130000/56,-30.,100.);
  TH1F* hcoin_bg  = new TH1F("hcoin_bg","",130000/56,-30.,100.);
  TH1F* hcoin_wobg  = new TH1F("hcoin_wobg","",130000/56,-30.,100.);
  TH1F* hcoin_new  = new TH1F("hcoin_new","",130000/56,-30.,100.);
  TH1F* hcoin_new_bg  = new TH1F("hcoin_new_bg","",130000/56,-30.,100.);
  TH1F* hcoin_newwobg  = new TH1F("hcoin_newwobg","",130000/56,-30.,100.);
  TH1F* hcoin_pi  = new TH1F("hcoin_pi","",130000/56,-30.,100.);
  TH1F* hcoin_pi_bg  = new TH1F("hcoin_pi_bg","",130000/56,-30.,100.);
  TH1F* hcoin_piwobg  = new TH1F("hcoin_piwobg","",130000/56,-30.,100.);
  TH1F* hcoin_p  = new TH1F("hcoin_p","",130000/56,-30.,100.);
  TH1F* hcoin_p_bg  = new TH1F("hcoin_p_bg","",130000/56,-30.,100.);
  TH1F* hcoin_pwobg  = new TH1F("hcoin_pwobg","",130000/56,-30.,100.);
  TH1F* hcoin_best  = new TH1F("hcoin_best","",130000/56,-30.,100.);
  TH1F* hcoin_best_bg  = new TH1F("hcoin_best_bg","",130000/56,-30.,100.);
  TH1F* hcoin_bestwobg  = new TH1F("hcoin_bestwobg","",130000/56,-30.,100.);
  TH1F* hcoin_strict  = new TH1F("hcoin_strict","",130000/56,-30.,100.);
  TH1F* hcoin_strict_bg  = new TH1F("hcoin_strict_bg","",130000/56,-30.,100.);
  TH1F* hcoin_strictwobg  = new TH1F("hcoin_strictwobg","",130000/56,-30.,100.);
  
  h1 ->SetLineColor(2);
  h1->SetLineWidth(2);

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

/***********************************/
/**	What parameter do you change? **/
/***********************************/
  bool zcut_flag  = true;
  bool accut_flag = false;
  bool fpcut_flag = false;
  bool trcut_flag = false;
/***********************************/
/***********************************/
/***********************************/



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


	

		bool ac_cut=false;
		bool ac_cut_new=false;
		bool ac_cut_pi=false;
		bool ac_cut_p=false;
		bool best_cut=false;
		bool strict_cut=false;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.)ac_cut=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.)ac_cut_new=true;
		if(ac1sum>3.75&&ac2sum>10.)ac_cut_pi=true;
		if(ac1sum<3.75&&ac2sum<0.001&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2)ac_cut_p=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2)best_cut=true;
		if(ac1sum<3.75&&ac2sum>3.&&ac2sum<10.&&fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2)strict_cut=true;

		if(fabs(ct)<1.006)ct_cut=true;
		else ct_cut=false;
		//if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection=true;
		else event_selection=false;

		event_selection_nocut=false;
		if(zcut_flag){
			if(ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(accut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&R_Tr&&R_FP&&L_Tr&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(fpcut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_Tr&&L_Tr)event_selection_nocut=true;
			else event_selection_nocut=false;
		}
		if(trcut_flag){
			if(fabs(L_tr_vz-R_tr_vz)<0.025&&fabs(R_tr_vz+L_tr_vz)<0.2&&ac1sum<3.75&&ac2sum>3.&&ac2sum<20.&&R_FP&&L_FP)event_selection_nocut=true;
			else event_selection_nocut=false;
		}

		if(ac_cut)hcoin->Fill(ct);
		if(ac_cut_new)hcoin_new->Fill(ct);
		if(ac_cut_pi)hcoin_pi->Fill(ct);
		if(ac_cut_p)hcoin_p->Fill(ct);
		//if(best_cut)hcoin_best->Fill(ct);
		if(strict_cut)hcoin_best->Fill(ct);
		if(strict_cut)hcoin_strict->Fill(ct);
		if(ac_cut&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_new&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_new_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_pi&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_pi_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(ac_cut_p&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_p_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(best_cut&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_best_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime
		if(strict_cut&&20.<ct && ct<60.){
		   double ct_ = ct;
		       while(1){
		     if(-20.<ct && ct<20.){
		   	 hcoin_strict_bg->Fill(ct);
		   	 break;}
		          else if(ct<-20.){ct=ct+40.;}
		          else if(20.<ct){ct=ct-40.;}
		    }
		   ct = ct_;
		   }//cointime

//	    //===== Right Hand Coordinate ====//
//	    //th and phi are originally meant tan(theta) and tan(phi),
//	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama
//	    double R_pz = R_mom/sqrt(1.0*1.0 + pow((R_tr_tg_th), 2.0) + pow(( R_tr_tg_ph),2.0) );
//	    double R_px = R_pz * (R_tr_tg_th );
//	    double R_py = R_pz * ( R_tr_tg_ph );
//
//	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(( L_tr_tg_th ), 2.0) + pow(( L_tr_tg_ph),2.0));
//	    double L_px = L_pz * ( L_tr_tg_th );
//	    double L_py = L_pz * ( L_tr_tg_ph );
//
//	    double B_E =sqrt(B_mom*B_mom + Me*Me);
//	    double R_E =sqrt(R_mom*R_mom + MK*MK);
//	    double L_E =sqrt(L_mom*L_mom + Me*Me);
//
//		TLorentzVector L_4vec;//Left
//		TLorentzVector R_4vec;//Right
//		TLorentzVector B_4vec;//Beam
//		TLorentzVector T_4vec;//Target
//		TLorentzVector G_4vec;//Gamma (Virtual Photon)
//		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
//        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
//        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
//        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);
//
//	    L_4vec.RotateX( -13.2/180.*PI );
//	    R_4vec.RotateX(  13.2/180.*PI );
//
//        double mass,mm;
//		TLorentzVector Missing;
//		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
//		mass = Missing.M();
//	    mm=mass - mh;//shift by ML
//		if(event_selection&&ct_cut)hmm_L_fom_best->Fill(mm);
//		if(event_selection_nocut&&ct_cut)hmm_L_fom_nocut->Fill(mm);



}//ENum

//BG subtraction
     //hcoin_bg->Scale(40./80.);
     //hcoin_new_bg->Scale(40./80.);
     //hcoin_pi_bg->Scale(40./80.);
     //hcoin_p_bg->Scale(40./80.);
     //hcoin_best_bg->Scale(40./80.);
     //hcoin_strict_bg->Scale(40./80.);
     hcoin_bg->Sumw2();
     hcoin_new_bg->Sumw2();
     hcoin_pi_bg->Sumw2();
     hcoin_p_bg->Sumw2();
     hcoin_best_bg->Sumw2();
     hcoin_strict_bg->Sumw2();
     hcoin_wobg->Add(hcoin,hcoin_bg,1.0,-1.0);
     hcoin_newwobg->Add(hcoin_new,hcoin_new_bg,1.0,-1.0);
     hcoin_piwobg->Add(hcoin_pi,hcoin_pi_bg,1.0,-1.0);
     hcoin_pwobg->Add(hcoin_p,hcoin_p_bg,1.0,-1.0);
     hcoin_bestwobg->Add(hcoin_best,hcoin_best_bg,1.0,-1.0);
     //hcoin_strictwobg->Add(hcoin_strict,hcoin_strict_bg,1.0,-1.0);
	double chisq, chisq2, chisq3;
	double dof, dof2, dof3;


/*---------------*/
/*   strict Cut  */
/*---------------*/
cout<<"%%%%%%%%%%%%%%%%%%%%"<<endl;
cout<<"%%% strict Cut %%%%%%%"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%"<<endl;


		TCanvas *c4 = new TCanvas("c4", "c4", 1000, 800);
//Pion Fitting
cout<<"Pion Fitting is performed as a first step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	 TF1 *fcoin_first_strict= new TF1("fcoin_first_strict",F_Voigt,1.,5.,4);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
	 fcoin_first_strict->SetNpx(20000);
	 fcoin_first_strict->SetParameter(0,95500.);
	 fcoin_first_strict->SetParameter(1,3.18);
	 fcoin_first_strict->SetParameter(2,0.2);
	 fcoin_first_strict->SetParameter(3,0.4);
	 //hcoin_piwobg->Fit("fcoin_first_strict","","",2.,4.);
	 //change
	 hcoin_pi->Fit("fcoin_first_strict","","",2.,4.);
	 double pion_strict_par0 = fcoin_first_strict->GetParameter(0);
	 double pion_strict_par1 = fcoin_first_strict->GetParameter(1);
	 double pion_strict_par2 = fcoin_first_strict->GetParameter(2);
	 double pion_strict_par3 = fcoin_first_strict->GetParameter(3);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Pion Fitting after BG subtraction (fcoin_first_strict) "<<endl;
	 cout<<"Pion Fitting w/o BG subtraction (fcoin_first_strict) "<<endl;
	 chisq = fcoin_first_strict->GetChisquare();
	 dof  = fcoin_first_strict->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Proton Fitting
cout<<"Proton Fitting is performed as a next step"<<endl;
cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	 TF1 *fcoin_firstp_strict= new TF1("fcoin_firstp_strict",Dgauss,-9,-7.,6);
    // par[0] : area
    // par[1] : location
    // par[2] : gaussian sigma
    // par[3] : lorentz fwhm
		fcoin_firstp_strict->SetNpx(20000);
	 	fcoin_firstp_strict->SetParameter(0,2000.);
	 	fcoin_firstp_strict->FixParameter(1,-8.303);
	 	//fcoin_firstp_strict->SetParLimits(1,-8.8,-8.1);
	 	fcoin_firstp_strict->SetParameter(2,0.4);
	 	fcoin_firstp_strict->SetParLimits(2,0.,0.8);
	 	fcoin_firstp_strict->SetParameter(3,2000.);
	 	fcoin_firstp_strict->FixParameter(4,-7.463);
	 	//fcoin_firstp_strict->SetParLimits(4,-8.1,-7.5);
	 	fcoin_firstp_strict->SetParameter(5,0.4);
	 	fcoin_firstp_strict->SetParLimits(5,0.,0.8);
	 //hcoin_pwobg->Fit("fcoin_firstp_strict","","",-10.,-6.);
	 //change
	 hcoin_p->Fit("fcoin_firstp_strict","","",-10.,-6.);
	 double proton_strict_par0 = fcoin_firstp_strict->GetParameter(0);
	 double proton_strict_par1 = fcoin_firstp_strict->GetParameter(1);
	 double proton_strict_par2 = fcoin_firstp_strict->GetParameter(2);
	 double proton_strict_par3 = fcoin_firstp_strict->GetParameter(5);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 //cout<<"Proton Fitting after BG subtraction (fcoin_firstp_strict) "<<endl;
	 cout<<"Proton Fitting w/o BG subtraction (fcoin_firstp_strict) "<<endl;
	 chisq = fcoin_firstp_strict->GetChisquare();
	 dof  = fcoin_firstp_strict->GetNDF();
	 cout<<"chisq="<<chisq<<endl;
	 cout<<"dof="<<dof<<endl;
	 cout<<"Reduced chi-square = "<<chisq/dof<<endl;

//Accidentals Fitting
	 reject_flag=true;//flag for ignoring pi,K,p
	 TF1 *fcoin_strict= new TF1("fcoin_strict",fcoin_total,-30.,100.,22);
	 fcoin_strict->SetNpx(20000);
	 fcoin_strict->SetParameter(0,20.);//acc
	 fcoin_strict->FixParameter(1,0.);//if par[1]=0, then no meaning. (Not used)
	 fcoin_strict->FixParameter(2,pion_strict_par2);//
	 fcoin_strict->FixParameter(3,pion_strict_par3);//
	 fcoin_strict->FixParameter(4,2.012);//shift
	 fcoin_strict->SetParameter(5,0.5);//pi vs p
	 fcoin_strict->SetParLimits(5,0.,1.);//
	 fcoin_strict->FixParameter(6,proton_strict_par2);//
	 fcoin_strict->FixParameter(7,proton_strict_par3);//
	 fcoin_strict->FixParameter(8,0.);//pi
	 fcoin_strict->FixParameter(9,0.);//pi
	 fcoin_strict->FixParameter(10,0.);//pi
	 fcoin_strict->FixParameter(11,0.);//pi
	 fcoin_strict->FixParameter(12,0.);//k
	 fcoin_strict->FixParameter(13,0.);//k
	 fcoin_strict->FixParameter(14,0.);//k
	 fcoin_strict->FixParameter(15,0.);//k
	 fcoin_strict->FixParameter(16,0.);//p
	 fcoin_strict->FixParameter(17,0.);//p
	 fcoin_strict->FixParameter(18,0.);//p
	 fcoin_strict->FixParameter(19,0.);//p
	 fcoin_strict->FixParameter(20,0.);//p
	 fcoin_strict->FixParameter(21,0.);//p
	 fcoin_strict->SetLineColor(kCyan);
	 hcoin_strict->Fit("fcoin_strict","","",-30.,30.);
	 cout<<"%%%%%Information%%%%%"<<endl;
	 cout<<"Accidentals Fitting [-30,30] (fcoin) "<<endl;
	 chisq2 = fcoin_strict->GetChisquare();
	 dof2  = fcoin_strict->GetNDF();
	 cout<<"chisq="<<chisq2<<endl;
	 cout<<"dof="<<dof2<<endl;
	 cout<<"Reduced chi-square = "<<chisq2/dof2<<endl;
	 hcoin_strict->Draw("");
	 fcoin_strict->Draw("same");
	reject_flag=false;
	 cout<<"MAX ACC. = "<<fcoin_strict->GetMaximum()<<endl;
	 cout<<"MIN ACC. = "<<fcoin_strict->GetMinimum()<<endl;

//Global Fitting
	 fcoin_strict->ReleaseParameter(0);
	 fcoin_strict->ReleaseParameter(1);
	 fcoin_strict->ReleaseParameter(2);
	 fcoin_strict->ReleaseParameter(3);
	 fcoin_strict->ReleaseParameter(4);
	 fcoin_strict->ReleaseParameter(5);
	 fcoin_strict->ReleaseParameter(6);
	 fcoin_strict->ReleaseParameter(7);
	 fcoin_strict->ReleaseParameter(8);
	 fcoin_strict->ReleaseParameter(9);
	 fcoin_strict->ReleaseParameter(10);
	 fcoin_strict->ReleaseParameter(11);
	 fcoin_strict->ReleaseParameter(12);
	 fcoin_strict->ReleaseParameter(13);
	 fcoin_strict->ReleaseParameter(14);
	 fcoin_strict->ReleaseParameter(15);
	 fcoin_strict->ReleaseParameter(16);
	 fcoin_strict->ReleaseParameter(17);
	 fcoin_strict->ReleaseParameter(18);
	 fcoin_strict->ReleaseParameter(19);
	 fcoin_strict->ReleaseParameter(20);
	 fcoin_strict->ReleaseParameter(21);
//ACC
	 fcoin_strict->FixParameter(0,fcoin_strict->GetParameter(0));
	 fcoin_strict->FixParameter(1,fcoin_strict->GetParameter(1));
	 fcoin_strict->FixParameter(2,fcoin_strict->GetParameter(2));
	 fcoin_strict->FixParameter(3,fcoin_strict->GetParameter(3));
	 fcoin_strict->FixParameter(4,fcoin_strict->GetParameter(4));
	 fcoin_strict->FixParameter(5,fcoin_strict->GetParameter(5));
	 fcoin_strict->FixParameter(6,fcoin_strict->GetParameter(6)); fcoin_strict->FixParameter(7,fcoin_strict->GetParameter(7));
//pion
	 fcoin_strict->SetParameter(8,16000.);//pi scale
	 fcoin_strict->SetParLimits(8,0.,100000.);//pi scale
	 fcoin_strict->SetParameter(9,3.18);
	 fcoin_strict->SetParameter(10,pion_strict_par2);
	 fcoin_strict->SetParLimits(10,0.,0.8);
	 fcoin_strict->SetParameter(11,pion_strict_par3);
	 fcoin_strict->SetParLimits(11,0.,10.);
//kaon
	 fcoin_strict->SetParameter(12,500.);//k scale
	 fcoin_strict->SetParLimits(12,0.,50000.);//k scale
	 fcoin_strict->SetParameter(13,0.);
	 fcoin_strict->SetParLimits(13,-0.5,0.5);
	 //fcoin_strict->SetParameter(14,0.308788-0.00);//change
	 //fcoin_strict->SetParameter(15,0.0976195+0.0);//contami. min
	 //fcoin_strict->FixParameter(15,0.);//contami. max
	 fcoin_strict->SetParameter(14,0.5);
	 fcoin_strict->SetParLimits(14,0.,0.7);
	 fcoin_strict->SetParameter(15,0.4);
	 fcoin_strict->SetParLimits(15,0.,0.7);
//proton
		fcoin_strict->SetNpx(20000);
	 	fcoin_strict->SetParameter(16,500.);
		fcoin_strict->SetParLimits(16,0.,50000.);
	 	fcoin_strict->FixParameter(17,-8.303);
	 	//fcoin_strict->FixParameter(18,proton_strict_par2);
	 	fcoin_strict->SetParameter(18,0.4);
	 	fcoin_strict->SetParLimits(18,0.,0.8);
	 	fcoin_strict->SetParameter(19,500.);
		fcoin_strict->SetParLimits(19,0.,50000.);
	 	fcoin_strict->FixParameter(20,-7.463);
	 	//fcoin_strict->FixParameter(21,proton_strict_par3);
	 	fcoin_strict->SetParameter(21,0.4);
	 	fcoin_strict->SetParLimits(21,0.,0.8);

	 hcoin_strict->Fit("fcoin_strict","","",-20.,20.);

	 TF1 *fcoin_bg_strict= new TF1("fcoin_bg_strict",fcoin_acc,-20.,20.,8);
	 TF1 *fcoin_pi_strict= new TF1("fcoin_pi_strict",F_Voigt,-20.,20.,4);
	 TF1 *fcoin_k_strict= new TF1("fcoin_k_strict",F_Voigt,-20.,20.,4);
	 TF1 *fcoin_p_strict= new TF1("fcoin_p_strict",Dgauss,-20.,20.,6);
	 fcoin_pi_strict->SetNpx(20000);
	 fcoin_k_strict->SetNpx(20000);
	 fcoin_p_strict->SetNpx(20000);
	 fcoin_bg_strict->SetParameter(0,fcoin_strict->GetParameter(0));
	 fcoin_bg_strict->SetParameter(1,fcoin_strict->GetParameter(1));
	 fcoin_bg_strict->SetParameter(2,fcoin_strict->GetParameter(2));
	 fcoin_bg_strict->SetParameter(3,fcoin_strict->GetParameter(3));
	 fcoin_bg_strict->SetParameter(4,fcoin_strict->GetParameter(4));
	 fcoin_bg_strict->SetParameter(5,fcoin_strict->GetParameter(5));
	 fcoin_bg_strict->SetParameter(6,fcoin_strict->GetParameter(6));
	 fcoin_bg_strict->SetParameter(7,fcoin_strict->GetParameter(7));
	 fcoin_pi_strict->SetParameter(0,fcoin_strict->GetParameter(8));
	 fcoin_pi_strict->SetParameter(1,fcoin_strict->GetParameter(9));
	 fcoin_pi_strict->SetParameter(2,fcoin_strict->GetParameter(10));
	 fcoin_pi_strict->SetParameter(3,fcoin_strict->GetParameter(11));
	 fcoin_k_strict->SetParameter(0,fcoin_strict->GetParameter(12));
	 fcoin_k_strict->SetParameter(1,fcoin_strict->GetParameter(13));
	 fcoin_k_strict->SetParameter(2,fcoin_strict->GetParameter(14));
	 fcoin_k_strict->SetParameter(3,fcoin_strict->GetParameter(15));
	 fcoin_p_strict->SetParameter(0,fcoin_strict->GetParameter(16));
	 fcoin_p_strict->SetParameter(1,fcoin_strict->GetParameter(17));
	 fcoin_p_strict->SetParameter(2,fcoin_strict->GetParameter(18));
	 fcoin_p_strict->SetParameter(3,fcoin_strict->GetParameter(19));
	 fcoin_p_strict->SetParameter(4,fcoin_strict->GetParameter(20));
	 fcoin_p_strict->SetParameter(5,fcoin_strict->GetParameter(21));
//
//cout<<"fcoin(4)="<<fcoin->GetParameter(4)<<endl;
//
		c4->SetLogy(0);
		TH1 *frame_strict = c4->DrawFrame(-20.,1.,20.,400.);
		frame_strict->GetXaxis()->SetTitle("Coincidence Time [ns]");
		frame_strict->GetYaxis()->SetTitle("Counts/0.056ns");
		frame_strict->Draw("");
		hcoin_strict->Draw("same");
		fcoin_strict->Draw("same");
		double ktegrated_strict = fcoin_k_strict->Integral(-1.006,1.006);
		double k_strict = hcoin_strict->Integral(hcoin_strict->FindBin(-1.006),hcoin_strict->FindBin(1.006));
		double bgtegrated_strict = fcoin_bg_strict->Integral(-1.006,1.006);
		bgtegrated_strict /= 0.056;
		k_strict = k_strict - bgtegrated_strict;
		double pitegrated_strict = fcoin_pi_strict->Integral(-1.006,1.006);
		cout<<"%%%%%Information%%%%%"<<endl;
	 	chisq3 = fcoin_strict->GetChisquare();
	 	dof3  = fcoin_strict->GetNDF();
	 	cout<<"chisq="<<chisq3<<endl;
	 	cout<<"dof="<<dof3<<endl;
	 	cout<<"Reduced chi-square = "<<chisq3/dof3<<endl;
		cout<<"True (-1.006ns<ct<1.006ns): "<<k_strict<<endl;
		cout<<"Kaon (-1.006ns<ct<1.006ns): "<<ktegrated_strict/0.056<<endl;
		cout<<"Pion (-1.006ns<ct<1.006ns): "<<pitegrated_strict/0.056<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated_strict*100./(k_strict*0.056)<<endl;
		cout<<"Pion Contamination (-1.006ns<ct<1.006ns): "<<pitegrated_strict*100./(ktegrated_strict+pitegrated_strict)<<endl;
		fcoin_pi_strict->SetLineColor(kOrange);
		fcoin_k_strict->SetLineColor(kGreen);
		fcoin_p_strict->SetLineColor(kRed);
		fcoin_strict->Draw("same");
		fcoin_pi_strict->Draw("same");
		fcoin_k_strict->Draw("same");
		fcoin_p_strict->Draw("same");


		TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
		c5->SetLogy(1);
		TH1 *frame2_strict = c5->DrawFrame(-20.,1.,20.,400.);
		//fcoin_bg_strict->Draw("same");
		fcoin_strict->Draw("same");
		fcoin_pi_strict->Draw("same");
		fcoin_k_strict->Draw("same");
		fcoin_p_strict->Draw("same");

		TCanvas *c7 = new TCanvas("c7", "c7", 800, 800);
		TH1 *frame7 = c7->DrawFrame(-20.,1.,20.,40000.);
		//hcoin_piwobg->Draw("");
		hcoin_pi->Draw("same");
		fcoin_first_strict->SetLineColor(kMagenta);
		fcoin_first_strict->Draw("same");
		TCanvas *c8 = new TCanvas("c8", "c8", 800, 800);
		TH1 *frame8 = c8->DrawFrame(-20.,1.,20.,1000.);
		//hcoin_pwobg->Draw("");
		hcoin_p->Draw("same");
		fcoin_firstp_strict->SetLineColor(kMagenta);
		fcoin_firstp_strict->Draw("same");

	//	TCanvas *c9 = new TCanvas("c9", "c9", 800, 800);
	//	hcoin_strict->Fit("gausn","","",-10.,-6.);
	//	hcoin_strict->Fit("gausn","","",-1.,1.);
	//	hcoin_strict->Fit("gausn","","",2.,4.);
	//	TCanvas *c10 = new TCanvas("c10", "c10", 800, 800);
	//	TF1 *fcoinp_dg_strict= new TF1("fcoinp_dg_strict",Dgauss,-20.,20.,6);
	//	fcoinp_dg_strict->SetNpx(20000);
	// 	fcoinp_dg_strict->SetParameter(0,500.);
	// 	fcoinp_dg_strict->SetParameter(1,-8.3);
	// 	fcoinp_dg_strict->SetParameter(2,0.4);
	// 	fcoinp_dg_strict->SetParameter(3,500.);
	// 	fcoinp_dg_strict->SetParameter(4,-7.9);
	// 	fcoinp_dg_strict->SetParameter(5,0.4);
	//	hcoin_strict->Fit("fcoinp_dg_strict","","",-10.,-6.);
		TCanvas *c88 = new TCanvas("c88", "c88", 800, 800);
		TH1 *frame_88 = c88->DrawFrame(0.,0.,1.25,1.);
		TH1F* h_den  = new TH1F("h_den","denominator",25,0.,1.25);
		TH1F* h_num  = new TH1F("h_num","numerator",25,0.,1.25);
		double den, num;
		h_den->SetBinContent(1,0.);
		h_num->SetBinContent(1,0.);
		for(int i=1;i<25;i++){
		//den = fcoin_k_strict->Integral(-1.25,1.25)/0.056;
		den = fcoin_k_strict->Integral(-20.,20.)/0.056;
		num = fcoin_k_strict->Integral((double)i*(-0.05),(double)i*0.05)/0.056;
		h_den->SetBinContent(i+1,den);
		h_num->SetBinContent(i+1,num);
		cout<<"ct="<<(double)i*0.05<<" ns, SR="<<num*100./den<<" %"<<endl;
		}
		cout << "TEfficiency! (Gas SR)" << endl;
		TEfficiency *pEff1;
		if(TEfficiency::CheckConsistency(*h_num,*h_den,"w")){
		pEff1 = new TEfficiency(*h_num,*h_den);
		}
		frame_88->GetXaxis()->SetTitle("X [ns]");
		frame_88->GetYaxis()->SetTitle("K^{+} Survival Ratio");
		frame_88->GetYaxis()->SetDecimals();
		frame_88->Draw("");
		pEff1->Draw("same");

		c4->Print("./pdf/hcoin_fit_free.pdf");
		c88->Print("./pdf/hcoin_kSR.pdf");
}//fit
