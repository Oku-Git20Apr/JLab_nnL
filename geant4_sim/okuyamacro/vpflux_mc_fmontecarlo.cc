#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <string>

#include "TStyle.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TLine.h"
#include "TRandom.h"
using namespace std;
void CanvasSet( TCanvas*, int, double, double, double, double);
void SetTitle( TH2D*, const char*, const char*, const char*);
void SetTitle( TH1D*, const char*, const char*, const char*);
void TreeBranch( TTree*, const char*, double*);

static const double PI = 4.*atan(1.0);

//---------Not used anymore----------//
//---------VP Flux----------//
double vpflux_lab(double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = par[1];	

	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//if(theta>0.230&&theta<0.231){
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;
//}

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	return vpflux;
}

//TF1* func3 = new TF1("func3",vpflux_lab, 0.001, 0.3,1);
//func3->SetNpx(600);
//func3->SetParameter(0,2.1);
//func3->SetLineColor(kGreen);
//func3->SetLineWidth(4);
//func3->Draw("same");

int main(int argc, char** argv){
	//////////////////////
	// command argument //
	//////////////////////
	int option;
	string filename = "test2";
	bool BatchFlag = true;
	bool PDFFlag = true;
	while((option=getopt(argc, argv, "f:bp"))!=-1){
		switch(option){
			case 'f':
				filename = optarg;
				break;
			case 'b':
				cout << "batch mode" << endl;
				BatchFlag = true;
				break;
			case 'p':
				PDFFlag = true;
				break;
			case 'h':
				cout<<"-f : input root filename"<<endl;
				cout<<"-b : execute in batch mode"<<endl;
				cout<<"-p : print pdf file"<<endl;
				return 0;
				break;
			case '?':
				cout << "unknown option: " << option << endl;
				return 0;
				break;
			default:
				cout<<"type -h to see help!!"<<endl;
				return 0;
				break;
		}
	}
  TApplication theApp("App", &argc, argv);
	ostringstream RootFile;
	ostringstream LogFile;
	ostringstream PDFFile;
	RootFile << "../ana/root/" << filename << ".root";	
	LogFile << "../ana/root/" << filename << ".root_Log";	
	PDFFile << "./pdf/vpflux.pdf";	

	///////////////////////
	// General condition //
	///////////////////////
    gROOT->SetStyle("Plain");
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");
    gStyle->SetTitleOffset(0.9, "X");
    gStyle->SetTitleOffset(0.9, "Y");
	gStyle->SetOptStat(0);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	FILE *fp;
	const int SizeOfBuffer = 32;
    char str[SizeOfBuffer];
	double centralmom, centraltheta, thetawidth;
	double centralphi, phiwidth;

	if((fp=fopen(LogFile.str().c_str(), "r"))==NULL){
		std::cerr << "file open fail" << std::endl;
		exit(1);
	}

	while(fgets( str, SizeOfBuffer, fp)){
		if(sscanf( str, "Central Momentum [GeV] = %lf", &centralmom )==1){
		}
		else if(sscanf( str, "Theta central [rad] = %lf", &centraltheta )==1){
		}
		else if(sscanf( str, "Theta Gen. Range [rad] = %lf", &thetawidth )==1){
		}
		else if(sscanf( str, "Phi Gen. Range [rad] = %lf", &phiwidth )==1){
		}
	}
	fclose(fp);
	cout<< "central theta = " << centraltheta*180/PI << " [deg]" << endl;
	cout<< "theta width= " << thetawidth*180/PI << " [deg]" << endl;
	centralphi = 0.;

	TFile *f = new TFile(RootFile.str().c_str());
	TTree *t = (TTree*)f->Get("tree");
	TTree *tree = (TTree*)f->Get("tree");

	double thetamin = 1.*(1.*centraltheta - thetawidth);
	double thetamax = 1.*(1.*centraltheta + thetawidth);
	double phimin = 1.*(centralphi - phiwidth);
	double phimax = 1.*(centralphi + phiwidth);
	double omega = 2.*PI* (1-cos(thetawidth))*1000.; // [msr]
	double MAX = 0.01;
	double volume = omega*MAX/1000.; 
	//double volume = omega/1000.; 
	cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
	cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
	cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
	cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
	cout<< "omega = " << omega << " [msr]" <<endl;
	cout<< "volume = " << volume <<endl;

//	// ======== Cut condition (Aida) =========
//	// =========== E12-15-008 ============
//	TCut Sepi = "0.<EXSepi && EXSepi<50.";
//	TCut Sepe = "0.<EXSepe && EXSepe<50.";
//	TCut Q1i = "-30.<EXQ1i && EXQ1i<30.";
//	TCut Q1e = "-30.<EXQ1e && EXQ1e<30.";
//	TCut Q2i = "-40.<EXQ2i && EXQ2i<40.";
//	TCut Q2e = "-40.<EXQ2e && EXQ2e<40.";
////	TCut Di = "-159.5/2.<EXDi && EXDi<159.5/2.";
//	TCut Di = "-999.<EXDi";
////	TCut De = "-159.5/2.<EXDe && EXDe<159.5/2.";
//	TCut De = "-999.<EXDe";
//	TCut Q3i = "-40.<EXQ3i && EXQ3i<40.";
//	TCut Q3e = "-40.<EXQ3e && EXQ3e<40.";
////	TCut RP = "EVDTrig && -132./2.<EXFP && EXFP<132./2."; // Reference Plane
//	TCut RP = "EVDTrig && EXFP1>-999 && EXFP2>-999"; // Reference Plane
////	TCut ETOF2X = "VDTrig && -150.<EXFP5 && EXFP5<150.";

	// ======== Cut condition (Suzuki) =========
	// "HRS/ana/getHist/src/input/scan.input"
//Q1 Collimator[cm]: 11.959969
//Q2 Collimator[cm]: 30.0
//Q3 Collimator[cm]: 30.0
//Q1in Cut xmin[cm]: -1000.
//Q1in Cut xmax[cm]: 1000.
//Q1in Cut ymin[cm]: -1000.
//Q1in Cut ymax[cm]: 1000
//Z center[cm]: 0.0;
//Z min[cm]: -12.5
//Z max[cm]: 12.5
//Mom max[GeV]: 100.
//Mom min[GeV]: 0.
//PCSM ymax[cm]: 10.0
//PCSM ymin[cm]: -10.0

	// ======== Cut condition =========
	TCut Q1 = "sqrt(EXQ1i*EXQ1i+EYQ1i*EYQ1i)<11.959969";
	TCut Q2 = "sqrt(EXQ2i*EXQ2i+EYQ2i*EYQ2i)<30.";
	TCut Q3 = "sqrt(EXQ3i*EXQ3i+EYQ3i*EYQ3i)<30.";
//====K40====//
	//TCut Q1 = "sqrt(EXQ1i*EXQ1i+EYQ1i*EYQ1i)<10.";
	//TCut Q2 = "sqrt(EXQ2i*EXQ2i+EYQ2i*EYQ2i)<12.";
	//TCut Q3 = "sqrt(EXQ3i*EXQ3i+EYQ3i*EYQ3i)<30.";
//====K40====//
	TCut RP = "EVDTrig && EDCTrig"; // Reference Plane
	TCut Zpos = "-12.5<EZt && EZt<12.5";
	TCut top = "fabs(th-0.225)<0.025&&fabs(ph)<0.25";
	int bin_mom = 150;//expanded?
	double min_mom = 1.8;//nnL
	double max_mom = 2.4;//nnL
	//double min_mom = 2.6;//K40
	//double max_mom = 3.4;//K40
	int bin_th = 150;
	double min_th = 0.10;
	double max_th = 0.35;
	double min_costh = 0.95;//cos
	double max_costh = 1.0;//cos
	double min_ph = -0.4292;//default=0.57
	double max_ph = 0.4292;
	int bin_2D_mom = 50;
	int bin_2D_th = 50;
	int bin_2D_costh = 50;
	int bin_2D_ph = 50;

	TH1D *h_mom_gen = new TH1D( "h_mom_gen", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q1 = new TH1D( "h_mom_q1", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q2 = new TH1D( "h_mom_q2", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q3 = new TH1D( "h_mom_q3", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_rp = new TH1D( "h_mom_rp", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_result = new TH1D( "h_mom_result", "", bin_mom, min_mom, max_mom);
	t->Project( "h_mom_gen", "EMom");
	t->Project( "h_mom_q1", "EMom", Q1);
	t->Project( "h_mom_q2", "EMom", Q2);
	t->Project( "h_mom_q3", "EMom", Q3);
	t->Project( "h_mom_rp", "EMom", RP);
	//t->Project( "h_mom_result", "EMom",RP&&Q1&&Q2&&Q3);
	TH1D *h_th_gen = new TH1D( "h_th_gen", "", bin_th, min_th, max_th);
	TH1D *h_th_q1 = new TH1D( "h_th_q1", "", bin_th, min_th, max_th);
	TH1D *h_th_q2 = new TH1D( "h_th_q2", "", bin_th, min_th, max_th);
	TH1D *h_th_q3 = new TH1D( "h_th_q3", "", bin_th, min_th, max_th);
	TH1D *h_th_rp = new TH1D( "h_th_rp", "", bin_th, min_th, max_th);
	TH1D *h_th_result = new TH1D( "h_th_result", "", bin_th, min_th, max_th);
	t->Project( "h_th_gen", "ETheta");
	t->Project( "h_th_q1", "ETheta", Q1);
	t->Project( "h_th_q2", "ETheta", Q2);
	t->Project( "h_th_q3", "ETheta", Q3);
	t->Project( "h_th_rp", "ETheta", RP);
	t->Project( "h_th_result", "ETheta", RP&&Q1&&Q2&&Q3);



	TH2D *h_mom_th_gen = new TH2D( "h_mom_th_gen", "", bin_2D_mom, min_mom, max_mom, bin_2D_th, min_th, max_th);
	TH2D *h_mom_th = new TH2D( "h_mom_th", "", bin_2D_mom, min_mom, max_mom, bin_2D_th, min_th, max_th);
	t->Project( "h_mom_th_gen", "ETheta:EMom");
	t->Project( "h_mom_th", "ETheta:EMom", RP&&Q1&&Q2&&Q3);

	TH2D *h_th_ph_gen = new TH2D( "h_th_ph_gen", "", bin_2D_th, min_th, max_th, bin_2D_ph, min_ph, max_ph);
	TH2D *h_costh_ph_gen = new TH2D( "h_costh_ph_gen", "", bin_2D_costh, min_costh, max_costh, bin_2D_ph, min_ph, max_ph);
	TH2D *h_th_ph = new TH2D( "h_th_ph", "", bin_2D_th, min_th, max_th, bin_2D_ph, min_ph, max_ph);
	TH2D *h_costh_ph = new TH2D( "h_costh_ph", "", bin_2D_costh, min_costh, max_costh, bin_2D_ph, min_ph, max_ph);
	//t->Project( "h_th_ph_gen", "EPhi:cos(ETheta)");
	//t->Project( "h_th_ph", "EPhi:cos(ETheta)", RP&&Q1&&Q2&&Q3);

	// ========= Solid angle vs. momentum ===========
	TH1D *h_sa_mom_q1 = new TH1D("h_sa_mom_q1", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q2 = new TH1D("h_sa_mom_q2", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q3 = new TH1D("h_sa_mom_q3", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_rp = new TH1D("h_sa_mom_rp", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_result = new TH1D("h_sa_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result = new TH1D("h_vp_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result2 = new TH1D("h_vp_mom_result2", "", bin_mom, min_mom, max_mom);//w/ RHRS Acceptance

	int n1, n2;
	double val;
	double err;
	const double min_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindFirstBinAbove()) - 0.05;
	const double max_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindLastBinAbove()) + 0.05;
	const double min_mom_bin = h_mom_gen->FindBin(min_mom_gen);
	const double max_mom_bin = h_mom_gen->FindBin(max_mom_gen);
	cout <<"mom range: "<< min_mom_gen << " " << max_mom_gen << endl;
	cout <<"mom_bin: "<< min_mom_bin << " " << max_mom_bin << endl;
//----------------------------//
//------------Fill------------//
//----------------------------//
	TH1D *h_vp_mom = new TH1D( "h_vp_mom", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom2 = new TH1D( "h_vp_mom2", "w/ HRS-R Acceptance", bin_mom, min_mom, max_mom);
	TH1D *h_vp_th = new TH1D( "h_vp_th", "", bin_th, min_th, max_th);
	int ENum=0;
	double mom=0.;
	double th=0.;
	double ph=0.;
	double x_q1i=0.;
	double x_q2i=0.;
	double x_q3i=0.;
	double y_q1i=0.;
	double y_q2i=0.;
	double y_q3i=0.;
	char evdtrig=0;
	char edctrig=0;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("EMom",1);tree->SetBranchAddress("EMom",&mom);
	tree->SetBranchStatus("ETheta",1);tree->SetBranchAddress("ETheta",&th);
    tree->SetBranchStatus("EPhi"    ,1);tree->SetBranchAddress("EPhi"    ,&ph     );
    tree->SetBranchStatus("EXQ1i"   ,1);tree->SetBranchAddress("EXQ1i"   ,&x_q1i  );
    tree->SetBranchStatus("EXQ2i"   ,1);tree->SetBranchAddress("EXQ2i"   ,&x_q2i  );
    tree->SetBranchStatus("EXQ3i"   ,1);tree->SetBranchAddress("EXQ3i"   ,&x_q3i  );
    tree->SetBranchStatus("EYQ1i"   ,1);tree->SetBranchAddress("EYQ1i"   ,&y_q1i  );
    tree->SetBranchStatus("EYQ2i"   ,1);tree->SetBranchAddress("EYQ2i"   ,&y_q2i  );
    tree->SetBranchStatus("EYQ3i"   ,1);tree->SetBranchAddress("EYQ3i"   ,&y_q3i  );
    tree->SetBranchStatus("EVDTrig" ,1);tree->SetBranchAddress("EVDTrig" ,&evdtrig);
    tree->SetBranchStatus("EDCTrig" ,1);tree->SetBranchAddress("EDCTrig" ,&edctrig);
    ENum = tree->GetEntries();
  for(int i=0;i<ENum;i++){
    tree->GetEntry(i);
    if(i%100000==0)cout<<i<<" / "<<ENum<<endl;
	double cosine=cos(th);
	h_th_ph_gen->Fill(th,ph);
	h_costh_ph_gen->Fill(cosine,ph);
	
	if((sqrt(x_q1i*x_q1i+y_q1i*y_q1i)<11.959969) &&
	   (sqrt(x_q2i*x_q2i+y_q2i*y_q2i)<30.) &&
	   (sqrt(x_q3i*x_q3i+y_q3i*y_q3i)<30.) &&	
		evdtrig && edctrig){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = mom;//[GeV]
	double theta = th;	

	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double deltaE = Einc - Escat;
	double q2=Qsq+deltaE*deltaE;
	double kg=deltaE-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	double k = MAX*gRandom->Uniform();
//cout<<"vpflux vs k = "<<vpflux<<" : "<<k<<endl;
		//if(fabs(cos(th)-0.9808)<0.003&&fabs(ph-0.25)<0.125){
	//	if(mom>2.1&&mom<2.15){
		h_mom_result->Fill(mom);
		h_th_ph->Fill(th,ph);
		h_costh_ph->Fill(cos(th),ph);
	 //if(vpflux>k&&fabs(th-0.225)<0.025&&fabs(ph)<0.25){//new top-quality (6msr)
	 if(vpflux>k){//Full
		h_vp_mom->Fill(mom);
		//h_vp_mom->SetBinContent(h_vp_mom->FindBin(mom),vpflux);
		//if(fabs(mom-2.125)<0.025){//new top-quality (Lambda)
		//if(fabs(mom-2.075)<0.025){//new top-quality (Sigma)
		if(mom>2.1){//Full (Lambda)
		//if(mom<2.12){//Full (Sigma)
		h_vp_mom2->Fill(mom);
		h_vp_th->Fill(th);
	 //}//Full_S
	 }//Full_L
	 //}//Top_S
	 //}//Top_L
	}//VP Flux
	 //}//6msr
	 //}//HRS-R acceptance
	 //}//6msr (w/ cos(th) cut)
	}//Acceptance
	}//ENum

double vpflux_tot_mom=0.;
double vpflux_tot_mom_err=0.;
double vpflux_tot_th=0.;
double n1_tot=0.;
double n2_tot=0.;
double n1n2_tot=0.;
	for(int i=0; i<bin_mom; i++){
		n1 = 0;
		n1 = (int)h_mom_gen->GetBinContent(i+1);

		// === Q1 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q1->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q1->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q1->SetBinError(i+1, err);

		// === Q2 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q2->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q2->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q2->SetBinError(i+1, err);

		// === Q3 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q3->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q3->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q3->SetBinError(i+1, err);

		// === RP ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_rp->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_rp->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_rp->SetBinError(i+1, err);

		// === all cuts ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_result->GetBinContent(i+1);
		cout<<"nbin_mom:"<<i<<endl;
		cout<<"n1="<<n1<<endl;
		cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_result->SetBinError(i+1, err);
		n1_tot+=n1;
		n2_tot+=n2;
		n1n2_tot+=n1*n2;

		// === VP Flux ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom->GetBinContent(i+1);
		cout<<"nbin_mom:"<<i<<endl;
		cout<<"n1="<<n1<<endl;
		cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1);
		h_vp_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result->SetBinError(i+1, err);

		// === VP Flux (w/ RHRS Acceptance)===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom2->GetBinContent(i+1);
		cout<<"nbin_mom:"<<i<<endl;
		cout<<"n1="<<n1<<endl;
		cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1);
		h_vp_mom_result2->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result2->SetBinError(i+1, err);

		double bin_width=(max_mom-min_mom)/bin_mom;//GeV/c
		vpflux_tot_mom+=val*bin_width;
		vpflux_tot_mom_err+=err*bin_width*err*bin_width;
	}
		double ave_n1=n1_tot/bin_mom;
		double ave_n2=n2_tot/bin_mom;
		double ave_n1n2=n1n2_tot/bin_mom;
		cout<<"ave_n1="<<ave_n1<<endl;
		cout<<"ave_n2="<<ave_n2<<endl;
		cout<<"ave_n1n2="<<ave_n1n2<<endl;
		cout<<"Cov(n1,n2)="<<ave_n1n2-ave_n1*ave_n2<<endl;
		cout<<"rho(n1,n2)="<<(ave_n1n2-ave_n1*ave_n2)/(sqrt(n1_tot*n2_tot))<<endl;

	// ========= Solid angle vs. theta ===========
	TH1D *h_sa_th_q1 = new TH1D("h_sa_th_q1", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q2 = new TH1D("h_sa_th_q2", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q3 = new TH1D("h_sa_th_q3", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_rp = new TH1D("h_sa_th_rp", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_result = new TH1D("h_sa_th_result", "", bin_th, min_th, max_th);
	TH1D *h_vp_th_result = new TH1D("h_vp_th_result", "", bin_th, min_th, max_th);

	const double min_th_gen = h_th_gen->GetBinCenter(h_th_gen->FindFirstBinAbove()) - 0.02;
	const double max_th_gen = h_th_gen->GetBinCenter(h_th_gen->FindLastBinAbove()) + 0.02;
	const double min_th_bin = h_th_gen->FindBin(min_th_gen);
	const double max_th_bin = h_th_gen->FindBin(max_th_gen);
	cout <<"th range: "<< min_th_gen << " " << max_th_gen << endl;
	cout <<"th bin: "<< min_th_bin << " " << max_th_bin << endl;

	for(int i=0; i<bin_th; i++){
		n1 = 0;
		n1 = (int)h_th_gen->GetBinContent(i+1);

		// === Q1 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q1->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q1->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q1->SetBinError(i+1, err);

		// === Q2 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q2->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q2->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q2->SetBinError(i+1, err);

		// === Q3 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q3->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q3->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q3->SetBinError(i+1, err);

		// === RP ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_rp->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_rp->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_rp->SetBinError(i+1, err);

		// === all cuts ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_result->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_result->SetBinError(i+1, err);

		// === VP Flux ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_th->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1);
		h_vp_th_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_th_result->SetBinError(i+1, err);

		double bin_width=(max_th-min_th)/bin_th;
		vpflux_tot_th+=val*bin_width*2*PI*sin(h_vp_th_result->GetXaxis()->GetBinCenter(i));
	}


	// === Momentum vs. Theta ===
	TH2D *h_sa_mom_th = new TH2D("h_sa_mom_th", "", bin_2D_mom, min_mom, max_mom, bin_2D_th, min_th, max_th);
	for(int i=0; i<bin_2D_mom; i++){
		for(int j=0; j<bin_2D_th; j++){
			n1 = 0;
			n1 = (int)h_mom_th_gen->GetBinContent(i+1, j+1);
			n2 = 0;
			val = 0.;
			n2 = (int)h_mom_th->GetBinContent(i+1, j+1);
			if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
//		if(val!=0)cout<< n1 << " " << n2 << " " << val << endl;
			h_sa_mom_th->SetBinContent(i+1, j+1, val);
		}
	}

	// === Theta vs. Phi ===
	TH2D *h_sa_costh_ph = new TH2D("h_sa_costh_ph", "", bin_2D_costh, min_costh, max_costh, bin_2D_ph, min_ph, max_ph);
	for(int i=0; i<bin_2D_th; i++){
		for(int j=0; j<bin_2D_ph; j++){
			n1 = 0;
			n1 = (int)h_costh_ph_gen->GetBinContent(i+1, j+1);
			n2 = 0;
			val = 0.;
			n2 = (int)h_costh_ph->GetBinContent(i+1, j+1);
			if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
//		if(val!=0)cout<< n1 << " " << n2 << " " << val << endl;
			h_sa_costh_ph->SetBinContent(i+1, j+1, val);
		}
	}
cout<<"vpflux_tot_mom="<<vpflux_tot_mom<<endl;
cout<<"vpflux_tot_mom_err="<<sqrt(vpflux_tot_mom_err)<<endl;
cout<<"vpflux_tot_th="<<vpflux_tot_th<<endl;
	
////////////////////
// Draw histgrams //
////////////////////
	TH1D *hframe;
	TLine *lmom = new TLine( centralmom, 0., centralmom, omega);
	TLine *lth = new TLine( centraltheta, 0., centraltheta, omega);
	lmom->SetLineColor(4);
	lth->SetLineColor(4);

	TCanvas *c1 = new TCanvas("c1", "Momentum Acceptance (Exit)");
	TCanvas *c2 = new TCanvas("c2", "Angular Acceptance (Exit)");
	TCanvas *c3 = new TCanvas("c3", "Momentum Acceptance Expansion");
	TCanvas *c4 = new TCanvas("c4", "Angular Acceptance Expansion");
	TCanvas *c5 = new TCanvas("c5", "Momentum vs. Theta");
	TCanvas *c6 = new TCanvas("c6", "Theta vs. Phi");
	TCanvas *c7 = new TCanvas("c7", "VP Flux");


	//======= Momentum Acceptance (Entrance) ======
	c1->Divide(2,2);
//	gPad->SetGrid();
	c1->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q1 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q1->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q2 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q2->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q3 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q3->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum (w/ EVDTrig, EDCTrig", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_rp->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	//======= Angular Acceptance (Entrance) ======
	c2->Divide(2,2);
//	gPad->SetGrid();
	c2->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Theta at Q1 Entrance", "Theta [rad]", "Solid Angle [msr]");
	h_sa_th_q1->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Theta at Q2 Entrance", "Theta [rad]", "Solid Angle [msr]");
	h_sa_th_q2->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Theta at Q3 Entrance", "Theta [rad]", "Solid Angle [msr]");
	h_sa_th_q3->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Theta (w/ EVDTrig, EDCTrig", "Theta [rad]", "Solid Angle [msr]");
	h_sa_th_rp->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	//======= Momentum Acceptance Expansion ======
	c3->Divide(2,2);
	c3->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	//SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum (w/ all Cuts)", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_result->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
//	h_sa_mom_result->Draw("same");
	h_sa_mom_result->Draw();
//	lmom->Draw("same");
	c3->cd(2);
	SetTitle(h_mom_gen, "Solid Angle vs. Momentum (generated)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_gen->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_gen->Draw();
	c3->cd(3);
	SetTitle(h_mom_result, "Solid Angle vs. Momentum (cut)", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_result->Draw();
	c3->cd(4);
	h_mom_gen->Draw();
	h_mom_result->SetLineColor(kAzure);
	h_mom_result->Draw("same");

	//======= Angular Acceptance Expansion ======
	c4->Divide(2,2);
	c4->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(h_sa_th_result, "Solid Angle vs. Theta (w/ all Cuts)", "Theta [rad]", "Solid Angle [msr]");
//	h_sa_th_result->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_th_result->GetXaxis()->SetNdivisions(505, kFALSE);
//	h_sa_th_result->Draw("same");
	h_sa_th_result->Draw();
//	lth->Draw("same");
	c4->cd(2);
	SetTitle(h_th_gen, "Solid Angle vs. Theta (generated)", "Theta [rad]", "Solid Angle [msr]");
	h_th_gen->GetXaxis()->SetNdivisions(505, kFALSE);
	h_th_gen->Draw();
	c4->cd(3);
	SetTitle(h_th_result, "Solid Angle vs. Theta (cut)", "Theta [rad]", "Solid Angle [msr]");
	h_th_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_th_result->Draw();

	//======= Theta vs. Momentum  ======
	c5->cd();
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, max_th_gen, max_mom_gen, max_th_gen);
	SetTitle(h_sa_mom_th, "#theta_{e} vs. Momentum (w/ all cuts)", "Momentum [GeV/c]", "#theta_{e} [rad]");
//	h_sa_mom_th->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_th->GetXaxis()->SetNdivisions(505, kFALSE);
	h_sa_mom_th->GetYaxis()->SetNdivisions(304, kFALSE);
//	h_sa_mom_th->Draw("samecolz");
	h_sa_mom_th->Draw("colz");

	//======= Theta vs. Momentum  ======
	c6->cd();
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, max_th_gen, max_mom_gen, max_th_gen);
	SetTitle(h_sa_costh_ph, "cos(#theta_{e}) vs. #phi_{e} (w/ all cuts)",  "cos(#theta_{e})", "#phi_{e} [rad]");
//	h_sa_costh_ph->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_costh_ph->GetXaxis()->SetNdivisions(505, kFALSE);
	h_sa_costh_ph->GetYaxis()->SetNdivisions(304, kFALSE);
//	h_sa_costh_ph->Draw("samecolz");
	h_sa_costh_ph->Draw("colz");

	//======= VP Flux  ======
	c7->Divide(2,2);
	c7->cd(1);
	//SetTitle(h_vp_mom_result, "Integrated VP Flux vs. Momentum (w/ all Cuts)", "Momentum [GeV/c]", "Integrated VP Flux [/GeV]");
	h_vp_mom_result->Sumw2();
	h_vp_mom_result2->Sumw2();
	h_vp_mom_result->Scale(bin_mom);
	h_vp_mom_result2->Scale(bin_mom);
	h_vp_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result->Draw();
	h_vp_mom_result2->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result2->SetLineColor(kAzure);
	h_vp_mom_result2->Draw("same");
	double ymax = (h_vp_mom_result->GetBinContent(h_vp_mom_result->GetMaximumBin()));
	TLine *RHRS_min = new TLine( 2.1, 0., 2.1, 1.1*ymax);
	TLine *RHRS_max = new TLine( 2.22, 0., 2.22, 1.1*ymax);
	RHRS_min->SetLineColor(kRed);
	RHRS_max->SetLineColor(kRed);
	RHRS_min->Draw("same");
	RHRS_max->Draw("same");
	c7->cd(2);
	//SetTitle(h_vp_th_result, "Integrated VP Flux vs. Theta (w/ all Cuts)", "Theta [rad]", "Integrated VP Flux [/GeV]");
	h_vp_th_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_th_result->Draw();
	c7->cd(3);
	h_vp_mom->Draw();
	c7->cd(4);
	h_vp_th->Draw();

	if(!BatchFlag){
		theApp.Run();
	}

	if(PDFFlag){
		cout << "Creating a PDF file ... " << endl;
		c1 ->Print(Form("%s[",PDFFile.str().c_str()) );
		c1 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c3 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c4 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c5 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c6 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c7 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c7 ->Print(Form("%s]" ,PDFFile.str().c_str()) );
	}
c1->Close();
c2->Close();
c3->Close();
c4->Close();
c5->Close();
c6->Close();
c7->Close();

	cout<< "Finish!" << endl;
	return 0;
}

void TreeBranch( TTree *tree, const char *name, double *branch){
	tree->SetBranchStatus(name);
	tree->SetBranchAddress( name, branch);
}

void SetTitle( TH2D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}

void SetTitle( TH1D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}

void CanvasSet( TCanvas *c, int pad, double x1, double y1, double x2, double y2){
	c->cd(pad);
	gPad->SetGrid();
	c->DrawFrame(x1, y1, x2, y2);
}
