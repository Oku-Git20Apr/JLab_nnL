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
using namespace std;
void CanvasSet( TCanvas*, int, double, double, double, double);
void SetTitle( TH2D*, const char*, const char*, const char*);
void SetTitle( TH1D*, const char*, const char*, const char*);
void TreeBranch( TTree*, const char*, double*);

static const double PI = 4.*atan(1.0);

int main(int argc, char** argv){
	//////////////////////
	// command argument //
	//////////////////////
	int option;
	string filename = "test";
	bool BatchFlag = false;
	bool PDFFlag = false;
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
	PDFFile << "./pdf/Acceptance/" << filename << ".pdf";	

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
//	cout<< "central theta = " << centraltheta*180/PI << " [deg]" << endl;
	cout<< "theta width= " << thetawidth*180/PI << " [deg]" << endl;
	centralphi = 0.;

	TFile *f = new TFile(RootFile.str().c_str());
	TTree *t = (TTree*)f->Get("tree");

	double thetamin = 1.*(1.*centraltheta - thetawidth);
	double thetamax = 1.*(1.*centraltheta + thetawidth);
//	double phimin = 1.*(centralphi - phiwidth);
//	double phimax = 1.*(centralphi + phiwidth);
//	double omega = 1.0*(phimax - phimin) * (cos(thetamin) - cos(thetamax))*1000.; // [msr]
//	double thetamin = 0.;
//	double thetamax = thetawidth;
	double phimin = 0.;
	double phimax = 2.*PI;
	//double omega = 2.*PI* (cos(thetamin) - cos(thetamax))*1000.; // [msr]
	double omega = 2.*PI* (1-cos(thetawidth))*1000.; // [msr]
	cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
	cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
	cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
	cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
	cout<< "omega = " << omega << " [msr]" <<endl;

	// ======== Cut condition =========
	TCut Sepi = "0.<EXSepi && EXSepi<50.";
	TCut Sepe = "0.<EXSepe && EXSepe<50.";
	TCut Q1i = "-30.<EXQ1i && EXQ1i<30.";
	TCut Q1e = "-30.<EXQ1e && EXQ1e<30.";
	TCut Q2i = "-40.<EXQ2i && EXQ2i<40.";
	TCut Q2e = "-40.<EXQ2e && EXQ2e<40.";
//	TCut Di = "-159.5/2.<EXDi && EXDi<159.5/2.";
	TCut Di = "-999.<EXDi";
//	TCut De = "-159.5/2.<EXDe && EXDe<159.5/2.";
	TCut De = "-999.<EXDe";
	TCut Q3i = "-40.<EXQ3i && EXQ3i<40.";
	TCut Q3e = "-40.<EXQ3e && EXQ3e<40.";
//	TCut RP = "EVDTrig && -132./2.<EXFP && EXFP<132./2."; // Reference Plane
	TCut RP = "EVDTrig && EXFP1>-999 && EXFP2>-999"; // Reference Plane
//	TCut ETOF2X = "VDTrig && -150.<EXFP5 && EXFP5<150.";

	// ========== Tree information ===========-
//default
	//int bin_mom = 100;
	//double min_mom = 2.8;
	//double max_mom = 3.3;
	//int bin_mom_ex = 150;
	//double min_mom_ex = 2.8;
	//double max_mom_ex = 3.3;
	//int bin_th = 100;
	//double min_th = 0.00;
	//double max_th = 0.20;
	//int bin_th_ex = 150;
	//double min_th_ex = 0.06;
	//double max_th_ex = 0.18;
	//int bin_2D_mom = 50;
	//int bin_2D_th = 50;

	int bin_mom = 100;
	double min_mom = 1.8;
	double max_mom = 2.8;
	int bin_mom_ex = 150;
	double min_mom_ex = 1.8;
	double max_mom_ex = 2.8;
	int bin_th = 100;
	double min_th = 0.10;
	double max_th = 0.30;
	int bin_th_ex = 150;
	double min_th_ex = 0.06;
	double max_th_ex = 0.18;
	int bin_2D_mom = 50;
	int bin_2D_th = 50;

	TH1D *h_mom_gen = new TH1D( "h_mom_gen", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_gen_ex = new TH1D( "h_mom_gen_ex", "", bin_mom_ex, min_mom_ex, max_mom_ex);
	TH1D *h_mom_sepi = new TH1D( "h_mom_sepi", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_sepe = new TH1D( "h_mom_sepe", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q1i = new TH1D( "h_mom_q1i", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q1e = new TH1D( "h_mom_q1e", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q2i = new TH1D( "h_mom_q2i", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q2e = new TH1D( "h_mom_q2e", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_di = new TH1D( "h_mom_di", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_de = new TH1D( "h_mom_de", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q3i = new TH1D( "h_mom_q3i", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_q3e = new TH1D( "h_mom_q3e", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_rp = new TH1D( "h_mom_rp", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_rp_ex = new TH1D( "h_mom_rp_ex", "", bin_mom_ex, min_mom_ex, max_mom_ex); // expansion
	TH1D *h_mom_tof = new TH1D( "h_mom_tof", "", bin_mom, min_mom, max_mom);
	t->Project( "h_mom_gen", "EMom");
	t->Project( "h_mom_gen_ex", "EMom");
	t->Project( "h_mom_sepi", "EMom", Sepi);
	t->Project( "h_mom_sepe", "EMom", Sepe);
	t->Project( "h_mom_q1i", "EMom", Q1i);
	t->Project( "h_mom_q1e", "EMom", Q1e);
	t->Project( "h_mom_q2i", "EMom", Q2i);
	t->Project( "h_mom_q2e", "EMom", Q2e);
	t->Project( "h_mom_di", "EMom", Di);
	t->Project( "h_mom_de", "EMom", De);
	t->Project( "h_mom_q3i", "EMom", Q3i);
	t->Project( "h_mom_q3e", "EMom", Q3e);
	t->Project( "h_mom_rp", "EMom", RP);
	t->Project( "h_mom_rp_ex", "EMom", RP);
//	t->Project( "h_mom_tof", "EMom", ETOF2X);
	TH1D *h_th_gen = new TH1D( "h_th_gen", "", bin_th, min_th, max_th);
	TH1D *h_th_gen_ex = new TH1D( "h_th_gen_ex", "", bin_th_ex, min_th_ex, max_th_ex);
	TH1D *h_th_sepi = new TH1D( "h_th_sepi", "", bin_th, min_th, max_th);
	TH1D *h_th_sepe = new TH1D( "h_th_sepe", "", bin_th, min_th, max_th);
	TH1D *h_th_q1i = new TH1D( "h_th_q1i", "", bin_th, min_th, max_th);
	TH1D *h_th_q1e = new TH1D( "h_th_q1e", "", bin_th, min_th, max_th);
	TH1D *h_th_q2i = new TH1D( "h_th_q2i", "", bin_th, min_th, max_th);
	TH1D *h_th_q2e = new TH1D( "h_th_q2e", "", bin_th, min_th, max_th);
	TH1D *h_th_di = new TH1D( "h_th_di", "", bin_th, min_th, max_th);
	TH1D *h_th_de = new TH1D( "h_th_de", "", bin_th, min_th, max_th);
	TH1D *h_th_q3i = new TH1D( "h_th_q3i", "", bin_th, min_th, max_th);
	TH1D *h_th_q3e = new TH1D( "h_th_q3e", "", bin_th, min_th, max_th);
	TH1D *h_th_rp = new TH1D( "h_th_rp", "", bin_th, min_th, max_th);
	TH1D *h_th_rp_ex = new TH1D( "h_th_rp_ex", "", bin_th_ex, min_th_ex, max_th_ex); // expansion
	TH1D *h_th_tof = new TH1D( "h_th_tof", "", bin_th, min_th, max_th);
	t->Project( "h_th_gen", "ETheta");
	t->Project( "h_th_gen_ex", "ETheta");
	t->Project( "h_th_sepi", "ETheta", Sepi);
	t->Project( "h_th_sepe", "ETheta", Sepe);
	t->Project( "h_th_q1i", "ETheta", Q1i);
	t->Project( "h_th_q1e", "ETheta", Q1e);
	t->Project( "h_th_q2i", "ETheta", Q2i);
	t->Project( "h_th_q2e", "ETheta", Q2e);
	t->Project( "h_th_di", "ETheta", Di);
	t->Project( "h_th_de", "ETheta", De);
	t->Project( "h_th_q3i", "ETheta", Q3i);
	t->Project( "h_th_q3e", "ETheta", Q3e);
	t->Project( "h_th_rp", "ETheta", RP);
	t->Project( "h_th_rp_ex", "ETheta", RP);
//	t->Project( "h_th_tof", "ETheta", ETOF2X);

	TH2D *h_mom_th_gen = new TH2D( "h_mom_th_gen", "", bin_2D_mom, min_mom_ex, max_mom_ex, bin_2D_th, min_th_ex, max_th_ex);
	TH2D *h_mom_th = new TH2D( "h_mom_th", "", bin_2D_mom, min_mom_ex, max_mom_ex, bin_2D_th, min_th_ex, max_th_ex);
	t->Project( "h_mom_th_gen", "ETheta:EMom");
	t->Project( "h_mom_th", "ETheta:EMom", RP);

	// ========= Solid angle vs. momentum ===========
	TH1D *h_sa_mom_sepi = new TH1D("h_sa_mom_sepi", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_sepe = new TH1D("h_sa_mom_sepe", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q1i = new TH1D("h_sa_mom_q1i", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q1e = new TH1D("h_sa_mom_q1e", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q2i = new TH1D("h_sa_mom_q2i", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q2e = new TH1D("h_sa_mom_q2e", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_di = new TH1D("h_sa_mom_di", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_de = new TH1D("h_sa_mom_de", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q3i = new TH1D("h_sa_mom_q3i", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_q3e = new TH1D("h_sa_mom_q3e", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_rp = new TH1D("h_sa_mom_rp", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_rp_ex = new TH1D("h_sa_mom_rp_ex", "", bin_mom_ex, min_mom_ex, max_mom_ex);
	TH1D *h_sa_mom_tof = new TH1D("h_sa_mom_tof", "", bin_mom, min_mom, max_mom);

	int n1, n2;
	double val;
	double err;
	const double min_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindFirstBinAbove()) - 0.05;
	const double max_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindLastBinAbove()) + 0.05;
	const double min_mom_bin = h_mom_gen->FindBin(min_mom_gen);
	const double max_mom_bin = h_mom_gen->FindBin(max_mom_gen);
	cout <<"mom range: "<< min_mom_gen << " " << max_mom_gen << endl;
	cout <<"mom_bin: "<< min_mom_bin << " " << max_mom_bin << endl;

	for(int i=0; i<bin_mom; i++){
		n1 = 0;
		n1 = (int)h_mom_gen->GetBinContent(i+1);

		// === Septum entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_sepi->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_sepi->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_sepi->SetBinError(i+1, err);

		// === Septum exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_sepe->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_sepe->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_sepe->SetBinError(i+1, err);

		// === Q1 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q1i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q1i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q1i->SetBinError(i+1, err);

		// === Q1 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q1e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q1e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q1e->SetBinError(i+1, err);

		// === Q2 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q2i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q2i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q2i->SetBinError(i+1, err);

		// === Q2 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q2e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q2e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q2e->SetBinError(i+1, err);

		// === Dipole entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_di->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_di->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_di->SetBinError(i+1, err);

		// === Dipole exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_de->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_de->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_de->SetBinError(i+1, err);

		// === Q3 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q3i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q3i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q3i->SetBinError(i+1, err);

		// === Q3 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_q3e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_q3e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_q3e->SetBinError(i+1, err);

		// === RP ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_rp->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_rp->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_rp->SetBinError(i+1, err);

		// === HTOF2X ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_tof->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_tof->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_tof->SetBinError(i+1, err);
	}

	for(int i=0; i<bin_mom_ex; i++){
		n1 = 0;
		n1 = (int)h_mom_gen_ex->GetBinContent(i+1);

		// === RP Expansion ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_rp_ex->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_mom_rp_ex->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_rp_ex->SetBinError(i+1, err);
	}

	// ========= Solid angle vs. theta ===========
	TH1D *h_sa_th_sepi = new TH1D("h_sa_th_sepi", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_sepe = new TH1D("h_sa_th_sepe", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q1i = new TH1D("h_sa_th_q1i", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q1e = new TH1D("h_sa_th_q1e", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q2i = new TH1D("h_sa_th_q2i", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q2e = new TH1D("h_sa_th_q2e", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_di = new TH1D("h_sa_th_di", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_de = new TH1D("h_sa_th_de", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q3i = new TH1D("h_sa_th_q3i", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_q3e = new TH1D("h_sa_th_q3e", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_rp = new TH1D("h_sa_th_rp", "", bin_th, min_th, max_th);
	TH1D *h_sa_th_rp_ex = new TH1D("h_sa_th_rp_ex", "", bin_th_ex, min_th_ex, max_th_ex);
	TH1D *h_sa_th_tof = new TH1D("h_sa_th_tof", "", bin_th, min_th, max_th);

	const double min_th_gen = h_th_gen->GetBinCenter(h_th_gen->FindFirstBinAbove()) - 0.02;
	const double max_th_gen = h_th_gen->GetBinCenter(h_th_gen->FindLastBinAbove()) + 0.02;
	const double min_th_bin = h_th_gen->FindBin(min_th_gen);
	const double max_th_bin = h_th_gen->FindBin(max_th_gen);
	cout <<"th range: "<< min_th_gen << " " << max_th_gen << endl;
	cout <<"th bin: "<< min_th_bin << " " << max_th_bin << endl;

	for(int i=0; i<bin_th; i++){
		n1 = 0;
		n1 = (int)h_th_gen->GetBinContent(i+1);

		// === Septum entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_sepi->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_sepi->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_sepi->SetBinError(i+1, err);

		// === Septum exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_sepe->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_sepe->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_sepe->SetBinError(i+1, err);

		// === Q1 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q1i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q1i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q1i->SetBinError(i+1, err);

		// === Q1 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q1e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q1e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q1e->SetBinError(i+1, err);

		// === Q2 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q2i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q2i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q2i->SetBinError(i+1, err);

		// === Q2 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q2e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q2e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q2e->SetBinError(i+1, err);

		// === Dipole entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_di->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_di->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_di->SetBinError(i+1, err);

		// === Dipole exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_de->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_de->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_de->SetBinError(i+1, err);

		// === Q3 entrance ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q3i->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q3i->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q3i->SetBinError(i+1, err);

		// === Q3 exit ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_q3e->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_q3e->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_q3e->SetBinError(i+1, err);

		// === RP ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_rp->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_rp->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_rp->SetBinError(i+1, err);

		// === HTOF2X ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_tof->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_tof->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_tof->SetBinError(i+1, err);
	}

	for(int i=0; i<bin_th_ex; i++){
		n1 = 0;
		n1 = (int)h_th_gen_ex->GetBinContent(i+1);

		// === RP Expansion ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_th_rp_ex->GetBinContent(i+1);
		if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		h_sa_th_rp_ex->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_th_rp_ex->SetBinError(i+1, err);
	}

	// === Momentum vs. Theta ===
	TH2D *h_sa_mom_th = new TH2D("h_sa_mom_th", "", bin_2D_mom, min_mom_ex, max_mom_ex, bin_2D_th, min_th_ex, max_th_ex);
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
	TCanvas *c6 = new TCanvas("c6", "Momentum Acceptance (Entrance)");
	TCanvas *c7 = new TCanvas("c7", "Angular Acceptance (Entrance)");

	c1->Divide(3,2);
	c2->Divide(3,2);
	c6->Divide(3,2);
	c7->Divide(3,2);

	//======= Momentum Acceptance (Exit) ======
//	gPad->SetGrid();
	c1->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Septum Exit", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_sepe->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q1 Exit", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q1e->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q2 Exit", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q2e->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Dipole Exit", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_de->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(5);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q3 Exit", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q3e->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	c1->cd(6);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
//	SetTitle(hframe, "Solid Angle vs. Momentum at HTOF-2X", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_tof->Draw("same");
	SetTitle(hframe, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_rp->Draw("same");
	lmom->Draw("same");
//	c1->RedrawAxis();

	//======= Angular Acceptance (Exit) ======
	c2->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Septum Exit", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_sepe->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q1 Exit", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q1e->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q2 Exit", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q2e->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Dipole Exit", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_de->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(5);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q3 Exit", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q3e->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	c2->cd(6);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
//	SetTitle(hframe, "Solid Angle vs. #theta_{e} at HTOF-2X", "#theta_{e} [rad]", "Solid Angle [msr]");
//	h_sa_th_tof->Draw("same");
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Reference Plane", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_rp->Draw("same");
	lth->Draw("same");
//	c2->RedrawAxis();

	//======= Momentum Acceptance Expansion ======
	c3->cd();
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(h_sa_mom_rp_ex, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
  h_sa_mom_rp_ex->SetTitleOffset(1.2, "X");
  h_sa_mom_rp_ex->SetTitleOffset(1.2, "Y");
  h_sa_mom_rp_ex->SetLabelOffset(0.02, "X");
  h_sa_mom_rp_ex->SetLabelOffset(0.02, "Y");
  h_sa_mom_rp_ex->SetLabelSize(0.04, "X");
  h_sa_mom_rp_ex->SetLabelSize(0.06, "Y");
//	h_sa_mom_rp_ex->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_rp_ex->GetXaxis()->SetNdivisions(505, kFALSE);
//	h_sa_mom_rp_ex->Draw("same");
	h_sa_mom_rp_ex->Draw();
//	lmom->Draw("same");

	//======= Angular Acceptance Expansion ======
	c4->cd();
//	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(h_sa_th_rp_ex, "Solid Angle vs. #theta_{e} at Reference Plane", "#theta_{e} [rad]", "Solid Angle [msr]");
  h_sa_th_rp_ex->SetTitleOffset(1.2, "X");
  h_sa_th_rp_ex->SetTitleOffset(1.2, "Y");
  h_sa_th_rp_ex->SetLabelOffset(0.02, "X");
  h_sa_th_rp_ex->SetLabelOffset(0.02, "Y");
  h_sa_th_rp_ex->SetLabelSize(0.04, "X");
  h_sa_th_rp_ex->SetLabelSize(0.06, "Y");
	h_sa_th_rp_ex->GetXaxis()->SetNdivisions(304, kFALSE);
//	h_sa_th_rp_ex->Draw("same");
	h_sa_th_rp_ex->Draw();

	//======= Theta vs. Momentum  ======
	c5->cd();
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, max_th_gen, max_mom_gen, max_th_gen);
	SetTitle(h_sa_mom_th, "#theta_{e} vs. Momentum at Reference Plane", "Momentum [GeV/c]", "#theta_{e} [rad]");
  h_sa_mom_th->SetTitleOffset(1.2, "X");
  h_sa_mom_th->SetTitleOffset(1.2, "Y");
  h_sa_mom_th->SetLabelOffset(0.01, "X");
  h_sa_mom_th->SetLabelOffset(0.01, "Y");
  h_sa_mom_th->SetLabelSize(0.06, "X");
  h_sa_mom_th->SetLabelSize(0.055, "Y");
//	h_sa_mom_th->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_th->GetXaxis()->SetNdivisions(505, kFALSE);
	h_sa_mom_th->GetYaxis()->SetNdivisions(304, kFALSE);
//	h_sa_mom_th->Draw("samecolz");
	h_sa_mom_th->Draw("colz");

	//======= Momentum Acceptance (Entrance) ======
	c6->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Septum Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_sepi->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	c6->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q1 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q1i->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	c6->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q2 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q2i->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	c6->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Dipole Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_di->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	c6->cd(5);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. Momentum at Q3 Entrance", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_q3i->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	c6->cd(6);
	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, omega + 0.5);
//	SetTitle(hframe, "Solid Angle vs. Momentum at HTOF-2X", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_tof->Draw("same");
	SetTitle(hframe, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	h_sa_mom_rp->Draw("same");
	lmom->Draw("same");
//	c6->RedrawAxis();

	//======= Angular Acceptance (Entrance) ======
	c7->cd(1);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Septum Entrance", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_sepi->Draw("same");
//	c7->RedrawAxis();

	c7->cd(2);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q1 Entrance", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q1i->Draw("same");
//	c7->RedrawAxis();

	c7->cd(3);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q2 Entrance", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q2i->Draw("same");
//	c7->RedrawAxis();

	c7->cd(4);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Dipole Entrance", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_di->Draw("same");
//	c7->RedrawAxis();

	c7->cd(5);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Q3 Entrance", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_q3i->Draw("same");
//	c7->RedrawAxis();

	c7->cd(6);
	hframe = (TH1D*)gPad->DrawFrame( min_th_gen, 0., max_th_gen, omega + 0.5);
//	SetTitle(hframe, "Solid Angle vs. #theta_{e} at HTOF-2X", "#theta_{e} [rad]", "Solid Angle [msr]");
//	h_sa_th_tof->Draw("same");
	SetTitle(hframe, "Solid Angle vs. #theta_{e} at Reference Plane", "#theta_{e} [rad]", "Solid Angle [msr]");
	h_sa_th_rp->Draw("same");
//	c7->RedrawAxis();

//	c1->Close();
//	c2->Close();
//	c3->Close();
//	c4->Close();
//	c5->Close();
//	c6->Close();
//	c7->Close();
	if(!BatchFlag){
		theApp.Run();
	}

	if(PDFFlag){
		cout << "Creating a PDF file ... " << endl;
		c6 ->Print(Form("%s[",PDFFile.str().c_str()) );
		c6 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c1 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c7 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c3 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c4 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c5 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c5 ->Print(Form("%s]" ,PDFFile.str().c_str()) );
	}

	cout<< "YEAH!!!" << endl;
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
