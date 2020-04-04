//#include <iostream>
//#include <fstream>
//using namespace std;
//#include "TApplication.h"
//#include "hrs_tuningAC.h"


////////////////////////////////////////////////
////////Copy from hrs_tuningAC.h ///////////////
////////////////////////////////////////////////
//
//const double c=299792458e-9;// [m/ns]
//const double mk=493.7e-3;// Kaon mass [GeV/c^2]
//const double me=0.511e-3;// electron mass [GeV/c^2] 
//const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
//const double mn=939.6e-3; // neutron mass [GeV/c^2]
//const double mpi=139.6e-3;// pion mass [GeV/c^2]
//#include <iostream>
//#include <fstream>
//#include <math.h>
//#include <string>
//#include <time.h>
//#include <stdio.h>
//#include <unistd.h>
//#include <sstream>
//using namespace std;
//#include "TApplication.h"
//#include "TH1F.h"
//#include "TH2F.h"
//#include "TF1.h"
//#include "TFile.h"
//#include "TLeaf.h"
//#include "TTree.h"
//#include "TCut.h"
//#include "TChain.h"
//#include "TCanvas.h"
//#include "TVector3.h"
//#include "TGraph.h"
//#include "TLine.h"
//#include "TLatex.h"
//#include "TText.h"
//#include "TStyle.h"
//#include "TROOT.h"
//#include "TGraphErrors.h"
//#include "TProfile.h"
//#include "TSystem.h"
//#include "TColor.h"
//#include "TPaveText.h"
//#include "TRandom.h"
//#include "Setting.h


void ac_tuning_oku(){

//  set = new Setting();
//  set -> Initialize();

  gStyle->SetOptFit(111111111);
  int ch;
  string ifname = "nothing";
  string ofname = "./test.pdf";
  string root_name;
  string print_name;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag  = true;
  bool coin_flag  = false;
  bool print_flag = false;
  bool root_flag  = false;
  //  bool ac2_min=true;
  string itel;  
  string pngname;
  //  bool single_flag=false;

  cout << " input file name is " << ifname << endl;
  cout << " output file name is " << ofname << endl;
  
  //extern char *optarg;
  //while((ch=getopt(argc,argv,"h:f:w:s:n:r:i:o:bcop:GHT12"))!=-1){
  //  switch(ch){
  //  case 'f':
  //    ifname = optarg;
  //    cout<<"input filename : "<<ifname<<endl;
  //    break;
  //  case 's':
  //    ifname = optarg;
  //    cout<<"input filename : "<<ifname<<endl;
  //    single_flag=true;
  //    root_flag = true;
  //    draw_flag = false;            
  //    break;

  //    
  //  case 'w':
  //    print_flag = true;
  //     draw_flag = false;
  //    print_name = optarg;
  //    cout<<"output PDF filename : "<<print_name<<endl;
  //    break;


  //  case 'r':
  //    root_flag = true;
  //    draw_flag = false;      
  //    root_name = optarg;
  //    cout<<"output root filename : "<<root_name<<endl;      
  //    break;

  //  case 'o':
  //    root_flag=true;
  //    draw_flag=false;
  //    print_flag=true;
  //    ofname = optarg;
  //    root_name="./../rootfiles/ACtuning/" + ofname + ".root";
  //    print_name="./../pdf/ACtuning/" +ofname + ".pdf";
  //    break;
  //    
  //  case 'i':
  //    itel= optarg;
  //    nth= atoi(itel.c_str());
  //    break;

  //    
  //  case 'b':
  //    draw_flag = false;
  //    cout<<"BACH MODE!"<<endl;
  //    break;
  //
  //  case 'G':
  //  mode="G";
  //    break;
  //
  //  case 'H':
  //  mode="H";
  //    break;

  //  case 'T':
  //    mode="T";    
  //  break;

  //  case 'U':
  //    ac2_min=false;    
  //  break;

  //case '1':
  //  tdc_time=56.23e-3;//[ns]
  //  kine=1;
  //    break;

  //case '2':
  //  tdc_time=58e-3;//[ns]
  //  kine=2;
  //    break;


  //  case 'h':
  //    cout<<"-f : input root filename"<<endl;
  //    cout<<"-w : output pdf filename"<<endl;
  //    cout<<"-r : output root filename"<<endl;      
  //    cout<<"-o : output pdf & root  filename"<<endl;
  //    cout<<"-1 : F1TDC resolution 56 ns"<<endl;
  //    cout<<"-2 : F1TDC resolution 58 ns"<<endl;      
  //    cout<<"-T or H or G : Mode of root"<<endl;      
  //    return 0;
  //    break;
  //  case '?':
  //    cout<<"unknown option...."<<endl;
  //    return 0;
  //    break;
  //  default:
  //    cout<<"type -h to see help!!"<<endl;
  //    return 0;
  //  }
  //}

	double Ra1a[100], Ra2a[100], Ra1a_c[100], Ra2a_c[100];
	double Ra1sum, Ra2sum;

	double adc1[100], adc2[100];

   TChain *tr = new TChain("T");
   tr->SetBranchStatus("*",0);
   tr->Add("~/rootfile_nnL/tritium_111157_1.root");//or runlist
   tr->Add("~/rootfile_nnL/tritium_111157_2.root");
   tr->Add("~/rootfile_nnL/tritium_111157_3.root");
   tr->Add("~/rootfile_nnL/tritium_111157_4.root");
   tr->Add("~/rootfile_nnL/tritium_111157.root");
   
	tr->SetBranchStatus("R.a1.a",1);
	tr->SetBranchAddress("R.a1.a",Ra1a);
	tr->SetBranchStatus("R.a2.a",1);
	tr->SetBranchAddress("R.a2.a",Ra2a);
	tr->SetBranchStatus("R.a1.a_c",1);
	tr->SetBranchAddress("R.a1.a_c",Ra1a_c);
	tr->SetBranchStatus("R.a2.a_c",1);
	tr->SetBranchAddress("R.a2.a_c",Ra2a_c);
	tr->SetBranchStatus("R.a1.asum_c",1);
	tr->SetBranchAddress("R.a1.asum_c",&Ra1sum);
	tr->SetBranchStatus("R.a2.asum_c",1);
	tr->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 // tr->SetBranchStatus("R.a1.asum_p",1);
 // tr->SetBranchAddress("R.a1.asum_p",&Ra1sum);
 // tr->SetBranchStatus("R.a2.asum_p",1);
 // tr->SetBranchAddress("R.a2.asum_p",&Ra2sum);
	
	double ENum = tr->GetEntries();
	cout << "ENum = " << ENum << endl;

	TH1D *h_a1[24], *h_a2[26];
	for(int i=0;i<24;i++){
	h_a1[i] = new TH1D(Form("h_a1[%d]",i),Form("h_a1[%d]",i),500,5000.,6000.);
	}
	for(int i=0;i<26;i++){
	h_a2[i] = new TH1D(Form("h_a2[%d]",i),Form("h_a2[%d]",i),1000,5000.,7000.);
	}

	for(int n=0;n<ENum;n++){
		tr->GetEntry(n);
	for(int i=0;i<24;i++){
		double adc1 = Ra1a_c[i];
		h_a1[i] -> Fill( adc1 );
	}
	for(int i=0;i<26;i++){
		double adc2 = Ra2a_c[i];
		h_a2[i] -> Fill( adc2 );
	}
	}
	TCanvas *c1;//AC1 viewer (Linear)
	TCanvas *c2;//AC1 viewer (Log)
	c1 = new TCanvas("c1","A1 Linear",800.,800.);
	c2 = new TCanvas("c2","A1 Log",800.,800.);
	c1->Divide(6,4);
	c2->Divide(6,4);
	for(int i=0;i<24;i++){
	c1->cd(i+1);
	h_a1[i]->Draw();
	h_a1[i]->SetStats(0);
	c2->cd(i+1)->SetLogy();
	h_a1[i]->Draw();
	h_a1[i]->SetStats(0);
	}
 
	TCanvas *c3;//AC2 viewer (Linear)
	TCanvas *c4;//AC2 viewer (Log)
	c3 = new TCanvas("c3","A2 Linear",800.,800.);
	c4 = new TCanvas("c4","A2 Log",800.,800.);
	c3->Divide(7,4);
	c4->Divide(7,4);
	for(int i=0;i<26;i++){
	c3->cd(i+1);
	h_a2[i]->Draw();
	h_a2[i]->SetStats(0);
	c4->cd(i+1)->SetLogy();
	h_a2[i]->Draw();
	h_a2[i]->SetStats(0);
	}

cout << "Print is starting" << endl;
	c1->Print(Form("%s[",ofname.c_str()));
	c1->Print(Form("%s",ofname.c_str()));
	c2->Print(Form("%s",ofname.c_str()));
	c3->Print(Form("%s",ofname.c_str()));
	c4->Print(Form("%s",ofname.c_str()));
	c4->Print(Form("%s]",ofname.c_str()));

	


//  tuningAC* AC=new tuningAC();
//  if(single_flag)AC->SetRun(ifname);
//  else AC->SetRunList(ifname);
//  if(root_flag)AC->SetRoot(root_name);
//  AC->SetBranch();
//  AC->SetParam();
//  AC->MakeHist();
//  AC->Fill();
//	cout << "AC->Fill() is done" << endl;
//  AC->Fitting();
//	cout << "AC->Fitting() is done" << endl;
//  AC->Tuning();
//	cout << "AC->Tuning() is done" << endl;
//  AC->Draw();
//  if(print_flag)AC->Print(print_name);
//  if(root_flag)AC->Write();
//  AC->Comment();
//  
//  cout<<"=========== Output files ============="<<endl;
//  cout<<"output rootfile: "<<root_name<<endl;
  
  
	
//  TApplication *theApp =new TApplication("App",&argc,argv);
// if(draw_flag==0)gROOT->SetBatch(1);
// if(draw_flag==0)gSystem->Exit(1);
// theApp->Run();
 return 0;

}


