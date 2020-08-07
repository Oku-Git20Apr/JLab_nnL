
using namespace std;
#include "TAanaHRS.hh"
#include "Setting.h"

const double PI = atan(1.)*4.;

  //------- Aida, 20180119, uni11-1-1-1 ----------------
  const double  XFPm = 2.63, XpFPm = -0.000761; 
  const double  YFPm = 0.460, YpFPm = 0.000510; 
  const double  Xtm = 0., Xptm = -0.247;
  const double  Ytm = 0., Yptm = -0.002119;
  const double  Momm = 1.209;
  const double  XFPr = 59., XpFPr = 0.111; 
  const double  YFPr = 4.6, YpFPr = 0.0144; 
  const double  Xtr = 0.1, Xptr = 0.0369;
  const double  Ytr = 0.1, Yptr = 0.019;
  const double  Momr = 0.15;
  const double  ztm=0.0 , ztr=0.2715;

//____________________________________________________________________________________________
TAanaHRS::TAanaHRS()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  set = new Setting();

} //TAanaHRS()
//____________________________________________________________________________________________
TAanaHRS::~TAanaHRS()
{
std::cout<<"TAanaHRS::deconstructor called"<<std::endl;
}
//____________________________________________________________________________________________
void TAanaHRS::AnaTest()
{
  OpenTestRootFile();
  MakeHistOfTest();
  FillHistOfTest();
  FillSolidAngleOfTest();
  DrawMomOfTest();
  DrawXOfTest();
  DrawThetaOfTest();
} //AnaTest()
//____________________________________________________________________________________________
void TAanaHRS::AnaAcceptance()
{
  OpenTestRootFile();
  OpenTestLogFile();
  MakeHistOfTest();
  FillHistOfTest();
  FillSolidAngleOfTest();
  DrawAcceptanceOfTest();
} //AnaAcceptance()
//____________________________________________________________________________________________
void TAanaHRS::OpenTestRootFile()
{
  cout<<"OpenTestRootFile"<<endl;
  tree = new TChain("tree","tree");
  //tree -> Add("/data/7c/akiyama/JLab_simulation/sim/HRS/ana/root/test.root");
  tree -> Add("/data/7c/akiyama/JLab_simulation/sim/HRS/ana/root/0002.root");

} //OpenTestRootFile()
//____________________________________________________________________________________________
void TAanaHRS::OpenTestLogFile()
{
//   cout<<"Merge runnum: "<< runnum <<endl;
  cout<<"OpenTestLogFile"<<endl;
  //std::string filename = "/data/7c/akiyama/JLab_simulation/sim/HRS/ana/root/test";
  std::string filename = "/data/7c/akiyama/JLab_simulation/sim/HRS/ana/root/0002";
  std::ostringstream LogFile;
  LogFile << filename << ".root_Log";	
  
  FILE *fp;
  const int SizeOfBuffer = 64;
  char str[SizeOfBuffer];

  centraltheta = thetawidth = 0.;
  centralphi   = phiwidth   = 0.;
  centraltheta = thetawidth = 0.;
  centralphi   = phiwidth   = 0.;
  thetamin     = thetamax   = 0.;
  phimin       = phimax     = 0.;
  omega        =              0.; // [msr]
  
  if((fp=fopen(LogFile.str().c_str(), "r"))==NULL){
    std::cerr << "file open fail" << std::endl;
    exit(1);
  }

  while(fgets( str, SizeOfBuffer, fp)){
    if(sscanf( str, "Theta central [rad] = %lf", &centraltheta )==1){
    }
    else if(sscanf( str, "Theta Gen. Range [rad] = %lf", &thetawidth )==1){
    }
    else if(sscanf( str, "Phi Gen. Range [rad] = %lf", &phiwidth )==1){
    }
  }
  fclose(fp);
  cout<< "central theta = " << centraltheta*180/PI << " [deg]" << endl;
  cout<< "theta width= " << thetawidth*180/PI << " [deg]" << endl;
  centralphi = 1.*PI;
  
  thetamin = 1.*(1.*centraltheta - thetawidth);
  thetamax = 1.*(1.*centraltheta + thetawidth);
  phimin = 1.*(centralphi - phiwidth);
  phimax = 1.*(centralphi + phiwidth);
  omega = 1.0*(phimax - phimin) * (cos(thetamin) - cos(thetamax))*1000.; // [msr] <-- wrong calculation
//  omega = 1.0*(phimax - phimin) * (1.0-cos(thetawidth))*1000.; // [msr]
  cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
  cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
  cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
  cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
  cout<< "omega = " << omega << " [msr]" <<endl;

} //OpenTestLogFile()
//____________________________________________________________________________________________
void TAanaHRS::MakeHistOfTest()
{
  cout<<"makehistoftest"<<endl;
  bin_mom = 200; min_mom =    2.5; max_mom =   3.5;
  bin_th  = 200; min_th  =    0.0; max_th  =   0.4;
  bin_x   = 200; min_x   =  -70.0; max_x   =  70.0;

  bin_2D_th = 200; bin_2D_mom = 200;

  bin_3D = 2000;

  h_test_mom_gen    = new TH1D("h_test_mom_gen"    ,"h_test_mom_gen"    ,bin_mom ,min_mom ,max_mom);
  h_test_mom_pcs    = new TH1D("h_test_mom_pcs"    ,"h_test_mom_pcs"    ,bin_mom ,min_mom ,max_mom);
  h_test_mom_q1     = new TH1D("h_test_mom_q1"     ,"h_test_mom_q1"     ,bin_mom ,min_mom ,max_mom);
  h_test_mom_q2     = new TH1D("h_test_mom_q2"     ,"h_test_mom_q2"     ,bin_mom ,min_mom ,max_mom);
  h_test_mom_q3     = new TH1D("h_test_mom_q3"     ,"h_test_mom_q3"     ,bin_mom ,min_mom ,max_mom);
  h_test_mom_di     = new TH1D("h_test_mom_di"     ,"h_test_mom_di"     ,bin_mom ,min_mom ,max_mom);
  h_test_mom_rp     = new TH1D("h_test_mom_rp"     ,"h_test_mom_rp"     ,bin_mom ,min_mom ,max_mom);
//
  h_test_th_gen     = new TH1D("h_test_th_gen"     ,"h_test_th_gen"     ,bin_th ,min_th ,max_th);
  h_test_th_pcs     = new TH1D("h_test_th_pcs"     ,"h_test_th_pcs"     ,bin_th ,min_th ,max_th);
  h_test_th_q1      = new TH1D("h_test_th_q1"      ,"h_test_th_q1"      ,bin_th ,min_th ,max_th);
  h_test_th_q2      = new TH1D("h_test_th_q2"      ,"h_test_th_q2"      ,bin_th ,min_th ,max_th);
  h_test_th_q3      = new TH1D("h_test_th_q3"      ,"h_test_th_q3"      ,bin_th ,min_th ,max_th);
  h_test_th_di      = new TH1D("h_test_th_di"      ,"h_test_th_di"      ,bin_th ,min_th ,max_th);
  h_test_th_rp      = new TH1D("h_test_th_rp"      ,"h_test_th_rp"      ,bin_th ,min_th ,max_th);
//
  h_test_x_pcs      = new TH1D("h_test_x_pcs"     ,"h_test_x_pcs"       ,bin_x   ,min_x   ,max_x  );
  h_test_x_q1       = new TH1D("h_test_x_q1"      ,"h_test_x_q1"        ,bin_x   ,min_x   ,max_x  );
  h_test_x_q2       = new TH1D("h_test_x_q2"      ,"h_test_x_q2"        ,bin_x   ,min_x   ,max_x  );
  h_test_x_q3       = new TH1D("h_test_x_q3"      ,"h_test_x_q3"        ,bin_x   ,min_x   ,max_x  );
  h_test_x_di       = new TH1D("h_test_x_di"      ,"h_test_x_di"        ,bin_x   ,min_x   ,max_x  );
  h_test_x_rp       = new TH1D("h_test_x_rp"      ,"h_test_x_rp"        ,bin_x   ,min_x   ,max_x  );
//
  h_test_x_pcs_hit  = new TH1D("h_test_x_pcs_hit"  ,"h_test_x_pcs_hit"  ,bin_x   ,min_x   ,max_x  );
  h_test_x_q1_hit   = new TH1D("h_test_x_q1_hit"   ,"h_test_x_q1_hit"   ,bin_x   ,min_x   ,max_x  );
  h_test_x_q2_hit   = new TH1D("h_test_x_q2_hit"   ,"h_test_x_q2_hit"   ,bin_x   ,min_x   ,max_x  );
  h_test_x_q3_hit   = new TH1D("h_test_x_q3_hit"   ,"h_test_x_q3_hit"   ,bin_x   ,min_x   ,max_x  );
  h_test_x_di_hit   = new TH1D("h_test_x_di_hit"   ,"h_test_x_di_hit"   ,bin_x   ,min_x   ,max_x  );
  h_test_x_rp_hit   = new TH1D("h_test_x_rp_hit"   ,"h_test_x_rp_hit"   ,bin_x   ,min_x   ,max_x  );
//
  h_test_sa_mom_rp  = new TH1D("h_test_sa_mom_rp"  ,"h_test_sa_mom_rp"  ,bin_mom ,min_mom ,max_mom);
//
  h_test_sa_th_pcs  = new TH1D("h_test_sa_th_pcs"  ,"h_test_sa_th_pcs"  ,bin_th  ,min_th  ,max_th);
  h_test_sa_th_rp   = new TH1D("h_test_sa_th_rp"   ,"h_test_sa_th_rp"   ,bin_th  ,min_th  ,max_th);
//
//
//  //TH2
  h2_test_mom_x_rp     = new TH2D("h2_test_mom_x_rp"    ,"h2_test_mom_x_rp"    ,bin_mom   ,min_mom,max_mom,bin_x    ,min_x ,max_x );
  h2_test_mom_th_gen   = new TH2D("h2_test_mom_th_gen"  ,"h2_test_mom_th_gen"  ,bin_2D_mom,min_mom,max_mom,bin_2D_th,min_th,max_th);
  h2_test_mom_th_rp    = new TH2D("h2_test_mom_th_rp"   ,"h2_test_mom_th_rp"   ,bin_2D_mom,min_mom,max_mom,bin_2D_th,min_th,max_th);
  h2_test_sa_mom_th_rp = new TH2D("h2_test_sa_mom_th_rp","h2_test_sa_mom_th_rp",bin_2D_mom,min_mom,max_mom,bin_2D_th,min_th,max_th);
//
//  //TH3
  h3_test_sa_gen   = new TH3D("h3_test_sa_gen"  ,"h3_test_sa_gen"  ,100,-1.,1.,100,-1.,1.,100,0.8,1.);
  h3_test_sa_rp    = new TH3D("h3_test_sa_rp"   ,"h3_test_sa_rp"   ,100,-1.,1.,100,-1.,1.,100,0.8,1.);

  set->SetTH1(h_test_mom_gen   ,"Momentum Distribution (Generated)"            ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_pcs   ,"Momentum Distribution (PCS Exit)"             ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_q1    ,"Momentum Distribution (Q1 Exit)"              ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_q2    ,"Momentum Distribution (Q2 Exit)"              ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_q3    ,"Momentum Distribution (Q3 Exit)"              ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_di    ,"Momentum Distribution (Dipole Exit)"          ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_mom_rp    ,"Momentum Distribution (Ref. Plane)"           ,"Momentum [GeV/c]" ,"Counts"           ,kPink +6,3005,kPink +6);

  set->SetTH1(h_test_th_gen    ,"Theta Distribution (Generated)"               ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_pcs    ,"Theta Distribution (PCS Exit)"                ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_q1     ,"Theta Distribution (Q1 Exit)"                 ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_q2     ,"Theta Distribution (Q2 Exit)"                 ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_q3     ,"Theta Distribution (Q3 Exit)"                 ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_di     ,"Theta Distribution (Dipole Exit)"             ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_th_rp     ,"Theta Distribution (Ref. Plane)"              ,"theta [rad]"      ,"Counts"           ,kPink +6,3005,kPink +6);

  set->SetTH1(h_test_x_pcs     ,"X Distribution (PCS Exit)"                    ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_x_q1      ,"X Distribution (Q1 Exit)"                     ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_x_q2      ,"X Distribution (Q2 Exit)"                     ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_x_q3      ,"X Distribution (Q3 Exit)"                     ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_x_di      ,"X Distribution (Dipole Exit)"                 ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);
  set->SetTH1(h_test_x_rp      ,"X Distribution (Ref. Plane)"                  ,"X [cm]"           ,"Counts"           ,kPink +6,3005,kPink +6);

  set->SetTH1(h_test_x_pcs_hit ,"X Distribution (PCS Exit)"                    ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);
  set->SetTH1(h_test_x_q1_hit  ,"X Distribution (Q1 Exit)"                     ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);
  set->SetTH1(h_test_x_q2_hit  ,"X Distribution (Q2 Exit)"                     ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);
  set->SetTH1(h_test_x_q3_hit  ,"X Distribution (Q3 Exit)"                     ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);
  set->SetTH1(h_test_x_di_hit  ,"X Distribution (Dipole Exit)"                 ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);
  set->SetTH1(h_test_x_rp_hit  ,"X Distribution (Ref. Plane)"                  ,"X [cm]"           ,"Counts"           ,kAzure+6,3005,kAzure+6);

  set->SetTH1(h_test_sa_mom_rp ,"Solid Angle:Momentum Dependence (Ref. Plane)" ,"Momentum [GeV/c]" ,"Solid Angle [msr]",kPink +6,3005,kPink +6);

  set->SetTH1(h_test_sa_th_pcs ,"Solid Angle:Theta Dependence (PCS Exit)"      ,"Theta [rad]"      ,"Solid Angle [msr]",kPink +6,3005,kPink +6);
  set->SetTH1(h_test_sa_th_rp  ,"Solid Angle:Theta Dependence (Ref. Plane)"    ,"Theta [rad]"      ,"Solid Angle [msr]",kPink +6,3005,kPink +6);


  //TH2
  set->SetTH2(h2_test_mom_x_rp     ,"X v.s Momentum (Ref. Plane)"    ,"Momentum [GeV/c]","X [cm]"     );
  set->SetTH2(h2_test_mom_th_gen   ,"Theta v.s Momentum (generated)" ,"Momentum [GeV/c]","Theta [rad]");
  set->SetTH2(h2_test_mom_th_rp    ,"Theta v.s Momentum (Ref. Plane)","Momentum [GeV/c]","Theta [rad]");
  set->SetTH2(h2_test_sa_mom_th_rp ,"Solid Angle (Ref. Plane)"       ,"Momentum [GeV/c]","Theta [rad]");

  //TH3
  set->SetTH3(h3_test_sa_gen,"Momentum Direction","X","Y","Z",0.8,1,1.0,kPink +6);
  set->SetTH3(h3_test_sa_rp ,"Momentum Direction","X","Y","Z",0.8,1,1.0,kAzure+6);
} //MakeHistOfTest()
//____________________________________________________________________________________________
void TAanaHRS::FillHistOfTest()
{
  cout<<"FillHistOfTest"<<endl;
  double mom = 0.;
  double th  = 0.; double ph  = 0.;
  double x_sepe = 0.; double x_q1e = 0.; double x_q2e   = 0.; double x_q3e   = 0.;
  double x_die  = 0.; double x_rp  = 0.;// double x_tof2x = 0.;
  double x_dc1  = 0.; double x_dc2 = 0.;

  // use only for trig
  double x_sepi = 0.; double x_q1i = 0.; double x_q2i   = 0.; double x_q3i   = 0.;
  double x_dii  = 0.;
  double y_sepe = 0.; double y_q1e = 0.; double y_q2e   = 0.; double y_q3e   = 0.;
  double y_die  = 0.; double y_rp  = 0.;// double y_tof2x = 0.;
  double y_dc1  = 0.; double y_dc2 = 0.;
  double y_sepi = 0.; double y_q1i = 0.; double y_q2i   = 0.; double y_q3i   = 0.;
  double y_dii  = 0.;

  int evdtrig = 0;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("EMom"    ,1);tree->SetBranchAddress("EMom"    ,&mom    );
  tree->SetBranchStatus("ETheta"  ,1);tree->SetBranchAddress("ETheta"  ,&th     );
  tree->SetBranchStatus("EPhi"    ,1);tree->SetBranchAddress("EPhi"    ,&ph     );
  tree->SetBranchStatus("EXSepe"  ,1);tree->SetBranchAddress("EXSepe"  ,&x_sepe );
  tree->SetBranchStatus("EXQ1e"   ,1);tree->SetBranchAddress("EXQ1e"   ,&x_q1e  );
  tree->SetBranchStatus("EXQ2e"   ,1);tree->SetBranchAddress("EXQ2e"   ,&x_q2e  );
  tree->SetBranchStatus("EXQ3e"   ,1);tree->SetBranchAddress("EXQ3e"   ,&x_q3e  );
  tree->SetBranchStatus("EXDe"    ,1);tree->SetBranchAddress("EXDe"    ,&x_die  );
  tree->SetBranchStatus("EXFP"    ,1);tree->SetBranchAddress("EXFP"    ,&x_rp   );
  tree->SetBranchStatus("EXFP1"   ,1);tree->SetBranchAddress("EXFP1"   ,&x_dc1  );
  tree->SetBranchStatus("EXFP2"   ,1);tree->SetBranchAddress("EXFP2"   ,&x_dc2  );
  tree->SetBranchStatus("EXSepi"  ,1);tree->SetBranchAddress("EXSepi"  ,&x_sepi );
  tree->SetBranchStatus("EXQ1i"   ,1);tree->SetBranchAddress("EXQ1i"   ,&x_q1i  );
  tree->SetBranchStatus("EXQ2i"   ,1);tree->SetBranchAddress("EXQ2i"   ,&x_q2i  );
  tree->SetBranchStatus("EXQ3i"   ,1);tree->SetBranchAddress("EXQ3i"   ,&x_q3i  );
  tree->SetBranchStatus("EXDi"    ,1);tree->SetBranchAddress("EXDi"    ,&x_dii  );
  tree->SetBranchStatus("EYSepe"  ,1);tree->SetBranchAddress("EYSepe"  ,&y_sepe );
  tree->SetBranchStatus("EYQ1e"   ,1);tree->SetBranchAddress("EYQ1e"   ,&y_q1e  );
  tree->SetBranchStatus("EYQ2e"   ,1);tree->SetBranchAddress("EYQ2e"   ,&y_q2e  );
  tree->SetBranchStatus("EYQ3e"   ,1);tree->SetBranchAddress("EYQ3e"   ,&y_q3e  );
  tree->SetBranchStatus("EYDe"    ,1);tree->SetBranchAddress("EYDe"    ,&y_die  );
  tree->SetBranchStatus("EYFP"    ,1);tree->SetBranchAddress("EYFP"    ,&y_rp   );
  tree->SetBranchStatus("EYFP1"   ,1);tree->SetBranchAddress("EYFP1"   ,&y_dc1  );
  tree->SetBranchStatus("EYFP2"   ,1);tree->SetBranchAddress("EYFP2"   ,&y_dc2  );
  tree->SetBranchStatus("EYSepi"  ,1);tree->SetBranchAddress("EYSepi"  ,&y_sepi );
  tree->SetBranchStatus("EYQ1i"   ,1);tree->SetBranchAddress("EYQ1i"   ,&y_q1i  );
  tree->SetBranchStatus("EYQ2i"   ,1);tree->SetBranchAddress("EYQ2i"   ,&y_q2i  );
  tree->SetBranchStatus("EYQ3i"   ,1);tree->SetBranchAddress("EYQ3i"   ,&y_q3i  );
  tree->SetBranchStatus("EYDi"    ,1);tree->SetBranchAddress("EYDi"    ,&y_dii  );
  tree->SetBranchStatus("EVDTrig" ,1);tree->SetBranchAddress("EVDTrig" ,&evdtrig);
  ENum = tree->GetEntries();

  for(int i=0;i<ENum;i++){
    //cout<<i<<endl;
    tree->GetEntry(i);
    if(i%100000==0){
      cout<<i<<" / "<<ENum<<endl;
      //tree->Show(i);
    }

    //////////////////////
    ///// X position /////
    //////////////////////

    // PCSM
    h_test_x_pcs->Fill(x_sepe);

    // Q1
    h_test_x_q1->Fill(x_q1e);

    // Q2
    h_test_x_q2->Fill(x_q2e);

    // Q3
    h_test_x_q3->Fill(x_q3e);


    // HRS Dipole
    h_test_x_di->Fill(x_die);

    // Ref. Plane
    h_test_x_rp->Fill(x_rp);

    // DC1

    // DC2

    // PCSM in

    // Q1 in

    // Q2 in

    // HRS Dipole in

    //////////////////////
    ///// Y position /////
    //////////////////////

    // PCSM

    // Q1

    // Q2

    // HRS Dipole

    // Ref. Plane

    // DC1

    // DC2

    // PCSM in

    // Q1 in

    // Q2 in

    // HRS Dipole in


    //////////////////
    ///// Target /////
    //////////////////
    h_test_mom_gen    ->Fill(mom);
    h_test_th_gen     ->Fill(th );
    h2_test_mom_th_gen->Fill(mom,th);
    double x_3d = sin(th)*cos(ph);
    double y_3d = sin(th)*sin(ph);
    double z_3d = cos(th);
    //cout << x_3d << " " << y_3d << " " << z_3d << endl;
    //h3_test_sa_gen->Fill(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
    h3_test_sa_gen->Fill(x_3d,y_3d,z_3d);

    //cout << "OK!" << endl;
    /////////////////////
    ///// PCSM Exit /////
    /////////////////////
    if(0.<x_sepe && x_sepe<50.){ //if Hit at the PCSM Exit
      h_test_mom_pcs   ->Fill(mom);
      h_test_th_pcs    ->Fill(th );
      h_test_x_pcs_hit ->Fill(x_sepe);
    } //if a Hit at PCSM Exit

    ///////////////////
    ///// Q1 Exit /////
    ///////////////////
    if(sqrt((x_q1e*x_q1e)+(y_q1e*y_q1e))<24.765/2.){ //if Hit at the Q1 Exit
      h_test_mom_q1   ->Fill(mom);
      h_test_th_q1    ->Fill(th );
      h_test_x_q1_hit ->Fill(x_q1e);
    } //if a Hit at Q1 Exit

    ///////////////////
    ///// Q2 Exit /////
    ///////////////////
    if(sqrt((x_q2e*x_q2e)+(y_q2e*y_q2e))<40.0  /2.){ //if Hit at the Q2 Exit
      h_test_mom_q2   ->Fill(mom);
      h_test_th_q2    ->Fill(th );
      h_test_x_q2_hit ->Fill(x_q2e);
    } //if a Hit at Q2 Exit

    ///////////////////
    ///// Q3 Exit /////
    ///////////////////
    if(sqrt((x_q3e*x_q3e)+(y_q3e*y_q3e))<40.0  /2.){ //if Hit at the Q3 Exit
      h_test_mom_q3   ->Fill(mom);
      h_test_th_q3    ->Fill(th );
      h_test_x_q3_hit ->Fill(x_q3e);
    } //if a Hit at Q3 Exit

    ///////////////////////
    ///// Dipole Exit /////
    ///////////////////////
    if(-159.5/2.<x_die && x_die<159.5/2.){ //if Hit at the Dipole Exit
      h_test_mom_di   ->Fill(mom);
      h_test_th_di    ->Fill(th );
      h_test_x_di_hit ->Fill(x_die);
    } //if a Hit at Dipole Exit


    //////////////////////
    ///// Ref. Plane /////
    //////////////////////
    bool f_x_sepi  = false; bool f_y_sepi  = false; bool f_x_sepe = false; bool f_y_sepe = false;
    bool f_x_q1i   = false; bool f_y_q1i   = false; bool f_x_q1e  = false; bool f_y_q1e  = false;
    bool f_x_q2i   = false; bool f_y_q2i   = false; bool f_x_q2e  = false; bool f_y_q2e  = false;
    bool f_x_q3i   = false; bool f_y_q3i   = false; bool f_x_q3e  = false; bool f_y_q3e  = false;
    bool f_x_dii   = false; bool f_y_dii   = false; bool f_x_die  = false; bool f_y_die  = false;
    bool f_x_dc1   = false; bool f_y_dc1   = false; bool f_x_dc2  = false; bool f_y_dc2  = false;
    bool f_x_rp    = false; bool f_y_rp    = false;
    if(evdtrig==1){ //Trig Condition
      // Cut Condition
      if(   0.<x_sepi && x_sepi<50.){f_x_sepi = true;}
      //if( -30.<y_sepi && y_sepi<30.){f_y_sepi = true;}
      if(-999.<y_sepi              ){f_y_sepi = true;}
      if(   0.<x_sepe && x_sepe<50.){f_x_sepe = true;}
      //if( -30.<y_sepe && y_sepe<30.){f_y_sepe = true;}
      if(-999.<y_sepe              ){f_y_sepe = true;}
      if(sqrt((x_q1i*x_q1i)+(y_q1i*y_q1i))<24.765/2.){f_x_q1i=true; f_y_q1i=true;}
      if(sqrt((x_q1e*x_q1e)+(y_q1e*y_q1e))<24.765/2.){f_x_q1e=true; f_y_q1e=true;}
      if(sqrt((x_q2i*x_q2i)+(y_q2i*y_q2i))<40.0  /2.){f_x_q2i=true; f_y_q2i=true;}
      if(sqrt((x_q2e*x_q2e)+(y_q2e*y_q2e))<40.0  /2.){f_x_q2e=true; f_y_q2e=true;}
      if(sqrt((x_q3i*x_q3i)+(y_q3i*y_q3i))<40.0  /2.){f_x_q3i=true; f_y_q3i=true;}
      if(sqrt((x_q3e*x_q3e)+(y_q3e*y_q3e))<40.0  /2.){f_x_q3e=true; f_y_q3e=true;}
      if(-159.5/2.<x_dii && x_dii<159.5/2.){f_x_dii   = true;}
      //if(- 20./2. <y_dii && y_dii<  20./2.){f_y_dii   = true;}
      if(-999.    <y_dii                  ){f_y_dii   = true;}
      if(-159.5/2.<x_die && x_die<159.5/2.){f_x_die   = true;}
      if(-999.    <y_die                  ){f_y_die   = true;}
      //if(- 28.8/2.<x_dc1 && x_dc1< 28.8/2.){f_x_dc1   = true;}
      if(-999.    <x_dc1                  ){f_x_dc1   = true;}
      //if(-211.8/2.<y_dc1 && y_dc1<211.8/2.){f_y_dc1   = true;}
      if(-999.    <y_dc1                  ){f_y_dc1   = true;}
      //if(- 28.8/2.<x_dc2 && x_dc2< 28.8/2.){f_x_dc2   = true;}
      if(-999.    <x_dc2                  ){f_x_dc2   = true;}
      //if(-211.8/2.<y_dc2 && y_dc2<211.8/2.){f_y_dc2   = true;}
      if(-999.    <y_dc2                  ){f_y_dc2   = true;}
      if(-211.8/2.<x_rp  && x_rp <211.8/2.){f_x_rp    = true;}
      if(- 28.8/2.<y_rp  && y_rp < 28.8/2.){f_y_rp    = true;}

      // event fill
      if( //f_x_sepi  && f_y_sepi  && f_x_sepe && f_y_sepe &&
          //f_x_q1i   && f_y_q1i   &&// f_x_q1e  && f_y_q1e  &&
          //f_x_q2i   && f_y_q2i   &&// f_x_q2e  && f_y_q2e  &&
          //f_x_q3i   && f_y_q3i   &&// f_x_q3e  && f_y_q3e  &&
          //f_x_dii   && f_y_dii  // && f_x_die  && f_y_die//  &&
          f_x_dc1   && f_y_dc1   && f_x_dc2  && f_y_dc2//  &&
          //f_x_rp    && f_y_rp
      ){ //if triggerd & Hit at the ref. plrane

        h_test_mom_rp    ->Fill(mom);
        h_test_th_rp     ->Fill(th );
        h_test_x_rp_hit  ->Fill(x_rp);
        h2_test_mom_x_rp ->Fill(mom,x_rp);
        h2_test_mom_th_rp->Fill(mom,th);
        h3_test_sa_rp->Fill(x_3d,y_3d,z_3d);
        //h3_test_sa_rp->Fill(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
      } //if a Hit at RP
    } //if Trig condition

  } //for i (event fill)

} //FillHistOfTest()
//____________________________________________________________________________________________
void TAanaHRS::FillSolidAngleOfTest()
{
  cout<<"FillSolidAngleOfTest"<<endl;
  int n1, n2;
  double val;
  double err;

  //
  // Solid Angle v.s Momentum
  //
  for(int i=0; i<bin_mom; i++){
    n1 = 0;
    n1 = h_test_mom_gen->GetBinContent(i+1);
    
    // Ref. Plane
    n2 = 0;
    val = 0.;
    err = 0.;
    n2 = h_test_mom_rp->GetBinContent(i+1);
    if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
//      if(val!=0)cout<< n1 << " " << n2 << " " << val << endl;
    h_test_sa_mom_rp->SetBinContent(i+1, val);
    if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
    h_test_sa_mom_rp->SetBinError(i+1, err);
    
  } //for bin(mom)

  //
  // Solid Angle v.s Momentum
  //
  for(int i=0; i<bin_th; i++){
    n1 = 0;
    n1 = h_test_th_gen->GetBinContent(i+1);
    
    //PCSM Exit
    n2 = 0;
    val = 0.;
    err = 0.;
    n2 = h_test_th_pcs->GetBinContent(i+1);
    if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
    h_test_sa_th_pcs->SetBinContent(i+1, val);
    if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
    h_test_sa_th_pcs->SetBinError(i+1, err);
    
    //Ref. Plane
    n2 = 0;
    val = 0.;
    err = 0.;
    n2 = h_test_th_rp->GetBinContent(i+1);
    if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
    h_test_sa_th_rp->SetBinContent(i+1, val);
    if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
    h_test_sa_th_rp->SetBinError(i+1, err);

  } //for bin(th)

  //
  // Solid Angle v.s Momentum v.s Theta
  //
  for(int i=0; i<bin_2D_mom; i++){
    for(int j=0; j<bin_2D_th; j++){
      n1 = 0;
      n1 = h2_test_mom_th_gen->GetBinContent(i+1, j+1);
      n2 = 0;
      val = 0.;
      n2 = h2_test_mom_th_rp->GetBinContent(i+1, j+1);
      if(n1!=0 && n2!=0)val = omega*(1.0*n2/n1);
		//if(val!=0)cout<< n1 << " " << n2 << " " << val << endl;
      h2_test_sa_mom_th_rp->SetBinContent(i+1, j+1, val);
    } //for j
  } //for i

} //FillSolidAngleOfTest()
//____________________________________________________________________________________________
void TAanaHRS::DrawMomOfTest()
{
  cout<<"DrawMomOfTest"<<endl;
  c00[0] = new TCanvas("c00_0","c00_0",800,600);
  c00[0]->Clear();c00[0]->Divide(1,1,1E-4,1E-4);
    c00[0]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_gen->Draw();

  c00[1] = new TCanvas("c00_1","c00_1",800,600);
  c00[1]->Clear();c00[1]->Divide(1,1,1E-4,1E-4);
    c00[1]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_pcs->Draw();

  c00[2] = new TCanvas("c00_2","c00_2",800,600);
  c00[2]->Clear();c00[2]->Divide(1,1,1E-4,1E-4);
    c00[2]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_q1->Draw();

  c00[3] = new TCanvas("c00_3","c00_3",800,600);
  c00[3]->Clear();c00[3]->Divide(1,1,1E-4,1E-4);
    c00[3]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_q2->Draw();

  c00[4] = new TCanvas("c00_4","c00_4",800,600);
  c00[4]->Clear();c00[4]->Divide(1,1,1E-4,1E-4);
    c00[4]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_q3->Draw();

  c00[5] = new TCanvas("c00_5","c00_5",800,600);
  c00[5]->Clear();c00[5]->Divide(1,1,1E-4,1E-4);
    c00[5]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_di->Draw();

  c00[6] = new TCanvas("c00_6","c00_6",800,600);
  c00[6]->Clear();c00[6]->Divide(1,1,1E-4,1E-4);
    c00[6]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_rp->Draw();

} // DrawMomOfTest()
//____________________________________________________________________________________________
void TAanaHRS::DrawXOfTest()
{
  cout<<"DrawXOfTest"<<endl;
  c01[0] = new TCanvas("c01_0","c01_0",800,600);
  c01[0]->Clear();c01[0]->Divide(1,1,1E-4,1E-4);
    c01[0]->cd(1);/*gPad->SetLogy(1);*/h_test_x_pcs->Draw();h_test_x_pcs_hit->Draw("same");

  c01[1] = new TCanvas("c01_1","c01_1",800,600);
  c01[1]->Clear();c01[1]->Divide(1,1,1E-4,1E-4);
    c01[1]->cd(1);/*gPad->SetLogy(1);*/h_test_x_q1->Draw();h_test_x_q1_hit->Draw("same");

  c01[2] = new TCanvas("c01_2","c01_2",800,600);
  c01[2]->Clear();c01[2]->Divide(1,1,1E-4,1E-4);
    c01[2]->cd(1);/*gPad->SetLogy(1);*/h_test_x_q2->Draw();h_test_x_q2_hit->Draw("same");

  c01[3] = new TCanvas("c01_3","c01_3",800,600);
  c01[3]->Clear();c01[3]->Divide(1,1,1E-4,1E-4);
    c01[3]->cd(1);/*gPad->SetLogy(1);*/h_test_x_q3->Draw();h_test_x_q3_hit->Draw("same");

  c01[4] = new TCanvas("c01_4","c01_4",800,600);
  c01[4]->Clear();c01[4]->Divide(1,1,1E-4,1E-4);
    c01[4]->cd(1);/*gPad->SetLogy(1);*/h_test_x_di->Draw();h_test_x_di_hit->Draw("same");

  c01[5] = new TCanvas("c01_5","c01_5",800,600);
  c01[5]->Clear();c01[5]->Divide(1,1,1E-4,1E-4);
    c01[5]->cd(1);/*gPad->SetLogy(1);*/h_test_x_rp->Draw();h_test_x_rp_hit->Draw("same");

} // DrawXOfTest()
//____________________________________________________________________________________________
void TAanaHRS::DrawThetaOfTest()
{
  cout<<"DrawThetaOfTest"<<endl;
  c02[0] = new TCanvas("c02_0","c02_0",800,600);
  c02[0]->Clear();c02[0]->Divide(1,1,1E-4,1E-4);
    c02[0]->cd(1);/*gPad->SetLogy(1);*/h_test_th_gen->Draw();

  c02[1] = new TCanvas("c02_1","c02_1",800,600);
  c02[1]->Clear();c02[1]->Divide(1,1,1E-4,1E-4);
    c02[1]->cd(1);/*gPad->SetLogy(1);*/h_test_th_pcs->Draw();

  c02[2] = new TCanvas("c02_2","c02_2",800,600);
  c02[2]->Clear();c02[2]->Divide(1,1,1E-4,1E-4);
    c02[2]->cd(1);/*gPad->SetLogy(1);*/h_test_th_q1->Draw();

  c02[3] = new TCanvas("c02_3","c02_3",800,600);
  c02[3]->Clear();c02[3]->Divide(1,1,1E-4,1E-4);
    c02[3]->cd(1);/*gPad->SetLogy(1);*/h_test_th_q2->Draw();

  c02[4] = new TCanvas("c02_4","c02_4",800,600);
  c02[4]->Clear();c02[4]->Divide(1,1,1E-4,1E-4);
    c02[4]->cd(1);/*gPad->SetLogy(1);*/h_test_th_q3->Draw();

  c02[5] = new TCanvas("c02_5","c02_5",800,600);
  c02[5]->Clear();c02[5]->Divide(1,1,1E-4,1E-4);
    c02[5]->cd(1);/*gPad->SetLogy(1);*/h_test_th_di->Draw();

  c02[6] = new TCanvas("c02_6","c02_6",800,600);
  c02[6]->Clear();c02[6]->Divide(1,1,1E-4,1E-4);
    c02[6]->cd(1);/*gPad->SetLogy(1);*/h_test_th_rp->Draw();

} // DrawThetaOfTest()
//____________________________________________________________________________________________
void TAanaHRS::DrawAcceptanceOfTest()
{
  cout<<"DrawAcceptanceOfTest"<<endl;
  c03[0] = new TCanvas("c03_0","c03_0",800,600);
  c03[0]->Clear();c03[0]->Divide(1,1,1E-4,1E-4);
    c03[0]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_gen->Draw();

  c03[1] = new TCanvas("c03_1","c03_1",800,600);
  c03[1]->Clear();c03[1]->Divide(1,1,1E-4,1E-4);
    c03[1]->cd(1);/*gPad->SetLogy(1);*/h_test_th_gen->Draw();

  c03[2] = new TCanvas("c03_2","c03_2",800,600);
  c03[2]->Clear();c03[2]->Divide(1,1,1E-4,1E-4);
    c03[2]->cd(1);/*gPad->SetLogy(1);*/h_test_x_rp_hit->Draw();

  c03[3] = new TCanvas("c03_3","c03_3",800,600);
  c03[3]->Clear();c03[3]->Divide(1,1,1E-4,1E-4);
    c03[3]->cd(1);/*gPad->SetLogy(1);*/h_test_mom_rp->Draw();

  c03[4] = new TCanvas("c03_4","c03_4",800,600);
  c03[4]->Clear();c03[4]->Divide(1,1,1E-4,1E-4);
    c03[4]->cd(1);/*gPad->SetLogy(1);*/h_test_th_rp->Draw();

  c03[5] = new TCanvas("c03_5","c03_5",800,600);
  c03[5]->Clear();c03[5]->Divide(1,1,1E-4,1E-4);
    c03[5]->cd(1);/*gPad->SetLogy(1);*/h2_test_mom_x_rp->Draw("colz");

  c03[6] = new TCanvas("c03_6","c03_6",800,600);
  c03[6]->Clear();c03[6]->Divide(1,1,1E-4,1E-4);
    c03[6]->cd(1);/*gPad->SetLogy(1);*/h2_test_mom_th_rp->Draw("colz");

  c03[7] = new TCanvas("c03_7","c03_7",800,600);
  c03[7]->Clear();c03[7]->Divide(1,1,1E-4,1E-4);
    c03[7]->cd(1);/*gPad->SetLogy(1);*/h_test_sa_mom_rp->Draw();

  c03[8] = new TCanvas("c03_8","c03_8",800,600);
  c03[8]->Clear();c03[8]->Divide(1,1,1E-4,1E-4);
    c03[8]->cd(1);/*gPad->SetLogy(1);*/h_test_sa_th_rp->Draw();

  c03[9] = new TCanvas("c03_9","c03_9",800,600);
  c03[9]->Clear();c03[9]->Divide(1,1,1E-4,1E-4);
    c03[9]->cd(1);/*gPad->SetLogy(1);*/h2_test_sa_mom_th_rp->Draw("colz");

  c03[10] = new TCanvas("c03_10","c03_10",800,600);
  c03[10]->Clear();c03[10]->Divide(1,1,1E-4,1E-4);
    c03[10]->cd(1);/*gPad->SetLogy(1);*/h3_test_sa_gen->Draw();h3_test_sa_rp->Draw("same");

} // DrawAcceptanceOfTest()
