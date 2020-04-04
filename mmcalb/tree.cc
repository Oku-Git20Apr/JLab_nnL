#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
                    
#include "tree.h"
//=====================================//
//=========== Tree   ==================//
//=====================================//

tree::tree(){   T=new TChain("T");};
tree::~tree(){};

void tree::SetRun(string ifname){
  T->Add(ifname.c_str());
  ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;
  
}


void tree::ChainTree(string ifname){

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //    cout<<buf<<endl;
  }
  ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

};


void tree::SetBranch(){

  
 T->SetBranchStatus("*",0);  
 T->SetBranchStatus("DR.evtype",1);
 T->SetBranchAddress("DR.evtype",&DRevtype); 
 // T->SetBranchStatus("fEvtHdr.fRun",1);
 // T->SetBranchAddress("fEvtHdr.fRun",&runnum); 
 T->SetBranchStatus("HALLA_p",1);
 T->SetBranchAddress("HALLA_p",&hallap); 
 //------------------------------// 
 //------ Right Arm -------------//
 //------------------------------//
 // Scintillation //
 T->SetBranchStatus("RTDC.F1FirstHit",1);
 T->SetBranchAddress("RTDC.F1FirstHit",RF1);
 T->SetBranchStatus("R.s2.t_pads",1);
 T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
 T->SetBranchStatus("R.s2.trpad",1);
 T->SetBranchAddress("R.s2.trpad",Rs2trpad);
 T->SetBranchStatus("R.tr.n",1);
 T->SetBranchAddress("R.tr.n",&Rtrn);
 // PID //
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);
 T->SetBranchStatus("R.cer.asum_c",1);
 T->SetBranchAddress("R.cer.asum_c",&Rgssum);

 // path length//
 T->SetBranchStatus("R.s2.trpath",1); 
 T->SetBranchAddress("R.s2.trpath",rs2pathl); 
 T->SetBranchStatus("R.tr.pathl",1);  
 T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // target positon information //
 T->SetBranchAddress("Rrb.Raster2.rawcur.x", &R_Ras_x); // raster current
 T->SetBranchAddress("Lrb.Raster2.rawcur.x", &L_Ras_x); // raster current 
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.px",1);
 T->SetBranchAddress("R.tr.px",Rpx);
 T->SetBranchStatus("R.tr.py",1);
 T->SetBranchAddress("R.tr.py",Rpy);
 T->SetBranchStatus("R.tr.pz",1);
 T->SetBranchAddress("R.tr.pz",Rpz); 
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rz); 
 T->SetBranchStatus("R.tr.vx",1);    
 T->SetBranchAddress("R.tr.vx",Rx);
 T->SetBranchStatus("R.tr.vy",1);    
 T->SetBranchAddress("R.tr.vy",Ry);
 T->SetBranchStatus("R.tr.tg_th",1);    
 T->SetBranchAddress("R.tr.tg_th",Rth);
 T->SetBranchStatus("R.tr.tg_ph",1);    
 T->SetBranchAddress("R.tr.tg_ph",Rph);
 // Focal Plane infomation //
 T->SetBranchStatus("R.tr.x",1);    
 T->SetBranchAddress("R.tr.x",Rx_fp);
 T->SetBranchStatus("R.tr.y",1);    
 T->SetBranchAddress("R.tr.y",Ry_fp);
 T->SetBranchStatus("R.tr.th",1);    
 T->SetBranchAddress("R.tr.th",Rth_fp); 
 T->SetBranchStatus("R.tr.ph",1);    
 T->SetBranchAddress("R.tr.ph",Rph_fp); 

 //-------------------------------// 
 //------ Left Arm ---------------//
 //-------------------------------//
 // Scintillation//
 T->SetBranchStatus("LTDC.F1FirstHit",1);
 T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
 T->SetBranchStatus("L.s2.t_pads",1);
 T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
 T->SetBranchStatus("L.s2.trpad",1);
 T->SetBranchAddress("L.s2.trpad",Ls2trpad);
 T->SetBranchStatus("L.tr.n",1);
 T->SetBranchAddress("L.tr.n",&Ltrn);
 // PID //
 // T->SetBranchStatus("L.sh.asum_c",1);
 // T->SetBranchAddress("L.sh.asum_c",&Lshsum);
 // T->SetBranchStatus("L.ps.asum_c",1);
 // T->SetBranchAddress("L.ps.asum_c",&Lpssum);
 T->SetBranchStatus("L.cer.asum_c",1);
 T->SetBranchAddress("L.cer.asum_c",&Lcersum);  
 // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 // target positon information //
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);
 T->SetBranchStatus("L.tr.px",1);
 T->SetBranchAddress("L.tr.px",Lpx);
 T->SetBranchStatus("L.tr.py",1);
 T->SetBranchAddress("L.tr.py",Lpy);
 T->SetBranchStatus("L.tr.pz",1);
 T->SetBranchAddress("L.tr.pz",Lpz); 
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lz); 
 T->SetBranchStatus("L.tr.vx",1);    
 T->SetBranchAddress("L.tr.vx",Lx);
 T->SetBranchStatus("L.tr.vy",1);    
 T->SetBranchAddress("L.tr.vy",Ly);
 T->SetBranchStatus("L.tr.tg_th",1);    
 T->SetBranchAddress("L.tr.tg_th",Lth);
 T->SetBranchStatus("L.tr.tg_ph",1);    
 T->SetBranchAddress("L.tr.tg_ph",Lph);
 // Focal Plane infomation //
 T->SetBranchStatus("L.tr.x",1);    
 T->SetBranchAddress("L.tr.x",Lx_fp);
 T->SetBranchStatus("L.tr.y",1);    
 T->SetBranchAddress("L.tr.y",Ly_fp);
 T->SetBranchStatus("L.tr.th",1);    
 T->SetBranchAddress("L.tr.th",Lth_fp); 
 T->SetBranchStatus("L.tr.ph",1);    
 T->SetBranchAddress("L.tr.ph",Lph_fp); 

};


/*
void tree::NewBranch(string ofname, bool rarm){

  fnew = new TFile(ofname.c_str(),"recreate");


  if(rarm==true) tnew = new TTree("T","For momentum calibration (RHRS)");
    else tnew = new TTree("T","For momentum calibration (LHRS)");


  tnew->Branch("DR.evtype",&DRevtype,"DR.evtype/D");
  tnew->Branch("fEvtHdr.fRun",&runnum,"fEvtHdr.frun/D"); 
  tnew->Branch("HALLA_p",&hallap,"HALLA_p/D"); 
 
 //------------------------------// 
 //------ Right Arm -------------//
 //------------------------------//
 // Scintillation //

  tnew->Branch("RTDC.F1FirstHit",RF1,"RTDC.F1FirstHit[100]/D");
  tnew->Branch("R.s2.t_pads",Rs2tpads,"R.s2.t_pads[100]/D");
  tnew->Branch("R.s2.trpad",Rs2trpad,"R.s2.trpad[100]/D");
  tnew->Branch("R.tr.n",&Rtrn,"R.tr.n");
  // PID //

  tnew->Branch("R.a1.asum_c",&Ra1sum,"R.a1.asum_c/D");
  tnew->Branch("R.a2.asum_c",&Ra2sum,"R.a2.asum_c/D");
  tnew->Branch("R.cer.asum_c",&Rgssum,"R.cer.asun_c/D");

 // path length//
  tnew->Branch("R.s2.trpath",rs2pathl,"R.s2.trpath[100]/D"); 
  tnew->Branch("R.tr.pathl",rtrpathl,"R.tr.pathl[100]/D");
 // target positon information //
  tnew->Branch("R.tr.p",Rp,"R.tr.p[100]/D");    
  tnew->Branch("R.tr.vz",Rz,"R.tr.vz[100]/D");     
  tnew->Branch("R.tr.tg_x",Rx,"R.tr.tg_x[100]/D");  
  tnew->Branch("R.tr.tg_y",Ry,"R.tr.tg_y[100]/D");   
  tnew->Branch("R.tr.tg_th",Rth,"R.tr.tg_th[100]/D");   
  tnew->Branch("R.tr.tg_ph",Rph,"R.tr.tg_ph[100]/D");
 // Focal Plane infomation //
  tnew->Branch("R.tr.x",Rx_fp,"R.tr.x[100]/D");    
  tnew->Branch("R.tr.y",Ry_fp,"R.tr.y[100]/D");
  tnew->Branch("R.tr.th",Rth_fp,"R.tr.th[100]/D"); 
  tnew->Branch("R.tr.ph",Rph_fp,"R.tr.ph[100]/D"); 

 //-------------------------------// 
 //------ Left Arm ---------------//
 //-------------------------------//
 // Scintillation//
  tnew->Branch("LTDC.F1FirstHit",LF1,"LTDC.F1FirstHit[100]/D"); 
  tnew->Branch("L.s2.t_pads",Ls2tpads,"L.s2.t_pads[100]/D");
  tnew->Branch("L.s2.trpad",Ls2trpad,"L.s2.trpad[100]/D");
  tnew->Branch("L.tr.n",&Ltrn,"L.tr.n");
  // PID //
 tnew->Branch("L.sh.asum_c",&Lshsum,"L.sh.asum_c[100]");
 tnew->Branch("L.ps.asum_c",&Lpssum,"L.ps.asum_c[100]/D"); 
 // path length//
 tnew->Branch("L.s2.trpath",ls2pathl,"L.s2.trpath[100]/D"); 
 tnew->Branch("L.tr.pathl",ltrpathl,"L.tr.pathl[100]/D");
 tnew->Branch("L.tr.p",Lp,"L.tr.p[100]/D");      
 tnew->Branch("L.tr.vz",Lz,"L.tr.vz[100]");
 // target positon information //
 tnew->Branch("L.tr.p",Lp,"L.tr.p[100]/D");
 tnew->Branch("L.tr.vz",Lz,"L.tr.vz[100]/D");   
 tnew->Branch("L.tr.tg_x",Lx,"L.tr.tg_x[100]/D");  
 tnew->Branch("L.tr.tg_y",Ly,"L.tr.tg_y[100]/D"); 
 tnew->Branch("L.tr.tg_th",Lth,"L.tr.tg_th[100]/D");
 tnew->Branch("L.tr.tg_ph",Lph,"L.tr.tg_ph[100]/D");
 // Focal Plane infomation //
 tnew->Branch("L.tr.x",Lx_fp,"L.tr.x[100]/D");
 tnew->Branch("L.tr.y",Ly_fp,"L.tr.y[100]/D");
 tnew->Branch("L.tr.th",Lth_fp,"L.tr.th[100]/D");
 tnew->Branch("L.tr.ph",Lph_fp,"L.tr.ph[100]/D"); 
};
*/
