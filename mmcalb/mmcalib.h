#ifndef mmcalib_h
#define mmcalib_h 1

#include <iostream>
#include <fstream>
#include "Tree.h"
#include "Tuning.h"
#include "Param.h"
#include "Setting.h"
#include "define.h"


  bool RHRS=false;
  bool LHRS=false;


//=================================================//
//============= Mmcalib Class ====================//
//=================================================//

class mmcalib :public Tree, public Tuning{

  
 public:
  mmcalib();
  ~mmcalib();

  Setting* set;
    
  void SetRunList(string ifname);
  void NewRoot(string ofroot, string ifname);
  void MTParam_R();
  void MTParam_L();
  void MTParam(string mtparam);
  void Mpt(string ofname);
  void Tuning(string ofname);
  void MakeHist();
  void Fill();
  void Close_tree();

  //------- MakeHist -----//
  TH1F* hRvz;
  TH1F* hRvz_c;
  TH1F* hRvz_Rc;
  TH1F* hRth;
  TH1F* hRth_c;
  TH1F* hRph;
  TH1F* hRph_c;
  TH1F* hLvz;
  TH1F* hLvz_c;
  TH1F* hLvz_Rc;
  TH1F* hLth;
  TH1F* hLth_c;
  TH1F* hLph;
  TH1F* hLph_c;
  
  //------- New Root ----//
  TFile* fnew;
  TTree* tnew;
  double R_tr_vx_c[MAX],R_tr_vy_c[MAX],R_tr_vz_c[MAX],R_tr_tg_th_c[MAX],R_tr_tg_ph_c[MAX];
  double L_tr_vx_c[MAX],L_tr_vy_c[MAX],L_tr_vz_c[MAX],L_tr_tg_th_c[MAX],L_tr_tg_ph_c[MAX];

  //----- MTParam ------//
  string param[10];
  double Pzt[nParamTz],Pzt_L[nParamTz],Pxpt[nParamT],Pxpt_L[nParamT],Pypt[nParamT],Pypt_L[nParamT];
  double parRaster[nParamT_ras],parRaster_L[nParamT_ras];
  double parRaster_L_0,parRaster_L_1,parRaster_L_2,parRaster_L_3;
  double parRaster_R_0,parRaster_R_1,parRaster_R_2,parRaster_R_3;
  //------ Fill -------//
  double XFP, XpFP;
  double YFP, YpFP;
  double XFP_L, XpFP_L;
  double YFP_L, YpFP_L;  


  
};

/////////////////////////////////////////////////////////////////////
mmcalib::mmcalib(){
 set = new Setting();
 set->Initialize();
 cout<<"RHRS: "<<RHRS<<endl;
 cout<<"LHRS: "<<LHRS<<endl;
};

mmcalib::~mmcalib(){};


///////////////////////////////////////////////////////////////////////
void mmcalib::SetRunList(string ifname){
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
    cout<<buf<<endl;
  }

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}

///////////////////////////////////////////////////////////////////////////
void mmcalib::MakeHist(){

  //======== RHRS Hist ==================//
  hRvz=new TH1F("hRvz","",1000,-0.2,0.2);
  set->SetTH1(hRvz,"RHRS z posi hist w/o Tuning ","[m]","Counts");  
  hRvz_c=new TH1F("hRvz_c","",1000,-0.2,0.2);
  set->SetTH1(hRvz_c,"RHRS z posi hist w/ Tuning ","[m]","Counts");  
  hRvz_Rc=new TH1F("hRvz_Rc","",1000,-0.2,0.2);
  set->SetTH1(hRvz_Rc,"RHRS z posi hist w/ Raster & Matrix Tuning ","[m]","Counts");  
  hRth=new TH1F("hRth","",1000,-0.2,0.2);
  set->SetTH1(hRth,"RHRS theta hist w/o Tuning","px/p","Counts");
  hRth_c=new TH1F("hRth_c","",1000,-0.2,0.2);
  set->SetTH1(hRth_c,"RHRS theta hist w/ Tuning","px/p","Counts");
  hRph=new TH1F("hRph","",1000,-0.2,0.2);
  set->SetTH1(hRph,"RHRS phi hist w/o Tuning","px/p","Counts");
  hRph_c=new TH1F("hRph_c","",1000,-0.2,0.2);
  set->SetTH1(hRph_c,"RHRS phi hist w/ Tuning","py/p","Counts");

  //======== LHRS Hist ==================//
  hLvz=new TH1F("hLvz","",1000,-0.2,0.2);
  set->SetTH1(hLvz,"LHRS z posi hist w/o Tuning ","[m]","Counts");  
  hLvz_c=new TH1F("hLvz_c","",1000,-0.2,0.2);
  set->SetTH1(hLvz_c,"LHRS z posi hist w/ Tuning ","[m]","Counts");  
  hLvz_Rc=new TH1F("hLvz_Rc","",1000,-0.2,0.2);
  set->SetTH1(hLvz_Rc,"LHRS z posi hist w/ Raster & Matrix Tuning ","[m]","Counts");  
  hLth=new TH1F("hLth","",1000,-0.2,0.2);
  set->SetTH1(hLth,"LHRS theta hist w/o Tuning","px/p","Counts");
  hLth_c=new TH1F("hLth_c","",1000,-0.2,0.2);
  set->SetTH1(hLth_c,"LHRS theta hist w/ Tuning","px/p","Counts");
  hLph=new TH1F("hLph","",1000,-0.2,0.2);
  set->SetTH1(hLph,"LHRS phi hist w/o Tuning","py/p","Counts");
  hLph_c=new TH1F("hLph_c","",1000,-0.2,0.2);
  set->SetTH1(hLph_c,"LHRS phi hist w/ Tuning","py/p","Counts");

  //======= Set Line Color ========//
  hRvz_c->SetLineColor(2);
  hRvz_Rc->SetLineColor(3);
  hRth_c->SetLineColor(2);
  hRph_c->SetLineColor(2);
  hLvz_c->SetLineColor(2);
  hLvz_Rc->SetLineColor(3);
  hLth_c->SetLineColor(2);
  hLph_c->SetLineColor(2);

  
  
  
}


////////////////////////////////////////////////////////////////////////////

void mmcalib::NewRoot(string ofroot, string ifname){

  fnew = new TFile(Form("%s",ofroot.c_str()),"recreate");

  //========== Get Time ==============//

  time_t t=time(nullptr);
  const tm* lt=localtime(&t);
  stringstream s;
  s<<"20";
  s<<lt->tm_year-100;
  s<<"-";
  s<<lt->tm_mon+1;
  s<<"-";
  s<<lt->tm_mday;
  s<<" ";
  s<<lt->tm_hour;
  s<<":";
  s<<lt->tm_min;
  string TIME=s.str();

  tnew=new TTree("T",Form("matrix tuning root run file: %s [%s]",ifname.c_str(),TIME.c_str()));

    tnew = tree->CloneTree(0);

  /*
  //------ RHRS ---------//
  tnew->Branch("R.tr.tg_th_c"           , R_tr_tg_th_c ,"R.tr.tg_th_c[10]/D");
  tnew->Branch("R.tr.tg_ph_c"           , R_tr_tg_ph_c ,"R.tr.tg_ph_c[10]/D");
  tnew->Branch("R.tr.vx_c"              , R_tr_vx_c    ,"R.tr.vx_c[10]/D"   );
  tnew->Branch("R.tr.vy_c"              , R_tr_vy_c    ,"R.tr.vy_c[10]/D"   );
  tnew->Branch("R.tr.vz_c"              , R_tr_vz_c    ,"R.tr.vz_c[10]/D"   );
  //------ LHRS --------//
  tnew->Branch("L.tr.tg_th_c"           , L_tr_tg_th_c ,"L.tr.tg_th_c[10]/D");
  tnew->Branch("L.tr.tg_ph_c"           , L_tr_tg_ph_c ,"L.tr.tg_ph_c[10]/D");
  tnew->Branch("L.tr.vx_c"              , L_tr_vx_c    ,"L.tr.vx_c[10]/D"   );
  tnew->Branch("L.tr.vy_c"              , L_tr_vy_c    ,"L.tr.vy_c[10]/D"   );
  tnew->Branch("L.tr.vz_c"              , L_tr_vz_c    ,"L.tr.vz_c[10]/D"   );
     */

}

/////////////////////////////////////////////////////////////////////////////

void mmcalib::MTParam(string mtparam){


  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }

  if(RHRS) {MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;}
  if(LHRS) {MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;}
 
  
}

/////////////////////////////////////////////////////////////////////////////
void mmcalib::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//

  

  
  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);

   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;

   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param[1].c_str()); // optimized
    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);

  for (int i=0;i<nParamT_ras;i++){

    Mras >> parRaster[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) parRaster[i] = 0.0;
  }

  parRaster_R_0=parRaster[0];
  parRaster_R_1=parRaster[1];
  parRaster_R_2=parRaster[2];
  parRaster_R_3=parRaster[3];
  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);

    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param[3].c_str()); // optimized  
  ifstream Mypt(name_Mypt);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
  }
  Mypt.close();    



  
};
//////////////////////////////////////////////////////////////

void mmcalib::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L); 
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param[5].c_str()); // optimized
    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);

  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> parRaster_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) parRaster_L[i] = 0.0;
  }

  parRaster_L_0=parRaster_L[0];
  parRaster_L_1=parRaster_L[1];
  parRaster_L_2=parRaster_L[2];
  parRaster_L_3=parRaster_L[3];
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param[7].c_str()); // optimized
    cout<<"LHRS phi parameters file: "<<name_Mypt_L<<endl;
  ifstream Mypt_L(name_Mypt_L);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;
  }
  Mypt_L.close();    





}


/////////////////////////////////////////////////////////////

void mmcalib::Fill(){
  cout<<"============================"<<endl;
  cout<<"======= start fill ========="<<endl;
  cout<<"============================"<<endl;
  ENum=tree->GetEntries();
  cout<<"Events: "<<ENum<<endl;
  div_t d=div(ENum,10000);

  for(int i=0;i<ENum;i++){

   tree->GetEntry(i);

   //======== w/o Tuning Hist =======//
	hRvz->Fill(R_tr_vz[0]);
	hRth->Fill(R_tr_tg_th[0]);
	hRph->Fill(R_tr_tg_ph[0]);	
	hLvz->Fill(L_tr_vz[0]);
	hLth->Fill(L_tr_tg_th[0]);
	hLph->Fill(L_tr_tg_ph[0]);
	
   // ===== Initialization ====== //    
      for(int j=0;j<MAX;j++){
	if(RHRS){
      //      R_tr_vx[j]  = -2222.0;
      //      R_tr_vy[j]  = -2222.0;
           R_tr_vz[j]  = -2222.0;
      //      R_tr_tg_th[j] = -2222.0;
      //      R_tr_tg_ph[j] = -2222.0;
	}
	if(LHRS){
      //      L_tr_vx[j]  = -2222.0;
      //      L_tr_vy[j]  = -2222.0;
           L_tr_vz[j]  = -2222.0;
      //      L_tr_tg_th[j] = -2222.0;
      //      L_tr_tg_ph[j] = -2222.0;
	}
          }



    
      // ----------------------------------------------------------- //
      // --------------------- RHRS ------------------------ ------- //
      // ----------------------------------------------------------- //
    
    if(RHRS){
      R_tr_x[0]  = (R_tr_x[0]-XFPm)/XFPr;
      R_tr_th[0] = (R_tr_th[0]-XpFPm)/XpFPr;
      R_tr_y[0]  = (R_tr_y[0]-YFPm)/YFPr;
      R_tr_ph[0] = (R_tr_ph[0]-YpFPm)/YpFPr;

      R_tr_vz[0] = calcf2t_zt(Pzt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0]);
      R_tr_vz[0] = R_tr_vz[0] * Ztr + Ztm;

      //-------w/ Tuning z calib ----------//
      hRvz_c->Fill(R_tr_vz[0]);

      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //


      double RasterCor = calcRasterCor(R_Ras_x, parRaster_R_2, parRaster_R_0);
      RasterCor = RasterCor/tan(hrs_ang);
      R_tr_vz[0]=R_tr_vz[0]+RasterCor;// w/ Raster correction
      hRvz_Rc->Fill(R_tr_vz[0]);


      
      double Zt=R_tr_vz[0];
      Zt= (Zt-Ztm)/Ztr;
      
      R_tr_tg_th[0]  =  calcf2t_ang(Pxpt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0],Zt);
      R_tr_tg_ph[0]  =  calcf2t_ang(Pypt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0],Zt);
      R_tr_tg_th[0]  =  R_tr_tg_th[0]*Xptr +Xptm;
      R_tr_tg_ph[0]  =  R_tr_tg_ph[0]*Yptr +Yptm;

      R_tr_x[0]  = R_tr_x[0]  * XFPr  + XFPm;
      R_tr_y[0]  = R_tr_y[0]  * YFPr  + YFPm;
      R_tr_th[0] = R_tr_th[0] * XpFPr + XpFPm;
      R_tr_ph[0] = R_tr_ph[0] * YpFPr + YpFPm;

  	hRth_c->Fill(R_tr_tg_th[0]);
	hRph_c->Fill(R_tr_tg_ph[0]);


    }//end RHRS Fill
    


	
      // ----------------------------------------------------------- //
      // --------------------- LHRS ------------------------ ------- //
      // ----------------------------------------------------------- //
	

    if(LHRS){
      L_tr_x[0]  = (L_tr_x[0]-XFPm)/XFPr;
      L_tr_th[0] = (L_tr_th[0]-XpFPm)/XpFPr;
      L_tr_y[0]  = (L_tr_y[0]-YFPm)/YFPr;
      L_tr_ph[0] = (L_tr_ph[0]-YpFPm)/YpFPr;
      L_tr_vz[0] =  calcf2t_zt(Pzt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0]);
      L_tr_vz[0] = L_tr_vz[0] * Ztr + Ztm;


      //-------w/ Tuning z calib ----------//
      hLvz_c->Fill(L_tr_vz[0]);

      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //

      double RasterCor_L = calcRasterCor(L_Ras_x, parRaster_L_2, parRaster_L_0);
      RasterCor_L = RasterCor_L/tan(hrs_ang);
      L_tr_vz[0]=L_tr_vz[0]+RasterCor;
      hLvz_Rc->Fill(L_tr_vz[0]);

      double Zt_L=L_tr_vz[0];
      Zt_L= (Zt_L-Ztm)/Ztr;
      L_tr_tg_th[0]  =  calcf2t_ang(Pxpt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0],Zt_L);
      L_tr_tg_ph[0] =   calcf2t_ang(Pypt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0],Zt_L);
      L_tr_tg_th[0]  = L_tr_tg_th[0]*Xptr +Xptm;
      L_tr_tg_ph[0]  = L_tr_tg_ph[0]*Yptr +Yptm;

      L_tr_x[0]  = L_tr_x[0]  * XFPr  + XFPm;
      L_tr_y[0]  = L_tr_y[0]  * YFPr  + YFPm;
      L_tr_th[0] = L_tr_th[0] * XpFPr + XpFPm;
      L_tr_ph[0] = L_tr_ph[0] * YpFPr + YpFPm;


	hLth_c->Fill(L_tr_tg_th[0]);
	hLph_c->Fill(L_tr_tg_ph[0]);	

    }//End LHRS Fill

    tnew->Fill();
    if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }

  cout<<"Fill is done !"<<endl;
}

//======================================================================//

void mmcalib::Close_tree(){


  hRvz->Write();
  hRth->Write();
  hRph->Write();
  hLvz->Write();
  hLth->Write();
  hLph->Write();
  if(RHRS){
  hRvz_c->Write();
  hRvz_Rc->Write();
  hRth_c->Write();
  hRph_c->Write();
  }
  if(LHRS){
  hLvz_c->Write();
  hLvz_Rc->Write();  
  hLth_c->Write();
  hLph_c->Write();
  }
  
  tnew->Write();
  fnew->Close();
  
};



#endif
