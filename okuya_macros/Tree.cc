#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Tree.h"

Tree::Tree()
{
  tree = new TChain("T");
  //  cout<<"start tree class"<<endl;
}
//////////////////////////////////////////////////
Tree::~Tree()
{
  //cout<<"end tree class"<<endl;
}
//////////////////////////////////////////////////
void Tree::pack_tree()
{
 
  ENum = tree->GetEntries();  
  tree->SetBranchStatus("*" ,0);
  //tree->SetBranchStatus("fEvtHdr.fRun"               ,1);  tree->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  //  tree->SetBranchStatus("fEvtHdr.fEvtNu"               ,1);   tree->SetBranchAddress("fEvtHdr.fEvtNum", &evid    );
  tree->SetBranchStatus("runnum"               ,1);  tree->SetBranchAddress("runnum", &runnum );
  tree->SetBranchStatus("HALLA_p"               ,1);  tree->SetBranchAddress("HALLA_p", &HALLA_p );
  tree->SetBranchStatus("L.tr.n"               ,1);  tree->SetBranchAddress("L.tr.n"               ,&L_tr_n              );
  tree->SetBranchStatus("R.tr.n"               ,1);  tree->SetBranchAddress("R.tr.n"               ,&R_tr_n              );
  //  tree->SetBranchStatus("DR.evtypebits"        ,1);  tree->SetBranchAddress("DR.evtypebits",&DR_evtypebits);
  //tree->SetBranchStatus("DR.T1"                ,1);  tree->SetBranchAddress("DR.T1",&DR_T1);
  //tree->SetBranchStatus("DR.T2"                ,1);  tree->SetBranchAddress("DR.T2",&DR_T2);
  //tree->SetBranchStatus("DR.T3"                ,1);  tree->SetBranchAddress("DR.T3",&DR_T3);
  //tree->SetBranchStatus("DR.T4"                ,1);  tree->SetBranchAddress("DR.T4",&DR_T4);
  //  tree->SetBranchStatus("DR.T5"                ,1);  tree->SetBranchAddress("DR.T5",&DR_T5);
  //tree->SetBranchStatus("DR.T6"                ,1);  tree->SetBranchAddress("DR.T6",&DR_T6);
  //tree->SetBranchStatus("DR.T7"                ,1);  tree->SetBranchAddress("DR.T7",&DR_T7);
  //tree->SetBranchStatus("DR.T8"                ,1);  tree->SetBranchAddress("DR.T8",&DR_T8);
  tree->SetBranchStatus("DL.evtype"            ,1);  tree->SetBranchAddress("DL.evtype"       ,&L_evtype        );
  tree->SetBranchStatus("DR.evtype"            ,1);  tree->SetBranchAddress("DR.evtype"       ,&R_evtype        );
}
//////////////////////////////////////////////////
void Tree::chain_tree(string ifname)
{
  tree->Add(Form("%s",ifname.c_str()));
  //  cout<<ifname<<endl;

  ENum = tree->GetEntries();
  tree->SetBranchStatus("*" ,0);
  tree->SetBranchStatus("L.tr.n"               ,1);  tree->SetBranchAddress("L.tr.n"               ,&L_tr_n              );
  tree->SetBranchStatus("R.tr.n"               ,1);  tree->SetBranchAddress("R.tr.n"               ,&R_tr_n              );
}
//////////////////////////////////////////////////
void Tree::readtreeHRSLR()
{
  readtreeHRSL();
  readtreeHRSR();
}
//////////////////////////////////////////////////
void Tree::readtreeHRSL()
{
  
  readtreeRasterL();
  readtreeS0L();
  readtreeS2L();
  readtreeTrackL();
  readtreeF1TDCL();
  readtreeGCL();
  //  readtreePSL();
  //  readtreeShL();
  //readtreeMiscL();
}
//////////////////////////////////////////////////
void Tree::readtreeHRSR()
{
  readtreeRasterR();
  readtreeS0R();
  readtreeS2R();
  readtreeA1R();
  readtreeA2R();
  readtreeTrackR();
  readtreeF1TDCR();
  //  readtreeGCR();
  //  readtreePSR();
  //  readtreeShR();

  //readtreeMiscR();
}
//////////////////////////////////////////////////
void Tree::readtreeRasterL(){
  tree->SetBranchStatus("Lrb.Raster2.rawcur.x",  1); tree->SetBranchAddress("Lrb.Raster2.rawcur.x"              , &L_Ras_x   );
  tree->SetBranchStatus("Lrb.Raster2.rawcur.y",  1); tree->SetBranchAddress("Lrb.Raster2.rawcur.y"              , &L_Ras_y   );
}

//////////////////////////////////////////////////
void Tree::readtreeS0L()
{
    tree->SetBranchStatus("L.s0.la"              ,1);  tree->SetBranchAddress("L.s0.la"              , L_s0_la             );
    tree->SetBranchStatus("L.s0.la_c"            ,1);  tree->SetBranchAddress("L.s0.la_c"            , L_s0_la_c           );
    tree->SetBranchStatus("L.s0.la_p"            ,1);  tree->SetBranchAddress("L.s0.la_p"            , L_s0_la_p           );
    tree->SetBranchStatus("L.s0.lt"              ,1);  tree->SetBranchAddress("L.s0.lt"              , L_s0_lt             );
    tree->SetBranchStatus("L.s0.lt_c"            ,1);  tree->SetBranchAddress("L.s0.lt_c"            , L_s0_lt_c           );
    tree->SetBranchStatus("L.s0.ra"              ,1);  tree->SetBranchAddress("L.s0.ra"              , L_s0_ra             );
    tree->SetBranchStatus("L.s0.ra_c"            ,1);  tree->SetBranchAddress("L.s0.ra_c"            , L_s0_ra_c           );
    tree->SetBranchStatus("L.s0.ra_p"            ,1);  tree->SetBranchAddress("L.s0.ra_p"            , L_s0_ra_p           );
    tree->SetBranchStatus("L.s0.rt"              ,1);  tree->SetBranchAddress("L.s0.rt"              , L_s0_rt             );
    tree->SetBranchStatus("L.s0.rt_c"            ,1);  tree->SetBranchAddress("L.s0.rt_c"            , L_s0_rt_c           );
  tree->SetBranchStatus("L.s0.t_pads"          ,1);  tree->SetBranchAddress("L.s0.t_pads"          , L_s0_t_pads         );
  tree->SetBranchStatus("L.s0.time"            ,1);  tree->SetBranchAddress("L.s0.time"            , L_s0_time           );
  tree->SetBranchStatus("L.s0.dedx"            ,1);  tree->SetBranchAddress("L.s0.dedx"            , L_s0_dedx           );
  tree->SetBranchStatus("L.s0.trdy"            ,1);  tree->SetBranchAddress("L.s0.trdy"            , L_s0_trdy           );
  //    tree->SetBranchStatus("L.s0.troff"           ,1);  tree->SetBranchAddress("L.s0.troff"           , L_s0_troff          );
    tree->SetBranchStatus("L.s0.trpad"           ,1);  tree->SetBranchAddress("L.s0.trpad"           , L_s0_trpad          );
  tree->SetBranchStatus("L.s0.trpath"          ,1);  tree->SetBranchAddress("L.s0.trpath"          , L_s0_trpath         );
  tree->SetBranchStatus("L.s0.trx"             ,1);  tree->SetBranchAddress("L.s0.trx"             , L_s0_trx            );
  tree->SetBranchStatus("L.s0.try"             ,1);  tree->SetBranchAddress("L.s0.try"             , L_s0_try            );
  tree->SetBranchStatus("L.s0.lnhits"          ,1);  tree->SetBranchAddress("L.s0.lnhits"          , L_s0_lnhits         );
    tree->SetBranchStatus("L.s0.lpeak"           ,1);  tree->SetBranchAddress("L.s0.lpeak"           , L_s0_lpeak          );
    tree->SetBranchStatus("L.s0.lt_fadc"         ,1);  tree->SetBranchAddress("L.s0.lt_fadc"         , L_s0_lt_fadc        );
    tree->SetBranchStatus("L.s0.ltc_fadc"        ,1);  tree->SetBranchAddress("L.s0.ltc_fadc"        , L_s0_ltc_fadc       );
    //  tree->SetBranchStatus("L.s0.lbadped"         ,1);  tree->SetBranchAddress("L.s0.lbadped"         , L_s0_lbadped        );
    //  tree->SetBranchStatus("L.s0.lunderflow"      ,1);  tree->SetBranchAddress("L.s0.lunderflow"      , L_s0_lunderflow     );
    //  tree->SetBranchStatus("L.s0.loverflow"       ,1);  tree->SetBranchAddress("L.s0.loverflow"       , L_s0_loverflow      );
    tree->SetBranchStatus("L.s0.rnhits"          ,1);  tree->SetBranchAddress("L.s0.rnhits"          , L_s0_rnhits         );
    tree->SetBranchStatus("L.s0.rpeak"           ,1);  tree->SetBranchAddress("L.s0.rpeak"           , L_s0_rpeak          );
    tree->SetBranchStatus("L.s0.rt_fadc"         ,1);  tree->SetBranchAddress("L.s0.rt_fadc"         , L_s0_rt_fadc        );
    tree->SetBranchStatus("L.s0.rtc_fadc"        ,1);  tree->SetBranchAddress("L.s0.rtc_fadc"        , L_s0_rtc_fadc       );
    //  tree->SetBranchStatus("L.s0.rbadped"         ,1);  tree->SetBranchAddress("L.s0.rbadped"         , L_s0_rbadped        );
    //  tree->SetBranchStatus("L.s0.runderflow"      ,1);  tree->SetBranchAddress("L.s0.runderflow"      , L_s0_runderflow     );
    //  tree->SetBranchStatus("L.s0.roverflow"       ,1);  tree->SetBranchAddress("L.s0.roverflow"       , L_s0_roverflow      );
}
/////////////////////////////////////////////////////////////////////////////////////////
void Tree::readtreeS2L()
{
    tree->SetBranchStatus("L.s2.la"              ,1);  tree->SetBranchAddress("L.s2.la"              , L_s2_la             );
    tree->SetBranchStatus("L.s2.la_c"            ,1);  tree->SetBranchAddress("L.s2.la_c"            , L_s2_la_c           );
    tree->SetBranchStatus("L.s2.la_p"            ,1);  tree->SetBranchAddress("L.s2.la_p"            , L_s2_la_p           );
    tree->SetBranchStatus("L.s2.lt"              ,1);  tree->SetBranchAddress("L.s2.lt"              , L_s2_lt             );
    tree->SetBranchStatus("L.s2.lt_c"            ,1);  tree->SetBranchAddress("L.s2.lt_c"            , L_s2_lt_c           );
    tree->SetBranchStatus("L.s2.ra"              ,1);  tree->SetBranchAddress("L.s2.ra"              , L_s2_ra             );
    tree->SetBranchStatus("L.s2.ra_c"            ,1);  tree->SetBranchAddress("L.s2.ra_c"            , L_s2_ra_c           );
    tree->SetBranchStatus("L.s2.ra_p"            ,1);  tree->SetBranchAddress("L.s2.ra_p"            , L_s2_ra_p           );
    tree->SetBranchStatus("L.s2.rt"              ,1);  tree->SetBranchAddress("L.s2.rt"              , L_s2_rt             );
    tree->SetBranchStatus("L.s2.rt_c"            ,1);  tree->SetBranchAddress("L.s2.rt_c"            , L_s2_rt_c           );
    tree->SetBranchStatus("L.s2.t_pads"          ,1);  tree->SetBranchAddress("L.s2.t_pads"          , L_s2_t_pads         );
    tree->SetBranchStatus("L.s2.time"            ,1);  tree->SetBranchAddress("L.s2.time"            , L_s2_time           );
  tree->SetBranchStatus("L.s2.dedx"            ,1);  tree->SetBranchAddress("L.s2.dedx"            , L_s2_dedx           );
  tree->SetBranchStatus("L.s2.trdx"            ,1);  tree->SetBranchAddress("L.s2.trdx"            , L_s2_trdx           );
  //    tree->SetBranchStatus("L.s2.troff"           ,1);  tree->SetBranchAddress("L.s2.troff"           , L_s2_troff          );
  tree->SetBranchStatus("L.s2.trpad"           ,1);  tree->SetBranchAddress("L.s2.trpad"           , L_s2_trpad          );
  tree->SetBranchStatus("L.s2.trpath"          ,1);  tree->SetBranchAddress("L.s2.trpath"          , L_s2_trpath         );
  tree->SetBranchStatus("L.s2.trx"             ,1);  tree->SetBranchAddress("L.s2.trx"             , L_s2_trx            );
  tree->SetBranchStatus("L.s2.try"             ,1);  tree->SetBranchAddress("L.s2.try"             , L_s2_try            );
    tree->SetBranchStatus("L.s2.lnhits"          ,1);  tree->SetBranchAddress("L.s2.lnhits"          , L_s2_lnhits         );
    tree->SetBranchStatus("L.s2.lpeak"           ,1);  tree->SetBranchAddress("L.s2.lpeak"           , L_s2_lpeak          );
    tree->SetBranchStatus("L.s2.lt_fadc"         ,1);  tree->SetBranchAddress("L.s2.lt_fadc"         , L_s2_lt_fadc        );
    tree->SetBranchStatus("L.s2.ltc_fadc"        ,1);  tree->SetBranchAddress("L.s2.ltc_fadc"        , L_s2_ltc_fadc       );
    //  tree->SetBranchStatus("L.s2.lbadped"         ,1);  tree->SetBranchAddress("L.s2.lbadped"         , L_s2_lbadped        );
    //  tree->SetBranchStatus("L.s2.lunderflow"      ,1);  tree->SetBranchAddress("L.s2.lunderflow"      , L_s2_lunderflow     );
    //  tree->SetBranchStatus("L.s2.loverflow"       ,1);  tree->SetBranchAddress("L.s2.loverflow"       , L_s2_loverflow      );
    tree->SetBranchStatus("L.s2.rnhits"          ,1);  tree->SetBranchAddress("L.s2.rnhits"          , L_s2_rnhits         );
    tree->SetBranchStatus("L.s2.rpeak"           ,1);  tree->SetBranchAddress("L.s2.rpeak"           , L_s2_rpeak          );
    tree->SetBranchStatus("L.s2.rt_fadc"         ,1);  tree->SetBranchAddress("L.s2.rt_fadc"         , L_s2_rt_fadc        );
    tree->SetBranchStatus("L.s2.rtc_fadc"        ,1);  tree->SetBranchAddress("L.s2.rtc_fadc"        , L_s2_rtc_fadc       );
    //  tree->SetBranchStatus("L.s2.rbadped"         ,1);  tree->SetBranchAddress("L.s2.rbadped"         , L_s2_rbadped        );
    //  tree->SetBranchStatus("L.s2.runderflow"      ,1);  tree->SetBranchAddress("L.s2.runderflow"      , L_s2_runderflow     );
    //  tree->SetBranchStatus("L.s2.roverflow"       ,1);  tree->SetBranchAddress("L.s2.roverflow"       , L_s2_roverflow      );
}                                 
/////////////////////////////////
void Tree::readtreeTrackL(){
  tree->SetBranchStatus("L.tr.n"               ,1);  tree->SetBranchAddress("L.tr.n"               ,&L_tr_n              );
    tree->SetBranchStatus("L.tr.flag"            ,1);  tree->SetBranchAddress("L.tr.flag"            , L_tr_flag           );
    tree->SetBranchStatus("L.tr.ndof"            ,1);  tree->SetBranchAddress("L.tr.ndof"            , L_tr_ndof           );
   tree->SetBranchStatus("L.tr.chi2"            ,1);  tree->SetBranchAddress("L.tr.chi2"            , L_tr_chi2           );
     tree->SetBranchStatus("L.tr.beta"            ,1);  tree->SetBranchAddress("L.tr.beta"            , L_tr_beta           );
     tree->SetBranchStatus("L.tr.d_x"             ,1);  tree->SetBranchAddress("L.tr.d_x"             , L_tr_d_x            );
     tree->SetBranchStatus("L.tr.d_y"             ,1);  tree->SetBranchAddress("L.tr.d_y"             , L_tr_d_y            );
     tree->SetBranchStatus("L.tr.d_th"            ,1);  tree->SetBranchAddress("L.tr.d_th"            , L_tr_d_th           );
     tree->SetBranchStatus("L.tr.d_ph"            ,1);  tree->SetBranchAddress("L.tr.d_ph"            , L_tr_d_ph           );
     tree->SetBranchStatus("L.tr.r_x"             ,1);  tree->SetBranchAddress("L.tr.r_x"             , L_tr_r_x            );
     tree->SetBranchStatus("L.tr.r_y"             ,1);  tree->SetBranchAddress("L.tr.r_y"             , L_tr_r_y            );
     tree->SetBranchStatus("L.tr.r_th"            ,1);  tree->SetBranchAddress("L.tr.r_th"            , L_tr_r_th           );
     tree->SetBranchStatus("L.tr.r_ph"            ,1);  tree->SetBranchAddress("L.tr.r_ph"            , L_tr_r_ph           );
  tree->SetBranchStatus("L.tr.x"               ,1);  tree->SetBranchAddress("L.tr.x"               , L_tr_x              );
  tree->SetBranchStatus("L.tr.y"               ,1);  tree->SetBranchAddress("L.tr.y"               , L_tr_y              );
  tree->SetBranchStatus("L.tr.th"              ,1);  tree->SetBranchAddress("L.tr.th"              , L_tr_th             );
  tree->SetBranchStatus("L.tr.ph"              ,1);  tree->SetBranchAddress("L.tr.ph"              , L_tr_ph             );
   tree->SetBranchStatus("L.tr.time"            ,1);  tree->SetBranchAddress("L.tr.time"            , L_tr_time           );
  tree->SetBranchStatus("L.tr.p"               ,1);  tree->SetBranchAddress("L.tr.p"               , L_tr_p              );
  tree->SetBranchStatus("L.tr.pathl"           ,1);  tree->SetBranchAddress("L.tr.pathl"           , L_tr_pathl          );
  tree->SetBranchStatus("L.tr.px"              ,1);  tree->SetBranchAddress("L.tr.px"              , L_tr_px             );
  tree->SetBranchStatus("L.tr.py"              ,1);  tree->SetBranchAddress("L.tr.py"              , L_tr_py             );
  tree->SetBranchStatus("L.tr.pz"              ,1);  tree->SetBranchAddress("L.tr.pz"              , L_tr_pz             );
  tree->SetBranchStatus("L.tr.tg_dp"           ,1);  tree->SetBranchAddress("L.tr.tg_dp"           , L_tr_tg_dp          );
  tree->SetBranchStatus("L.tr.tg_y"            ,1);  tree->SetBranchAddress("L.tr.tg_y"            , L_tr_tg_y           );
  tree->SetBranchStatus("L.tr.tg_th"           ,1);  tree->SetBranchAddress("L.tr.tg_th"           , L_tr_tg_th          );
  tree->SetBranchStatus("L.tr.tg_ph"           ,1);  tree->SetBranchAddress("L.tr.tg_ph"           , L_tr_tg_ph          );
  tree->SetBranchStatus("L.tr.vx"              ,1);  tree->SetBranchAddress("L.tr.vx"              , L_tr_vx             );
  tree->SetBranchStatus("L.tr.vy"              ,1);  tree->SetBranchAddress("L.tr.vy"              , L_tr_vy             );
  tree->SetBranchStatus("L.tr.vz"              ,1);  tree->SetBranchAddress("L.tr.vz"              , L_tr_vz             );
}
//////////////////////////////////////////////////
void Tree::readtreeF1TDCL()
{
  tree->SetBranchStatus("LTDC.F1FirstHit"      ,1);  tree->SetBranchAddress("LTDC.F1FirstHit"      , LTDC_F1FirstHit);
}
//////////////////////////////////////////////////
void Tree::convertF1TDCL(ParamMan *param)
{

  bool F1Shift_L=true;

//  double F1res=  param->GetF1reso();
  //  doulbe F1shift= param->GetF1ShiftOffset();
  double LF1shift_off[LS2];

  //===== Reference check ======//

  LF1Ref[0] = LTDC_F1FirstHit[30];
  LF1Ref[1] = LTDC_F1FirstHit[37];
  if(LF1Ref[0]==0)LF1Ref[0]=LF1Ref[1];
  if(LF1Ref[1]==0)LF1Ref[1]=LF1Ref[0];
  
  //=============================//  
  
  //S2
  for(int i=0;i<LS2;i++){
    LS2T_F1TDC[i]  =  LTDC_F1FirstHit[i];
    LS2B_F1TDC[i]  =  LTDC_F1FirstHit[i+48];
      
    if(LF1Ref[0]>1. && LF1Ref[1]>1. && LS2T_F1TDC[i]>1. && LS2B_F1TDC[i]>1.){

      
      
      LS2T_F1[i] = LS2T_F1TDC[i] - LF1Ref[0] ;
      LS2B_F1[i] = LS2B_F1TDC[i] - LF1Ref[1] ;

      LS2T_F1time[i] = param->time(1,i,0,0, LS2T_F1[i]  , L_s2_ra_p[i]);
      LS2B_F1time[i] = param->time(1,i,0,1, LS2B_F1[i]  , L_s2_la_p[i]);
      LS2_F1time[i]  = 0.5*(LS2T_F1time[i] + LS2B_F1time[i]);
      ltof[i] = LS2_F1time[i];

      
      //====== modified by itabashi ========//
      if(F1Shift_L){


	LS2T_F1_b[i] = - LS2T_F1[i];
	LS2B_F1_b[i] = - LS2B_F1[i];

	
	LF1shift_off[i]=param->GetF1Shift(i,0,0);
	if(LS2T_F1[i]>pow(2,15))LS2T_F1[i] -= LF1shift_off[i];
	if(LS2B_F1[i]>pow(2,15))LS2B_F1[i] -= LF1shift_off[i];

	LS2T_F1_c[i]=LS2T_F1[i];
	LS2B_F1_c[i]=LS2B_F1[i];

	LS2T_F1time[i] = param->time(1,i,0,0, LS2T_F1[i]  , L_s2_ra_p[i]);
	LS2B_F1time[i] = param->time(1,i,0,1, LS2B_F1[i]  , L_s2_la_p[i]);
	LS2_F1time_c[i]  = 0.5*(LS2T_F1time[i] + LS2B_F1time[i]);
	ltof_c[i]=LS2_F1time[i];


	//LS2T_F1_g[i] =  - LF1Ref[0];
	//	LS2B_F1_g[i] =  - LF1Ref[1];


	//		LS2T_F1_g[i] = LF1Ref[0] - LS2T_F1[40];
	//	LS2B_F1_g[i] = LF1Ref[1] - LS2T_F1[40];

	//	LS2T_F1_g[i] = LS2T_F1[40]  - LF1Ref[0] ;
	//	LS2B_F1_g[i] = LS2T_F1[40]  - LF1Ref[1] ;

	LS2T_F1_g[i] = 0.0 ;
	LS2B_F1_g[i] = 0.0 ;
	LS2T_F1time_g[i] = 0.0;
	LS2B_F1time_g[i] = 0.0;

	//	LS2T_F1time_g[i] = param->time(1,i,0,0, LS2T_F1_g[i]  , L_s2_ra_p[i]);
	//	LS2B_F1time_g[i] = param->time(1,i,0,1, LS2B_F1_g[i]  , L_s2_la_p[i]);
	LS2_F1time_g[i]  = 0.5*(LS2T_F1time_g[i] + LS2B_F1time_g[i]);
	ltof_g[i] = LS2_F1time_g[i];
	
	
      }
      //=====================================//




      
    }

    else {

      LS2T_F1time[i] = -9999.;
      LS2B_F1time[i] = -9999.;
      LS2_F1time[i]  = -9999.;

      
    }

    
  }

  //S0
  for(int i=0;i<LS0;i++){
    LS0T_F1TDC[i]  =  LTDC_F1FirstHit[i+27];
    LS0B_F1TDC[i]  =  LTDC_F1FirstHit[i+28];
    LS0T_F1time[i] = param->time(2,i,0,0, LS0T_F1TDC[i] - LTDC_F1FirstHit[37]  );
    LS0B_F1time[i] = param->time(2,i,0,1, LS0B_F1TDC[i] - LTDC_F1FirstHit[37] );
    LS0_F1time[i]  = 0.5*(LS0T_F1time[i] + LS0B_F1time[i]);
  }

  //RF
  LRF_F1time = param->time(3,0,0,0, LTDC_F1FirstHit[37] - LTDC_F1FirstHit[37]  );

  

}

//////////////////////////////////////////////////

void Tree::convertF1TDCL_rhit(ParamMan *param, int seg)
{



  bool F1Shift_L=true;

//  double F1res=  param->GetF1reso();
  double F1shift= param->GetF1ShiftOffset();
  double LF1shift_off[LS2];

  
  //S2
  LF1Ref[0] = LTDC_F1FirstHit[30];
  LF1Ref[1] = LTDC_F1FirstHit[37];
  
  LS2T_F1TDC[seg]  =  LTDC_F1FirstHit[seg];
  LS2B_F1TDC[seg]  =  LTDC_F1FirstHit[seg+48];

  //  cout<<" LS2T "<<LS2T_F1TDC[seg]<<" LS2B "<<LS2B_F1TDC[seg]<<" LS2T_ref "<<LF1Ref[0]<<" LS2B_ref "<<LF1Ref[1]<<endl;
  if(LF1Ref[0]>1. && LF1Ref[1]>1. && LS2T_F1TDC[seg]<1. && LS2B_F1TDC[seg]<1.){  

      LS2T_F1[seg] = LS2T_F1TDC[seg] - LF1Ref[0] - F1shift;
      LS2B_F1[seg] = LS2B_F1TDC[seg] - LF1Ref[1] - F1shift;
      
      LS2T_F1time[seg] = param->time(1,seg,0,0, LS2T_F1[seg]  , L_s2_ra_p[seg]);
      LS2B_F1time[seg] = param->time(1,seg,0,1, LS2B_F1[seg]  , L_s2_la_p[seg]);
      LS2_F1time[seg]  = 0.5*(LS2T_F1time[seg] + LS2B_F1time[seg]);
      ltof[seg]=LS2_F1time[seg];

      
      //====== modsegfseged by segtabashseg ========//
      if(F1Shift_L){
	
	LF1shift_off[seg]=param->GetF1Shift(seg,0,0);
	if(LS2T_F1[seg]>pow(2,15))LS2T_F1[seg] -= LF1shift_off[seg];
	if(LS2B_F1[seg]>pow(2,15))LS2B_F1[seg] -= LF1shift_off[seg];
	LS2T_F1time[seg] = param->time(1,seg,0,0, LS2T_F1[seg]  , L_s2_ra_p[seg]);
	LS2B_F1time[seg] = param->time(1,seg,0,1, LS2B_F1[seg]  , L_s2_la_p[seg]);
	LS2_F1time_c[seg]  = 0.5*(LS2T_F1time[seg] + LS2B_F1time[seg]);
	ltof_c[seg]=LS2_F1time[seg];

	
      }
      //=====================================//

      
  }
    
}




//////////////////////////////////////////////////
void Tree::readtreeRasterR(){
  tree->SetBranchStatus("Rrb.Raster2.rawcur.x",  1); tree->SetBranchAddress("Rrb.Raster2.rawcur.x"              , &R_Ras_x   );
  tree->SetBranchStatus("Rrb.Raster2.rawcur.y",  1); tree->SetBranchAddress("Rrb.Raster2.rawcur.y"              , &R_Ras_y   );
}
//////////////////////////////////////////////////
void Tree::readtreeS0R()
{
    tree->SetBranchStatus("R.s0.la"              ,1);  tree->SetBranchAddress("R.s0.la"              , R_s0_la             );
    tree->SetBranchStatus("R.s0.la_c"            ,1);  tree->SetBranchAddress("R.s0.la_c"            , R_s0_la_c           );
   tree->SetBranchStatus("R.s0.la_p"            ,1);  tree->SetBranchAddress("R.s0.la_p"            , R_s0_la_p           );
    tree->SetBranchStatus("R.s0.lt"              ,1);  tree->SetBranchAddress("R.s0.lt"              , R_s0_lt             );
     tree->SetBranchStatus("R.s0.lt_c"            ,1);  tree->SetBranchAddress("R.s0.lt_c"            , R_s0_lt_c           );
    tree->SetBranchStatus("R.s0.ra"              ,1);  tree->SetBranchAddress("R.s0.ra"              , R_s0_ra             );
    tree->SetBranchStatus("R.s0.ra_c"            ,1);  tree->SetBranchAddress("R.s0.ra_c"            , R_s0_ra_c           );
  tree->SetBranchStatus("R.s0.ra_p"            ,1);  tree->SetBranchAddress("R.s0.ra_p"            , R_s0_ra_p           );
    tree->SetBranchStatus("R.s0.rt"              ,1);  tree->SetBranchAddress("R.s0.rt"              , R_s0_rt             );
    tree->SetBranchStatus("R.s0.rt_c"            ,1);  tree->SetBranchAddress("R.s0.rt_c"            , R_s0_rt_c           );
    tree->SetBranchStatus("R.s0.t_pads"          ,1);  tree->SetBranchAddress("R.s0.t_pads"          , R_s0_t_pads         );
    tree->SetBranchStatus("R.s0.time"            ,1);  tree->SetBranchAddress("R.s0.time"            , R_s0_time           );
  tree->SetBranchStatus("R.s0.dedx"            ,1);  tree->SetBranchAddress("R.s0.dedx"            , R_s0_dedx           );
  tree->SetBranchStatus("R.s0.trdy"            ,1);  tree->SetBranchAddress("R.s0.trdy"            , R_s0_trdy           );
    tree->SetBranchStatus("R.s0.troff"           ,1);  tree->SetBranchAddress("R.s0.troff"           , R_s0_troff          );
    tree->SetBranchStatus("R.s0.trpad"           ,1);  tree->SetBranchAddress("R.s0.trpad"           , R_s0_trpad          );
  tree->SetBranchStatus("R.s0.trpath"          ,1);  tree->SetBranchAddress("R.s0.trpath"          , R_s0_trpath         );
  tree->SetBranchStatus("R.s0.trx"             ,1);  tree->SetBranchAddress("R.s0.trx"             , R_s0_trx            );
  tree->SetBranchStatus("R.s0.try"             ,1);  tree->SetBranchAddress("R.s0.try"             , R_s0_try            );
  tree->SetBranchStatus("R.s0.lnhits"          ,1);  tree->SetBranchAddress("R.s0.lnhits"          , R_s0_lnhits         );
    tree->SetBranchStatus("R.s0.lpeak"           ,1);  tree->SetBranchAddress("R.s0.lpeak"           , R_s0_lpeak          );
    tree->SetBranchStatus("R.s0.lt_fadc"         ,1);  tree->SetBranchAddress("R.s0.lt_fadc"         , R_s0_lt_fadc        );
    tree->SetBranchStatus("R.s0.ltc_fadc"        ,1);  tree->SetBranchAddress("R.s0.ltc_fadc"        , R_s0_ltc_fadc       );
    //  tree->SetBranchStatus("R.s0.lbadped"         ,1);  tree->SetBranchAddress("R.s0.lbadped"         , R_s0_lbadped        );
    //  tree->SetBranchStatus("R.s0.lunderflow"      ,1);  tree->SetBranchAddress("R.s0.lunderflow"      , R_s0_lunderflow     );
    //  tree->SetBranchStatus("R.s0.loverflow"       ,1);  tree->SetBranchAddress("R.s0.loverflow"       , R_s0_loverflow      );
  tree->SetBranchStatus("R.s0.rnhits"          ,1);  tree->SetBranchAddress("R.s0.rnhits"          , R_s0_rnhits         );
    tree->SetBranchStatus("R.s0.rpeak"           ,1);  tree->SetBranchAddress("R.s0.rpeak"           , R_s0_rpeak          );
    tree->SetBranchStatus("R.s0.rt_fadc"         ,1);  tree->SetBranchAddress("R.s0.rt_fadc"         , R_s0_rt_fadc        );
    tree->SetBranchStatus("R.s0.rtc_fadc"        ,1);  tree->SetBranchAddress("R.s0.rtc_fadc"        , R_s0_rtc_fadc       );
    //  tree->SetBranchStatus("R.s0.rbadped"         ,1);  tree->SetBranchAddress("R.s0.rbadped"         , R_s0_rbadped        );
  //  tree->SetBranchStatus("R.s0.runderflow"      ,1);  tree->SetBranchAddress("R.s0.runderflow"      , R_s0_runderflow     );
  //  tree->SetBranchStatus("R.s0.roverflow"       ,1);  tree->SetBranchAddress("R.s0.roverflow"       , R_s0_roverflow      );
}                                 
//////////////////////////////////////////////////
void Tree::readtreeS2R()
{
    tree->SetBranchStatus("R.s2.la"              ,1);  tree->SetBranchAddress("R.s2.la"              , R_s2_la             );
    tree->SetBranchStatus("R.s2.la_c"            ,1);  tree->SetBranchAddress("R.s2.la_c"            , R_s2_la_c           );
  tree->SetBranchStatus("R.s2.la_p"            ,1);  tree->SetBranchAddress("R.s2.la_p"            , R_s2_la_p           );
    tree->SetBranchStatus("R.s2.lt"              ,1);  tree->SetBranchAddress("R.s2.lt"              , R_s2_lt             );
    tree->SetBranchStatus("R.s2.lt_c"            ,1);  tree->SetBranchAddress("R.s2.lt_c"            , R_s2_lt_c           );
    tree->SetBranchStatus("R.s2.ra"              ,1);  tree->SetBranchAddress("R.s2.ra"              , R_s2_ra             );
    tree->SetBranchStatus("R.s2.ra_c"            ,1);  tree->SetBranchAddress("R.s2.ra_c"            , R_s2_ra_c           );
  tree->SetBranchStatus("R.s2.ra_p"            ,1);  tree->SetBranchAddress("R.s2.ra_p"            , R_s2_ra_p           );
    tree->SetBranchStatus("R.s2.rt"              ,1);  tree->SetBranchAddress("R.s2.rt"              , R_s2_rt             );
    tree->SetBranchStatus("R.s2.rt_c"            ,1);  tree->SetBranchAddress("R.s2.rt_c"            , R_s2_rt_c           );
    tree->SetBranchStatus("R.s2.t_pads"          ,1);  tree->SetBranchAddress("R.s2.t_pads"          , R_s2_t_pads         );
    tree->SetBranchStatus("R.s2.time"            ,1);  tree->SetBranchAddress("R.s2.time"            , R_s2_time           );
  tree->SetBranchStatus("R.s2.dedx"            ,1);  tree->SetBranchAddress("R.s2.dedx"            , R_s2_dedx           );
  tree->SetBranchStatus("R.s2.trdx"            ,1);  tree->SetBranchAddress("R.s2.trdx"            , R_s2_trdx           );
  //    tree->SetBranchStatus("R.s2.troff"           ,1);  tree->SetBranchAddress("R.s2.troff"           , R_s2_troff          );
    tree->SetBranchStatus("R.s2.trpad"           ,1);  tree->SetBranchAddress("R.s2.trpad"           , R_s2_trpad          );
  tree->SetBranchStatus("R.s2.trpath"          ,1);  tree->SetBranchAddress("R.s2.trpath"          , R_s2_trpath         );
  tree->SetBranchStatus("R.s2.trx"             ,1);  tree->SetBranchAddress("R.s2.trx"             , R_s2_trx            );
  tree->SetBranchStatus("R.s2.try"             ,1);  tree->SetBranchAddress("R.s2.try"             , R_s2_try            );
  tree->SetBranchStatus("R.s2.lnhits"          ,1);  tree->SetBranchAddress("R.s2.lnhits"          , R_s2_lnhits         );
  tree->SetBranchStatus("R.s2.lpeak"           ,1);  tree->SetBranchAddress("R.s2.lpeak"           , R_s2_lpeak          );
    tree->SetBranchStatus("R.s2.lt_fadc"         ,1);  tree->SetBranchAddress("R.s2.lt_fadc"         , R_s2_lt_fadc        );
    tree->SetBranchStatus("R.s2.ltc_fadc"        ,1);  tree->SetBranchAddress("R.s2.ltc_fadc"        , R_s2_ltc_fadc       );
    //  tree->SetBranchStatus("R.s2.lbadped"         ,1);  tree->SetBranchAddress("R.s2.lbadped"         , R_s2_lbadped        );
    //  tree->SetBranchStatus("R.s2.lunderflow"      ,1);  tree->SetBranchAddress("R.s2.lunderflow"      , R_s2_lunderflow     );
    //  tree->SetBranchStatus("R.s2.loverflow"       ,1);  tree->SetBranchAddress("R.s2.loverflow"       , R_s2_loverflow      );
  tree->SetBranchStatus("R.s2.rnhits"          ,1);  tree->SetBranchAddress("R.s2.rnhits"          , R_s2_rnhits         );
  tree->SetBranchStatus("R.s2.rpeak"           ,1);  tree->SetBranchAddress("R.s2.rpeak"           , R_s2_rpeak          );
  tree->SetBranchStatus("R.s2.rt_fadc"         ,1);  tree->SetBranchAddress("R.s2.rt_fadc"         , R_s2_rt_fadc        );
    tree->SetBranchStatus("R.s2.rtc_fadc"        ,1);  tree->SetBranchAddress("R.s2.rtc_fadc"        , R_s2_rtc_fadc       );
    //  tree->SetBranchStatus("R.s2.rbadped"         ,1);  tree->SetBranchAddress("R.s2.rbadped"         , R_s2_rbadped        );
    //  tree->SetBranchStatus("R.s2.runderflow"      ,1);  tree->SetBranchAddress("R.s2.runderflow"      , R_s2_runderflow     );
    //  tree->SetBranchStatus("R.s2.roverflow"       ,1);  tree->SetBranchAddress("R.s2.roverflow"       , R_s2_roverflow      );
    tree->SetBranchStatus("DR.rtRFtime"          ,1);  tree->SetBranchAddress("DR.rtRFtime"          , R_rtRFtime          );
}                                 
//////////////////////////////////////////////////
void Tree::readtreeA1R(){
  tree->SetBranchStatus("R.a1.a"               ,1);  tree->SetBranchAddress("R.a1.a"               , R_a1_a              ); 
  //  tree->SetBranchStatus("R.a1.a_c"             ,1);  tree->SetBranchAddress("R.a1.a_c"             , R_a1_a_c            ); 
  tree->SetBranchStatus("R.a1.a_p"             ,1);  tree->SetBranchAddress("R.a1.a_p"             , R_a1_a_p            ); 
  //  tree->SetBranchStatus("R.a1.t"               ,1);  tree->SetBranchAddress("R.a1.t"               , R_a1_t              ); 
  //  tree->SetBranchStatus("R.a1.t_c"             ,1);  tree->SetBranchAddress("R.a1.t_c"             , R_a1_t_c            ); 
  tree->SetBranchStatus("R.a1.asum_p"          ,1);  tree->SetBranchAddress("R.a1.asum_p"          ,&R_a1_asum_p         ); 
  tree->SetBranchStatus("R.a1.asum_c"          ,1);  tree->SetBranchAddress("R.a1.asum_c"          ,&R_a1_asum_c         ); 
  //  tree->SetBranchStatus("R.a1.trpath"          ,1);  tree->SetBranchAddress("R.a1.trpath"          , R_a1_trpath         ); 
  tree->SetBranchStatus("R.a1.trx"             ,1);  tree->SetBranchAddress("R.a1.trx"             , R_a1_trx            ); 
  tree->SetBranchStatus("R.a1.try"             ,1);  tree->SetBranchAddress("R.a1.try"             , R_a1_try            ); 
  tree->SetBranchStatus("R.a1.nhits"           ,1);  tree->SetBranchAddress("R.a1.nhits"           , R_a1_nhits          ); 
  //  tree->SetBranchStatus("R.a1.peak"            ,1);  tree->SetBranchAddress("R.a1.peak"            , R_a1_peak           ); 
  //  tree->SetBranchStatus("R.a1.t_fadc"          ,1);  tree->SetBranchAddress("R.a1.t_fadc"          , R_a1_t_fadc         ); 
  //  tree->SetBranchStatus("R.a1.tc_fadc"         ,1);  tree->SetBranchAddress("R.a1.tc_fadc"         , R_a1_tc_fadc        ); 
  //tree->SetBranchStatus("R.a1.nbadped"         ,1);  tree->SetBranchAddress("R.a1.nbadped"         , R_a1_nbadped        ); 
  //tree->SetBranchStatus("R.a1.noverflow"       ,1);  tree->SetBranchAddress("R.a1.noverflow"       , R_a1_noverflow      ); 
  //tree->SetBranchStatus("R.a1.nunderflow"      ,1);  tree->SetBranchAddress("R.a1.nunderflow"      , R_a1_nunderflow     ); 

}
//////////////////////////////////////////////////
void Tree::readtreeA2R(){
  tree->SetBranchStatus("R.a2.a"               ,1);  tree->SetBranchAddress("R.a2.a"               , R_a2_a              ); 
  //  tree->SetBranchStatus("R.a2.a_c"             ,1);  tree->SetBranchAddress("R.a2.a_c"             , R_a2_a_c            ); 
  tree->SetBranchStatus("R.a2.a_p"             ,1);  tree->SetBranchAddress("R.a2.a_p"             , R_a2_a_p            ); 
  //  tree->SetBranchStatus("R.a2.t"               ,1);  tree->SetBranchAddress("R.a2.t"               , R_a2_t              ); 
  //  tree->SetBranchStatus("R.a2.t_c"             ,1);  tree->SetBranchAddress("R.a2.t_c"             , R_a2_t_c            ); 
  tree->SetBranchStatus("R.a2.asum_p"          ,1);  tree->SetBranchAddress("R.a2.asum_p"          ,&R_a2_asum_p         ); 
  tree->SetBranchStatus("R.a2.asum_c"          ,1);  tree->SetBranchAddress("R.a2.asum_c"          ,&R_a2_asum_c         ); 
  //  tree->SetBranchStatus("R.a2.trpath"          ,1);  tree->SetBranchAddress("R.a2.trpath"          , R_a2_trpath         ); 
  tree->SetBranchStatus("R.a2.trx"             ,1);  tree->SetBranchAddress("R.a2.trx"             , R_a2_trx            ); 
  tree->SetBranchStatus("R.a2.try"             ,1);  tree->SetBranchAddress("R.a2.try"             , R_a2_try            ); 
  tree->SetBranchStatus("R.a2.nhits"           ,1);  tree->SetBranchAddress("R.a2.nhits"           , R_a2_nhits          ); 
  //  tree->SetBranchStatus("R.a2.peak"            ,1);  tree->SetBranchAddress("R.a2.peak"            , R_a2_peak           ); 
  //  tree->SetBranchStatus("R.a2.t_fadc"          ,1);  tree->SetBranchAddress("R.a2.t_fadc"          , R_a2_t_fadc         ); 
  //  tree->SetBranchStatus("R.a2.tc_fadc"         ,1);  tree->SetBranchAddress("R.a2.tc_fadc"         , R_a2_tc_fadc        ); 
  //tree->SetBranchStatus("R.a2.nbadped"         ,1);  tree->SetBranchAddress("R.a2.nbadped"         , R_a2_nbadped        ); 
  //tree->SetBranchStatus("R.a2.noverflow"       ,1);  tree->SetBranchAddress("R.a2.noverflow"       , R_a2_noverflow      ); 
  //tree->SetBranchStatus("R.a2.nunderflow"      ,1);  tree->SetBranchAddress("R.a2.nunderflow"      , R_a2_nunderflow     ); 

}
//////////////////////////////////////////////////


/////////////////////////////////////////////////

void Tree::readtreeGCR(){

  tree->SetBranchStatus("R.cer.asum_c"           ,1);  tree->SetBranchAddress("R.cer.asum_c"           ,& R_cer_asum_c          );

}
//////////////////////////////////////////////////
void Tree::readtreeShR(){
    tree->SetBranchStatus("R.sh.asum_c"           ,1);  tree->SetBranchAddress("R.sh.asum_c"           ,& R_sh_asum_c          ); 
}
//////////////////////////////////////////////////
void Tree::readtreePSR(){
     tree->SetBranchStatus("R.ps.asum_c"           ,1);  tree->SetBranchAddress("R.ps.asum_c"           ,& R_ps_asum_c          ); 
}

//////////////////////////////////////////////////
void Tree::readtreeGCL(){

  tree->SetBranchStatus("L.cer.asum_c"           ,1);  tree->SetBranchAddress("L.cer.asum_c"           ,& L_cer_asum_c          );

}
//////////////////////////////////////////////////
void Tree::readtreeShL(){
    tree->SetBranchStatus("L.sh.asum_c"           ,1);  tree->SetBranchAddress("L.sh.asum_c"           ,& L_sh_asum_c          ); 
}
//////////////////////////////////////////////////
void Tree::readtreePSL(){
    tree->SetBranchStatus("L.ps.asum_c"           ,1);  tree->SetBranchAddress("L.ps.asum_c"           ,& L_ps_asum_c          ); 
}



//////////////////////////////////////////////////
void Tree::readtreeTrackR(){
  //  tree->SetBranchStatus("R.tr.flag"            ,1);  tree->SetBranchAddress("R.tr.flag"            , R_tr_flag           );
  //  tree->SetBranchStatus("R.tr.ndof"            ,1);  tree->SetBranchAddress("R.tr.ndof"            , R_tr_ndof           );
  tree->SetBranchStatus("R.tr.chi2"            ,1);  tree->SetBranchAddress("R.tr.chi2"            , R_tr_chi2           );
  //  tree->SetBranchStatus("R.tr.beta"            ,1);  tree->SetBranchAddress("R.tr.beta"            , R_tr_beta           );
  //  tree->SetBranchStatus("R.tr.d_x"             ,1);  tree->SetBranchAddress("R.tr.d_x"             , R_tr_d_x            );
  //  tree->SetBranchStatus("R.tr.d_y"             ,1);  tree->SetBranchAddress("R.tr.d_y"             , R_tr_d_y            );
  //  tree->SetBranchStatus("R.tr.d_th"            ,1);  tree->SetBranchAddress("R.tr.d_th"            , R_tr_d_th           );
  //  tree->SetBranchStatus("R.tr.d_ph"            ,1);  tree->SetBranchAddress("R.tr.d_ph"            , R_tr_d_ph           );
  //  tree->SetBranchStatus("R.tr.r_x"             ,1);  tree->SetBranchAddress("R.tr.r_x"             , R_tr_r_x            );
  //  tree->SetBranchStatus("R.tr.r_y"             ,1);  tree->SetBranchAddress("R.tr.r_y"             , R_tr_r_y            );
  //  tree->SetBranchStatus("R.tr.r_th"            ,1);  tree->SetBranchAddress("R.tr.r_th"            , R_tr_r_th           );
  //  tree->SetBranchStatus("R.tr.r_ph"            ,1);  tree->SetBranchAddress("R.tr.r_ph"            , R_tr_r_ph           );
  tree->SetBranchStatus("R.tr.x"               ,1);  tree->SetBranchAddress("R.tr.x"               , R_tr_x              );
  tree->SetBranchStatus("R.tr.y"               ,1);  tree->SetBranchAddress("R.tr.y"               , R_tr_y              );
  tree->SetBranchStatus("R.tr.th"              ,1);  tree->SetBranchAddress("R.tr.th"              , R_tr_th             );
  tree->SetBranchStatus("R.tr.ph"              ,1);  tree->SetBranchAddress("R.tr.ph"              , R_tr_ph             );
  //  tree->SetBranchStatus("R.tr.time"            ,1);  tree->SetBranchAddress("R.tr.time"            , R_tr_time           );
  tree->SetBranchStatus("R.tr.p"               ,1);  tree->SetBranchAddress("R.tr.p"               , R_tr_p              );
  tree->SetBranchStatus("R.tr.pathl"           ,1);  tree->SetBranchAddress("R.tr.pathl"           , R_tr_pathl          );
  tree->SetBranchStatus("R.tr.px"              ,1);  tree->SetBranchAddress("R.tr.px"              , R_tr_px             );
  tree->SetBranchStatus("R.tr.py"              ,1);  tree->SetBranchAddress("R.tr.py"              , R_tr_py             );
  tree->SetBranchStatus("R.tr.pz"              ,1);  tree->SetBranchAddress("R.tr.pz"              , R_tr_pz             );
  tree->SetBranchStatus("R.tr.tg_dp"           ,1);  tree->SetBranchAddress("R.tr.tg_dp"           , R_tr_tg_dp          );
  tree->SetBranchStatus("R.tr.tg_y"            ,1);  tree->SetBranchAddress("R.tr.tg_y"            , R_tr_tg_y           );
  tree->SetBranchStatus("R.tr.tg_th"           ,1);  tree->SetBranchAddress("R.tr.tg_th"           , R_tr_tg_th          );
  tree->SetBranchStatus("R.tr.tg_ph"           ,1);  tree->SetBranchAddress("R.tr.tg_ph"           , R_tr_tg_ph          );
  tree->SetBranchStatus("R.tr.vx"              ,1);  tree->SetBranchAddress("R.tr.vx"              , R_tr_vx             );
  tree->SetBranchStatus("R.tr.vy"              ,1);  tree->SetBranchAddress("R.tr.vy"              , R_tr_vy             );
  tree->SetBranchStatus("R.tr.vz"              ,1);  tree->SetBranchAddress("R.tr.vz"              , R_tr_vz             );

}
//////////////////////////////////////////////////
void Tree::readtreeF1TDCR()
{
  tree->SetBranchStatus("RTDC.F1FirstHit"      ,1);  tree->SetBranchAddress("RTDC.F1FirstHit"      , RTDC_F1FirstHit);
}
//////////////////////////////////////////////////
void Tree::convertF1TDCR(ParamMan *param)
{
  
  


  
  bool F1Shift_R=true;
  //bool F1Shift_R=false;
  double RF1shift_off[RS2];

  //===== Reference check ======//
  RF1Ref[0]=RTDC_F1FirstHit[9];
  RF1Ref[1]=RTDC_F1FirstHit[46];
  
  if(RF1Ref[0]==0)RF1Ref[0]=RF1Ref[1];
  if(RF1Ref[1]==0)RF1Ref[1]=RF1Ref[0];

  //=============================//  
  
  //S2
  for(int i=0;i<RS2;i++){
    RS2T_F1TDC[i]  =  RTDC_F1FirstHit[i+16];
    RS2B_F1TDC[i]  =  RTDC_F1FirstHit[i+48];

    
    if(RF1Ref[0]>1. && RF1Ref[1]>1. && RS2T_F1TDC[i]>1. && RS2B_F1TDC[i]>1.){
      RS2T_F1[i] = RS2T_F1TDC[i] - RF1Ref[0] ;
      RS2B_F1[i] = RS2B_F1TDC[i] - RF1Ref[1] ;

      RS2T_F1time[i] = param->time(1,i,1,0, RS2T_F1[i]  , R_s2_ra_p[i]);
      RS2B_F1time[i] = param->time(1,i,1,1, RS2B_F1[i]  , R_s2_la_p[i]);
      RS2_F1time[i]  = 0.5*(RS2T_F1time[i] + RS2B_F1time[i]);      

      rtof[i]=RS2_F1time[i];

      //====== modified by itabashi ========//
      
      if(F1Shift_R){


	RF1shift_off[i]=param->GetF1Shift(i,1,0);
	if(RS2T_F1[i]>pow(2,15))RS2T_F1[i] -= RF1shift_off[i];
	if(RS2B_F1[i]>pow(2,15))RS2B_F1[i] -= RF1shift_off[i];

	RS2T_F1_c[i]=RS2T_F1[i];
	RS2B_F1_c[i]=RS2B_F1[i];

	
	RS2T_F1time[i] = param->time(1,i,1,0, RS2T_F1[i]  , R_s2_ra_p[i]);
	RS2B_F1time[i] = param->time(1,i,1,1, RS2B_F1[i]  , R_s2_la_p[i]);
	RS2_F1time_c[i]  = 0.5*(RS2T_F1time[i] + RS2B_F1time[i]);
	rtof_c[i]=RS2_F1time_c[i];
      }
      

      //=====================================//
    }
    else{
      //      RS2T_F1TDC[i]  = -9999.;
      //      RS2B_F1TDC[i]  = -9999.;
      RS2T_F1time[i] = -9999.;
      RS2B_F1time[i] = -9999.;
      RS2_F1time[i]  = -9999.;
    }
  }

  //S0
  for(int i=0;i<RS0;i++){
    RS0T_F1TDC[i]  =  RTDC_F1FirstHit[i+43];
    RS0B_F1TDC[i]  =  RTDC_F1FirstHit[i+44];
    RS0T_F1time[i] = param->time(2,i,1,0, RS0T_F1TDC[i] - RTDC_F1FirstHit[46] );
    RS0B_F1time[i] = param->time(2,i,1,1, RS0B_F1TDC[i] - RTDC_F1FirstHit[46] );
    RS0_F1time[i]  = 0.5*(RS0T_F1time[i] + RS0B_F1time[i]);
  }

  //RF
  RRF_F1time = param->time(3,0,1,0, RTDC_F1FirstHit[15] - RTDC_F1FirstHit[9]  );

}
//////////////////////////////////////////////////
double Tree::GetBeta_S0S2wF1TDCR(int itrack)
{
  double betaR = -999.;
  betaR = (R_s2_trpath[itrack] - R_s0_trpath[itrack])/(RS2_F1time[(int)R_s2_trpad[itrack]] - RS0_F1time[(int)R_s0_trpad[itrack]])/LightVelocity;

  return betaR;
}

void Tree::F1Shift(bool RHRS, bool F1TOP, int seg){

  //======= RHRS =================//

  RF1Ref[0]=RTDC_F1FirstHit[9] ;
  RF1Ref[1]=RTDC_F1FirstHit[46];
  RS2T_F1TDC[seg]  =  RTDC_F1FirstHit[seg+16];
  RS2B_F1TDC[seg]  =  RTDC_F1FirstHit[seg+48];

  
  
}

