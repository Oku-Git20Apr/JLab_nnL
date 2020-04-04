#ifndef Param_h
#define Param_h 1


// =================================================== //
// ==== Offset and scaling factors for matrices ====== //
// =================================================== //
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  Ztm = -0.15,Ztr=0.35;
//==== momentum scaled  parameters =====//
const double  PLm = 2.0, PLr=0.22; 
const double  PRm =1.74   ,PRr=0.2;


//const double  Ztm = 0.0,Ztr=1.0; //NO scale  


const int nn = 4; // 4th order matrix using xf, xpf, yf, ypf, and zt
const int nParamT = 126;  // Number of parameters
const int nnz =3;
const int nParamTz=35;
const int nParamT_ras=4;

//const int nnp=4;// 4th order matrix using xf, xpf, yf, ypf, and zt
//const int nParamTp=126;
const int nnp=3;// 3th order matrix using xf, xpf, yf, ypf, and zt
const int nParamTp=56;
//const int nnp=2;// 2nd order matrix using xf, xpf, yf, ypf, and zt
//const int nParamTp=21;
///==== MTtuning =====//
double Pras[nParamT_ras];
double Pras_L[nParamT_ras];
double Pzt[nParamTz],Pzt_L[nParamTz];
double Pxpt[nParamT],Pxpt_L[nParamT];
double Pypt[nParamT],Pypt_L[nParamT];
double Prp[nParamTp],Plp[nParamTp];
double Opt_par_R[nParamTp],Opt_par_L[nParamTp];
double Opt_par[nParamTp*2+nParamT*2];

//const int nmax = 1000; // Number of events used for tuning
const int nmax=10000;// Number of events used for tuning;
//const int nmax = 50000; // Number of events used for tuning
//const int nite =0;
const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
double selection_width = 0.0125; 
//double selection_width = 0.008; // event selection width for z

 double RasterCor, RasterCor_L;
const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double refx[nsshole];
double refy[nsshole];
double selec_widthx = 0.60; // selection width in x (dispersive plane)
double selec_widthy = 0.45; // selection width in y 
int ntune_event;





double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double z_recon[nmax]; // reconstructed z position
double Rras_x[nmax],Lras_x[nmax];
int foil_flag[nmax];  
int holegroup[nmax];
int mass_flag[nmax];
double mass_ref[nmax];
double MM[nmax],rx_fp[nmax],rth_fp[nmax],ry_fp[nmax],rph_fp[nmax],lx_fp[nmax],lth_fp[nmax],ly_fp[nmax],lph_fp[nmax];
double rx[nmax],ry[nmax],rth[nmax],rph[nmax],rz[nmax],lx[nmax],ly[nmax],lth[nmax],lph[nmax],lz[nmax];
double beam_p[nmax];
double bp[nmax],rp[nmax],lp[nmax];
double dth[nfoil];
double l[nfoil];
double projectf[nfoil];
double OptPar1[nParamT];
double OptPar2[nParamT];


//const int nParamT2 = 4; 
//double parRaster[nParamT2];
//double Opt_Par[nParamT2];
//double Ras_curx[nmax];
//double Ras_cury[nmax];

const double hrs_ang = 13.2 * 3.14159 / 180.;
//const double tdc_time =56.0e-3;//[ns]

//======= Coincidence Offset Paramters ===========//

const double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
const double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};

const double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
const double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};


const  double pathl_off=-470.5;
const  double s2_offset=-499.75;
//const  double coin_offset=456.65;
const  double coin_offset=464.65;
#endif
