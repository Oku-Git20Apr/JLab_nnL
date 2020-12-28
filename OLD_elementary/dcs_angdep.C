{
//   double W = 1670; // MeV
//   double Q2 = 0.050e6; // (MeV/c)^2
//   double epsilon = 0.54;
//   double phi = 0; // degree

  //double W = 1750; // MeV
  //double W = 1835; // MeV
  double W = 2130; // MeV
  //double W = 2240; // MeV
  //double Q2 = 0.036e6; // (MeV/c)^2
  //double Q2 = 0.; // (MeV/c)^2
  double Q2 = 0.5e6; // (MeV/c)^2
  //double epsilon = 0.395;
  double epsilon = 0.7;
  double phi = 0; // degree
  int isospin = 1;

  TStrangeModel model("test_pion","pionHighE");
  
//RPR3
  //model.SetStrangeModel("/usr/local/apps/strangecalc/share/strangecalc-wrapper/models/rpr-2007/iso1+4/init/fit_specification");
//RPR2011
  model.SetStrangeModel("/usr/local/apps/strangecalc/share/strangecalc-wrapper/models/rpr-2011/iso1+4/init/fit_specification");

  TKinematics kin("kin","example kinematics",isospin,"w:costhkcm:qsquared",W,0.25,Q2);
  TCalcInfo cs_l(TCalcInfo::kElectro,isospin,"diff_l"); cs_l->SetElectroEpsilon(epsilon);
  TCalcInfo cs_t(TCalcInfo::kElectro,isospin,"diff_t"); cs_t->SetElectroEpsilon(epsilon);
  TCalcInfo cs_tl(TCalcInfo::kElectro,isospin,"diff_t+l"); cs_tl->SetElectroEpsilon(epsilon);
  TCalcInfo cs_tt_unpol(TCalcInfo::kElectro,isospin,"diff_tt_unpol"); cs_tt_unpol->SetElectroEpsilon(epsilon);
  TCalcInfo cs_tl_unpol(TCalcInfo::kElectro,isospin,"diff_tl_unpol"); cs_tl_unpol->SetElectroEpsilon(epsilon);

  std::cout << " epsilon = " << cs_l->GetElectroEpsilon() << std::endl;
  std::cout << model.GetCalcpoint(kin,&cs_l) << " "
	    << model.GetCalcpoint(kin,&cs_t) << " "
	    << model.GetCalcpoint(kin,&cs_tl) << " "
	    << std::endl;

  ofstream output;
  output.open("dcs.dat");

  int N = 37;
  double costcm, costcmmin=-1, costcmmax=1;
  double tcm, tcmmin=0.,tcmmax=4.*atan(1.);
  double csl,cst,cstl,cstt_unpol,cstl_unpol;
  double cs;
  for(int i = 0; i <= N; i++)
    {
      //costcm = costcmmin + (costcmmax-costcmmin)*double(i)/double(N);
      tcm = tcmmin + (tcmmax-tcmmin)*double(i)/double(N);
      
      //kin.SetVar(2,costcm);
      kin.SetVar(2,cos(tcm));
      // unit : ub/sr
      csl = model.GetCalcpoint(kin,&cs_l);
      cst = model.GetCalcpoint(kin,&cs_t);
      cstl = model.GetCalcpoint(kin,&cs_tl);
      cstt_unpol = model.GetCalcpoint(kin,&cs_tt_unpol);
      cstl_unpol = model.GetCalcpoint(kin,&cs_tl_unpol);
      cs = cstl + sqrt(2*epsilon*(1+epsilon))*cstl_unpol*cos(phi*TMath::DegToRad()) + epsilon*cstt_unpol*cos(2*phi*TMath::DegToRad());
      output //<< cos(tcm) << "\t\t\t\t\t\t"
	//     << csl << "\t"
	//     << cst << "\t"
	//     << cstl << "\t"
	//     << cstt_unpol << "\t"
	//     << cstl_unpol << "\t"
	     << cs //<< "\t"
	     << endl;
    }

  output.close();
}
