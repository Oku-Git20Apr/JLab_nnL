// This macro can make histograms of all segments of AC1,2
// If you want to change input files, you just have to change "Nrun".
// However, it starts to count from "runnum".
// 
// -------a1shift & a2shift are no longer used.-------
// Their fit ranges have been roughly tuned by Okuyama (2020.3.21)
// (fit range width was fixed at 200[ch])
// -------a1shift & a2shift are no longer used.-------
// 
//
//
//-------Modification Log----------
//
//--------------------------------
//---------Kazuki Okuyama---------
//--------------------------------
//
// I modified the fitting method. (2020.4.5)
// The center value wes changed to be obtained from first fit result.
// (fit range width is still fixed at 200[ch])
//
// I modified again. (2020.4.21)
// Fitting range became much larger during a first fit.
// If it is a wrong fitting, a maximum value is used as a gain.
// a1shift & a2shift are deleted.
// Using small root
//
// I added that we can monitor RMS value. (2020.4.30)
// Run111132~Run111156 were added.
// Using root/, not small/
// 2ch/bin -> 6ch/bin
// but this didn't work well. (2020.5.4)
//
// 6ch/bin -> 2ch/bin
// Using root/ (2020.5.5) 

void gain_all_segment(){


//	const int Nrun = 64;//Run111157 ~ Run111220
	const int Nrun = 115;

	for(int loop=0;loop<Nrun;loop++){
	//int runnum = 157;
	int runnum = 738;
		runnum += loop;


	string ofname = Form("gain_run111%d.dat",runnum);
	string pdfname = Form("gain_run111%d.pdf",runnum);
cout << " output file name is " << ofname << endl;
cout << " output pdf file name is " << pdfname << endl;
  
	double Ra1a[100], Ra2a[100], Ra1a_c[100], Ra2a_c[100], Ra1a_p[100], Ra2a_p[100];
	double Ra1sum, Ra2sum;

	double adc1[100], adc2[100];

    TChain *tr = new TChain("T");

	string ifname = "runlist.all";//Run List including all run
    ifstream ifp(Form("%s",ifname.c_str()),ios::in);
    if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
    string buf, runname;
    while(1){
      getline(ifp,buf);
      if( buf[0]=='#' ){ continue; }
      if( ifp.eof() ) break;
	  if( buf.find(Form("tritium_111%d",runnum))!= -1){
      istringstream sbuf(buf);
      sbuf >> runname;
	  tr->Add(runname.c_str());
      cout<<buf<<endl;
      }
    }
	int ENum;
    ENum=tr->GetEntries();
    cout<<"Events: "<<ENum<<endl; 

	const int nofsegment1 = 24;//a1 segment number
	const int nofsegment2 = 26;//a2 segment number
	double count1[nofsegment1], count2[nofsegment2];
	double count_err1[nofsegment1], count_err2[nofsegment2];
	double mean1[nofsegment1], mean2[nofsegment2];
	double mean_err1[nofsegment1], mean_err2[nofsegment2];
	double sigma1[nofsegment1], sigma2[nofsegment2];
	double sigma_err1[nofsegment1], sigma_err2[nofsegment2];
	double RMS1[nofsegment1], RMS2[nofsegment2];

	double main_part1[nofsegment1], main_part2[nofsegment2];//1p.e. center position
	int index = runnum - 157;
	TH1D *H1_temp[nofsegment1], *H2_temp[nofsegment2];
	TH1D *H1[nofsegment1], *H2[nofsegment2];
	TF1 *f1[nofsegment1], *f2[nofsegment2];


///////////////////////////////////////////////////
///////////////////////////////////////////////////
///I've really wanted to know the maximum bin!!////
///////////////////////////////////////////////////
	for(int i=0;i<nofsegment1;i++){
	H1_temp[i] = new TH1D(Form("H1_temp[%d]",i),Form("R.a1.a_p[%d]",i),999,2.,2000.);
	tr->Project(Form("H1_temp[%d]",i),Form("R.a1.a_p[%d]",i),Form("R.a1.a_p[%d]>200.",i),"");
	main_part1[i] = (H1_temp[i]->GetBinCenter(H1_temp[i]->GetMaximumBin()));
cout << "a1_center[" << i << "] = " << main_part1[i] << endl;
	}
	for(int i=0;i<nofsegment2;i++){
	H2_temp[i] = new TH1D(Form("H2_temp[%d]",i),Form("R.a2.a_p[%d]",i),999,2.,2000.);
	tr->Project(Form("H2_temp[%d]",i),Form("R.a2.a_p[%d]",i),Form("R.a2.a_p[%d]>100.",i),"");
	main_part2[i] = (H2_temp[i]->GetBinCenter(H2_temp[i]->GetMaximumBin()));
cout << "a2_center[" << i << "] = " << main_part2[i] << endl;
	}
///////////////////////////////////////////////////
///////////////////////////////////////////////////
	

	TCanvas *c1 = new TCanvas("c1","A1",800.,600.);
	TCanvas *c2 = new TCanvas("c2","A2",800.,600.);


	for(int i=0;i<nofsegment1;i++){
	H1[i] = new TH1D(Form("H1[%d]",i),Form("R.a1.a_p[%d]",i),999,2.,2000.);
	tr->Project(Form("H1[%d]",i),Form("R.a1.a_p[%d]",i),"","");
	f1[i] = new TF1(Form("f1[%d]",i),"gausn",0.,2000.);
	f1[i]->SetParameter(1,main_part1[i]);
	}
	for(int i=0;i<nofsegment2;i++){
	H2[i] = new TH1D(Form("H2[%d]",i),Form("R.a2.a_p[%d]",i),999,2.,2000.);
	tr->Project(Form("H2[%d]",i),Form("R.a2.a_p[%d]",i),"","");
	f2[i] = new TF1(Form("f2[%d]",i),"gausn",0.,2000.);
	f2[i]->SetParameter(1,main_part2[i]);
	}
	

	c1->Divide(6,4);
	////////////////A1 Fit Range
	for(int i=0;i<nofsegment1;i++){
	c1->cd(i+1)->SetLogy(1);
	H1[i]->Draw("");
	H1[i]->Fit(Form("f1[%d]",i),"","",main_part1[i]-300.,main_part1[i]+300.);//First fit (width:600ch)
	mean1[i]	  = f1[i]->GetParameter(1);
	RMS1[i]=H1[i]->GetRMS();
	H1[i]->Fit(Form("f1[%d]",i),"","",mean1[i]-100.,mean1[i]+100.);//fitting again (Width is fixed at 200ch)
	count1[i]	  = f1[i]->GetParameter(0);
	count_err1[i] = f1[i]->GetParError(0);
	mean1[i]	  = f1[i]->GetParameter(1);
	mean_err1[i]  = f1[i]->GetParError(1);
	sigma1[i]	  = f1[i]->GetParameter(2);
	sigma_err1[i] = f1[i]->GetParError(2);
	f1[i]->SetLineStyle(2);
	f1[i]->Draw("same");
	H1[i]->SetStats(0);
	cout << "ac1[" << i <<"]: " << count1[i] <<"/"<< count_err1[i] <<"/"<<mean1[i]<<"/"<<mean_err1[i]<<"/"<<sigma1[i]<<"/"<<sigma_err1[i]<<"/"<<main_part1[i]<<"/"<<RMS1[i]<<"/"<<endl;
	}

	c2->Divide(7,4);
	////////////////A2 Fit Range
	for(int i=0;i<nofsegment2;i++){
	c2->cd(i+1)->SetLogy(1);
	H2[i]->Draw("");
	H2[i]->Fit(Form("f2[%d]",i),"","",main_part2[i]-300.,main_part2[i]+300.);
	mean2[i]	  = f2[i]->GetParameter(1);
	RMS2[i]=H2[i]->GetRMS();
	H2[i]->Fit(Form("f2[%d]",i),"","",mean2[i]-100.,mean2[i]+100.);
	count2[i]	  = f2[i]->GetParameter(0);
	count_err2[i] = f2[i]->GetParError(0);
	mean2[i]	  = f2[i]->GetParameter(1);
	mean_err2[i]  = f2[i]->GetParError(1);
	sigma2[i]	  = f2[i]->GetParameter(2);
	sigma_err2[i] = f2[i]->GetParError(2);
	f2[i]->SetLineStyle(2);
	f2[i]->Draw("same");
	H2[i]->SetStats(0);
	cout << "ac2[" << i <<"]: " << count2[i] <<"/"<< count_err2[i] <<"/"<<mean2[i]<<"/"<<mean_err2[i]<<"/"<<sigma2[i]<<"/"<<sigma_err2[i]<<"/"<<main_part2[i]<<"/"<<RMS2[i]<<"/"<<endl;
	}


//	ofstream writing_file;
//	writing_file.open(ofname,ios::out);
	ofstream fout(Form("gain_run111%d.dat",runnum));
	cout << "Fit data is filled in " << ofname << endl;
		for(int i=0;i<nofsegment1;i++){ 
	fout << "1 " <<i <<" "<< count1[i] <<" "<< count_err1[i] <<" "<<mean1[i]<<" "<<mean_err1[i]<<" "<<sigma1[i]<<" "<<sigma_err1[i]<<" "<<main_part1[i]<<" "<<RMS1[i]<<" "<<endl;
		} 
		for(int i=0;i<nofsegment2;i++){ 
	fout << "2 " <<i << " "<< count2[i] <<" "<< count_err2[i] <<" "<<mean2[i]<<" "<<mean_err2[i]<<" "<<sigma2[i]<<" "<<sigma_err2[i]<<" "<<main_part2[i]<<" "<<RMS2[i]<<" "<<endl;
		}	

cout << "Print is starting" << endl;
	c1->Print(Form("%s[",pdfname.c_str()));
	c1->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s]",pdfname.c_str()));

	c1->Close();
	c2->Close();
	for(int i=0;i<nofsegment1;i++){
delete gROOT->Get(Form("H1[%d]",i));
delete gROOT->Get(Form("H1_temp[%d]",i));
	}
	for(int i=0;i<nofsegment2;i++){
delete gROOT->Get(Form("H2[%d]",i));
delete gROOT->Get(Form("H2_temp[%d]",i));
	}
	
	}//loop	

 return 0;

}


