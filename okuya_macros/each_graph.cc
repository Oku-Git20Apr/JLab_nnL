void each_graph(){

	string pname[900];
	static const int total_run = 709;
	int n_run = 0;
	static const int nofsegment1 = 24;
	static const int nofsegment2 = 26;
//	int nofcanvas = 4;
	double gain1[nofsegment1][total_run];
	double err1[nofsegment1][total_run];
	double segment1[nofsegment1][total_run];
	double gain2[nofsegment2][total_run];
	double err2[nofsegment2][total_run];
	double segment2[nofsegment2][total_run];
	double null1[nofsegment2][total_run];
	double null2[nofsegment2][total_run];
	double peak1[nofsegment1][total_run];
	double peak2[nofsegment2][total_run];
	double rms1[nofsegment1][total_run];
	double rms2[nofsegment2][total_run];
	double runnum[total_run];

	int n;
	for(int i=132;i-132<total_run;i++){
	if(i==193 || i==199) continue;
	n = i - 132;
	pname[n] = Form("/data/41a/ELS/okuyama/nnL_output/ac_gain_20200516/gain_run111%d.dat",i); 
	//pname[n] = Form("./ac_gain_20200501/gain_run111%d.dat",i); 
	//pname[n] = Form("./ac_gain_20200406/gain_run111%d.dat",i); 
	//pname[n] = Form("./ac_gain_20200421/gain_run111%d.dat",i); 
	ifstream ifp(pname[n].c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname[n].c_str() << endl;

	string buf;
	int ac, seg;
	int n_seg1 = 0;
	int n_seg2 = 0;
	double err_val;
	double count, count_err, mean, mean_err, sigma, sigma_err, peak, rms;



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> ac >> seg >> count >> count_err >> mean >> mean_err >> sigma >> sigma_err >> peak >> rms; //0501,0516
		//sbuf >> ac >> seg >> count >> count_err >> mean >> mean_err >> sigma >> sigma_err; //0406
		//sbuf >> ac >> seg >> count >> count_err >> mean >> mean_err >> sigma >> sigma_err >> peak; //0421
		//cout << ac << seg << count << count_err << mean << mean_err << sigma << sigma_err << endl;

		if (ac == 1){
		//err_val  = mean_err;
		err_val  = sigma_err;
		//err_val  = sigma;
		//gain1[n_seg1][n_run] = mean;//y_val
		gain1[n_seg1][n_run] = sigma;//y_val
		err1[n_seg1][n_run] = err_val;//y_error
		segment1[n_seg1][n_run] = (double)seg;
		null1[n_seg1][n_run] = 0.;//x_error
		peak1[n_seg1][n_run] = peak;
		rms1[n_seg1][n_run] = rms;
		//if (rms == 0) gain1[n_seg1][n_run]=err1[n_seg1][n_run]=0.;
		n_seg1++;
		}
		else if (ac == 2){
		//err_val  = mean_err;
		err_val  = sigma_err;
		//err_val  = sigma;
		//gain2[n_seg2][n_run] = mean;//y_val
		gain2[n_seg2][n_run] = sigma;//y_val
		err2[n_seg2][n_run] = err_val;//y_error
		segment2[n_seg2][n_run] = (double)seg;
		null2[n_seg2][n_run] = 0.;//x_error
		peak2[n_seg2][n_run] = peak;
		rms2[n_seg2][n_run] = rms;
		//if (rms == 0) gain2[n_seg2][n_run]=err2[n_seg2][n_run]=0.;
		n_seg2++;
		}
	}
//	cout << "n_seg1= " << n_seg1 << endl;
//	cout << "n_seg2= " << n_seg2 << endl;

		ifp.close();

		runnum[n_run] = (double)i;//x_val
		n_run++;
	}

	cout << "n_run= " << n_run << endl;

		TGraphErrors *g1[nofsegment1], *g2[nofsegment2];
		for(int j=0;j<nofsegment1;j++){
		g1[j] = new TGraphErrors(n_run, runnum, gain1[j], null1[j], err1[j] );
		//g1[j] = new TGraphErrors(n_run, runnum, peak1[j], null1[j], null1[j] );
		g1[j]->SetMarkerStyle(21);
		g1[j]->SetMarkerColor(kAzure+j);
		g1[j]->SetMarkerSize(1.0);
		g1[j]->SetTitle(Form("a1[%d] gain",j));
//	cout << "g1 is created now:  " << j << endl;
		}
		for(int j=0;j<nofsegment2;j++){
		g2[j] = new TGraphErrors(n_run, runnum, gain2[j], null2[j], err2[j] );
		//g2[j] = new TGraphErrors(n_run, runnum, peak2[j], null2[j], null2[j] );
		g2[j]->SetMarkerStyle(21);
		g2[j]->SetMarkerColor(kAzure+j);
		g2[j]->SetMarkerSize(1.0);
		g2[j]->SetTitle(Form("a2[%d] gain",j));
//	cout << "g2 is created now:  " << j << endl;
		}

//	string paraname = "/data/41a/ELS/okuyama/Itabashi_20200310/ac/param/offset_ac.dat"; 
	string paraname = "/home/kazuki/Git_projects/JLab_nnL/ac/param/offset_ac.dat"; 
	cout << "gain.dat" << endl;
	ifstream ifp(paraname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << paraname.c_str() << endl;

	string pdfname = "each_graph.pdf~~~";
	cout << "output pdf file name is " << pdfname << endl;
	
	string buf;
	double line1[nofsegment1], line2[nofsegment2];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		int ac=0;
		int seg=0;
		double ped=0.;
		double ope=0.;
		double gain_val=0.;
		sbuf >> ac >> seg >> ped >> ope;
//		cout << ac << "/" << seg << "/" << ped << "/"  << ope <<  endl;
		gain_val = ope - ped;
		if(ac == 1){
		line1[seg] = gain_val;
//cout << "line1[" << seg << "] = " << line1[seg] << endl;
		}
		if(ac == 2){
		line2[seg] = gain_val;
//cout << "line2[" << seg << "] = " << line2[seg] << endl;
		}
	}
	

/*----Gain----------*/
	TCanvas *c[6];
	TLine *tl1[nofsegment1];
	for(int k=0;k<6;k++){
	c[k] =  new TCanvas(Form("c[%d]",k),Form("a1[%d]~a1[%d] gain",4*k,4*k+3),600.,600.);
	c[k]->Divide(2,2);
	for(int l=0;l<4;l++){
	int seg = 4*k+l;
	c[k]->cd(l+1)->SetGrid();
	if(seg==14){
	TH1 *frame = c[k]->cd(l+1)->DrawFrame(132.,gain1[seg][157]-500.,132.+(double)total_run,gain1[seg][157]+500.);
	frame->SetTitle(Form("a1[%d] gain",seg));
	}else{
	TH1 *frame = c[k]->cd(l+1)->DrawFrame(132.,gain1[seg][157]-100.,132.+(double)total_run,gain1[seg][157]+100.);
	frame->SetTitle(Form("a1[%d] gain",seg));
	}
	tl1[seg] = new TLine(132.,gain1[seg][157],132.+(double)total_run,gain1[seg][157]);
	tl1[seg]->SetLineWidth(3);
	tl1[seg]->SetLineColor(kBlack);
	tl1[seg]->Draw("same");
	g1[seg]->Draw("psame");
//cout << "ac1 " << seg << endl;
	}
}
	TCanvas *c2[7];
	TLine *tl2[nofsegment2];
	for(int k=0;k<7;k++){
	c2[k] =  new TCanvas(Form("c2[%d]",k),Form("a2[%d]~a2[%d] gain",4*k,4*k+3),600.,600.);
	c2[k]->Divide(2,2);
	for(int l=0;l<4;l++){
	int seg = 4*k+l;
	if(seg >= 26)break;
	c2[k]->cd(l+1)->SetGrid();
	if(seg==3 || seg==5){//range change
	TH1 *frame = c2[k]->cd(l+1)->DrawFrame(132.,gain2[seg][157]-500.,132.+(double)total_run,gain2[seg][157]+500.);
	frame->SetTitle(Form("a2[%d] gain",seg));
	}else{
	TH1 *frame = c2[k]->cd(l+1)->DrawFrame(132.,gain2[seg][157]-100.,132.+(double)total_run,gain2[seg][157]+100.);
	frame->SetTitle(Form("a2[%d] gain",seg));
	}
	tl2[seg] = new TLine(132.,gain2[seg][157],132.+(double)total_run,gain2[seg][157]);
	tl2[seg]->SetLineWidth(3);
	tl2[seg]->SetLineColor(kBlack);
	tl2[seg]->Draw("same");
	g2[seg]->Draw("psame");
//cout << "ac2 " << seg << endl;
	}
}
/*----Gain----------*/
//	TCanvas *c[6];
//	TLine *tl1[nofsegment1];
//	for(int k=0;k<6;k++){
//	c[k] =  new TCanvas(Form("c[%d]",k),Form("a1[%d]~a1[%d] gain",4*k,4*k+3),600.,600.);
//	c[k]->Divide(2,2);
//	for(int l=0;l<4;l++){
//	int seg = 4*k+l;
//	c[k]->cd(l+1)->SetGrid();
//	if(seg==14){
//	TH1 *frame = c[k]->cd(l+1)->DrawFrame(132.,line1[seg]-500.,132.+(double)total_run,line1[seg]+500.);
//	frame->SetTitle(Form("a1[%d] gain",seg));
//	}else{
//	TH1 *frame = c[k]->cd(l+1)->DrawFrame(132.,line1[seg]-100.,132.+(double)total_run,line1[seg]+100.);
//	frame->SetTitle(Form("a1[%d] gain",seg));
//	}
//	tl1[seg] = new TLine(132.,line1[seg],132.+(double)total_run,line1[seg]);
//	tl1[seg]->SetLineWidth(3);
//	tl1[seg]->SetLineColor(kBlack);
//	tl1[seg]->Draw("same");
//	g1[seg]->Draw("psame");
////cout << "ac1 " << seg << endl;
//	}
//}
//	TCanvas *c2[7];
//	TLine *tl2[nofsegment2];
//	for(int k=0;k<7;k++){
//	c2[k] =  new TCanvas(Form("c2[%d]",k),Form("a2[%d]~a2[%d] gain",4*k,4*k+3),600.,600.);
//	c2[k]->Divide(2,2);
//	for(int l=0;l<4;l++){
//	int seg = 4*k+l;
//	if(seg >= 26)break;
//	c2[k]->cd(l+1)->SetGrid();
//	if(seg==3 || seg==5){//range change
//	TH1 *frame = c2[k]->cd(l+1)->DrawFrame(132.,line2[seg]-500.,132.+(double)total_run,line2[seg]+500.);
//	frame->SetTitle(Form("a2[%d] gain",seg));
//	}else{
//	TH1 *frame = c2[k]->cd(l+1)->DrawFrame(132.,line2[seg]-100.,132.+(double)total_run,line2[seg]+100.);
//	frame->SetTitle(Form("a2[%d] gain",seg));
//	}
//	tl2[seg] = new TLine(132.,line2[seg],132.+(double)total_run,line2[seg]);
//	tl2[seg]->SetLineWidth(3);
//	tl2[seg]->SetLineColor(kBlack);
//	tl2[seg]->Draw("same");
//	g2[seg]->Draw("psame");
////cout << "ac2 " << seg << endl;
//	}
//}
cout << "Print is starting" << endl;
	c[0]->Print(Form("%s[",pdfname.c_str()));
	c[0]->Print(Form("%s",pdfname.c_str()));
	c[1]->Print(Form("%s",pdfname.c_str()));
	c[2]->Print(Form("%s",pdfname.c_str()));
	c[3]->Print(Form("%s",pdfname.c_str()));
	c[4]->Print(Form("%s",pdfname.c_str()));
	c[5]->Print(Form("%s",pdfname.c_str()));
	c2[0]->Print(Form("%s",pdfname.c_str()));
	c2[1]->Print(Form("%s",pdfname.c_str()));
	c2[2]->Print(Form("%s",pdfname.c_str()));
	c2[3]->Print(Form("%s",pdfname.c_str()));
	c2[4]->Print(Form("%s",pdfname.c_str()));
	c2[5]->Print(Form("%s",pdfname.c_str()));
	c2[6]->Print(Form("%s",pdfname.c_str()));
	c2[6]->Print(Form("%s]",pdfname.c_str()));

}
