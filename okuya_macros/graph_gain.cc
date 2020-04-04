void graph_gain(){

	string pname[100];
	int total_run = 64;
	int n_run = 0;
	int nofsegment1 = 24;
	int nofsegment2 = 26;
	int nofcanvas = 4;

	double gain1[nofsegment1][total_run], err1[nofsegment1][total_run], segment1[nofsegment1][total_run];
	double gain2[nofsegment2][total_run], err2[nofsegment2][total_run], segment2[nofsegment2][total_run];
	double null1[nofsegment2][total_run], null2[nofsegment2][total_run], runnum[total_run];

	int n;
	for(int i=157;i-157<total_run;i++){
	if(i==193 || i==199) continue;
	n = i - 157;
	pname[n] = Form("./gain_run111%d.dat",i); 
	ifstream ifp(pname[n].c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname[n].c_str() << endl;

	string buf;
	int ac, seg;
	int n_seg1 = 0;
	int n_seg2 = 0;
	double err_val;
	double count, count_err, mean, mean_err, sigma, sigma_err;



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> ac >> seg >> count >> count_err >> mean >> mean_err >> sigma >> sigma_err;
		//cout << ac << seg << count << count_err << mean << mean_err << sigma << sigma_err << endl;

		if (ac == 1){
		err_val  = mean_err;
		//err_val  = sigma;
		gain1[n_seg1][n_run] = mean;//y_val
		err1[n_seg1][n_run] = err_val;//y_error
		segment1[n_seg1][n_run] = (double)seg;
		null1[n_seg1][n_run] = 0.;//x_error
		n_seg1++;
		}
		else if (ac == 2){
		err_val  = mean_err;
		//err_val  = sigma;
		gain2[n_seg2][n_run] = mean;//y_val
		err2[n_seg2][n_run] = err_val;//y_error
		segment2[n_seg2][n_run] = (double)seg;
		null2[n_seg2][n_run] = 0.;//x_error
		n_seg2++;
		}
	}

		ifp.close();

		runnum[n_run] = (double)i;//x_val
		n_run++;
	}

	TF1 *acg1[nofsegment1], *acg2[nofsegment2];

		TGraphErrors *g1[nofsegment1], *g2[nofsegment2];
		for(int j=0;j<nofsegment1;j++){
		acg1[j] = new TF1(Form("acg1[%d]",j),"",157.,157+(double)total_run);
		g1[j] = new TGraphErrors(n_run, runnum, gain1[j], null1[j], err1[j] );
		g1[j]->SetMarkerStyle(21);
		g1[j]->SetMarkerColor(kAzure+j);
		g1[j]->SetMarkerSize(1.0);
		g1[j]->SetTitle(Form("a1[%d] gain",j));
		g1[j]->Fit(acg1[j],"",157.,157+(double)total_run);
		}
		for(int j=0;j<nofsegment2;j++){
		g2[j] = new TGraphErrors(n_run, runnum, gain2[j], null2[j], err2[j] );
		g2[j]->SetMarkerStyle(21);
		g2[j]->SetMarkerColor(kAzure+j);
		g2[j]->SetMarkerSize(1.0);
		g2[j]->SetTitle(Form("a2[%d] gain",j));
		}

	string paraname = "/data/40a/okuyama/Itabashi_20200310/ac/param/offset_ac.dat"; 
	cout << "gain.dat" << endl;
	ifstream ifp(paraname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << paraname.c_str() << endl;

	string buf;
	int ac, seg;
	double ped, ope;
	double gain_val;
	double line1[nofsegment1], line2[nofsegment2];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> ac >> seg >> ped >> ope;
//		cout << ac << "/" << seg << "/" << ped << "/" << ped_er << "/" << ope << "/" << ope_er << endl;

		gain_val = ope - ped;
		if(ac == 1){
		line1[seg] = gain_val;
cout << "line1[" << seg << "] = " << line1[seg] << endl;
		}
		if(ac == 2){
		line2[seg] = gain_val;
cout << "line2[" << seg << "] = " << line2[seg] << endl;
		}
	}
	

	TCanvas *c[nofcanvas];
	for(int k=0;k<nofcanvas;k++){
	c[k] =  new TCanvas(Form("c[%d]",k),"a1 gain",600.,600.);
//	c[k]->Divide(2,2);
//	for(int l=0;l<4;l++){
//	c[k]->cd(l+1)->SetGrid();
	TH1 *frame = c[k]->DrawFrame(157.,100.,157.+(double)total_run,550.);
	TLine *tl1 = new TLine(157.,line1[6*k+0],157.+(double)total_run,line1[6*k+0]);
	tl1->SetLineWidth(5);
	tl1->SetLineColor(kAzure);
	//tl1->Draw("same");
	g1[6*k+0]->Draw("plsame");
	acg1[6*k+0]->Draw("plsame");
	g1[6*k+1]->Draw("plsame");
	g1[6*k+2]->Draw("plsame");
	g1[6*k+3]->Draw("plsame");
	g1[6*k+4]->Draw("plsame");
	g1[6*k+5]->Draw("plsame");
//	}
}
	TCanvas *c2[nofcanvas];
	for(int k=0;k<nofcanvas;k++){
	c2[k] =  new TCanvas(Form("c2[%d]",k),"a2 gain",600.,600.);
//	c2[k]->Divide(2,2);
//	for(int l=0;l<4;l++){
//	c2[k]->cd(l+1)->SetGrid();
	TH1 *frame = c2[k]->DrawFrame(157.,100.,157.+(double)total_run,550.);
	g2[7*k+0]->Draw("plsame");
	g2[7*k+1]->Draw("plsame");
	g2[7*k+2]->Draw("plsame");
	g2[7*k+3]->Draw("plsame");
	g2[7*k+4]->Draw("plsame");
	if(k==3)continue;
	g2[7*k+5]->Draw("plsame");
	g2[7*k+6]->Draw("plsame");
//	}
}

}
