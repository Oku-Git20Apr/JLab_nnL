// This macro can make histograms of all segments of AC1,2
// If you want to change its run number, you just have to change "runnum" variable and inside "tr->Add(...)".
// Only two things.
// 
// Their fit ranges have been roughly tuned by Okuyama (2020.3.21)
// (fit range width was fixed at 200[ch])
//


void gain_all_segment(){

	int runnum = 219;
	string ofname = Form("gain_run111%d.dat",runnum);
	string pdfname = Form("gain_run111%d.pdf",runnum);
cout << " output file name is " << ofname << endl;
cout << " output pdf file name is " << pdfname << endl;
  
	double Ra1a[100], Ra2a[100], Ra1a_c[100], Ra2a_c[100], Ra1a_p[100], Ra2a_p[100];
	double Ra1sum, Ra2sum;

	double adc1[100], adc2[100];

    TChain *tr = new TChain("T");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111157_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111157_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111157_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111157_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111157.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111158_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111158_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111158_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111158.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111159_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111159_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111159_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111159.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111160.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111161.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111162.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111163.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111164.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111165.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111166.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111167.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111168.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111169.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111170.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111171.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111172.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111173.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111174_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111174.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111175_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111175_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111175_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111175_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111175.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111176.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111177.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111178.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111179.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111180_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111180_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111180_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111180_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111180.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111181.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111182.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111183_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111183_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111183_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111183.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111184.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111185.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111186.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111187_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111187_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111187_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111187.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111188.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111189.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111190.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111191.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111192_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111194_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111194_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111194_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111194_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111194.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111195.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111196.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111197.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111198.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111199.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111200.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111201.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_7.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202_8.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111202.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111203.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111204.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111205.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111206.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207_6.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111207.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111208.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111209.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111210.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111211.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111212.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111213.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111214.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111215.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111216.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111217.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218_3.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218_4.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218_5.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111218.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219_1.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219_2.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219_3.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219_4.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219_5.root");
tr->Add("/data/11b/itabashi/root_ole/tritium_111219.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111220_1.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111220_2.root");
//tr->Add("/data/11b/itabashi/root_ole/tritium_111220.root");

	const int nofsegment1 = 24;//a1 segment number
	const int nofsegment2 = 26;//a2 segment number
	const int Nrun = 64;//Run111157 ~ Run111220
	double count1[nofsegment1], count2[nofsegment2];
	double count_err1[nofsegment1], count_err2[nofsegment2];
	double mean1[nofsegment1], mean2[nofsegment2];
	double mean_err1[nofsegment1], mean_err2[nofsegment2];
	double sigma1[nofsegment1], sigma2[nofsegment2];
	double sigma_err1[nofsegment1], sigma_err2[nofsegment2];

	double main_part1[nofsegment1], main_part2[nofsegment2];//1p.e. center position
	int index = runnum - 157;
	//a1 fit range shift
	double a1shift[Nrun][nofsegment1]={ 
	//   0  1  2    3  4  5  6  7  8  9   10 11 12 13 14 15 16 17 18 19 20 21 22 23
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111157
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111158
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111159
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111160
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111161
	//   0   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21  22 23
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111162
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111163
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111164
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111165
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111166
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111167
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111168
		{50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111169
		{50, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 50, 0, 0, 50, 50, 50, 0, 0, 0, 0, 50, 0},//Run111170
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111171
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111172
		{0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111173
		{50, 0, 50, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0},//Run111174
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111175
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111176
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111177
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111178
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111179
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111180
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111181
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111182
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111183
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111184
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111185
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111186
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111187
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111188
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111189
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111190
	//   0   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21  22 23
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111191
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 50, 0},//Run111192
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111193 No data
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111194
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 50},//Run111195
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111196
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111197
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 50, 0, 0, 0, 0, 50, 0, 0, 50, 0},//Run111198
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111199 No data
		{50, 50, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111200
	//   0   1  2  3  4  5  6  7  8  9 10 11  12 13 14 15 16 17 18 19 20 21  22 23
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111201
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111202
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111203
		{50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 50, 0},//Run111204
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111205
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 50, 0},//Run111206
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111207
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111208
		{50, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111209
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111210
		{50, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111211
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111212
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111213
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111214
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 50, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111215
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 50, 0},//Run111216
		{50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 50, 0, 0, 50, 0},//Run111217
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 50},//Run111218
		{50, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0},//Run111219
		{50, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 50, 0}};//Run111220
	//a2 fit range shift
	double a2shift[Nrun][nofsegment2]={
	//   0  1  2    3  4  5  6  7  8  9   10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
		{0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111157
		{0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111158
		{50, 0, 50, -50, 0, 0, 0, 0, 0, -50, -50, 0, 0, 0, 0, 0, -50, 0, 0, 50, -50, 0, 0, 0, 0, 100},//Run111159
		{50, 0, 0, -50, 0, 0, 0, 0, 0, -50, -100, 0, 0, 0, 0, 0, -50, 0, 0, 50, -50, 0, 0, 0, 0, 100},//Run111160
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, -100, 0, 0, 0, 0, 0, -50, 0, -100, 50, -50, 0, 0, 0, 0, 0},//Run111161
	//   0   1  2    3  4    5  6  7  8    9 10 11 12 13 14 15  16 17    18  19   20 21 22 23 24  25
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 50, 0, 0, 0, 0, 0, 50, 0, -100, 50, -50, 0, 0, 0, 0, 0},//Run111162
	//   0   1  2    3  4    5  6  7  8     9 10 11 12 13 14 15  16 17 18 19 20 21 22 23 24  25
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111163
		{-50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, -100},//Run111164
		{-50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -100},//Run111165
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111166
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111167
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -100},//Run111168
		{50, 0, 0, 0, 0, -50, 0, 0, 50, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111169
	//   0   1  2    3  4    5  6  7  8    9 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111170
		{50, 0, 0, -50, 0, 0, 0, 0, 0, -50, 0, 0, 0, 50, 0, 0, -50, 0, 0, 50, 0, 0, 0, 0, -50, 0},//Run111171
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111172
		{50, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111173
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, -50, 0, 0, 0, 0, 0, 50},//Run111174
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 50, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111175
		{50, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 50},//Run111176
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111177
		{0, 0, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111178
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, -50, 0},//Run111179
		{50, 0, 0, -50, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 50},//Run111180
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111181
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 50, 0, 0, 0, 0, 0, 0},//Run111182
		{-50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111183
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111184
		{0, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 50},//Run111185
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111186
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111187
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, -50, 0},//Run111188
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111189
		{0, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111190
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111191
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111192
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},         //Run111193 No data
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111194
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111195
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, -50},//Run111196
	//   0   1  2    3  4    5  6  7  8    9 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111197
		{50, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111198
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},         //Run111199 No data
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111200
		{50, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111201
		{0, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111202
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 50, 0, 0, 0, 0, -50, 0, 0, 50, 0, 0, 0, 0, 0, 0},//Run111203
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111204
		{50, 0, 0, 50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111205
		{50, 0, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111206
		{50, 0, 50, 50, 0, -50, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111207
		{50, 0, 0, 50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111208
		{50, 0, 0, 50, 0, -50, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111209
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 50, -50, 0, 0, 0, 0, 0, 0},//Run111210
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111211
		{50, 0, 0, 50, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111212
		{50, 0, 0, -50, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, -50, 0, 0, 50, 0, 0, 0},//Run111213
		{50, 0, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111214
		{50, 0, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111215
		{0, 0, 0, 0, 0, -50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111216
		{50, 0, 0, -50, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111217
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111218
		{50, 0, 0, -50, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0},//Run111219
		{50, 0, 0, -50, 0, -50, 0, 0, 0, -50, 0, 0, 50, 0, 0, 0, -50, 0, 0, 0, 0, 0, 0, 0, 0, 0}};//Run111220
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
	H1[i]->Fit(Form("f1[%d]",i),"","",main_part1[i]-100.+a1shift[index][i],main_part1[i]+100.+a1shift[index][i]);
	count1[i]	  = f1[i]->GetParameter(0);
	count_err1[i] = f1[i]->GetParError(0);
	mean1[i]	  = f1[i]->GetParameter(1);
	mean_err1[i]  = f1[i]->GetParError(1);
	sigma1[i]	  = f1[i]->GetParameter(2);
	sigma_err1[i] = f1[i]->GetParError(2);
	f1[i]->SetLineStyle(2);
	f1[i]->Draw("same");
	H1[i]->SetStats(0);
	cout << "ac1[" << i <<"]: " << count1[i] <<"/"<< count_err1[i] <<"/"<<mean1[i]<<"/"<<mean_err1[i]<<"/"<<sigma1[i]<<"/"<<sigma_err1[i]<<"/"<<endl;
	}

	c2->Divide(7,4);
	////////////////A2 Fit Range
	for(int i=0;i<nofsegment2;i++){
	c2->cd(i+1)->SetLogy(1);
	H2[i]->Draw("");
	H2[i]->Fit(Form("f2[%d]",i),"","",main_part2[i]-100+a2shift[index][i],main_part2[i]+100.+a2shift[index][i]);
	count2[i]	  = f2[i]->GetParameter(0);
	count_err2[i] = f2[i]->GetParError(0);
	mean2[i]	  = f2[i]->GetParameter(1);
	mean_err2[i]  = f2[i]->GetParError(1);
	sigma2[i]	  = f2[i]->GetParameter(2);
	sigma_err2[i] = f2[i]->GetParError(2);
	f2[i]->SetLineStyle(2);
	f2[i]->Draw("same");
	H2[i]->SetStats(0);
	cout << "ac2[" << i <<"]: " << count2[i] <<"/"<< count_err2[i] <<"/"<<mean2[i]<<"/"<<mean_err2[i]<<"/"<<sigma2[i]<<"/"<<sigma_err2[i]<<"/"<<endl;
	}


	ofstream writing_file;
	writing_file.open(ofname,ios::out);
	cout << "Fit data is filled in " << ofname << endl;

	for(int i=0;i<nofsegment1;i++){
	writing_file << "1 " <<i <<" "<< count1[i] <<" "<< count_err1[i] <<" "<<mean1[i]<<" "<<mean_err1[i]<<" "<<sigma1[i]<<" "<<sigma_err1[i]<<" "<<endl;
	}
	for(int i=0;i<nofsegment2;i++){
	writing_file << "2 " <<i << " "<< count2[i] <<" "<< count_err2[i] <<" "<<mean2[i]<<" "<<mean_err2[i]<<" "<<sigma2[i]<<" "<<sigma_err2[i]<<" "<<endl;
	}	

cout << "Print is starting" << endl;
	c1->Print(Form("%s[",pdfname.c_str()));
	c1->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s",pdfname.c_str()));
	c2->Print(Form("%s]",pdfname.c_str()));

	

 return 0;

}


