void LE()
{
	TChain *tree = new TChain("nudir","nudir");
	tree->Add("out_1112.root");

	double nuE, vz, wgt;
	int nu_pdg;
	tree->SetBranchAddress("nuE", &nuE);
	tree->SetBranchAddress("nu_pdg", &nu_pdg);
	tree->SetBranchAddress("vz", &vz);
	tree->SetBranchAddress("wght", &wgt);

        const Int_t nbinsX = 29;
        Double_t xEdges[nbinsX+1];
  	xEdges[0]=0.;
  	for(int i=0; i<nbinsX+1; i++)
  	{
    	if(i<20)               xEdges[i+1] = xEdges[i] + 0.2;
    	else if(i>=20 && i<24) xEdges[i+1] = xEdges[i] + 1.0;
    	else if(i>=24 && i<28) xEdges[i+1] = xEdges[i] + 8.0;
    	//else if(i>=30 && i<34) xEdges[i+1] = xEdges[i] + 5.0;
    	else 		       xEdges[i+1] = xEdges[i] + 40.0;
  	}

	TH2D *h_e  = new TH2D("h_e","",nbinsX,xEdges,90,0.3,0.6);
	TH2D *hz_e = new TH2D("hz_e","",200,0,8,90,0.3,0.6);
	TH2D *h_m  = new TH2D("h_m","",nbinsX,xEdges,90,0.3,0.6);
	TH2D *hz_m = new TH2D("hz_m","",200,0,8,90,0.3,0.6);

	double nuE_max=0.;

	const int N = tree->GetEntries();
	for(int ii=0; ii<N; ii++){
		tree->GetEntry(ii);
                if(ii%100000 == 0) printf("%f percent of %d events...\n", ii*100./N, N);

		if(nuE>nuE_max) nuE_max = nuE;
		//hE->Fill(nuE);
		if(nu_pdg == 14 ) {
 			h_m->Fill(nuE,(57400-vz)/1E5,wgt);
 			hz_m->Fill(nuE,(57400-vz)/1E5,wgt);
		}
		if(nu_pdg == 12 ) {
 			h_e->Fill(nuE,(57400-vz)/1E5,wgt);
 			hz_e->Fill(nuE,(57400-vz)/1E5,wgt);
		}
	}

	//std::cout << nuE_max << "\n";

	gStyle->SetPalette(kColorPrintableOnGrey);
	TColor::InvertPalette();
	gStyle->SetNumberContours(999);

	h_e->SetTitle("L-E distribution (nue)");
	h_e->GetXaxis()->SetTitle("E_nu (GeV)");
	h_e->GetYaxis()->SetTitle("L (km)");
	hz_e->SetTitle("L-E distribution (nue)");
	hz_e->GetXaxis()->SetTitle("E_nu (GeV)");
	hz_e->GetYaxis()->SetTitle("L (km)");

	h_m->SetTitle("L-E distribution (numu)");
	h_m->GetXaxis()->SetTitle("E_nu (GeV)");
	h_m->GetYaxis()->SetTitle("L (km)");
	hz_m->SetTitle("L-E distribution (numu)");
	hz_m->GetXaxis()->SetTitle("E_nu (GeV)");
	hz_m->GetYaxis()->SetTitle("L (km)");

	TCanvas *c = new TCanvas("c","",800,600);
	h_e->Draw("colz");
	c->SaveAs("LE_e.png");
	hz_e->Draw("colz");
	c->SaveAs("LEz_e.png");

	h_m->Draw("colz");
	c->SaveAs("LE_m.png");
	hz_m->Draw("colz");
	c->SaveAs("LEz_m.png");


	TCanvas *c_log = new TCanvas("c_log","",800,600);
	c_log->SetLogx();
	h_e->Draw("colz");
	c_log->SaveAs("LE_e_log.png");

	h_m->Draw("colz");
	c_log->SaveAs("LE_m_log.png");


        TFile *f = new TFile("LE_1112_2.root","RECREATE");
        h_e->Write();
        h_m->Write();
        f->Close();

        //for(int i=410;i<430;i++){
        //std:cout << i << "\t" << h->ProjectionY("",i+1,i+1)->GetRandom() << "\n";}
}
