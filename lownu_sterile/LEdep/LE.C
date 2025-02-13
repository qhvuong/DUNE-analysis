#include <DUNEStyle.h>
#include <TLegend.h>

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void plot1D(TCanvas *c, TH1D *h, const char *name)
{
  c->Clear();
  c->cd();
  h->Draw();
  h->GetYaxis()->SetRangeUser(h->GetYaxis()->GetXmin(), h->GetMaximum()*1.25);  // make room for watermark
  dunestyle::CenterTitles(h);
  dunestyle::Simulation();
  c->SaveAs(Form("%s.png",name));
}

void plot2D(TCanvas *c, TH2D *h, const char *name)
{
  c->Clear();
  c->cd();
  //h->GetYaxis()->SetRangeUser(h->GetYaxis()->GetXmin(), h->GetMaximum()*2.5);  // make room for watermark
  dunestyle::CenterTitles(h);
  h->Draw("colz");
  dunestyle::Simulation();
  c->SaveAs(Form("%s.png",name));
}


void LE()
{
	TChain *tree = new TChain("nudir","nudir");
	tree->Add("out_1007.root");

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

	//TH2D *h_e  = new TH2D("h_e","",nbinsX,xEdges,90,0.3,0.6);
	TH2D *h_e  = new TH2D("h_e","",200,0,80,80,0.25,0.65);
	TH2D *hz_e = new TH2D("hz_e","",200,0,16,80,0.25,0.65);
	//TH2D *h_m  = new TH2D("h_m","",nbinsX,xEdges,90,0.3,0.6);
	TH2D *h_m  = new TH2D("h_m","",200,0,80,80,0.25,0.65);
	TH2D *hz_m = new TH2D("hz_m","",200,0,16,80,0.25,0.65);

	double nuE_max=0.;

	const int N = tree->GetEntries();
	//const int N = 1000000;
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

	//gStyle->SetPalette(kColorPrintableOnGrey);
	TColor::InvertPalette();
	//gStyle->SetNumberContours(999);

	h_e->SetTitle("L-E distribution (#nu_{e})");
	h_e->GetXaxis()->SetTitle("E_{#nu} (GeV)");
	h_e->GetYaxis()->SetTitle("L (km)");
	hz_e->SetTitle("L-E distribution (#nu_{e})");
	hz_e->GetXaxis()->SetTitle("E_{#nu} (GeV)");
	hz_e->GetYaxis()->SetTitle("L (km)");

	h_m->SetTitle("L-E distribution (#nu_{#mu})");
	h_m->GetXaxis()->SetTitle("E_{#nu} (GeV)");
	h_m->GetYaxis()->SetTitle("L (km)");
	hz_m->SetTitle("L-E distribution (#nu_{#mu})");
	hz_m->GetXaxis()->SetTitle("E_{#nu} (GeV)");
	hz_m->GetYaxis()->SetTitle("L (km)");

	TCanvas *c = new TCanvas("c","",800,600);

	plot2D(c, h_e, "LEe");
	plot2D(c, h_m, "LEm");

	c->SetLogx();
	plot2D(c, h_e, "LEe_log");
	plot2D(c, h_m, "LEm_log");
	plot2D(c, hz_e, "LEe_zlog");
	plot2D(c, hz_m, "LEm_zlog");

	c->Close();

	TFile *f = new TFile("LE_1007.root","RECREATE");
	h_e->Write();
	h_m->Write();
	f->Close();

        //for(int i=410;i<430;i++){
        //std:cout << i << "\t" << h->ProjectionY("",i+1,i+1)->GetRandom() << "\n";}
}
