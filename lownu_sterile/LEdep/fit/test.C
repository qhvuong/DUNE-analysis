#include <iostream>

void test()
{
	
  TTree *pardir = new TTree("pardir", "fitParameters");

  Double_t para[29][8], para0, para1, para2, para3, para4, para5, norm0, norm1, norm2, norm;

  pardir->Branch("para0",&para0,"para0/D");
  pardir->Branch("para1",&para1,"para1/D");
  pardir->Branch("para2",&para2,"para2/D");
  pardir->Branch("para3",&para3,"para3/D");
  pardir->Branch("para4",&para4,"para4/D");
  pardir->Branch("para5",&para5,"para5/D");
  pardir->Branch("norm",&norm,"norm/D");

  TFile *f = new TFile("/dune/app/users/qvuong/lownu/LEdep/LE_1112_2.root","READ");
  TH2D *h = (TH2D*)f->Get("h_m");
  TH1D *hL[29];

  TF1 *e = new TF1("e","exp([0]+[1]*x)",0.5,0.57);		//totally expo

  TF1 *p1 = new TF1("p1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.35,0.55);	//totally poly

  TF1 *cs1 = new TF1("cs1","[0]",0.34,0.35);
  TF1 *cs2 = new TF1("cs2","[0]",0.55,0.6);

  gStyle->SetOptFit();

  for(int i = 0; i < 29; i++) {

    TCanvas *c = new TCanvas("c","",800,600);
 
    hL[i] = (TH1D*)h->ProjectionY("",i+1,i+1);
    hL[i]->Scale(1./hL[i]->Integral());
    hL[i]->SetStats(0);
    hL[i]->SetMinimum(0.);
    hL[i]->SetTitle("L distribution (numu)");

    hL[i]->Fit(p1,"R");
    p1->GetParameters(&para[i][0]);
    norm0 = p1->Integral(0.35,0.55);
    hL[i]->Draw();

    double p4 = 0.;
    double p5 = 0.;
    int count4 = 0;
    int count5 = 0;
    for(int j=12; j<15; j++){
      if(hL[i]->GetBinContent(j+1) != 0) {
        count4 = count4 + 1;
        p4 = p4 + hL[i]->GetBinContent(j+1);
	//if(i==28) std::cout << hL[i]->GetBinContent(j+1) << "\n";
      }
    }
    for(int j=75; j<90; j++){
      if(hL[i]->GetBinContent(j+1) != 0) {
        count5 = count5 + 1;
        p5 = p5 + hL[i]->GetBinContent(j+1);
      }
    }
    if(count4 == 0) para[i][4] = 0;
    else            para[i][4] = p4/count4;

    if(count5 == 0) para[i][5] = 0;
    else            para[i][5] = p5/count5;

    cs1->SetParameter(0,para[i][4]);
    cs2->SetParameter(0,para[i][5]);
      
    norm1 = cs1->Integral(0.34,0.35);
    norm2 = cs2->Integral(0.55,0.6);

    norm = 1./(norm0 + norm1 + norm2);
    
    para0 = para[i][0];
    para1 = para[i][1];
    para2 = para[i][2];
    para3 = para[i][3];
    para4 = para[i][4];
    para5 = para[i][5];
    para[i][6] = norm;
    para[i][7] = 1./norm0;

    cs1->Draw("same");
    cs2->Draw("same");
    c->SaveAs(Form("hL_m_fit_%d.png",i));
    pardir->Fill();

  }

  ofstream myfile;
  myfile.open("/dune/app/users/qvuong/lownu/LEdep/fitPara_m.txt");
  for(int i = 0; i < 29; i++) {
    //if(i==1) para[i][2] = para[i][3] = 0.;
    myfile << para[i][0] << "\t" << para[i][1] << "\t" << para[i][2] << "\t" << para[i][3] << "\t" << para[i][4] << "\t" << para[i][5] << "\t" << para[i][6] << "\t" << para[i][7] << "\n";
  }
  myfile.close();


  TFile *fout = new TFile("/dune/app/users/qvuong/lownu/LEdep/fitPara_m.root","RECREATE");
  pardir->Write();
  fout->Close();

  
/*
  TF1 *e = new TF1("e","exp([0]+[1]*x) * pow(sin([2]*x),2)",0.3,0.6);
  TF1 *p = new TF1("p","([0]+[1]*x+[2]*x*x+[3]*x*x*x) * pow(sin([4]*x),2)",0.3,0.6);

  double par[3];
  par[0] = -29.2169;
  par[1] = 0.0473503;
  par[2] = 1.27*6/4;

  e->SetParameter(0,par[0]);
  e->SetParameter(1,par[1]);
  e->SetParameter(2,par[2]);

  p->SetParameter(0,-0.559594);
  p->SetParameter(1,0.00422911);
  p->SetParameter(2,-1.0622e-5);
  p->SetParameter(3,8.91259e-9);

  std::cout << e->Integral(0.3,0.6) << "\t" << p->Integral(0.3,0.6) << "\n";
*/

}

/*
  for(int i = 0; i < 30; i++) {
    hL_e[i] = (TH1D*)h_e->ProjectionY("",i+1,i+1);
    hL_e[i]->Scale(1./hL_e[i]->Integral());
    hL_e[i]->SetStats(0);
    hL_e[i]->Fit(p1,"R");
    p1->GetParameters(&para_e[i][0]);
    norm0_e = p1->Integral(0.35,0.55);
    //hL_e[i]->Draw();
    //c->SaveAs(Form("hL_e_fit_%d.png",i));
      
    para_e[i][4] = 0.;
    para_e[i][5] = 0.;
    for(int j=0; j<15; j++){
      para_e[i][4] += hL_e[i]->GetBinContent(j);
    }
    for(int j=75; j<90; j++){
      para_e[i][5] += hL_e[i]->GetBinContent(j);
    }
    para_e[i][4] = para_e[i][4]/15;
    para_e[i][5] = para_e[i][5]/15;

    cs1->SetParameter(0,para_e[i][4]);
    cs2->SetParameter(0,para_e[i][5]);
      
    norm1_e = cs1->Integral(0.3,0.35);
    norm2_e = cs2->Integral(0.55,0.6);

    norm_e = 1./(norm0_e + norm1_e + norm2_e);

    para0_e = para_e[i][0];
    para1_e = para_e[i][1];
    para2_e = para_e[i][2];
    para3_e = para_e[i][3];
    para4_e = para_e[i][4];
    para5_e = para_e[i][5];
    para_e[i][6] = norm_e;
    para_e[i][7] = 1./norm0_e;

    //pardir->Fill();
  }
*/  
