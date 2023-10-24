void plot()
{
  double Uee2, Umm2, dm2, sen;
  int N=5000;
  TH1D *h = new TH1D("h","",100,0,80);

  for(int i=0; i<N; i++) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/FC/output_%d.txt",i));

    if(i%20==0) std::cout << i*100./N << " percent\n";

    if (f) {
      f >> Uee2 >> Umm2 >> dm2 >> sen;

      h->Fill(sen);    
    }
    
    else std::cout << i << "\n"; 

    f.close();
  }


  //h->SetStats(0);
  //h->SetMaximum(20);
  //h->SetMinimum(0.1);
  h->SetTitle("Sensitivity");
  h->GetXaxis()->SetTitle("sensitivity");
  h->GetYaxis()->SetTitle("Number of throws");

  TCanvas *c = new TCanvas("c","",800,800);
  c->SetGrid();
  //cL->SetLogx();
  //cL->SetLogy();
  h->Draw();
  c->SaveAs("sensitivity.png");
/*
  TFile *out = new TFile(Form("noNuCut_dm%d.root",j),"RECREATE");
  h->Write();
  out->Close();
*/
}
