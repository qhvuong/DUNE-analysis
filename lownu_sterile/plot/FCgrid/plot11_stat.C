void plot11_stat()
{
  char seed[10]="11";
  
  double Uee2, Umm2, dm2, chi2, nochi2, dchi2;
  int N=10000;


  double xMin = 1e-5;
  double xMax = 1.0;
  double mMin = 1e-4;
  double mMax = 1e+3;

  int nBins = 100;
  double logXMin = TMath::Log10(xMin);
  double logXMax = TMath::Log10(xMax);
  double binWidth = (logXMax - logXMin) / nBins;

  double logmMin = TMath::Log10(mMin);
  double logmMax = TMath::Log10(mMax);
  double mbinWidth = (logmMax - logmMin) / nBins;

  std::vector<double> binEdges(nBins + 1, 0);
  std::vector<double> mbinEdges(nBins + 1, 0);
  for (int i = 0; i <= nBins; ++i) {
    binEdges[i] = TMath::Power(10, logXMin + i * binWidth);
    mbinEdges[i] = TMath::Power(10, logmMin + i * mbinWidth);
  }

  TH2D *chi2Cor = new TH2D("chi2Cor","",100,0,500, 100,0,500);
  TH2D *dChi2Vsdm2 = new TH2D("dChi2Vsdm2","", nBins, &mbinEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUmm2_s = new TH2D("dChi2VsUmm2_s","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUmm2_m = new TH2D("dChi2VsUmm2_m","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUmm2_l = new TH2D("dChi2VsUmm2_l","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUee2_s = new TH2D("dChi2VsUee2_s","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUee2_m = new TH2D("dChi2VsUee2_m","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dChi2VsUee2_l = new TH2D("dChi2VsUee2_l","", nBins, &binEdges[0], 100,-0.5,20.);
  TH2D *dm2VsUmm2_s = new TH2D("dm2VsUmm2_s","", nBins, &binEdges[0], nBins, &mbinEdges[0]);
  TH2D *dm2VsUmm2_l = new TH2D("dm2VsUmm2_l","", nBins, &binEdges[0], nBins, &mbinEdges[0]);
  TH2D *dm2VsUee2_s = new TH2D("dm2VsUee2_s","", nBins, &binEdges[0], nBins, &mbinEdges[0]);
  TH2D *dm2VsUee2_l = new TH2D("dm2VsUee2_l","", nBins, &binEdges[0], nBins, &mbinEdges[0]);
  TH2D *Uee2VsUmm2_s = new TH2D("Uee2VsUmm2_s","", nBins, &binEdges[0], nBins, &binEdges[0]);
  TH2D *Uee2VsUmm2_l = new TH2D("Uee2VsUmm2_l","", nBins, &binEdges[0], nBins, &binEdges[0]);

  int tot = 0, countNO = 0, count = 0;
  for(int i=0; i<N; i++) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/stat_new/s%s/output_%d.txt",seed,i));

    if(i%100==0) std::cout << i*100./N << " percent\n";

    if(!f) continue;

    if (f) {
      tot += 1;
      f >> Uee2 >> Umm2 >> dm2 >> chi2 >> nochi2;
      dchi2 = nochi2 - chi2;
 
      if(nochi2<0) {
        std::cout << i << "\t" << chi2 << "\t" << nochi2 << "\n";
        countNO += 1;
      }
      
      if(chi2<0) {
        std::cout << i << "\t" << chi2 << "\t" << nochi2 << "\n";
        count += 1;
      }


      chi2Cor->Fill(nochi2, chi2);   
      dChi2Vsdm2->Fill(dm2, dchi2);

      if(dm2<1.) {
        dChi2VsUmm2_s->Fill(Umm2, dchi2);
        dChi2VsUee2_s->Fill(Uee2, dchi2);}
      if(dm2>=1. && dm2<=10.) {
        dChi2VsUmm2_m->Fill(Umm2, dchi2);
        dChi2VsUee2_m->Fill(Uee2, dchi2);}
      if(dm2>10.) {
        dChi2VsUmm2_l->Fill(Umm2, dchi2);
        dChi2VsUee2_l->Fill(Uee2, dchi2);}

      if(dchi2<5.) {
        dm2VsUmm2_s->Fill(Umm2, dm2);
        dm2VsUee2_s->Fill(Uee2, dm2);
        Uee2VsUmm2_s->Fill(Umm2, Uee2);}
      if(dchi2>=5.) {
        dm2VsUmm2_l->Fill(Umm2, dm2);
        dm2VsUee2_l->Fill(Uee2, dm2);
        Uee2VsUmm2_l->Fill(Umm2, Uee2);}
    }
    

    f.close();
  }

  std::cout << tot << "\t" << countNO << "\t" << count << "\n";


  chi2Cor->SetStats(0);
  chi2Cor->SetTitle("#chi^{2} vs #chi^{2} at no osc");
  chi2Cor->GetXaxis()->SetTitle("#chi^{2} NO");
  chi2Cor->GetYaxis()->SetTitle("#chi^{2}");
  
  dChi2Vsdm2->SetStats(0);
  dChi2Vsdm2->SetTitle("#Delta#chi^{2} vs #Deltam^{2}");
  dChi2Vsdm2->GetXaxis()->SetTitle("#Deltam^{2}");
  dChi2Vsdm2->GetYaxis()->SetTitle("#Delta#chi^{2}");
  
  dChi2VsUmm2_s->SetStats(0);
  dChi2VsUmm2_s->SetTitle("#Delta#chi^{2} vs U_{#mu4}^{2} (#Deltam^{2} < 1.0)");
  dChi2VsUmm2_s->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  dChi2VsUmm2_s->GetYaxis()->SetTitle("#Delta#chi^{2}");
  dChi2VsUmm2_m->SetStats(0);
  dChi2VsUmm2_m->SetTitle("#Delta#chi^{2} vs U_{#mu4}^{2} (1.0 <= #Deltam^{2} <= 10.)");
  dChi2VsUmm2_m->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  dChi2VsUmm2_m->GetYaxis()->SetTitle("#Delta#chi^{2}");
  dChi2VsUmm2_l->SetStats(0);
  dChi2VsUmm2_l->SetTitle("#Delta#chi^{2} vs U_{#mu4}^{2} (#Deltam^{2} > 10.)");
  dChi2VsUmm2_l->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  dChi2VsUmm2_l->GetYaxis()->SetTitle("#Delta#chi^{2}");
  
  dChi2VsUee2_s->SetStats(0);
  dChi2VsUee2_s->SetTitle("#Delta#chi^{2} vs U_{e4}^{2} (#Deltam^{2} < 1.0)");
  dChi2VsUee2_s->GetXaxis()->SetTitle("U_{e4}^{2}");
  dChi2VsUee2_s->GetYaxis()->SetTitle("#Delta#chi^{2}");
  dChi2VsUee2_m->SetStats(0);
  dChi2VsUee2_m->SetTitle("#Delta#chi^{2} vs U_{e4}^{2} (1.0 <= #Deltam^{2} <= 10.)");
  dChi2VsUee2_m->GetXaxis()->SetTitle("U_{e4}^{2}");
  dChi2VsUee2_m->GetYaxis()->SetTitle("#Delta#chi^{2}");
  dChi2VsUee2_l->SetStats(0);
  dChi2VsUee2_l->SetTitle("#Delta#chi^{2} vs U_{e4}^{2} (#Deltam^{2} > 10.)");
  dChi2VsUee2_l->GetXaxis()->SetTitle("U_{e4}^{2}");
  dChi2VsUee2_l->GetYaxis()->SetTitle("#Delta#chi^{2}");
  
  dm2VsUmm2_s->SetStats(0);
  dm2VsUmm2_s->SetTitle("#Deltam^{2} vs U_{#mu4}^{2} (#Delta#chi^{2} < 5.)");
  dm2VsUmm2_s->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  dm2VsUmm2_s->GetYaxis()->SetTitle("#Deltam^{2}");
  dm2VsUmm2_l->SetStats(0);
  dm2VsUmm2_l->SetTitle("#Deltam^{2} vs U_{#mu4}^{2} (#Delta#chi^{2} >= 5.)");
  dm2VsUmm2_l->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  dm2VsUmm2_l->GetYaxis()->SetTitle("#Deltam^{2}");
  
  dm2VsUee2_s->SetStats(0);
  dm2VsUee2_s->SetTitle("#Deltam^{2} vs U_{e4}^{2} (#Delta#chi^{2} < 5.)");
  dm2VsUee2_s->GetXaxis()->SetTitle("U_{e4}^{2}");
  dm2VsUee2_s->GetYaxis()->SetTitle("#Deltam^{2}");
  dm2VsUee2_l->SetStats(0);
  dm2VsUee2_l->SetTitle("#Deltam^{2} vs U_{e4}^{2} (#Delta#chi^{2} >= 5.)");
  dm2VsUee2_l->GetXaxis()->SetTitle("U_{e4}^{2}");
  dm2VsUee2_l->GetYaxis()->SetTitle("#Deltam^{2}");

  Uee2VsUmm2_s->SetStats(0);
  Uee2VsUmm2_s->SetTitle("U_{e4}^{2} vs U_{#mu4}^{2} (#Delta#chi^{2} < 5.)");
  Uee2VsUmm2_s->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  Uee2VsUmm2_s->GetYaxis()->SetTitle("U_{e4}^{2}");
  Uee2VsUmm2_l->SetStats(0);
  Uee2VsUmm2_l->SetTitle("U_{e4}^{2} vs U_{#mu4}^{2} (#Delta#chi^{2} >= 5.)");
  Uee2VsUmm2_l->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  Uee2VsUmm2_l->GetYaxis()->SetTitle("U_{e4}^{2}");

  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,800);
  c->SetGrid();
  chi2Cor->SetMaximum(450);
  chi2Cor->Draw("colz");
  c->SaveAs(Form("new/chi2Cor_%s_stat.png",seed));
 
  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->SetGrid();
  c1->SetLogx();
  dChi2Vsdm2->SetMaximum(200);
  dChi2Vsdm2->Draw("colz");
  c1->SaveAs(Form("new/dChi2Vsdm2_%s_stat.png",seed));
  
  TCanvas *c2 = new TCanvas("c2","",1800,500);
  c2->Divide(3,1);
  c2->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUmm2_s->SetMaximum(100);
  dChi2VsUmm2_s->Draw("colz");
  c2->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUmm2_m->SetMaximum(20);
  dChi2VsUmm2_m->Draw("colz");
  c2->cd(3);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUmm2_l->SetMaximum(2);
  dChi2VsUmm2_l->Draw("colz");
  c2->SaveAs(Form("new/dChi2VsUmm2_%s_stat.png",seed));

  TCanvas *c3 = new TCanvas("c3","",1800,500);
  c3->Divide(3,1);
  c3->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUee2_s->SetMaximum(100);
  dChi2VsUee2_s->Draw("colz");
  c3->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUee2_m->SetMaximum(12);
  dChi2VsUee2_m->Draw("colz");
  c3->cd(3);
  gPad->SetGrid();
  gPad->SetLogx();
  dChi2VsUee2_l->SetMaximum(4);
  dChi2VsUee2_l->Draw("colz");
  c3->SaveAs(Form("new/dChi2VsUee2_%s_stat.png",seed));
  
  TCanvas *c4 = new TCanvas("c4","",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  dm2VsUmm2_s->SetMaximum(100);
  dm2VsUmm2_s->Draw("colz");
  c4->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  dm2VsUmm2_l->SetMaximum(40);
  dm2VsUmm2_l->Draw("colz");
  c4->SaveAs(Form("new/dm2VsUmm2_%s_stat.png",seed));
  
  TCanvas *c5 = new TCanvas("c5","",1200,600);
  c5->Divide(2,1);
  c5->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  dm2VsUee2_s->SetMaximum(240);
  dm2VsUee2_s->Draw("colz");
  c5->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  dm2VsUee2_l->SetMaximum(50);
  dm2VsUee2_l->Draw("colz");
  c5->SaveAs(Form("new/dm2VsUee2_%s_stat.png",seed));

  TCanvas *c6 = new TCanvas("c6","",1200,600);
  c6->Divide(2,1);
  c6->cd(1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  Uee2VsUmm2_s->SetMaximum(350);
  Uee2VsUmm2_s->Draw("colz");
  c6->cd(2);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  Uee2VsUmm2_l->SetMaximum(400);
  Uee2VsUmm2_l->Draw("colz");
  c6->SaveAs(Form("new/Uee2VsUmm2_%s_stat.png",seed));


  TFile *out = new TFile(Form("new/seed%s_stat.root",seed),"RECREATE");
  chi2Cor->Write();
  dChi2Vsdm2->Write();
  dChi2VsUmm2_s->Write();
  dChi2VsUmm2_m->Write();
  dChi2VsUmm2_l->Write();
  dChi2VsUee2_s->Write();
  dChi2VsUee2_m->Write();
  dChi2VsUee2_l->Write();
  dm2VsUmm2_s->Write();
  dm2VsUmm2_l->Write();
  dm2VsUee2_s->Write();
  dm2VsUee2_l->Write();
  Uee2VsUmm2_s->Write();
  Uee2VsUmm2_l->Write();
  out->Close();

}
