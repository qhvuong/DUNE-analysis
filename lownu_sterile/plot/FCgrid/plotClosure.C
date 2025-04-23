#include <DUNEStyle.h>
#include <TLegend.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}

void drawBestfit(TCanvas *c, TH1D*h, double trueVal, const char *outName){
  c->cd();
  c->Clear();
  h->Scale(1.0 / h->Integral()); 
  //h->SetMaximum(h->GetMaximum()*10.);
  auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  TLine *line = new TLine(trueVal, 0, trueVal, h->GetMaximum());
  line->SetLineColor(color);
  h->Draw("HIST");
  line->Draw("same");
  TLegend *legend = MakeLegend(0.6, 0.6, 0.85, 0.8);
  legend->AddEntry(line, Form("true value = %.4f", trueVal), "l");
  legend->Draw();
  dunestyle::CenterTitles(h);
  dunestyle::WIP();
  //c->SaveAs("FCdchi2.png");
  c->SaveAs(Form("%s.png",outName));
}

void drawDchi2(TCanvas *c, TH1D*h, double trueVal[4], const char *outName){
  c->cd();
  c->Clear();
  h->Scale(1.0 / h->Integral()); 
  //h->SetMaximum(h->GetMaximum()*10.);
  h->SetTitle(Form("U_{e4}^{2} = %.5f, U_{#mu4}^{2} = %.5f, U_{#tau4}^{2} = %.5f, #Deltam^{2} = %.1f", trueVal[0], trueVal[1], trueVal[2], trueVal[3]));
  //auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  // TLine *line = new TLine(trueVal, 0, trueVal, h->GetMaximum());
  // line->SetLineColor(color);
  h->Draw("HIST");
  // line->Draw("same");
  // TLegend *legend = MakeLegend(0.6, 0.6, 0.85, 0.8);
  // legend->AddEntry(line, Form("true value = %.3f", trueVal), "l");
  // legend->Draw();
  dunestyle::CenterTitles(h);
  dunestyle::WIP();
  c->RedrawAxis();
  //c->SaveAs("FCdchi2.png");
  c->SaveAs(Form("%s.png",outName));
}

void plotClosure()
{
  
  double ue42, um42, ut42, dm2, bfchi2, tgtchi2, nochi2, dchi2, ratio;
  double ue42Err, um42Err, ut42Err, dm2Err;

  double trueVals[3][4] = {
    {0.01, 0.02, 0., 5.},
    {0.2, 0.2, 0., 0.1},
    {1E-4, 1E-4, 0., 100.}
  };


  int N=100;
  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output/FCclosure";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  TH1D *hUe42[3], *hUm42[3], *hDm2[3];
  TH1D *h[3];

  TCanvas *c = new TCanvas("c","",800,600);

  for(int run=0; run<3; run++){
    TFile *out = new TFile(Form("FCclosure_%d.root",run),"RECREATE");
    h[run] = new TH1D(Form("h%d",run), "", 20, 0, 100);
    h[run]->GetXaxis()->SetTitle("FC #Delta#chi^{2}");
    h[run]->GetYaxis()->SetTitle("arbitrary unit");
    
    hUe42[run] = new TH1D(Form("hUe42%d",run), "", 20, 0, 2*trueVals[run][0]);
    hUe42[run]->GetXaxis()->SetTitle("U_{e4}^{2}");
    hUe42[run]->GetYaxis()->SetTitle("arbitrary unit");

    hUm42[run] = new TH1D(Form("hUm42%d",run), "", 20, 0, 2*trueVals[run][1]);
    hUm42[run]->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    hUm42[run]->GetYaxis()->SetTitle("arbitrary unit");

    hDm2[run] = new TH1D(Form("hDm2%d",run), "", 20, 0, 2*trueVals[run][3]);
    hDm2[run]->GetXaxis()->SetTitle("#Deltam^{2}");
    hDm2[run]->GetYaxis()->SetTitle("arbitrary unit");

    for(int i=0; i<N; i++) {
      if(i%10==0) std::cout << i*100./N << " percent...\n";
      TString filepath = Form("%s/%02d/output_%d.txt", data_path, run, i);
      ifstream f(filepath.Data());
      if(f) {
        //std::cout << filepath << "\n";
        f >> ue42 >> um42 >> ut42 >> dm2 >> ue42Err >> um42Err >> ut42Err >> dm2Err >> bfchi2 >> tgtchi2 >> nochi2 >> ratio;
        
        dchi2 = tgtchi2 - bfchi2;
        if(dchi2 < 0.) std::cout << i << "\t" << dchi2 << "\t" << tgtchi2 << "\t" << bfchi2 << "\n";
        h[run]->Fill(dchi2);
        hUe42[run]->Fill(ue42);
        hUm42[run]->Fill(um42);
        hDm2[run]->Fill(dm2);
      }
    }

    drawBestfit(c, hUe42[run], trueVals[run][0], Form("ue42_%d",run));
    drawBestfit(c, hUm42[run], trueVals[run][1], Form("um42_%d",run));
    drawBestfit(c, hDm2[run],  trueVals[run][3], Form("dm2_%d",run));
    drawDchi2(c, h[run], trueVals[run], Form("dchi2_%d",run));

    h[run]->Write();
    hUe42[run]->Write();
    hUm42[run]->Write();
    hDm2[run]->Write();


    out->Close();

  }
  


  //h[0]->Draw();
  

  //std::cout << GetCriticalChi2(h[0], 68.2/100.);

  /*
  auto color1 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
  auto color2 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  auto color3 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 3);



  h->Scale(1.0 / h->Integral()); 
  h->SetMaximum(h->GetMaximum()*10.);


  // Compute Confidence Levels (percentiles of a Gaussian)
  double mean = h->GetMean();
  double sigma = h->GetStdDev();
  
  double level1 = mean + 1 * sigma; // 68.2% (1σ)
  double level2 = mean + 2 * sigma; // 95.4% (2σ)
  double level3 = mean + 3 * sigma; // 99.7% (3σ)

  // Draw Confidence Level Lines
  TLine *line1 = new TLine(level1, 0, level1, h->GetMaximum());
  TLine *line2 = new TLine(level2, 0, level2, h->GetMaximum());
  TLine *line3 = new TLine(level3, 0, level3, h->GetMaximum());

  line1->SetLineColor(color1);
  line2->SetLineColor(color2);
  line3->SetLineColor(color3);

  line1->SetLineStyle(kSolid);
  line2->SetLineStyle(kSolid);
  line3->SetLineStyle(kSolid);

  line1->SetLineWidth(2);
  line2->SetLineWidth(2);
  line3->SetLineWidth(2);


  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogy();
  h->Draw("HIST");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  TLegend *legend = MakeLegend(0.6, 0.6, 0.85, 0.8);
  legend->AddEntry(line1, "1#sigma (68.2%)", "l");
  legend->AddEntry(line2, "2#sigma (95.4%)", "l");
  legend->AddEntry(line3, "3#sigma (99.7%)", "l");
  legend->Draw();
  dunestyle::CenterTitles(h);
  dunestyle::WIP();
  //c->SaveAs("FCdchi2.png");
  c->SaveAs("FCdchi2_tot.png");



  TFile *out = new TFile("FCtot.root","RECREATE");
  h->Write();
  out->Close();
  */
}

/*
  TH1D *hd[4];
  TH2D *hme[4], *hmt[4], *hdme[4], *hdmm[4], *hdmt[4];

  for(int j=0; j<4; j++){
  //h[j]  = new TH1D(Form("h%d",j), "", nbins, 0, 500);
  hd[j] = new TH1D(Form("hd%d",j), "", nbins, 0, 50);
  hme[j] = new TH2D(Form("hme%d",j), "", nbins, binEdges, nbins, binEdges);
  hmt[j] = new TH2D(Form("hmt%d",j), "", nbins, binEdges, nbins, binEdges);
  hdme[j] = new TH2D(Form("hdme%d",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmm[j] = new TH2D(Form("hdmm%d",j), "", nbins, binEdges, nbins, binEdgesM);
  hdmt[j] = new TH2D(Form("hdmt%d",j), "", nbins, binEdges, nbins, binEdgesM);
  }

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  for(int i=0; i<N; i++) {
    if(i%100==0) std::cout << i*100./N << " percent\n";

    //ifstream f0(Form("%s/FC_3pars/ut42_0/output_%d.txt",data_path,i));
    //ifstream f1(Form("%s/FC_3pars/ut42_1/output_%d.txt",data_path,i));
    //ifstream f2(Form("%s/FC_3pars/ut42_2/output_%d.txt",data_path,i));
    //ifstream f3(Form("%s/FC_3pars/test/output_%d.txt",data_path,i));
    ifstream f3(Form("%s/FC_4pars/0827/output_%d.txt",data_path,i));
    if(f3) {
      count+=1;
      //f0 >> ue42[0] >> um42[0] >> ut42[0] >> dm2[0] >> chi2[0] >> nochi2[0] >> ratio[0];
      //f1 >> ue42[1] >> um42[1] >> ut42[1] >> dm2[1] >> chi2[1] >> nochi2[1] >> ratio[1];
      //f2 >> ue42[2] >> um42[2] >> ut42[2] >> dm2[2] >> chi2[2] >> nochi2[2] >> ratio[2];
      f3 >> ue42[3] >> um42[3] >> ut42[3] >> dm2[3] >> chi2[3] >> nochi2[3] >> ratio[3];
      //double dchi2 = nochi2[0]-chi2[0];
      //if(dchi2 > 20)
      //std::cout << i << "\t" << nochi2[0] << "\t" << chi2[0] << "\t" << ue42[0] << "\t" << um42[0] << "\t" << ut42[0] << "\t" << dm2[0] << "\n";
      
      for(int j=0; j<4; j++){
      dchi2[j] = nochi2[j] - chi2[j];
      //if(dchi2[j] > 40) std::cout << i << "\t" << j << "\t" << dchi2[j] << "\t" << ue42[j] << "\t" << um42[j] << "\t" << ut42[j] << "\t" << dm2[j] << "\n";
      //h[j]->Fill(chi2[j]);
      hd[j]->Fill(dchi2[j]);
      hme[j]->Fill(um42[j], ue42[j]);
      hmt[j]->Fill(um42[j], ut42[j]);
      hdme[j]->Fill(ue42[j], dm2[j]);
      hdmm[j]->Fill(um42[j], dm2[j]);
      hdmt[j]->Fill(ut42[j], dm2[j]);
      }
      
      //f0.close();
      //f1.close();
      //f2.close();
      f3.close();
    }

  }

  double chi2_crit = GetCriticalChi2(hd[3], 99.7/100.);

  std::cout << "count = " << count << "\n";
  std::cout << "critical chi2 = " << chi2_crit << "\n";

  TFile *fout = new TFile("FC_chi2_500bins.root", "RECREATE");
  hd[3]->Write();
  fout->Close();
*/

/*

  hd[3]->SetTitle("Feldman-Cousins #Delta#chi^{2} distribution");
  hd[3]->GetXaxis()->SetTitle("#Delta#chi^{2}");
  hd[3]->GetYaxis()->SetTitle("Number of universes");

  TCanvas c;
  c.Clear();
  c.cd();
  c.SetGrid();
  gPad->SetLogy();
  //TLegend * leg = MakeLegend(0.50, 0.55, 0.90, 0.80);
  //TH1 * hFirst = nullptr;
  //for (std::size_t histIdx = 0; histIdx < 4; histIdx++)
  int histIdx = 3;
  
    //if(histIdx != 3) continue;
    TH1D *h = hd[histIdx];
    //auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto);
    h->SetLineColor(color);
    h->SetFillStyle(0);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy();  // need to leak it so it doesn't disappear
    //if (!hFirst)
      //hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("true #nu < %.1f GeV", nu_cuts[histIdx]), "l"));
  
  //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[0], Form("fixed U_{#tau4}^{2} = 0.0"), "l"));
  //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[1], Form("fixed U_{#tau4}^{2} = 0.3"), "l"));
  //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[2], Form("fixed U_{#tau4}^{2} = 0.7"), "l"));
  //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(hd[3], Form("free U_{#tau4}^{2}"), "l"));
  c.RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  //leg->Draw();
  h->SetMaximum(h->GetMaximum()*3.); // make some space for the watermark
  dunestyle::WIP();
  c.SaveAs("dchi2_500bins_0909_free.png");

*/



