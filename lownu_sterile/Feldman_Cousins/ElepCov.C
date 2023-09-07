static const int N = 500; // number of universes
static const int nbins = 100;
static const int nbins_E = 100;


TMatrixD covmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );

TMatrixD E_m(N, nbins_E);
TMatrixD E_e(N, nbins_E);
TMatrixD E_nue(N, nbins_E);
  
TMatrixD ECovars_mm    ( nbins_E, nbins_E );
TMatrixD ECovars_ee    ( nbins_E, nbins_E );
TMatrixD ECovars_nuenue  ( nbins_E, nbins_E );

TMatrixD ECovars_me  ( nbins_E, nbins_E );
TMatrixD ECovars_mnue( nbins_E, nbins_E );
TMatrixD ECovars_em  ( nbins_E, nbins_E );
TMatrixD ECovars_enue( nbins_E, nbins_E );
TMatrixD ECovars_nuem( nbins_E, nbins_E );
TMatrixD ECovars_nuee( nbins_E, nbins_E );

/*
TMatrixD ECorrel_m  ( nbins_E, nbins_E );
TMatrixD ECorrel_e  ( nbins_E, nbins_E );
TMatrixD ECorrel_nue( nbins_E, nbins_E );

TMatrixD ECorrel_me  ( nbins_E, nbins_E );
TMatrixD ECorrel_mnue( nbins_E, nbins_E );
TMatrixD ECorrel_em  ( nbins_E, nbins_E );
TMatrixD ECorrel_enue( nbins_E, nbins_E );
TMatrixD ECorrel_nuem( nbins_E, nbins_E );
TMatrixD ECorrel_nuee( nbins_E, nbins_E );
*/

TMatrixD ECovars( 3*nbins_E, 3*nbins_E );
//TMatrixD ECorrel( 3*nbins_E, 3*nbins_E );

void ElepCov()
{

  // covariance matrix
  TFile * covfile = new TFile( "/dune/app/users/qvuong/lownu/sigma_mtr/total_sigma_mtr.root", "OLD" );
  TH2D * hcovmx = (TH2D*) covfile->Get( "cv" );
  //TFile * covfile = new TFile( "/dune/app/users/qvuong/lownu/sigma_mtr/new_mtr/wgt_MaCCQE_Covmtr32.root", "OLD" );
  //TH2D * hcovmx = (TH2D*) covfile->Get( "hcv" );

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      covmx[x][y] = hcovmx->GetBinContent( x+1, y+1 );

      //std::cout << x << "\t" << y << "\t" << covmx[x][y] << "\n";
    }
  }

    // cholesky decomposition
  TDecompChol decomp( covmx );
  if( !decomp.Decompose() ) {
    printf( "Main covariance matrix failed Cholesky decomposition\n" );
    return;
  }
  const TMatrixD chol = decomp.GetU();

  TRandom3 * rand = new TRandom3(12345);

  // make random number matrix
  for( int i = 0; i < N; ++i ) {
    //if( i>0 && i%50==0 ) std::cout << "Progress: " << i/N*100. << "..." << "\n";
    double mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
    }

  }

  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns

      // compute column covariance
      double covar = 0.;
      for( int k = 0; k < N; ++k ) { // rows
        covar += (scales[k][i] * scales[k][j]); // means are 0 already by construction
      }
      scaleCovars[i][j] = covar/N;

      //std::cout << i << "\t" << j << "\t" << scaleCovars[i][j] << "\n";
    }
  }

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;

  char name[20] = "wgt_MaCCQE";
  int cutNu, cutEv;


  //for(para = 1; para <3; para++) {
  TFile *f     = new TFile("/dune/app/users/qvuong/lownu/gen_data/new_tgt/total_output_050223.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/lownu/gen_data/nue_output.root");
  //for(cutNu = 0; cutNu < 4; cutNu++) {
  //if(cutNu != 0 && cutNu != 3 ) continue;
  cutNu = 3;
  //std::cout << cutNu << "\n";
  //for(cutEv = 0; cutEv < 3; cutEv++) {
  cutEv = 2;

  TH1D *CC_m_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CC_e_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get(Form("hElep%d",cutEv));

  TH1D *m     = new TH1D("m","",100,0,16);
  TH1D *e     = new TH1D("e","",100,0,16);
  TH1D *nue   = new TH1D("nue","",100,0,16);

  
  for( int u = 0; u < N; ++u ) {
    m->Reset();
    e->Reset();
    nue->Reset();

    for(int CCmb=0; CCmb<nbins_E; CCmb++) {
      int fluxbin = CCmb+1;
      double evtwgt = scales[u][fluxbin];

      m->Fill(CC_m_nom->GetBinContent(CCmb+1),1.+evtwgt);
    }
    for(int CCeb=0; CCeb<nbins_E; CCeb++) {
      int fluxbin = CCeb+101;
      double evtwgt = scales[u][fluxbin];

      e->Fill(CC_e_nom->GetBinContent(CCeb+1),1.+evtwgt);
    }
    for(int nueb=0; nueb<nbins_E; nueb++) {
      int fluxbin = nueb+201;
      double evtwgt = scales[u][fluxbin];

      nue->Fill(CC_e_nom->GetBinContent(nueb+1),1.+evtwgt);
    }


    for(int i=0; i<nbins_E; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
      //E_m_nue[u][i] = m_nue->GetBinContent(i+1);
      //E_e_nue[u][i] = e_nue->GetBinContent(i+1);
      E_nue[u][i]   = nue->GetBinContent(i+1);
    }
   }

  //c->SaveAs("ratio_nue.png");
  //c_m3->SaveAs("ratio_CC_m3.png");


  for( int i = 0; i < nbins_E; ++i ) { // columns
    for( int j = 0; j < nbins_E; ++j ) { // columns
      // compute column covariance
      double covar_mm = 0.;
      double covar_ee = 0.;
      double covar_nuenue = 0.;
      double covar_me = 0.;
      double covar_mnue = 0.;
      double covar_em = 0.;
      double covar_enue = 0.;
      double covar_nuem = 0.;
      double covar_nuee = 0.;

      double var_m_i = 0.;
      double var_m_j = 0.;
      double var_e_i = 0.;
      double var_e_j = 0.;
      double var_nue_i = 0.;
      double var_nue_j = 0.;

      for( int k = 0; k < N; ++k ) { // rows
        covar_me     += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1));
        covar_ee     += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 
        covar_nuenue += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 

        covar_me     += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1));
        covar_mnue   += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_em     += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_enue   += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem   += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_nuee   += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 

/*
        var_m_i += (E_m[k][i] - CC_m_nom->GetBinContent(i+1)) * (E_m[k][i] - CC_m_nom->GetBinContent(i+1));
        var_m_j += (E_m[k][j] - CC_m_nom->GetBinContent(j+1)) * (E_m[k][j] - CC_m_nom->GetBinContent(j+1));

        var_e_i += (E_e[k][i] - CC_e_nom->GetBinContent(i+1)) * (E_e[k][i] - CC_e_nom->GetBinContent(i+1));
        var_e_j += (E_e[k][j] - CC_e_nom->GetBinContent(j+1)) * (E_e[k][j] - CC_e_nom->GetBinContent(j+1));

        var_nue_i += (E_nue[k][i] - nue_nom->GetBinContent(i+1)) * (E_nue[k][i] - nue_nom->GetBinContent(i+1));
        var_nue_j += (E_nue[k][j] - nue_nom->GetBinContent(j+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
*/
      }

      ECovars_mm[i][j]     = covar_mm    /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      ECovars_ee[i][j]     = covar_ee    /(N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      ECovars_nuenue[i][j] = covar_nuenue/(N * nue_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j+1));

      ECovars_me[i][j]   = covar_me  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      ECovars_mnue[i][j] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_em[i][j]   = covar_em  /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      ECovars_enue[i][j] = covar_enue/(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      ECovars_nuem[i][j] = covar_nuem/(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      ECovars_nuee[i][j] = covar_nuee/(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
/*
      ECorrel_m[i][j]   = covar_m  /sqrt(var_m_i   * var_m_j);
      ECorrel_e[i][j]   = covar_e  /sqrt(var_e_i   * var_e_j);
      ECorrel_nue[i][j] = covar_nue/sqrt(var_nue_i * var_nue_j);

      ECorrel_me[i][j]   = covar_me  /sqrt(var_m_i   * var_e_j);
      ECorrel_mnue[i][j] = covar_mnue/sqrt(var_m_i   * var_nue_j);
      ECorrel_em[i][j]   = covar_em  /sqrt(var_e_i   * var_m_j);
      ECorrel_enue[i][j] = covar_enue/sqrt(var_e_i   * var_nue_j);
      ECorrel_nuem[i][j] = covar_nuem/sqrt(var_nue_i * var_m_j);
      ECorrel_nuee[i][j] = covar_nuee/sqrt(var_nue_i * var_e_j);
*/
    }
  }

  for(int i = 0; i < nbins_E; i++) {
    for(int j = 0; j < nbins_E; j++) {
      ECovars[i][j]           = ECovars_mm[i][j];
      ECovars[i][j+nbins_E]   = ECovars_me[i][j];
      ECovars[i][j+2*nbins_E] = ECovars_mnue[i][j];

      ECovars[i+nbins_E][j]           = ECovars_em[i][j];
      ECovars[i+nbins_E][j+nbins_E]   = ECovars_ee[i][j];
      ECovars[i+nbins_E][j+2*nbins_E] = ECovars_enue[i][j];

      ECovars[i+2*nbins_E][j]           = ECovars_nuem[i][j];
      ECovars[i+2*nbins_E][j+nbins_E]   = ECovars_nuee[i][j];
      ECovars[i+2*nbins_E][j+2*nbins_E] = ECovars_nuenue[i][j];

/*
      ECorrel[i][j]           = ECorrel_m[i][j];
      ECorrel[i][j+nbins_E]   = ECorrel_me[i][j];
      ECorrel[i][j+2*nbins_E] = ECorrel_mnue[i][j];

      ECorrel[i+nbins_E][j]           = ECorrel_em[i][j];
      ECorrel[i+nbins_E][j+nbins_E]   = ECorrel_e[i][j];
      ECorrel[i+nbins_E][j+2*nbins_E] = ECorrel_enue[i][j];

      ECorrel[i+2*nbins_E][j]           = ECorrel_nuem[i][j];
      ECorrel[i+2*nbins_E][j+nbins_E]   = ECorrel_nuee[i][j];
      ECorrel[i+2*nbins_E][j+2*nbins_E] = ECorrel_nue[i][j];
*/
    }
  }

  TH2D *hcv = new TH2D("hcv","",300,0,300,300,0,300);
  TH2D *hcr = new TH2D("hcr","",300,0,300,300,0,300);
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
      //hcr->SetBinContent(i+1, j+1, ECorrel[i][j]);
    }
  }
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;
  if(cutEv == 0)      Ev = 3.0;
  else if(cutEv == 1) Ev = 0.8;
  else if(cutEv == 2) Ev = 0.5;

  hcv->SetStats(0);
  hcv->SetTitle(Form("%s Covariance (nu<%.1fGeV && Etheta2<%.1fMeV)",name,nu,Ev));

  //hcr->SetStats(0);
  //hcr->SetTitle(Form("%s Correlation (nu<%.1fGeV && Etheta2<%.1fMeV)",name,nu,Ev));

  TCanvas *ccv = new TCanvas("ccv","",900,800);
  hcv->Draw("colz");
  ccv->SaveAs(Form("FC_sig_%d.png",N));

/*
  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 1.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  TCanvas *ccr = new TCanvas("ccr","",900,800);
  hcr->GetZaxis()->SetRangeUser(-1., 1.);
  hcr->Draw("colz");
  ccr->SaveAs(Form("%s_Cor%d%d%d_%d.png",name,para,cutNu,cutEv,N));
*/

  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/Feldman_Cousins/sig_%d.root",N),"RECREATE");
  hcv->Write();
  //hcr->Write();
  out->Close();
}

/*
  TDecompSVD svd0(ElepCovars0);
  //svd0.Decompose();
  TMatrixD inv0 = svd0.Invert();
  TDecompSVD svd3(ElepCovars3);
  //svd3.Decompose();
  TMatrixD inv3 = svd3.Invert();

  TH2D *test0 = new TH2D(inv0);
  TH2D *test3 = new TH2D(inv3);
  TCanvas *c2 = new TCanvas("c2","",1800,1350);
  c2->Divide(2,1);
  c2->cd(1);
  test0->Draw("colz");
  c2->cd(2);
  test3->Draw("colz");
  c2->SaveAs(Form("test1_total_Cov_%d.png",N));
*/
/*
    if(u<20) {

      c->cd();
      pad1->cd();
      TH1D *CC_m0 = (TH1D*)nue->Clone();
      CC_m0->SetTitle("Nue Elep in 20 universes (nu<10.0)");
      CC_m0->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m0->SetLineColor(u+1);
      CC_m0->SetStats(0);
      CC_m0->Draw("same");
      
      c->cd();
      pad2->cd();
      TH1D *ratio_m0 = (TH1D*)nue->Clone("ratio");
      ratio_m0->SetLineColor(u+1);
      ratio_m0->SetMaximum(1.7);
      ratio_m0->SetMinimum(0.3);
      ratio_m0->Sumw2(); 
      ratio_m0->SetStats(0); 
      ratio_m0->Divide(nue_nom);
      ratio_m0->SetTitle("");
      ratio_m0->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m0->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m0->Draw("same");

      c_m3->cd();
      pad1->cd();
      TH1D *CC_m3 = (TH1D*)m3->Clone();
      CC_m3->SetTitle("CC_m Elep in 20 universes (nu<0.3)");
      CC_m3->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m3->SetLineColor(u+1);
      CC_m3->SetStats(0);
      CC_m3->Draw("same");
      
      c_m3->cd();
      pad2->cd();
      TH1D *ratio_m3 = (TH1D*)m3->Clone("ratio");
      ratio_m3->SetLineColor(u+1);
      ratio_m3->SetMaximum(1.7);
      ratio_m3->SetMinimum(0.3);
      ratio_m3->Sumw2(); 
      ratio_m3->SetStats(0); 
      ratio_m3->Divide(CC_m3_nom);
      ratio_m3->SetTitle("");
      ratio_m3->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m3->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m3->Draw("same");

    }
    else continue;
*/
