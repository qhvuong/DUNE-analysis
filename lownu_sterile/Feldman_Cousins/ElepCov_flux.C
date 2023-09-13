static const int N = 500; // number of universes
static const int nbins = 100;


TMatrixD covmx( 3*nbins, 3*nbins );
TMatrixD scaleCovars( 3*nbins, 3*nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, 3*nbins );

TMatrixD E_m(N, nbins);
TMatrixD E_e(N, nbins);
TMatrixD E_nue(N, nbins);

 
TMatrixD ECovars_m    ( nbins, nbins );
TMatrixD ECovars_e    ( nbins, nbins );
TMatrixD ECovars_nue  ( nbins, nbins );

TMatrixD ECovars_me  ( nbins, nbins );
TMatrixD ECovars_mnue( nbins, nbins );
TMatrixD ECovars_em  ( nbins, nbins );
TMatrixD ECovars_enue( nbins, nbins );
TMatrixD ECovars_nuem( nbins, nbins );
TMatrixD ECovars_nuee( nbins, nbins );

TMatrixD ECovars( 3*nbins, 3*nbins );

/*
TMatrixD ECorrel_m  ( 3*nbins, 3*nbins );
TMatrixD ECorrel_e  ( 3*nbins, 3*nbins );
TMatrixD ECorrel_nue( 3*nbins, 3*nbins );

TMatrixD ECorrel_me  ( 3*nbins, 3*nbins );
TMatrixD ECorrel_mnue( 3*nbins, 3*nbins );
TMatrixD ECorrel_em  ( 3*nbins, 3*nbins );
TMatrixD ECorrel_enue( 3*nbins, 3*nbins );
TMatrixD ECorrel_nuem( 3*nbins, 3*nbins );
TMatrixD ECorrel_nuee( 3*nbins, 3*nbins );

TMatrixD ECorrel( 3*3*nbins, 3*3*nbins );
*/

void ElepCov_flux()
{

  // covariance matrix
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root");

  int cutNu=3, cutEv=2;
  const char name1[20] = "wgt_MaCCQE";

  TH1D *CC_m_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name1,cutNu));
  TH1D *CC_e_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name1,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get(Form("hElep%d",cutEv));

  std::cout << CC_m_nom->GetNbinsX() << "\n";
  std::cout << CC_e_nom->GetNbinsX() << "\n";
  std::cout << nue_nom->GetNbinsX() << "\n";

  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      if(x==y){
        covmx[x][y]                 = CC_m_nom->GetBinContent( x+1, y+1 );
        covmx[x+nbins][y+nbins]     = CC_e_nom->GetBinContent( x+1, y+1 );
        covmx[x+2*nbins][y+2*nbins] = nue_nom->GetBinContent( x+1, y+1 );
        //std::cout << covmx[x][y] << "\t" << covmx[x+nbins][y+nbins] << "\t" << covmx[x+2*nbins][y+2*nbins] << "\n";
      }

      else{
        covmx[x][y]                     = 0.;
        covmx[x+nbins][y+nbins]         = 0.;
        covmx[x+2*nbins][y+2*nbins]     = 0.;}
    }
  }

  TH2D *hcovmx = new TH2D("hcovmx","",300,0,300,300,0,300);
  for( int i = 0; i < 3*nbins; ++i ) {
    for( int j = 0; j < 3*nbins; ++j ) {
      hcovmx->SetBinContent(i+1, j+1, covmx[i][j]);
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
  TH1D *hscales = new TH1D("hscales", "", 300,0,300);
  //TCanvas *c = new TCanvas("c","",800,800);

  // make random number matrix
  for( int i = 0; i < N; ++i ) {
    //if( i>0 && i%50==0 ) std::cout << "Progress: " << i/N*100. << "..." << "\n";
    double mean = 0.;
    for( int j = 0; j < 3*nbins; ++j ) {
      double val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < 3*nbins; ++j ) {
      double old = scales[i][j];
      scales[i][j] = old - mean;
      hscales->SetBinContent(j+1, scales[i][j]);
    }
  
    //hscales->Draw();
    //c->SaveAs("test.png");

  }

  

  for( int i = 0; i < 3*nbins; ++i ) { // columns
    for( int j = 0; j < 3*nbins; ++j ) { // columns

      // compute column covariance
      double covar = 0.;
      for( int k = 0; k < N; ++k ) { // rows
        covar += (scales[k][i] * scales[k][j]); // means are 0 already by construction
      }
      scaleCovars[i][j] = covar/N;
    }
  }



  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales *= inverse;

  scales *= chol;



  TH1D *m     = new TH1D("m","",100,0,16);
  TH1D *e     = new TH1D("e","",100,0,16);
  TH1D *nue   = new TH1D("nue","",100,0,16);

  for( int u = 0; u < N; ++u ) {
    m->Reset();
    e->Reset();
    nue->Reset();
    
    for(int b=0; b<3*nbins; b++) {
      double evtwgt = scales[u][b];

      if(b<nbins)                 m->Fill(CC_m_nom->GetBinContent(b+1),1.+evtwgt);
      if(b>=nbins && b<2*nbins)   e->Fill(CC_e_nom->GetBinContent(b-nbins+1),1.+evtwgt);
      if(b>=2*nbins && b<3*nbins) nue->Fill(nue_nom->GetBinContent(b-2*nbins+1),1.+evtwgt);
    }
    
    for(int i=0; i<nbins; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
      E_nue[u][i]   = nue->GetBinContent(i+1);
    }
   }


  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
      // compute column covariance
      double covar_m = 0.;
      double covar_e = 0.;
      double covar_nue = 0.;
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
        covar_m     += (E_m[k][i]     - CC_m_nom->GetBinContent(i+1))  * (E_m[k][j]     - CC_m_nom->GetBinContent(j+1));
        covar_e     += (E_e[k][i]     - CC_e_nom->GetBinContent(i+1))  * (E_e[k][j]     - CC_e_nom->GetBinContent(j+1)); 
        covar_nue   += (E_nue[k][i]   - nue_nom->GetBinContent(i+1))   * (E_nue[k][j]   - nue_nom->GetBinContent(j+1)); 

        covar_me   += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1));
        covar_mnue += (E_m[k][i]   - CC_m_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
        covar_em   += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_enue += (E_e[k][i]   - CC_e_nom->GetBinContent(i+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1)); 
        covar_nuem += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_m[k][j]   - CC_m_nom->GetBinContent(j+1)); 
        covar_nuee += (E_nue[k][i] - nue_nom->GetBinContent(i+1))  * (E_e[k][j]   - CC_e_nom->GetBinContent(j+1)); 
/*
        var_m_i += (E_m[k][i] - CC_m_nom->GetBinContent(i+1)) * (E_m[k][i] - CC_m_nom->GetBinContent(i+1));
        var_m_j += (E_m[k][j] - CC_m_nom->GetBinContent(j+1)) * (E_m[k][j] - CC_m_nom->GetBinContent(j+1));

        var_e_i += (E_e[k][i] - CC_e_nom->GetBinContent(i+1)) * (E_e[k][i] - CC_e_nom->GetBinContent(i+1));
        var_e_j += (E_e[k][j] - CC_e_nom->GetBinContent(j+1)) * (E_e[k][j] - CC_e_nom->GetBinContent(j+1));

        var_nue_i += (E_nue[k][i] - nue_nom->GetBinContent(i+1)) * (E_nue[k][i] - nue_nom->GetBinContent(i+1));
        var_nue_j += (E_nue[k][j] - nue_nom->GetBinContent(j+1)) * (E_nue[k][j] - nue_nom->GetBinContent(j+1));
*/
      }

      ECovars_m[i][j]     = covar_m    /(N * CC_m_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      ECovars_e[i][j]     = covar_e    /(N * CC_e_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
      ECovars_nue[i][j]   = covar_nue  /(N * nue_nom->GetBinContent(i+1)   * nue_nom->GetBinContent(j+1));

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



  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {
      ECovars[i][j]         = ECovars_m[i][j];
      ECovars[i][j+nbins]   = ECovars_me[i][j];
      ECovars[i][j+2*nbins] = ECovars_mnue[i][j];

      ECovars[i+nbins][j]         = ECovars_em[i][j];
      ECovars[i+nbins][j+nbins]   = ECovars_e[i][j];
      ECovars[i+nbins][j+2*nbins] = ECovars_enue[i][j];

      ECovars[i+2*nbins][j]         = ECovars_nuem[i][j];
      ECovars[i+2*nbins][j+nbins]   = ECovars_nuee[i][j];
      ECovars[i+2*nbins][j+2*nbins] = ECovars_nue[i][j];

/*
      ECorrel[i][j]           = ECorrel_m[i][j];
      ECorrel[i][j+3*nbins]   = ECorrel_me[i][j];
      ECorrel[i][j+2*3*nbins] = ECorrel_mnue[i][j];

      ECorrel[i+3*nbins][j]           = ECorrel_em[i][j];
      ECorrel[i+3*nbins][j+3*nbins]   = ECorrel_e[i][j];
      ECorrel[i+3*nbins][j+2*3*nbins] = ECorrel_enue[i][j];

      ECorrel[i+2*3*nbins][j]           = ECorrel_nuem[i][j];
      ECorrel[i+2*3*nbins][j+3*nbins]   = ECorrel_nuee[i][j];
      ECorrel[i+2*3*nbins][j+2*3*nbins] = ECorrel_nue[i][j];
*/
    }
  }



  TH2D *hcv = new TH2D("hcv","",300,0,300,300,0,300);
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
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
  hcv->SetTitle(Form("Covariance (nu<%.1fGeV && Etheta2<%.1fMeV)",nu,Ev));

  TCanvas *ccv = new TCanvas("ccv","",900,800);
  hcv->Draw("colz");
  ccv->SaveAs(Form("Cov%d%d_%d.png",cutNu,cutEv,N));

/*
  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/cov_matrix/%s_covmtr%d%d_%d.root",name,cutNu,cutEv,N),"RECREATE");
  hcv->Write();
  hcr->Write();
  out->Close();
  }
  //}
  //}
*/
}

