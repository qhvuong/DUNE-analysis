static const int N = 1000000; 
static const int nbins_CC = 100;
static const int nbins_nue = 50;
static const int nbins = 2*nbins_CC + nbins_nue;


TMatrixD covmx( nbins, nbins );
TMatrixD sysmx( nbins, nbins );
TMatrixD statmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );
  
//TH2D *hcv = new TH2D("hcv","",N,0,N,nbins,0,nbins);

void test()
{

  int cutNu = 3;
	
  const char name[20] = "wgt_MaCCQE";

  //for(para = 1; para <3; para++) {
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root", "READ");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root", "READ");
  TFile *f_sys = new TFile(Form("../xS_covmtr/total_sigmtr%d_5sig.root",cutNu), "READ");
  //TFile *f_sys = new TFile(Form("../xS_covmtr/%s_sigmtr%d.root",name,cutNu), "READ");
  TFile *f_fl  = new TFile(Form("../flux_covmtr/flux_covmtr%d_10000.root",cutNu),"READ");

  TH2D *hsys   = (TH2D*)f_sys->Get( "hcv" );
  TH2D *hfl    = (TH2D*)f_fl->Get("hcv");

  TH1D *CCm_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CCe_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep0");
  

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      sysmx[x][y]  = hfl->GetBinContent(x+1, y+1) + hsys->GetBinContent(x+1, y+1);
      //sysmx[x][y]  = hfl->GetBinContent(x+1, y+1);
      statmx[x][y] = 0.;

      if(x==y){
        if(x<nbins_CC)                  statmx[x][y] = CCm_nom->GetBinContent( x+1 );
        if(x>=nbins_CC && x<2*nbins_CC) statmx[x][y] = CCe_nom->GetBinContent( x-nbins_CC+1 );
        if(x>=2*nbins_CC)               statmx[x][y] = nue_nom->GetBinContent( x-2*nbins_CC+1 );
      }

    }
  }
  

  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      double sum = statmx[i][j] + sysmx[i][j];
        if(i<nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(CCm_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
      

        if(i>=nbins_CC && i<2*nbins_CC){ 
          if(j<nbins_CC)                  covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(CCe_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}


        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = sum/(nue_nom->GetBinContent(i-2*nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
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
  
  
  double mean;
  double val, old;
  for( int i = 0; i < N; ++i ) {
    mean = 0.;
    for( int j = 0; j < nbins; ++j ) {
      val = rand->Gaus( 0., 1. );
      scales[i][j] = val;
      mean += val;
    }

    mean /= N;
    // force each column mean to be exactly 0, this just eliminates tiny statistical fluctuations in the mean weight
    for( int j = 0; j < nbins; ++j ) {
      old = scales[i][j];
      scales[i][j] = old - mean;
    }
  }


  double covar;
  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
      // compute column covariance
      covar = 0.;
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
  scales = scales*inverse;
  scales *= chol;

/*
  for(int u=0; u<N; u++){
    for(int i=0; i<nbins; i++){
      hcv->SetBinContent(u+1, i+1, scales[u][i]);
    }
  }
*/

  TFile *out = new TFile(Form("FC%d_%d.root",cutNu,N),"RECREATE");
  scales.Write("scales");
  //hcv->Write();
  out->Close();

/*
  TH1D *m     = new TH1D("m","",nbins_CC,0,16);
  TH1D *e     = new TH1D("e","",nbins_CC,0,16);
  TH1D *nue   = new TH1D("nue","",nbins_nue,0,16);


  for( int u = 0; u < 1; ++u ) {
    //if(u%100==0) std::cout << u << "\n";

    m->Reset();
    e->Reset();
    nue->Reset();

    for(int b=0; b<nbins_CC; b++) {
      int fluxbin = b;
      double evtwgt = scales[u][fluxbin];
      double bincontent = CCm_nom->GetBinContent(b+1);

      m->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }
    for(int b=0; b<nbins_CC; b++) {
      int fluxbin = b+nbins_CC;
      double evtwgt = scales[u][fluxbin];
      double bincontent = CCe_nom->GetBinContent(b+1);

      e->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }
    for(int b=0; b<nbins_nue; b++) {
      int fluxbin = b+2*nbins_CC;
      double evtwgt = scales[u][fluxbin];
      double bincontent = nue_nom->GetBinContent(b+1);

      nue->SetBinContent(b+1, bincontent*(1.+evtwgt));
    }


    for(int i=0; i<nbins_CC; i++){
      E_m[u][i]     = m->GetBinContent(i+1);
      E_e[u][i]     = e->GetBinContent(i+1);
    }
    for(int i=0; i<nbins_nue; i++){
      E_nue[u][i]   = nue->GetBinContent(i+1);
    }

    TCanvas *c = new TCanvas("c", "", 1500,800);
    c->Divide(3,2);
    c->cd(1);
    CCm_nom->Draw();
    c->cd(2);
    CCe_nom->Draw();
    c->cd(3);
    nue_nom->Draw();
    c->cd(4);
    m->Draw();
    c->cd(5);
    e->Draw();
    c->cd(6);
    nue->Draw();
    c->SaveAs("test.png");
  }
*/




 

}

