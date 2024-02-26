#include <fstream>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TAxis.h>

static const int N = 150; 
static const int nbins_CC = 50;
static const int nbins_nue = 20;
static const int nbins = 2*nbins_CC + nbins_nue;

using namespace std;


TMatrixD covmx( nbins, nbins );
TMatrixD sysmx( nbins, nbins );
TMatrixD statmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );
  

void ElepCov()
{

  int cutNu = 3;
  const char name[20] = "wgt_MaCCQE";

  // Get data files
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output_test.root", "READ");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output_test.root", "READ");

  // Get systematic uncertainties
  TFile *f_sys = new TFile(Form("../xS_covmtr/total_sigmtr%d_5sig_test.root",cutNu), "READ");
  TFile *f_fl  = new TFile(Form("../flux_covmtr/flux_covmtr%d_test.root",cutNu),"READ");

  TH2D *hsys   = (TH2D*)f_sys->Get( "hcv" );
  TH2D *hfl    = (TH2D*)f_fl->Get("hcv");

  TH1D *CCm_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name,cutNu));
  TH1D *CCe_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name,cutNu));
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep0");
  
  TH1D *hstat = new TH1D("hstat","",nbins,0,nbins);

  // only need the ND FHC part, which is the first 52 bins probably
  for( int x = 0; x < nbins; ++x ) {
    for( int y = 0; y < nbins; ++y ) {
      sysmx[x][y]  = hfl->GetBinContent(x+1, y+1) + hsys->GetBinContent(x+1, y+1);		//total sys = flux unc + cross section unc
      //sysmx[x][y]  = hfl->GetBinContent(x+1, y+1);
      //sysmx[x][y]  = 0.;
      statmx[x][y] = 0.;

      if(x==y){
        if(x<nbins_CC)                  statmx[x][y] = CCm_nom->GetBinContent( x+1 );
        if(x>=nbins_CC && x<2*nbins_CC) statmx[x][y] = CCe_nom->GetBinContent( x-nbins_CC+1 );
        if(x>=2*nbins_CC)               statmx[x][y] = nue_nom->GetBinContent( x-2*nbins_CC+1 );
        hstat->SetBinContent(x+1, statmx[x][y]);
      }
    }
  }
 
  covmx = statmx + sysmx; 



/*
  gStyle->SetNumberContours(999);
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  hstat->SetStats(0);
  hstat->SetTitle("Statistical Uncertainty");
  hfl->SetStats(0);
  hfl->SetTitle("Flux Uncertainty");
  hsys->SetStats(0);
  hsys->SetTitle("Cross section Uncertainty");

  TAxis *a = hstat->GetXaxis();
  //a->LabelsOption("h");
  a->ChangeLabel(1,-1,0.);
  a->ChangeLabel(2,-1,-1,-1,-1,-1,"numuCC");
  a->ChangeLabel(3,-1,0.);
  a->ChangeLabel(4,-1,0.);
  a->ChangeLabel(3,-1,0.);
  //a->SetBinLabel(25,"numuCC");
  //a->SetBinLabel(75,"nueCC");
  //a->SetBinLabel(110,"nu+e");
  hstat->Draw();
*/


/*
  TCanvas *c = new TCanvas("c","",1800,400);
  c->Divide(3,1);
  c->cd(1);
  gPad->SetLogy();
  hstat->SetMaximum(1e9);
  hstat->Draw();
  c->cd(2);
  gPad->SetLogz();
  hfl->SetMaximum(1e9);
  //hfl->Draw("colz");
  c->cd(3);
  gPad->SetLogz();
  hsys->SetMaximum(1e9);
  //hsys->Draw("colz");
  c->SaveAs("unc.png");
*/


  // MAKE FRACTIONAL COVARIANCE MATRIX
  for( int i = 0; i < nbins; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      //double sum = statmx[i][j] + sysmx[i][j];
        if(i<nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = covmx[i][j]/(CCm_nom->GetBinContent(i+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = covmx[i][j]/(CCm_nom->GetBinContent(i+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = covmx[i][j]/(CCm_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
      

        if(i>=nbins_CC && i<2*nbins_CC){ 
          if(j<nbins_CC)                  covmx[i][j] = covmx[i][j]/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = covmx[i][j]/(CCe_nom->GetBinContent(i-nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = covmx[i][j]/(CCe_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}


        if(i>=2*nbins_CC){
          if(j<nbins_CC)                  covmx[i][j] = covmx[i][j]/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCm_nom->GetBinContent(j+1));
          if(j>=nbins_CC && j<2*nbins_CC) covmx[i][j] = covmx[i][j]/(nue_nom->GetBinContent(i-2*nbins_CC+1) * CCe_nom->GetBinContent(j-nbins_CC+1));
          if(j>=2*nbins_CC)               covmx[i][j] = covmx[i][j]/(nue_nom->GetBinContent(i-2*nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
  }
  }

  ofstream fout("fractional_covmx.txt");
  if(fout){
  for(int i=0; i<nbins; i++){
  for(int j=0; j<nbins; j++){
    fout << covmx[i][j];
    if(j<nbins-1) fout << ", ";}
    fout << "\n";}}
  fout.close();  

/*
  TMatrixDEigen eigencov(covmx);
  TVectorD eiV = eigencov.GetEigenValuesRe();
  eiV.Print();

  double eiVmin=1e-8, eiVmax=6.0;
  int NBins = 50;

  double logMin = TMath::Log10(eiVmin);
  double logMax = TMath::Log10(eiVmax);
  double binWidth = (logMax - logMin) / NBins;

  std::vector<double> binEdges(NBins + 1, 0);
  for (int i = 0; i <= NBins; ++i) {
    binEdges[i] = TMath::Power(10, logMin + i * binWidth);
  }

  TH1D *heiV = new TH1D("heiV","",NBins,&binEdges[0]);

  for(int i=0;i<nbins;i++){
  heiV->Fill(eiV[i]);}
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogx();
  heiV->Draw();


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
  
  TMatrixDEigen eigen(scaleCovars);
  TVectorD eigenValues = eigen.GetEigenValuesRe();
  eigenValues.Print();

  TDecompChol scaleDecomp( scaleCovars );
  if( !scaleDecomp.Decompose() ) printf( "Scale matrix didn't decompolse\n" );
  TMatrixD toInvert = scaleDecomp.GetU();
  TMatrixD inverse = toInvert.Invert();
  scales = scales*inverse;
  scales *= chol;

*/

/*
  TFile *out = new TFile(Form("FC%d_stat_100binsCC.root",cutNu),"RECREATE");
  scales.Write("hscales");
  out->Close();
*/



 

}

