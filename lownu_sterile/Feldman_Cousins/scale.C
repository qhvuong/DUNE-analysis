#include <fstream>

static const int N = 1000000; 
static const int nbins_CC = 100;
static const int nbins_nue = 50;
static const int nbins = 2*nbins_CC + nbins_nue;

using namespace std;


TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );
  

void scale()
{
  ofstream fout("scale.txt");

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


  for( int i = 0; i < N; ++i ) {
    for( int j = 0; j < nbins; ++j ) {
      fout << scales[i][j] << "\t";
    }
    fout << "\n";
  }

  fout.close();
}

