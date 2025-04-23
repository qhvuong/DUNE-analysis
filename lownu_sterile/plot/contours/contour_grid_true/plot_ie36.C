#include "DUNEStyle.h"

const int Nbins = 100;

void CalculateBinEdges(double val_max, double val_min, double binEdges[Nbins])
{
  double log_min = log10(val_min);
  double log_max = log10(val_max);
  double binWidth = (log_max - log_min) / Nbins;
  for(int i=0; i<Nbins+1; i++){
    binEdges[i] = pow(10, (log_min + i * binWidth));
  }    
}

void plot_ie36()
{
  const int N_samples = 6;
  double ue42, um42, ut42, dm2, bf_ue42, bf_um42, bf_ut42, bf_dm2, chi2, nochi2, dchi2;
  const char out_path[] = "/pnfs/dune/scratch/users/qvuong/output/SenContour/Ut42_0";
  double U_vals[N_samples] = {0.00104563, 0.0104497, 0.0194204, 0.0514304, 0.104431, 0.194082};
  double M_vals[N_samples] = {0.103576, 0.507293, 1.01218, 4.95746, 9.89143, 90.2109};
  int M_bins[N_samples] = {0, 23, 33, 56, 66, 98};
  int U_bins[N_samples] = {26, 52, 59, 70, 78, 85};
  
  double U_max = 0.7;
  double U_min = 0.0001;
  double M_max = 100.;
  double M_min = 0.1;

  double U_binEdges[Nbins+1], M_binEdges[Nbins+1];
  CalculateBinEdges(U_max, U_min, U_binEdges);
  CalculateBinEdges(M_max, M_min, M_binEdges);

  /*
  for(int i=0; i<Nbins; i++){
    for(int j=0; j<6; j++){
      if( ((M_binEdges[i]+M_binEdges[i+1])/2) == M_vals[j] ) M_bins[j] = i;
      if( ((U_binEdges[i]+U_binEdges[i+1])/2) == U_vals[j] ) U_bins[j] = i;
    }
  }

  for(int j=0; j<6; j++){
    std::cout << U_bins[j] << "\t" << M_bins[j] << "\n";}

  
  for(int i=0; i<Nbins; i++){
    std::cout << i << "\t" << (M_binEdges[i]+M_binEdges[i+1])/2 << "\n";
    std::cout << i << "\t" << (U_binEdges[i]+U_binEdges[i+1])/2 << "\n";
  }
  for(int i=0; i<10; i++){
    std::cout << int( i % Nbins ) << "\n";
    std::cout << int( (i % (Nbins * Nbins)) / Nbins ) << "\n";
    std::cout << int( i / (Nbins * Nbins) ) << "\n";
  }
  */

  TH2D *h_me[N_samples];
  TH2D *h_Me[N_samples], *h_Mm[N_samples];

  
  TFile *out = new TFile(Form("contours_ue42.root"), "UPDATE");
  for(int i=3; i<6; i++){
    std::cout << "ie\t" << U_vals[i] << "\n";
    int x = U_bins[i];
    h_Mm[i] = new TH2D(Form("h_Mm_%d",i), Form("U_{e4}^{2} = %.3f", U_vals[i]), Nbins, U_binEdges, Nbins, M_binEdges);
    h_Mm[i]->SetStats(0);

    for(int y=0; y<Nbins; y++){
      if(y%10==0) std::cout << y*1. << " percent...\n";

      for(int z=0; z<Nbins; z++){
        int ifile = x + Nbins * y + Nbins * Nbins * z;
        TString filepath = Form("%s/output_%d.txt",out_path,ifile);
        ifstream f(filepath.Data());
        //if(i%100==0) std::cout << M_vals[iM] << ": \t" << i*100./(Nbins*Nbins) << " percent...\n";
        if(!f) {
          std::cout << "failed: " << x << "\t" << y << "\t" << z << "\t" << ifile << "\n";
          continue;
        }
        else {
          f >> ue42 >> um42 >> ut42 >> dm2 >> nochi2;
          if(ue42 == U_vals[i]){
            //std::cout << "good: " << x << "\t" << y << "\t" << z << "\t" << ifile << "\n";
            //std::cout << "good: " << ue42 << "\t" << um42 << "\t" << dm2 << "\n";
            dchi2 = sqrt(nochi2 - chi2);
            h_Mm[i]->Fill(um42, dm2, dchi2);
          }
          f.close();
        }
      }
    }
    h_Mm[i]->Write();
  }
  out->Close(); 

}

