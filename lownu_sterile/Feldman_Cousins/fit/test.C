#include <iostream>
#include <fstream>

void test()
{
  int para=2, cutNu=3;

  int n = 17;
  int N[n];
  double c2, chi2[n];
  N[0]=300; 

  //TGraph *g = new TGraph(n,N,chi2);

  for(int i=0; i<n; i++){
  if(i<8){
  N[i]=N[0]+i*100;
  std::ifstream file (Form("FC_chi2_%d%d_%d.txt",para,cutNu,N[i]));
  file >> c2;
  std::cout << i << "\t" << N[i] << "\t" << c2 << "\n";

  chi2[i] = c2;
  }
  else{
  N[i]=N[i-1]+1000;
  std::ifstream file (Form("FC_chi2_%d%d_%d.txt",para,cutNu,N[i]));
  file >> c2;
  std::cout << i << "\t" << N[i] << "\t" << c2 << "\n";

  chi2[i] = c2;
  }



  }
}
