void failed()
{
  const char out_path[] = "/pnfs/dune/scratch/users/qvuong/output/SenContour/Ut42_0";
  double dm2[5] = {0.5, 1, 5, 10, 100};
  int N = 10000;
  for(int j=0; j<1; j++){
  std::cout << "file = " << dm2[j] << "\n";
  for(int i=0; i<N; i++){
    ifstream f(Form("%s/fixed23_0%d/output_%d.txt",out_path,j,i));
    if(i%50==0) std::cout << i*100./N << " percent...\n";
    if(!f) std::cout << "failed: " << i << "\n";
  }
  }
}
