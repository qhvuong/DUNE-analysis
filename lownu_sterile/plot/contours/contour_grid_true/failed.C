void failed()
{
  const char out_path[] = "/pnfs/dune/scratch/users/qvuong/output/contour_4pars";
  int dm2[4] = {1, 5, 10, 100};
  int N = 10000;
  for(int j=3; j<4; j++){
  std::cout << "file = " << dm2[j] << "\n";
  for(int i=0; i<N; i++){
    ifstream f(Form("%s/dm%d_true/output_%d.txt",out_path,dm2[j],i));
    if(i%50==0) std::cout << i*100./N << " percent...\n";
    if(!f) std::cout << "failed: " << i << "\n";
  }
  }
}
