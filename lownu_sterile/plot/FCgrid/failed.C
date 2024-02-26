void failed()
{
  int N = 1000;
  int seed = 0;
  ofstream out (Form("FC%d.txt",seed));

  int fail=0;
  for(int i=0; i<N; i+=1) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/FC/s%d/output_%d.txt",seed,i));

    if(i%100==0) std::cout << i*100./N << " percent" << "\n";

    if(!f){
      fail+=1;
      out << i << "\n";}

    f.close();
  }

  std::cout << fail << "\n";
  out.close();
}
