void failed()
{
  int N = 5000;
  ofstream out (Form("FC_%d.txt",N));

  for(int i=0; i<N; i+=1) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/FC/output_%d.txt",i));

    if(i%100==0) std::cout << i*100./N << " percent" << "\n";

    if(!f) out << i << "\n";

    f.close();
  }

  out.close();
}
