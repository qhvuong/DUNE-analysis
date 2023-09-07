void failed()
{
  int j=2;

  ofstream out (Form("dm2_%d.txt",j));

  for(int i=0; i<10000; i+=1) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/LvsE/nCov/dm2_%d/output_%d.txt",j,i));

    if(i%100==0) std::cout << i*1./100 << " percent" << "\n";

    if(!f) out << i << "\n";

    f.close();
  }

  out.close();
}
