void failed()
{
  int j=1;

  ofstream out (Form("NuCut_dm%d.txt",j));

  int N = 100;
  double min = 0.0001;
  double max = 0.1;
  double step = (max-min)/N; 


  for(int i=0; i<N*N; i+=1) {
    ifstream f(Form("/pnfs/dune/scratch/users/qvuong/output/NuCut/100dm%d/output_%d.txt",j,i));

    if(i%100==0) std::cout << i*1./100 << " percent" << "\n";
  
    int bx = int(i%N);
    int by = int(i/N);

    double x = min + bx*step + step/2;
    double y = min + by*step + step/2;
    double z = 6.0;

    if(!f) out << i << "\t" << x << "\t" << y << "\t" << z << "\n";
    //if(!f) out << x << "\t" << y << "\t" << z << "\n";

    f.close();
  }

  out.close();
}
