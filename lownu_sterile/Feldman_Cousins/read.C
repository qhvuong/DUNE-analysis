void read()
{
  TFile *f = new TFile("FC3_5000.root","READ");

  f->Print();

  //std::cout << f->GetListOfKeys();

  //TMatrixD mtr = (TMatrixD<double>)f->Get("scales");
  //TMatrixD mtr;

  TMatrixD* mtr = dynamic_cast<TMatrixD*>(f->Get("scales"));

  std::cout << mtr->GetNrows() << "\t" << mtr->GetNcols() << "\n";

  //if(f->IsOpen()) f->GetObject("scales", mtr);

}
