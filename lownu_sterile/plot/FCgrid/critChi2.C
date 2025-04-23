#include <DUNEStyle.h>
#include <TLegend.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}




double GetCriticalChi2(TH1D* h, double confidenceLevel) {
    if (!h) {
        std::cerr << "Histogram is null!" << std::endl;
        return -1;
    }

    // Total number of entries in the histogram
    double totalEntries = h->GetEntries();

    // Desired cumulative number of entries corresponding to the confidence level
    double targetEntries = confidenceLevel * totalEntries;

    // Initialize the cumulative sum
    double cumulativeEntries = 0.0;

    // Loop over the bins to find the bin where the cumulative sum exceeds targetEntries
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        cumulativeEntries += h->GetBinContent(bin);
        
        if (cumulativeEntries >= targetEntries) {
            // Found the critical bin
            double criticalChi2 = h->GetBinLowEdge(bin) + h->GetBinWidth(bin) / 2;
            return criticalChi2;
        }
    }

    // If we didn't find it, return the highest bin center as a fallback
    return h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()) / 2;
}





void critChi2()
{
  TFile *f = new TFile("FC_chi2_500bins.root", "READ");
  TH1D *h = (TH1D*)f->Get("hd3");



  h->Scale(1./h->Integral());

    // Compute the CDF from the PDF using GetCumulative()
    TH1D *hCDF = (TH1D*)h->GetCumulative();

      // Define sigma levels as cumulative probabilities
    const double sigma1 = 0.6827;  // 1 sigma (68.27%)
    const double sigma2 = 0.9;  // 2 sigma (95.45%)
    const double sigma3 = 0.9973;  // 3 sigma (99.73%)

    // Variables to store the x-value where CDF reaches the sigma levels
    double x1sigma = 0, x2sigma = 0, x3sigma = 0;

    // Loop over the bins to find where the CDF crosses the sigma levels
    for (int bin = 1; bin <= hCDF->GetNbinsX(); ++bin) {
        double cumulativeProb = hCDF->GetBinContent(bin);
        double xValue = hCDF->GetBinCenter(bin);

        if (cumulativeProb >= sigma1 && x1sigma == 0) {
            x1sigma = xValue;
        }
        if (cumulativeProb >= sigma2 && x2sigma == 0) {
            x2sigma = xValue;
        }
        if (cumulativeProb >= sigma3 && x3sigma == 0) {
            x3sigma = xValue;
            break;  // Exit loop once 3 sigma is found
        }
    }

    // Output the results
    std::cout << "1 sigma (68.27%) corresponds to x = " << x1sigma << std::endl;
    std::cout << "2 sigma (95.45%) corresponds to x = " << x2sigma << std::endl;
    std::cout << "3 sigma (99.73%) corresponds to x = " << x3sigma << std::endl;
/*
  TCanvas *c = new TCanvas("c", "PDF and CDF", 800, 600);
    c->Divide(1, 2);

    // Draw the PDF
    c->cd(1);
    gPad->SetLogy();
    h->Draw();
    h->SetTitle("PDF");
    // Draw the CDF
    c->cd(2);
    //hCDF->SetMinimum(0.1);
    gPad->SetLogy();
    hCDF->Draw();
    hCDF->SetTitle("CDF");

    // Save the canvas as an image
    c->SaveAs("pdf_and_cdf_500bins.png");
*/
}
