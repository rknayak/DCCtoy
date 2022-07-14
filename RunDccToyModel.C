/********************************************************************
 This is a Toy Model programme for Disoriented Chiral Condensate (DCC)
 analysis. Nudyn is the statistical observable taken to calculate the
 amount of flutuation in various centrality. The programme aims to study
 the nudyn for both binomial and DCC case.
 *********************************************************************
 Written by - Ranjit Nayak & Claude A. Pruneau
 email - ranjit.nayak@cern.ch, capruneau@gmail.com
 dated - 1st May 2018
 version 2 - Aug 1, 2018
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <TMinuit.h>
#include <fstream>
#include <TLatex.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include "TGraphErrors.h"
#include <TRandom3.h>
#include "THnSparse.h"
#include <TCanvas.h>
#endif

////////////////////////////////////////////
// Helper method to create a legend
////////////////////////////////////////////
TLegend * createLegend(double xLow, double xHigh, double yLow, double yHigh, double textSize, int nCols)
{
  TLegend * l = new TLegend(xLow, yLow, xHigh, yHigh,"","NBNDC");
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->SetTextSize(textSize);
  l->SetNColumns(nCols);
  return l;
}

class KaonEvent
{
public:
  
  double centrality;
  int  k;
  int  k0;
  int  k0s;
  int  kc;
  
  int  k_Binomial;
  int  k0_Binomial;
  int  k0s_Binomial;
  int  kc_Binomial;
  
  int  k_DCC;
  int  k0_DCC;
  int  k0s_DCC;
  int  kc_DCC;
  
  double k0Fraction;
  double k0sFraction;
  double k0Fraction_DCC;
  double k0sFraction_DCC;
  double k0Fraction_Binomial;
  double k0sFraction_Binomial;
  
  KaonEvent()
  {
  cout << "-I- KaonEvent() Creating event instance." << endl;
  }
  
  ~KaonEvent()
  {}
  
  void print(ostream & os)
  {
  cout << "----------------------------------------------------" << endl;
  cout << "        centrality: " << centrality << endl;
  cout << "                 k: " << k << endl;
  cout << "                k0: " << k0 << endl;
  cout << "               k0s: " << k0s << endl;
  cout << "                kc: " << kc << endl;
  cout << "        k_Binomial: " << k_Binomial << endl;
  cout << "       k0_Binomial: " << k0_Binomial  << endl;
  cout << "      k0s_Binomial: " << k0s_Binomial  << endl;
  cout << "       kc_Binomial: " << kc_Binomial  << endl;
  cout << "             k_DCC: " << k_DCC << endl;
  cout << "            k0_DCC: " << k0_DCC  << endl;
  cout << "           k0s_DCC: " << k0s_DCC  << endl;
  cout << "            kc_DCC: " << kc_DCC  << endl;
  cout << "        k0Fraction: " << k0Fraction << endl;
  cout << "       k0sFraction: " << k0sFraction << endl;
  cout << "    k0Fraction_DCC: " << k0Fraction_DCC << endl;
  cout << "   k0sFraction_DCC: " << k0sFraction_DCC << endl;
  }
  
};


class KaonMoments
{
public:
  
  ////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////
  
  TString baseName;
  int    nBinsKmult;
  double minKmult;
  double maxKmult;
  
  int    nEvents;
 
  
  // total kaons sums used to calculate moments
  double sum_k;
  double sum_kSq;
  double sum_k0;
  double sum_k0Sq;
  double sum_k0s;
  double sum_k0sSq;
  double sum_kc;
  double sum_kcSq;
  
  double sum_kck0;
  double sum_kck0s;
  
  // kaons moments - binomial part
  double sum_k_Binomial;
  double sum_kSq_Binomial;
  double sum_k0_Binomial;
  double sum_k0Sq_Binomial;
  double sum_k0s_Binomial;
  double sum_k0sSq_Binomial;
  double sum_kc_Binomial;
  double sum_kcSq_Binomial;
  
  // kaons moments - DCC part
  double sum_k_DCC;
  double sum_kSq_DCC;
  double sum_k0_DCC;
  double sum_k0Sq_DCC;
  double sum_k0s_DCC;
  double sum_k0sSq_DCC;
  double sum_kc_DCC;
  double sum_kcSq_DCC;
  
  
  // moments
  
  // total kaons moments
  double k_avg;
  double kSq_avg;
  double kF_avg;
  double k0_avg;
  double k0Sq_avg;
  double k0F_avg;
  double k0s_avg;
  double k0sSq_avg;
  double k0sF_avg;
  double kc_avg;
  double kcSq_avg;
  double kcF_avg;
  
  double kck0_avg;
  double kck0s_avg;
  
  // kaons moments - binomial part
  double k_Binomial_avg;
  double kSq_Binomial_avg;
  double kF_Binomial_avg;
  double k0_Binomial_avg;
  double k0Sq_Binomial_avg;
  double k0F_Binomial_avg;
  double k0s_Binomial_avg;
  double k0sSq_Binomial_avg;
  double k0sF_Binomial_avg;
  double kc_Binomial_avg;
  double kcSq_Binomial_avg;
  double kcF_Binomial_avg;
  
  // kaons moments - DCC part
  double k_DCC_avg;
  double kSq_DCC_avg;
  double kF_DCC_avg;
  double k0_DCC_avg;
  double k0Sq_DCC_avg;
  double k0F_DCC_avg;
  double k0s_DCC_avg;
  double k0sSq_DCC_avg;
  double k0sF_DCC_avg;
  double kc_DCC_avg;
  double kcSq_DCC_avg;
  double kcF_DCC_avg;
  
  // ratios and nuDyn
  double r00;
  double r0s0s;
  double rchch;
  double rch0;
  double rch0s;
  double nuDyn_ch0;
  double nuDyn_ch0s;
  
  TH1D * h_centrality;
  
  TH1D * h_k;
  TH1D * h_kc;
  TH1D * h_k0;
  TH1D * h_k0s;
  
  TH1D * h_k_DCC;
  TH1D * h_kc_DCC;
  TH1D * h_k0_DCC;
  TH1D * h_k0s_DCC;
  
  TH1D * h_k_Binomial;
  TH1D * h_kc_Binomial;
  TH1D * h_k0_Binomial;
  TH1D * h_k0s_Binomial;
  
  TH2D * h_kck0;
  TH2D * h_kck0s;
  
  ////////////////////////////////////////////
  // CTOR
  ////////////////////////////////////////////
  KaonMoments(const TString & aName,
              int    nBins,
              double min,
              double max)
  {
  baseName = aName;
  nBinsKmult = nBins;
  minKmult   = min;
  maxKmult   = max;
  initialize();
  }
  
  ////////////////////////////////////////////
  // DTOR
  ////////////////////////////////////////////
  ~KaonMoments()
  {
  
  }
  
  ////////////////////////////////////////////
  // Initialize all counters and histograms
  ////////////////////////////////////////////
  void initialize()
  {
  cout << "-I- initialize() Creating KaonMoments instance." << endl;
  nEvents   = 0;
  
  sum_k     = 0.0;
  sum_kSq   = 0.0;
  sum_k0    = 0.0;
  sum_k0Sq  = 0.0;
  sum_k0s   = 0.0;
  sum_k0sSq = 0.0;
  sum_kc    = 0.0;
  sum_kcSq  = 0.0;
  
  sum_kck0   = 0.0;
  sum_kck0s  = 0.0;
  
  // kaons moments - binomial part
  sum_k_Binomial     = 0.0;
  sum_kSq_Binomial   = 0.0;
  sum_k0_Binomial    = 0.0;
  sum_k0Sq_Binomial  = 0.0;
  sum_k0s_Binomial   = 0.0;
  sum_k0sSq_Binomial = 0.0;
  sum_kc_Binomial    = 0.0;
  sum_kcSq_Binomial  = 0.0;
  
  // kaons moments - DCC part
  sum_k_DCC     = 0.0;
  sum_kSq_DCC   = 0.0;
  sum_k0_DCC    = 0.0;
  sum_k0Sq_DCC  = 0.0;
  sum_k0s_DCC   = 0.0;
  sum_k0sSq_DCC = 0.0;
  sum_kc_DCC    = 0.0;
  sum_kcSq_DCC  = 0.0;
  
  cout << "-I- initialize() Instantiating histograms." << endl;

  TH1D::SetDefaultSumw2(true);
  
  TString name;
  
  name = baseName + "centrality";
  h_centrality = new TH1D(name, name, 1000, 0.0, 1000.0);
  
  name  = baseName + "k";
  h_k   = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name  = baseName + "kc";
  h_kc  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name  = baseName + "k0";
  h_k0  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name  = baseName + "k0s";
  h_k0s = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  
  name      = baseName + "k_DCC";
  h_k_DCC   = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name      = baseName + "kc_DCC";
  h_kc_DCC  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name      = baseName + "k0_DCC";
  h_k0_DCC  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name      = baseName + "k0s_DCC";
  h_k0s_DCC = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  
  name           = baseName + "k_Binomial";
  h_k_Binomial   = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name           = baseName + "kc_Binomial";
  h_kc_Binomial  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name           = baseName + "k0_Binomial";
  h_k0_Binomial  = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  name           = baseName + "k0s_Binomial";
  h_k0s_Binomial = new TH1D(name, name, nBinsKmult, minKmult, maxKmult);
  
  name = baseName + "kck0";
  h_kck0 = new TH2D(name, name, nBinsKmult, minKmult, maxKmult, nBinsKmult, minKmult, maxKmult);
  name = baseName + "kck0s";
  h_kck0s = new TH2D(name, name, nBinsKmult, minKmult, maxKmult, nBinsKmult, minKmult, maxKmult);;
  
  cout << "-I- initialize() Done." << endl;

  }
  
  ////////////////////////////////////////////
  // increment the counters and
  // fill histograms
  ////////////////////////////////////////////
  void fill(KaonEvent& kaonEvent)
  {
  nEvents++;
  
  // kaons moments - total
  sum_k      += kaonEvent.k;
  sum_kSq    += kaonEvent.k * kaonEvent.k;
  sum_k0     += kaonEvent.k0;
  sum_k0Sq   += kaonEvent.k0 * kaonEvent.k0;
  sum_k0s    += kaonEvent.k0s;
  sum_k0sSq  += kaonEvent.k0s * kaonEvent.k0s;
  sum_kc     += kaonEvent.kc;
  sum_kcSq   += kaonEvent.kc * kaonEvent.kc;
  
  sum_kck0   += kaonEvent.kc * kaonEvent.k0;
  sum_kck0s  += kaonEvent.kc * kaonEvent.k0s;
  
  // kaons moments - binomial part
  sum_k_Binomial     += kaonEvent.k_Binomial;
  sum_kSq_Binomial   += kaonEvent.k_Binomial * kaonEvent.k_Binomial;
  sum_k0_Binomial    += kaonEvent.k0_Binomial;
  sum_k0Sq_Binomial  += kaonEvent.k0_Binomial * kaonEvent.k0_Binomial;
  sum_k0s_Binomial   += kaonEvent.k0s_Binomial;
  sum_k0sSq_Binomial += kaonEvent.k0s_Binomial * kaonEvent.k0s_Binomial;
  sum_kc_Binomial   += kaonEvent.kc_Binomial;
  sum_kcSq_Binomial += kaonEvent.kc_Binomial * kaonEvent.kc_Binomial;
  
  // kaons moments - DCC part
  sum_k_DCC     += kaonEvent.k_DCC;
  sum_kSq_DCC   += kaonEvent.k_DCC * kaonEvent.k_DCC;
  sum_k0_DCC    += kaonEvent.k0_DCC;
  sum_k0Sq_DCC  += kaonEvent.k0_DCC * kaonEvent.k0_DCC;
  sum_k0s_DCC   += kaonEvent.k0s_DCC;
  sum_k0sSq_DCC += kaonEvent.k0s_DCC * kaonEvent.k0s_DCC;
  sum_kc_DCC    += kaonEvent.kc_DCC;
  sum_kcSq_DCC  += kaonEvent.kc_DCC * kaonEvent.kc_DCC;
  
  h_centrality->Fill(kaonEvent.centrality);
  
  h_k   ->Fill( kaonEvent.k   );
  h_kc  ->Fill( kaonEvent.kc  );
  h_k0  ->Fill( kaonEvent.k0  );
  h_k0s ->Fill( kaonEvent.k0s );
  h_k_DCC   ->Fill( kaonEvent.k_DCC   );
  h_kc_DCC  ->Fill( kaonEvent.kc_DCC  );
  h_k0_DCC  ->Fill( kaonEvent.k0_DCC  );
  h_k0s_DCC ->Fill( kaonEvent.k0s_DCC );
  h_k_Binomial   ->Fill( kaonEvent.k_Binomial   );
  h_kc_Binomial  ->Fill( kaonEvent.kc_Binomial  );
  h_k0_Binomial  ->Fill( kaonEvent.k0_Binomial  );
  h_k0s_Binomial ->Fill( kaonEvent.k0s_Binomial );
  
  h_kck0  ->Fill( kaonEvent.kc,  kaonEvent.k0  );
  h_kck0s ->Fill( kaonEvent.kc,  kaonEvent.k0s );
  
  }
  
  ////////////////////////////////////////////
  // Calculate moments and scale
  // histograms to per event values
  ////////////////////////////////////////////
  void calculateMoments()
  {
  if (nEvents<2)
    {
    cout << "-E- calculateMoments() nEvents<2 -- Cannot calculate moments" << endl;
    return;
    }
  double nev   = nEvents;
  double scale = 1.0/nev;
  
  k_avg      = sum_k     / nev;
  kSq_avg    = sum_kSq   / nev;
  k0_avg     = sum_k0    / nev;
  k0Sq_avg   = sum_k0Sq  / nev;
  k0s_avg    = sum_k0s   / nev;
  k0sSq_avg  = sum_k0sSq / nev;
  kc_avg     = sum_kc    / nev;
  kcSq_avg   = sum_kcSq  / nev;
  
  k_Binomial_avg      = sum_k_Binomial     / nev;
  kSq_Binomial_avg    = sum_kSq_Binomial   / nev;
  k0_Binomial_avg     = sum_k0_Binomial    / nev;
  k0Sq_Binomial_avg   = sum_k0Sq_Binomial  / nev;
  k0s_Binomial_avg    = sum_k0s_Binomial   / nev;
  k0sSq_Binomial_avg  = sum_k0sSq_Binomial / nev;
  kc_Binomial_avg     = sum_kc_Binomial    / nev;
  kcSq_Binomial_avg   = sum_kcSq_Binomial  / nev;
  
  k_DCC_avg      = sum_k_DCC     / nev;
  kSq_DCC_avg    = sum_kSq_DCC   / nev;
  k0_DCC_avg     = sum_k0_DCC    / nev;
  k0Sq_DCC_avg   = sum_k0Sq_DCC  / nev;
  k0s_DCC_avg    = sum_k0s_DCC   / nev;
  k0sSq_DCC_avg  = sum_k0sSq_DCC / nev;
  kc_DCC_avg     = sum_kc_DCC    / nev;
  kcSq_DCC_avg   = sum_kcSq_DCC  / nev;
  
  // n(n-1) factors
  kF_avg  = kSq_avg - k_avg;
  k0F_avg  = k0Sq_avg - k0_avg;
  k0sF_avg  = k0sSq_avg - k0s_avg;
  kcF_avg  = kcSq_avg - kc_avg;
  
  // n(n-1) factors
  kF_Binomial_avg    = kSq_Binomial_avg   - k_Binomial_avg;
  k0F_Binomial_avg   = k0Sq_Binomial_avg  - k0_Binomial_avg;
  k0sF_Binomial_avg  = k0sSq_Binomial_avg - k0s_Binomial_avg;
  kcF_Binomial_avg   = kcSq_Binomial_avg  - kc_Binomial_avg;
  
  kF_DCC_avg    = kSq_DCC_avg   - k_DCC_avg;
  k0F_DCC_avg   = k0Sq_DCC_avg  - k0_DCC_avg;
  k0sF_DCC_avg  = k0sSq_DCC_avg - k0s_DCC_avg;
  kcF_DCC_avg   = kcSq_DCC_avg  - kc_DCC_avg;
  
  kck0_avg    = sum_kck0  / nev;
  kck0s_avg   = sum_kck0s / nev;
  
  r00   = k0F_avg/k0_avg/k0_avg - 1.0;
  r0s0s = k0sF_avg/k0s_avg/k0s_avg - 1.0;
  rchch = kcF_avg/kc_avg/kc_avg - 1.0;
  rch0  = kck0_avg/kc_avg/k0_avg - 1.0;
  rch0s = kck0s_avg/kc_avg/k0s_avg - 1.0;
  nuDyn_ch0  = rchch + r00 - 2.0*rch0;
  nuDyn_ch0s = rchch + r0s0s - 2.0*rch0s;
  
  h_centrality->Scale(scale);
  
  h_k   ->Scale(scale);
  h_kc  ->Scale(scale);
  h_k0  ->Scale(scale);
  h_k0s ->Scale(scale);
  h_k_DCC   ->Scale(scale);
  h_kc_DCC  ->Scale(scale);
  h_k0_DCC  ->Scale(scale);
  h_k0s_DCC ->Scale(scale);
  h_k_Binomial   ->Scale(scale);
  h_kc_Binomial  ->Scale(scale);
  h_k0_Binomial  ->Scale(scale);
  h_k0s_Binomial ->Scale(scale);
  h_kck0  ->Scale(scale);
  h_kck0s ->Scale(scale);
  }
  
  void setHisto(TH1 * h, int markerStyle, int markerColor, float markerSize, int lineStyle, int lineColor, const TString & xTitle, const TString & yTitle)
  {
  h->SetStats(0);
  
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.8);
  
  h->GetYaxis()->SetTitle(yTitle);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.95);
  
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerColor(markerColor);
  h->SetMarkerSize(markerSize);
  h->SetLineStyle(lineStyle);
  h->SetLineColor(lineColor);
  }
  

  
  ////////////////////////////////////////////
  // Print the moments
  ////////////////////////////////////////////
  void print(ostream & os)
  {
  cout << "----------------------------------------------------" << endl;
  cout << "                 k: " << k_avg << endl;
  cout << "                k0: " << k0_avg << endl;
  cout << "               k0s: " << k0s_avg << endl;
  cout << "                kc: " << kc_avg << endl;
  cout << "        k_Binomial: " << k_Binomial_avg << endl;
  cout << "       k0_Binomial: " << k0_Binomial_avg  << endl;
  cout << "      k0s_Binomial: " << k0s_Binomial_avg  << endl;
  cout << "       kc_Binomial: " << kc_Binomial_avg  << endl;
  cout << "             k_DCC: " << k_DCC_avg << endl;
  cout << "            k0_DCC: " << k0_DCC_avg  << endl;
  cout << "           k0s_DCC: " << k0s_DCC_avg  << endl;
  cout << "            kc_DCC: " << kc_DCC_avg  << endl;
  }
  
  ////////////////////////////////////////////
  // Calculate moments and scale
  // histograms to per event values
  ////////////////////////////////////////////
  void plotHistograms()
  {
  TString plotOption2D = "colz";
  
  setHisto(h_k,  20,1,0.5,1,1, "k", "counts");
  setHisto(h_kc, 21,2,0.5,2,1, "kc", "counts");
  setHisto(h_k0, 22,4,0.5,3,1, "k0", "counts");
  setHisto(h_k0s,23,6,0.5,4,1, "k0s", "counts");
  
  setHisto(h_k_DCC,  20,2,0.5,1,2, "k (DCC)", "counts");
  setHisto(h_kc_DCC, 21,2,0.5,2,2, "kc (DCC)", "counts");
  setHisto(h_k0_DCC, 22,2,0.5,3,2, "k0 (DCC)", "counts");
  setHisto(h_k0s_DCC,23,2,0.5,4,2, "ks (DCC)", "counts");
  
  setHisto(h_k_Binomial,  21,4,0.5,1,4, "k (Binomial)", "counts");
  setHisto(h_kc_Binomial, 22,4,0.5,2,4, "kc (Binomial)", "counts");
  setHisto(h_k0_Binomial, 23,4,0.5,3,4, "k0 (Binomial)", "counts");
  setHisto(h_k0s_Binomial,24,4,0.5,4,4, "k0s (Binomial)", "counts");
  
  TString canvasName;
  TString pdfName;

  canvasName = baseName;
  canvasName += "KaonMultDist";
  TCanvas * c1 = new TCanvas(canvasName,canvasName,5,5,1000,800);
  c1->Divide(1,3);
  c1->cd(1);
  gPad-> SetLogy();
  
  h_k->SetMinimum(1e-7);
  h_k->SetMaximum(5e-1);
  h_k->Draw();
  h_kc->Draw("SAME");
  h_k0->Draw("SAME");
  h_k0s->Draw("SAME");
  
  TLegend * l1 = createLegend(0.5, 0.9, 0.78, 0.85, 0.1, 4);
  l1->AddEntry(h_k,   "K",        "p");
  l1->AddEntry(h_kc,  "K^{c}",    "p");
  l1->AddEntry(h_k0,  "K^{0}",    "p");
  l1->AddEntry(h_k0s, "K_{s}^{0}","p");
  l1->Draw();
  
  c1->cd(2);
  gPad-> SetLogy();
  
  h_k_DCC->SetMinimum(1e-7);
  h_k_DCC->SetMaximum(5e-1);
  h_k_DCC->Draw();
  h_kc_DCC->Draw("SAME");
  h_k0_DCC->Draw("SAME");
  h_k0s_DCC->Draw("SAME");
  TLegend * l2 = createLegend(0.5, 0.9, 0.78, 0.85, 0.1, 4);
  l2->AddEntry(h_k_DCC,   "K_{DCC}",        "p");
  l2->AddEntry(h_kc_DCC,  "K_{DCC}^{c}",    "p");
  l2->AddEntry(h_k0_DCC,  "K_{DCC}^{0} ",   "p");
  l2->AddEntry(h_k0s_DCC, "K_{s,DCC}^{0}",  "p");
  l2->Draw();
  
  c1->cd(3);
  gPad-> SetLogy();
  
  h_k_Binomial->SetMinimum(1e-7);
  h_k_Binomial->SetMaximum(5e-1);
  h_k_Binomial->Draw();
  h_kc_Binomial->Draw("SAME");
  h_k0_Binomial->Draw("SAME");
  h_k0s_Binomial->Draw("SAME");
  TLegend * l3 = createLegend(0.5, 0.9, 0.78, 0.85, 0.1, 4);
  l3->AddEntry(h_k_Binomial,   "K_{b}",        "p");
  l3->AddEntry(h_kc_Binomial,  "K_{b}^{c}",    "p");
  l3->AddEntry(h_k0_Binomial,  "K_{b}^{0} ",   "p");
  l3->AddEntry(h_k0s_Binomial, "K_{s,b}^{0}",  "p");
  l3->Draw();
 
  pdfName = canvasName;
  pdfName += ".pdf";
  c1->Print(pdfName);
  
  canvasName = baseName;
  canvasName += "k0Vskc";
  TCanvas * c2 = new TCanvas(canvasName,canvasName,50,5,800,800);
  gPad-> SetLogz(true);
  h_kck0->Draw(plotOption2D);
  pdfName = canvasName;
  pdfName += ".pdf";
  c2->Print(pdfName);

  canvasName = baseName;
  canvasName += "k0sVskc";
  TCanvas * c3 = new TCanvas(canvasName,canvasName,150,5,800,800);
  gPad-> SetLogz(true);
  h_kck0s->Draw(plotOption2D);
  pdfName = canvasName;
  pdfName += ".pdf";
  c3->Print(pdfName);

  }
  
  
};

class KaonGenerator
{
public:
  
  ////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////
  int    minMult;
  int    maxMult;
  double kaonFraction;
  double dccFraction;
  TRandom *  random;
  
  ////////////////////////////////////////////
  // CTOR
  ////////////////////////////////////////////
  KaonGenerator(int min_Mult,
                int max_Mult,
                double frac_kaon,
                double frac_dcc)
  {
  minMult      = min_Mult;
  maxMult      = max_Mult;
  kaonFraction = frac_kaon;
  dccFraction  = frac_dcc;
  random       = new TRandom3();
  }
  
  ////////////////////////////////////////////
  // DTOR
  ////////////////////////////////////////////
  ~KaonGenerator()
  {
  }
  
  void generate(KaonEvent & event)
  {
  event.centrality = random->Uniform(minMult,maxMult);
  
  // Total kaons in event
  event.k  = random->Binomial(event.centrality,kaonFraction);
  
  if (dccFraction==0.0) // pure binomial case
    {
    // no dcc kaons
    event.k_DCC   = 0.0;
    event.k0_DCC  = 0.0;
    event.k0s_DCC = 0.0;
    event.kc_DCC  = 0.0;
    // binomial generation
    event.k_Binomial   = event.k;
    event.k0_Binomial  = random->Binomial(event.k_Binomial,  0.5); // half the kaons are neutral on average
    event.k0s_Binomial = random->Binomial(event.k0_Binomial, 0.5); // half the neutral kaons are k-short on average
    event.kc_Binomial  = event.k_Binomial - event.k0_Binomial;
    }
  else if (dccFraction==1.0)  // pure DCC case
    {
    // dcc generation
    event.k_DCC   = event.k;
    event.k0Fraction_DCC = random->Uniform(0.0, 1.0); // uniform fraction of neutral kaons
    event.k0_DCC  =  event.k_DCC * event.k0Fraction_DCC;
    event.k0s_DCC =  random->Binomial(event.k0_DCC, 0.5);; // half the neutral kaons are k-short on average
    event.kc_DCC  =  event.k_DCC - event.k0_DCC;
    // no binomial kaons
    event.k_Binomial   = 0.0;
    event.k0_Binomial  = 0.0;
    event.k0s_Binomial = 0.0;
    event.kc_Binomial  = 0.0;
    
    }
  else // dccFraction is between 0 and 1
    {
    // dcc generation
    event.k_DCC    = dccFraction * event.k;
    event.k0Fraction_DCC = random->Uniform(0.0, 1.0); // uniform fraction of neutral kaons
    event.k0_DCC   =  event.k_DCC * event.k0Fraction_DCC;
    event.k0s_DCC  =  random->Binomial(event.k0_DCC, 0.5);; // half the neutral kaons are k-short on average
    event.kc_DCC   =  event.k_DCC - event.k0_DCC;
    // binomial generation
    event.k_Binomial   = event.k - event.k_DCC;
    event.k0_Binomial  = random->Binomial(event.k_Binomial,  0.5); // half the kaons are neutral on average
    event.k0s_Binomial = random->Binomial(event.k0_Binomial, 0.5); // half the neutral kaons are k-short on average
    event.kc_Binomial  = event.k_Binomial - event.k0_Binomial;
    }
  

    event.k0  = event.k0_Binomial  + event.k0_DCC;
    event.k0s = event.k0s_Binomial + event.k0s_DCC;
    event.kc  = event.kc_Binomial  + event.kc_DCC;

  
  event.k0Fraction  = double(event.k0)/double(event.k);
  event.k0sFraction = double(event.k0s)/double(event.k);
  event.k0Fraction_DCC  = double(event.k0_DCC)/double(event.k_DCC);
  event.k0sFraction_DCC = double(event.k0s_DCC)/double(event.k_DCC);
  event.k0Fraction_Binomial  = double(event.k0_Binomial)/double(event.k_Binomial);
  event.k0sFraction_Binomial = double(event.k0s_Binomial)/double(event.k_Binomial);
  
  }
};

class KaonDccSimulator
{
public:
  
  ////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////
  KaonEvent     * kaonEvent;
  KaonGenerator * kaonGenerator;
  KaonMoments   * kaonMoments;
  bool debug;
  
  ////////////////////////////////////////////
  // CTOR
  ////////////////////////////////////////////
  KaonDccSimulator(double kaonFraction = 0.3,
                   double dccFraction = 0.9,
                   int minMult=1,
                   int maxMult=200,
                   bool debugOption=false)
  {
  cout << "-I- KaonDccSimulator(...) Setting up simulator for ..." << endl;
  cout << "    kaonFraction : " << kaonFraction << endl;
  cout << "    dccFraction : " << dccFraction << endl;
  cout << "        minMult : " << minMult << endl;
  cout << "        maxMult : " << maxMult << endl;
  cout << "    debugOption : " << debugOption << endl;

  debug = debugOption;
  
  TString baseName;
  
  baseName += "Kaonf=0.";
  baseName += int(100*kaonFraction);
  baseName += "_";
  
  baseName += "DCCf=0.";
  baseName += int(100*dccFraction);
  baseName += "_";
  baseName += minMult;
  baseName += "M";
  baseName += maxMult;
  baseName += "_";
  
  kaonEvent     = new KaonEvent();
  kaonGenerator = new KaonGenerator(minMult,maxMult, kaonFraction, dccFraction);
  kaonMoments   = new KaonMoments(baseName,400, 0.0, 400.0);
 
  }
  
  ////////////////////////////////////////////
  // DTOR
  ////////////////////////////////////////////
  ~KaonDccSimulator()
  {
  
  }
  
  ////////////////////////////////////////////
  // run simulator
  ////////////////////////////////////////////
  void run(int nEventsRequested)
  {
  cout << "-I- run(int nEventsRequested) Generating nEvents:" << nEventsRequested << endl;
  for (int iEvent=0; iEvent<nEventsRequested; iEvent++)
    {
    kaonGenerator->generate(*kaonEvent);
    kaonMoments->fill(*kaonEvent);
    if (debug) kaonEvent->print(cout);
    }
  
  kaonMoments->calculateMoments();
  kaonMoments->print(cout);
  kaonMoments->plotHistograms();
  cout << "-I- run(int nEventsRequested) Done" << endl;
  }
  
  
};


////////////////////////////////////////////
// Run the toy model
////////////////////////////////////////////
void RunDccToyModel(int nEventsRequested = 1000)
{
  double kaonFraction = 0.3;

  TString baseName = "DccToyModel_";
  TString canvasName;
  TString pdfName;

  KaonDccSimulator ** modelVsMult = new KaonDccSimulator*[5];

  double minMult[5];
  double maxMult[5];
  double centerMult[5];
  double nudync0VsMult[5];
  double nudync0sVsMult[5];

  double rcc[5];
  double r0s0s[5];
  double rc0s[5];

  for (int k=0;k<5;k++)
  {
    minMult[k] = 1 + 200*k;
    maxMult[k] = minMult[k] + 200;
    centerMult[k] = 0.5*(minMult[k]+maxMult[k]);
    modelVsMult[k] = new KaonDccSimulator(kaonFraction, 0.0,  minMult[k], maxMult[k], false);
    modelVsMult[k]->run(nEventsRequested);
    nudync0VsMult[k] = modelVsMult[k]->kaonMoments->nuDyn_ch0;
    nudync0sVsMult[k] = modelVsMult[k]->kaonMoments->nuDyn_ch0s;
    rcc[k]   = modelVsMult[k]->kaonMoments->rchch;
    r0s0s[k] = modelVsMult[k]->kaonMoments->r0s0s;
    rc0s[k]  = modelVsMult[k]->kaonMoments->rch0s;
    cout<<nudync0VsMult[k]<<endl;
  }
  
  baseName = "DccToyModel_";

  canvasName = baseName;
  canvasName += "nudync0VsMult";
  pdfName = canvasName + ".pdf";
  TCanvas * c1 = new TCanvas(canvasName,canvasName,5,5,1000,800);
  TGraph * g_nudync0VsMult = new TGraph(5, centerMult, nudync0VsMult);
  TGraph * g_nudync0sVsMult = new TGraph(5, centerMult, nudync0sVsMult);

  g_nudync0VsMult->SetTitle("#nu_{0,c,dyn} vs. Multiplicity");
  g_nudync0VsMult->SetMinimum(-2.0);
  g_nudync0VsMult->SetMaximum(2.0);
  g_nudync0VsMult->GetXaxis()->SetTitle("Multiplicity");
  g_nudync0VsMult->GetXaxis()->SetNdivisions(5);
  g_nudync0VsMult->GetYaxis()->SetTitle("#nu_{dyn}");
  g_nudync0VsMult->GetYaxis()->SetNdivisions(5);
  g_nudync0VsMult->SetLineColor(2);
  g_nudync0VsMult->SetMarkerStyle(20);
  g_nudync0VsMult->SetMarkerColor(2);
  g_nudync0VsMult->SetMarkerSize(0.99);

  g_nudync0sVsMult->SetLineColor(4);
  g_nudync0sVsMult->SetMarkerStyle(21);
  g_nudync0sVsMult->SetMarkerColor(4);
  g_nudync0sVsMult->SetMarkerSize(0.99);

  g_nudync0VsMult->Draw("APL");
  g_nudync0sVsMult->Draw("PL");
  TLegend * l = createLegend(0.5, 0.9, 0.78, 0.85, 0.05, 2);
  l->AddEntry(g_nudync0VsMult,  "#nu_{c0,dyn}",  "p");
  l->AddEntry(g_nudync0sVsMult, "#nu_{c0s,dyn}", "p");
  l->Draw();
  c1->Print(pdfName);

  canvasName = baseName;
  canvasName += "RVsMult";
  pdfName = canvasName + ".pdf";
  TCanvas * c1a = new TCanvas(canvasName,canvasName,5,5,1000,800);
  TGraph * g_rccVsMult   = new TGraph(5, centerMult, rcc);
  TGraph * g_r0s0sVsMult = new TGraph(5, centerMult, r0s0s);
  TGraph * g_rc0sVsMult  = new TGraph(5, centerMult, rc0s);

  g_rccVsMult->SetTitle("R vs. Multiplicity");
  g_rccVsMult->SetMinimum(-2.0);
  g_rccVsMult->SetMaximum( 2.0);
  g_rccVsMult->GetXaxis()->SetTitle("Multiplicity");
  g_rccVsMult->GetXaxis()->SetNdivisions(5);
  g_rccVsMult->GetYaxis()->SetTitle("R");
  g_rccVsMult->GetYaxis()->SetNdivisions(5);
  g_rccVsMult->SetLineColor(2);
  g_rccVsMult->SetMarkerStyle(22);
  g_rccVsMult->SetMarkerColor(2);
  g_rccVsMult->SetMarkerSize(0.99);

  g_r0s0sVsMult->SetLineColor(4);
  g_r0s0sVsMult->SetMarkerStyle(23);
  g_r0s0sVsMult->SetMarkerColor(4);
  g_r0s0sVsMult->SetMarkerSize(0.99);

  g_rc0sVsMult->SetLineColor(1);
  g_rc0sVsMult->SetMarkerStyle(24);
  g_rc0sVsMult->SetMarkerColor(1);
  g_rc0sVsMult->SetMarkerSize(0.99);

  g_rccVsMult->Draw("APL");
  g_r0s0sVsMult->Draw("PL");
  g_rc0sVsMult->Draw("PL");

  TLegend * l1 = createLegend(0.5, 0.9, 0.78, 0.85, 0.05, 3);
  l1->AddEntry(g_rccVsMult,   "R_{cc}",    "p");
  l1->AddEntry(g_r0s0sVsMult, "R_{0s0s}",  "p");
  l1->AddEntry(g_rc0sVsMult,  "R_{c0s}",   "p");
  l1->Draw();

  c1a->Print(pdfName);


  // plots vs dcc fraction


  int nFraction = 11;
  KaonDccSimulator ** modelVsDccF = new KaonDccSimulator*[nFraction];
  double * dccFraction    = new double[nFraction];
  double * nudync0VsDccF  = new double[nFraction];
  double * nudync0sVsDccF = new double[nFraction];
  for (int k=0;k<nFraction;k++)
  {
    dccFraction[k] = 0.1*k;
    modelVsDccF[k] = new KaonDccSimulator(kaonFraction, dccFraction[k],  800, 1000, false);
    modelVsDccF[k]->run(nEventsRequested);
    nudync0VsDccF[k]  = modelVsDccF[k]->kaonMoments->nuDyn_ch0;
    nudync0sVsDccF[k] = modelVsDccF[k]->kaonMoments->nuDyn_ch0s;
    cout << " fraction: " << dccFraction[k] << "  nudync0VsMult:" << nudync0VsDccF[k] << "  nudync0sVsMult:" << nudync0sVsDccF[k] << endl;
  }

  canvasName = baseName;
  canvasName += "nudync0VsDccFraction";
  pdfName = canvasName + ".pdf";
  TCanvas * c2 = new TCanvas(canvasName,canvasName,5,5,1000,800);
  TGraph * g_nudync0VsDccF  = new TGraph(10, dccFraction, nudync0VsDccF);
  TGraph * g_nudync0sVsDccF = new TGraph(10, dccFraction, nudync0sVsDccF);

  g_nudync0VsDccF->SetTitle("#nu_{0,c,dyn} vs. DCC Fraction");
   g_nudync0VsDccF->SetMinimum(0.0);
   g_nudync0VsDccF->SetMaximum(2.0);
   g_nudync0VsDccF->GetXaxis()->SetTitle("DCC Fraction");
  g_nudync0VsDccF->GetXaxis()->SetNdivisions(5);
  g_nudync0VsDccF->GetYaxis()->SetTitle("#nu_{dyn}");
  g_nudync0VsDccF->GetYaxis()->SetNdivisions(5);
  g_nudync0VsDccF->SetLineColor(2);
  g_nudync0VsDccF->SetMarkerStyle(20);
  g_nudync0VsDccF->SetMarkerColor(2);
  g_nudync0VsDccF->SetMarkerSize(0.99);

  g_nudync0sVsDccF->SetLineColor(4);
  g_nudync0sVsDccF->SetMarkerStyle(21);
  g_nudync0sVsDccF->SetMarkerColor(4);
  g_nudync0sVsDccF->SetMarkerSize(0.99);

  g_nudync0VsDccF->Draw("APL");
  g_nudync0sVsDccF->Draw("PL");

  TLegend * l2 = createLegend(0.5, 0.9, 0.78, 0.85, 0.05, 2);
  l2->AddEntry(g_nudync0VsDccF,  "#nu_{c0,dyn}",  "p");
  l2->AddEntry(g_nudync0sVsDccF, "#nu_{c0s,dyn}", "p");
  l2->Draw();
  c2->Print(pdfName);


}




