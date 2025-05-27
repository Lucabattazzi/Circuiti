#include <TCanvas.h>
#include <TComplex.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <cmath>
#include <iostream>

constexpr double R_1{9.99E+02};
constexpr double L_1{5.34E-02};
constexpr double C_1{1.01E-08};
constexpr double R_L2{2.05E+02};
constexpr double R_2{1.50E+03};
constexpr double C_2{1.01E-08};
constexpr double L_2{4.64E-02};
constexpr double V_0{2.5};
constexpr double R_gen{50};
constexpr double delta_R_1{0.55};
constexpr double delta_L_1{5.34E-4};
constexpr double delta_C_1{1.01e-10};
constexpr double delta_R_L2{0.19};
constexpr double delta_R_2{1.67};
constexpr double delta_C_2{1.01e-10};
constexpr double delta_L_2{0.19};

double V_S(double *f, double *par) {
  double V_0{par[0]};

  double W = 2 * TMath::Pi() * f[0];
  TComplex Z_L1{0, W * par[1]};
  TComplex Z_C1{0, -1 / (W * par[2])};
  TComplex Z_R1{par[3], 0};
  TComplex Z_L2{0, W * par[4]};
  TComplex Z_C2{0, -1 / (W * par[5])};
  TComplex Z_RM{par[6] + par[7], 0}; // resistenza totale del midrange
  TComplex Z_gen{par[8], 0.};

  TComplex Z_W = Z_R1 + Z_L1;
  TComplex Z_T = Z_R1 + Z_C1;
  TComplex Z_M = Z_RM + Z_L2 + Z_C2;

  TComplex Z_load = 1. / (1. / Z_W + 1. / Z_M + 1. / Z_T);

  TComplex Z_tot = Z_load + Z_gen;

  return V_0 * (TComplex::Abs(Z_load / (Z_gen + Z_load)));
}

double V_W(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];
  double L = par[1];
  double R = par[3];

  return V_S(f, par) * R / std::sqrt((R * R + W * W * L * L));
}

double V_T(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];
  double C = par[2];
  double R = par[3];

  return V_S(f, par) * R / std::sqrt(R * R + 1 / (W * W * C * C));
}

double V_M(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];

  double L = par[4];
  double C = par[5];
  double RL = par[6];
  double R = par[7];

  return V_S(f, par) * R /
         std::sqrt((R + RL) * (R + RL) +
                   (W * L - 1 / (W * C)) * (W * L - 1 / (W * C)));
}

// Funzioni di fase
double deg(double x) { return x * 180. / TMath::Pi(); }

double p_W(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];
  double R = par[0];
  double L = par[1];

  return deg(-atan2(L * W, R));
}

double p_T(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];
  double R = par[0];
  double C = par[1];

  return deg(atan2(1 / (C * W), R));
}

double p_M(double *f, double *par) {
  double W = 2 * TMath::Pi() * f[0];
  double R = par[0];
  double RL = par[1];
  double L = par[2];
  double C = par[3];

  return deg(-atan2(W * L - 1 / (W * C), R + RL));
}

void plotAmplitude() {
  // canvas
  TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

  // leggi i file di dati con errori su x e y
  TGraphErrors *g_source = new TGraphErrors("V_source.txt", "%lg %lg %lg %lg");
  TGraphErrors *g_woofer = new TGraphErrors("V_woofer.txt", "%lg %lg %lg %lg");
  TGraphErrors *g_tweeter =
      new TGraphErrors("V_tweeter.txt", "%lg %lg %lg %lg");
  TGraphErrors *g_mid = new TGraphErrors("V_mid.txt", "%lg %lg %lg %lg");

  // stile e colori dei marker
  g_source->SetMarkerStyle(1);
  g_woofer->SetMarkerStyle(1);
  g_tweeter->SetMarkerStyle(1);
  g_mid->SetMarkerStyle(1);
  g_source->SetMarkerColor(kRed);
  g_woofer->SetMarkerColor(kBlue);
  g_tweeter->SetMarkerColor(kGreen);
  g_mid->SetMarkerColor(kMagenta);

  // titoli e assi
  g_source->SetTitle("Filtro CrossOver;Frequenza [Hz];Ampiezza [V]");
  g_source->SetMinimum(0);
  g_source->Draw("APE"); // primo grafico in canvas
  g_woofer->Draw("PE same");
  g_tweeter->Draw("PE same");
  g_mid->Draw("PE same");

  // recupera estremo di frequenza per disegnare le funzioni
  // double fmin = g_source->GetXaxis()->GetXmin();
  // double fmax = g_source->GetXaxis()->GetXmax();
  double fmin = 4000.; // ristretto range del fit
  double fmax = 10500.;

  // TF1 per le funzioni (stile linea, senza marker)
  TF1 *f_S = new TF1("f_S", V_S, fmin, fmax, 9);
  TF1 *f_W = new TF1("f_W", V_W, fmin, fmax, 9);
  TF1 *f_T = new TF1("f_T", V_T, fmin, fmax, 9);
  TF1 *f_M = new TF1("f_M", V_M, fmin, fmax, 9);

  // colori e stile linea
  f_S->SetLineColor(kRed + 1);
  f_S->SetLineWidth(3);
  f_W->SetLineColor(kBlue + 1);
  f_W->SetLineWidth(3);
  f_T->SetLineColor(kGreen + 1);
  f_T->SetLineWidth(3);
  f_M->SetLineColor(kMagenta + 1);
  f_M->SetLineWidth(3);

  // parametri delle funzioni
  // TF1 *f[4] = {f_S, f_T};
  // for (int i = 0; i < 2; ++i)
  // {
  //     f[i]->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
  //     f[i]->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2",
  //     "R_2", "R_gen");

  //     f[i]->FixParameter(0, V_0);
  //     f[i]->SetParLimits(1, L_1 - 6 * delta_L_1, L_1 + 6 * delta_L_1);
  //     f[i]->SetParLimits(2, C_1 - 6 * delta_C_1, C_1 + 6 * delta_C_1);
  //     f[i]->SetParLimits(3, R_1 - 6 * delta_R_1, R_1 + 6 * delta_R_1);
  //     f[i]->SetParLimits(4, L_2 - 6 * delta_L_2, L_2 + 6 * delta_L_2);
  //     f[i]->SetParLimits(5, C_2 - 6 * delta_C_2, C_2 + 6 * delta_C_2);
  //     f[i]->SetParLimits(6, R_L2- 6 * delta_R_L2, R_L2+6 * delta_R_L2);
  //     f[i]->SetParLimits(7, R_2 - 6 * delta_R_2, R_2 + 6 * delta_R_2);
  //     f[i]->FixParameter(8, 56.5);
  // }

  TF1 *f[2]= {f_S, f_T};
  for (int i = 0; i < 2; ++i) {
    f[i]->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
    f[i]->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2", "R_2",
                      "R_gen");

    f[i]->FixParameter(0, V_0);
    f[i]->SetParLimits(1, L_1 - 6 * delta_L_1, L_1 + 6 * delta_L_1);
    f[i]->SetParLimits(2, C_1 - 6 * delta_C_1, C_1 + 6 * delta_C_1);
    f[i]->SetParLimits(3, R_1 - 6 * delta_R_1, R_1 + 6 * delta_R_1);
    f[i]->SetParLimits(4, L_2 - 6 * delta_L_2, L_2 + 6 * delta_L_2);
    f[i]->SetParLimits(5, C_2 - 6 * delta_C_2, C_2 + 6 * delta_C_2);
    f[i]->SetParLimits(6, R_L2 - 6 * delta_R_L2, R_L2 + 6 * delta_R_L2);
    f[i]->SetParLimits(7, R_2 - 6 * delta_R_2, R_2 + 6 * delta_R_2);
    f[i]->FixParameter(8, 56.5); 
  }

  // midrange
  f_M->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
  f_M->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2", "R_2",
                   "R_gen");
  f_M->FixParameter(0, V_0);
  f_M->SetParLimits(1, L_1 - 30 * delta_L_1, L_1 + 30 * delta_L_1);
  f_M->SetParLimits(2, C_1 - 30 * delta_C_1, C_1 + 30 * delta_C_1);
  f_M->SetParLimits(3, R_1 - 30 * delta_R_1, R_1 + 30 * delta_R_1);
  f_M->SetParLimits(4, L_2 - 30 * delta_L_2, L_2 + 30 * delta_L_2);
  f_M->SetParLimits(5, C_2 - 30 * delta_C_2, C_2 + 30 * delta_C_2);
  f_M->SetParLimits(6, R_L2 -30 * delta_R_L2, R_L2+30 * delta_R_L2);
  f_M->SetParLimits(7, R_2 - 30 * delta_R_2, R_2 + 30 * delta_R_2);
  f_M->SetParLimits(8, 56.5+10, 56.5-10);

  // woofer
  f_W->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
  f_W->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2", "R_2",
                   "R_gen");
  f_W->FixParameter(0, V_0);
  f_W->SetParLimits(1, L_1 - 30 * delta_L_1, L_1 + 30 * delta_L_1);
  f_W->SetParLimits(2, C_1 - 30 * delta_C_1, C_1 + 30 * delta_C_1);
  f_W->SetParLimits(3, R_1 - 30 * delta_R_1, R_1 + 30 * delta_R_1);
  f_W->SetParLimits(4, L_2 - 30 * delta_L_2, L_2 + 30 * delta_L_2);
  f_W->SetParLimits(5, C_2 - 30 * delta_C_2, C_2 + 30 * delta_C_2);
  f_W->SetParLimits(6, R_L2 -30 * delta_R_L2, R_L2+30 * delta_R_L2);
  f_W->SetParLimits(7, R_2 - 30 * delta_R_2, R_2 + 30 * delta_R_2);
  f_W->SetParameter(8, 56.5); // escludere midrange

  //     TF1 *f[4] = {f_S, f_W, f_T, f_M};
  // for (int i = 0; i < 4; ++i)
  // {
  //     f[i]->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
  //     f[i]->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2",
  //     "R_2", "R_gen");

  //     f[i]->FixParameter(0, V_0);
  //     f[i]->SetParLimits(1, L_1 - 6 * delta_L_1, L_1 + 6 * delta_L_1);
  //     f[i]->SetParLimits(2, C_1 - 6 * delta_C_1, C_1 + 6 * delta_C_1);
  //     f[i]->SetParLimits(3, R_1 - 6 * delta_R_1, R_1 + 6 * delta_R_1);
  //     f[i]->SetParLimits(4, L_2 - 6 * delta_L_2, L_2 + 6 * delta_L_2);
  //     f[i]->SetParLimits(5, C_2 - 6 * delta_C_2, C_2 + 6 * delta_C_2);
  //     f[i]->SetParLimits(6, R_L2- 6 * delta_R_L2, R_L2+6 * delta_R_L2);
  //     f[i]->SetParLimits(7, R_2 - 6 * delta_R_2, R_2 + 6 * delta_R_2);
  //     f[i]->FixParameter(8, 56.5);  // escludere midrange
  // }

  // fit
  std::cout << "Source fit";
  g_source->Fit(f_S, "", "", fmin, fmax);
  std::cout << "Woofer fit";
  g_woofer->Fit(f_W, "", "", fmin, fmax);
  std::cout << "Tweeter fit";
  g_tweeter->Fit(f_T, "", "", fmin, fmax);
  std::cout << "Midrange fit";
  g_mid->Fit(f_M, "", "", fmin, fmax);

  // disegna le funzioni sullo stesso grafico
  f_S->Draw("same");
  f_W->Draw("same");
  f_T->Draw("same");
  f_M->Draw("same");

  // legenda
  TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->AddEntry(g_source, "source", "P");
  leg->AddEntry(g_woofer, "woofer", "P");
  leg->AddEntry(g_tweeter, "tweeter", "P");
  leg->AddEntry(g_mid, "mid", "P");
  leg->AddEntry(f_S, "source TF", "L");
  leg->AddEntry(f_W, "woofer TF", "L");
  leg->AddEntry(f_T, "tweeter TF", "L");
  leg->AddEntry(f_M, "mid TF", "L");
  leg->Draw();

  c->Update();
}

void plotPhase() {
  // canvas
  TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

  // leggi i file di dati con errori su x e y
  TGraphErrors *p_source = new TGraphErrors("P_source.txt", "%lg %lg %lg %lg");
  TGraphErrors *p_woofer = new TGraphErrors("P_woofer.txt", "%lg %lg %lg %lg");
  TGraphErrors *p_tweeter =
      new TGraphErrors("P_tweeter.txt", "%lg %lg %lg %lg");
  TGraphErrors *p_mid = new TGraphErrors("P_mid.txt", "%lg %lg %lg %lg");

  // stile e colori dei marker
  p_source->SetMarkerStyle(1);
  p_woofer->SetMarkerStyle(1);
  p_tweeter->SetMarkerStyle(1);
  p_mid->SetMarkerStyle(1);
  p_source->SetMarkerColor(kRed);
  p_woofer->SetMarkerColor(kBlue);
  p_tweeter->SetMarkerColor(kGreen);
  p_mid->SetMarkerColor(kMagenta);

  // titoli e assi
  p_source->SetTitle("Filtro CrossOver;Frequency [Hz];Phase [deg]");
  p_source->SetMinimum(-90);
  p_source->SetMaximum(90);
  p_source->Draw("APE"); // primo grafico in canvas
  p_woofer->Draw("PE same");
  p_tweeter->Draw("PE same");
  p_mid->Draw("PE same");

  // recupera estremo di frequenza per disegnare le funzioni
  // double fmin = p_source->GetXaxis()->GetXmin();
  // double fmax = p_source->GetXaxis()->GetXmax();
  double fmin = 4000.; // ristretto range del fit
  double fmax = 10500.;

  // TF1 per le funzioni (stile linea, senza marker)
  TF1 *phase_S = new TF1("phase_S", "pol0", fmin, fmax);
  TF1 *phase_W = new TF1("phase_W", p_W, fmin, fmax, 2);
  TF1 *phase_T = new TF1("phase_T", p_T, fmin, fmax, 2);
  TF1 *phase_M = new TF1("phase_M", p_M, fmin, fmax, 4);

  // colori e stile linea
  phase_S->SetLineColor(kRed + 1);
  phase_S->SetLineWidth(3);
  phase_W->SetLineColor(kBlue + 1);
  phase_W->SetLineWidth(3);
  phase_T->SetLineColor(kGreen + 1);
  phase_T->SetLineWidth(3);
  phase_M->SetLineColor(kMagenta + 1);
  phase_M->SetLineWidth(3);

  // parametri delle funzioni
  phase_S->SetParameters(0.);
  phase_W->SetParameters(R_1, L_1);
  phase_T->SetParameters(R_1, C_1);
  phase_M->SetParameters(R_2, R_L2, L_2, C_2);
  phase_S->SetParNames("offset");
  phase_W->SetParNames("R_1", "L_1");
  phase_T->SetParNames("R_1", "C_1");
  phase_M->SetParNames("R_2", "R_L2", "L_2", "C_2");
  phase_W->SetParLimits(0, R_1 - 2 * delta_R_1, R_1 + 2 * delta_R_1);
  phase_W->SetParLimits(1, L_1 - 30 * delta_L_1, L_1 + 30 * delta_L_1);
  phase_T->SetParLimits(0, R_1 - 2 * delta_R_1, R_1 + 2 * delta_R_1);
  phase_T->SetParLimits(1, C_1 - 1 * delta_C_1, C_1 + 1 * delta_C_1);
  phase_M->SetParLimits(0, R_2 - 8 * delta_R_2, R_2 + 8 * delta_R_2);
  phase_M->SetParLimits(1, R_L2 - 8 * delta_R_L2, R_L2 + 8 * delta_R_L2);
  phase_M->SetParLimits(2, L_2 - 8 * delta_L_2, L_2 + 8 * delta_L_2);
  phase_M->SetParLimits(3, C_2 - 8 * delta_C_2, C_2 + 8 * delta_C_2);
  std::cout << "Source fit";
  p_source->Fit(phase_S, "", "", fmin, fmax);
  std::cout << "Woofer fit";
  p_woofer->Fit(phase_W, "", "", fmin, fmax);
  std::cout << "Tweeter fit fit";
  p_tweeter->Fit(phase_T, "", "", fmin, fmax);
  std::cout << "Midrange fit";
  p_mid->Fit(phase_M, "", "", fmin, fmax);

  // disegna le funzioni sullo stesso grafico
  phase_S->Draw("same");
  phase_W->Draw("same");
  phase_T->Draw("same");
  phase_M->Draw("same");

  // legenda
  TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->AddEntry(p_source, "source", "P");
  leg->AddEntry(p_woofer, "woofer", "P");
  leg->AddEntry(p_tweeter, "tweeter", "P");
  leg->AddEntry(p_mid, "mid", "P");
  leg->AddEntry(phase_S, "source fit", "L");
  leg->AddEntry(phase_W, "woofer fit", "L");
  leg->AddEntry(phase_T, "tweeter fit", "L");
  leg->AddEntry(phase_M, "mid fit", "L");
  leg->Draw();

  c->Update();
}
