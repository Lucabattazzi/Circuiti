#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>
#include <iostream>
#include <TComplex.h>
#include <TMath.h>

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

double V_S(double *f, double *par)
{
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

double V_W(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double L = par[1];
    double R = par[3];

    return V_S(f, par) * R / std::sqrt((R * R + W * W * L * L));
}

double V_T(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double C = par[2];
    double R = par[3];

    return V_S(f, par) * R / std::sqrt(R * R + 1 / (W * W * C * C));
}

double V_M(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];

    double L = par[4];
    double C = par[5];
    double RL = par[6];
    double R = par[7];

    return V_S(f, par) * R / std::sqrt((R + RL) * (R + RL) + (W * L - 1 / (W * C)) * (W * L - 1 / (W * C)));
}

// Funzioni di fase
double deg(double x)
{
    return x * 180. / TMath::Pi();
}

double p_W(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double R = par[0];
    double L = par[1];

    return deg(-atan2(L * W, R));
}

double p_T(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double R = par[0];
    double C = par[1];

    return deg(atan2(1 / (C * W), R));
}

double p_M(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double R = par[0];
    double RL = par[1];
    double L = par[2];
    double C = par[3];

    return deg(-atan2(W * L - 1 / (W * C), R + RL));
}

void plotAmplitude()
{
    // canvas
    TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

    // leggi i file di dati con errori su x e y
    TGraphErrors *g_source = new TGraphErrors("data/V_source.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_woofer = new TGraphErrors("data/V_woofer.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_tweeter = new TGraphErrors("data/V_tweeter.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_mid = new TGraphErrors("data/V_mid.txt", "%lg %lg %lg %lg");

    // stile e colori dei marker
    g_source->SetMarkerStyle(6);
    g_woofer->SetMarkerStyle(6);
    g_tweeter->SetMarkerStyle(6);
    g_mid->SetMarkerStyle(6);
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
    double fmax = 15000.;

    // TF1 per le funzioni (stile linea, senza marker)
    TF1 *f_S = new TF1("f_S", V_S, fmin, fmax, 9);
    TF1 *f_W = new TF1("f_W", V_W, fmin, fmax, 9);
    TF1 *f_T = new TF1("f_T", V_T, fmin, fmax, 9);
    TF1 *f_M = new TF1("f_M", V_M, fmin, fmax, 9);

    // colori e stile linea
    f_S->SetLineColor(kRed + 1);
    f_W->SetLineColor(kBlue + 1);
    f_T->SetLineColor(kGreen + 1);
    f_M->SetLineColor(kMagenta + 1);
    f_M->SetLineWidth(3);


    TF1 *f[4] = {f_S, f_W, f_T, f_M};
    for (int i = 0; i < 4; ++i)
    {
        f[i]->SetLineWidth(3);
        f[i]->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
        f[i]->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2", "R_2", "R_gen");

        f[i]->FixParameter(0, V_0);
        f[i]->SetParLimits(1, L_1 - 6 * delta_L_1, L_1 + 6 * delta_L_1);
        f[i]->SetParLimits(2, C_1 - 6 * delta_C_1, C_1 + 6 * delta_C_1);
        f[i]->SetParLimits(3, R_1 - 6 * delta_R_1, R_1 + 6 * delta_R_1);
        f[i]->SetParLimits(4, L_2 - 6 * delta_L_2, L_2 + 6 * delta_L_2);
        f[i]->SetParLimits(5, C_2 - 6 * delta_C_2, C_2 + 6 * delta_C_2);
        f[i]->SetParLimits(6, R_L2 - 6 * delta_R_L2, R_L2 + 6 * delta_R_L2);
        f[i]->SetParLimits(7, R_2 - 6 * delta_R_2, R_2 + 6 * delta_R_2);
        // f[i]->SetParLimits(8, 0, 100);
        f[i]->FixParameter(8, 55.2); // escludere midrange
    }

    // fit
    g_source->Fit(f_S, "", "", fmin, fmax);
    g_woofer->Fit(f_W, "", "", fmin, fmax);
    g_tweeter->Fit(f_T, "", "", fmin, fmax);
    g_tweeter->Fit(f_T, "", "", fmin, fmax);
    g_mid->Fit(f_M, "", "", fmin, fmax);

  for (int i = 0; i < 4; ++i)
    {
        std::cout << "Function " << f[i]->GetName() << " reduced chi2: " <<
        f[i]->GetChisquare() / f[i]->GetNDF() << '\n';
    }

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
    leg->AddEntry(f_S, "source fit", "L");
    leg->AddEntry(f_W, "woofer fit", "L");
    leg->AddEntry(f_T, "tweeter fit", "L");
    leg->AddEntry(f_M, "mid fit", "L");
    leg->Draw();

    c->Update();
}

void plotPhase()
{
    // canvas
    TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

    // leggi i file di dati con errori su x e y
    TGraphErrors *p_source = new TGraphErrors("data/P_source.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_woofer = new TGraphErrors("data/P_woofer.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_tweeter = new TGraphErrors("data/P_tweeter.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_mid = new TGraphErrors("data/P_mid.txt", "%lg %lg %lg %lg");

    // stile e colori dei marker
    p_source->SetMarkerStyle(6);
    p_woofer->SetMarkerStyle(6);
    p_tweeter->SetMarkerStyle(6);
    p_mid->SetMarkerStyle(6);
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
    double fmin = 4000.;
    double fmax = 15000.;
    double fmin = 4000.;
    double fmax = 15000.;

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
    phase_W->SetParLimits(0, R_1 - 20 * delta_R_1, R_1 + 20 * delta_R_1);
    phase_W->SetParLimits(1, L_1 - 20 * delta_L_1, L_1 + 20 * delta_L_1);
    phase_T->SetParLimits(0, R_1 - 8 * delta_R_1, R_1 + 8 * delta_R_1);
    phase_T->SetParLimits(1, C_1 - 8 * delta_C_1, C_1 + 8 * delta_C_1);
    phase_M->SetParLimits(0, R_2 - 8 * delta_R_2, R_2 + 8 * delta_R_2);
    phase_M->SetParLimits(1, R_L2 - 8 * delta_R_L2, R_L2 + 8 * delta_R_L2);
    phase_M->SetParLimits(2, L_2 - 8 * delta_L_2, L_2 + 8 * delta_L_2);
    phase_M->SetParLimits(3, C_2 - 8 * delta_C_2, C_2 + 8 * delta_C_2);

    p_source->Fit(phase_S, "", "", fmin, fmax);
    p_woofer->Fit(phase_W, "", "", fmin, fmax);
    p_tweeter->Fit(phase_T, "", "", fmin, fmax);
    p_mid->Fit(phase_M, "", "", fmin, fmax);

    // disegna le funzioni sullo stesso grafico
    phase_S->Draw("same");
    phase_W->Draw("same");
    phase_T->Draw("same");
    phase_M->Draw("same");

    TF1 *phase[4] = {phase_S, phase_W, phase_T, phase_M};
    for (int i = 0; i < 4; ++i)
    {
        std::cout << "Function " << phase[i]->GetName() << " reduced chi2: " <<
        phase[i]->GetChisquare() / phase[i]->GetNDF() << '\n';
    }

    // legenda
    TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->AddEntry(p_source, "source", "lep");
    leg->AddEntry(p_woofer, "woofer", "lep");
    leg->AddEntry(p_tweeter, "tweeter", "lep");
    leg->AddEntry(p_mid, "mid", "lep");
    leg->AddEntry(phase_S, "source fit", "L");
    leg->AddEntry(phase_W, "woofer fit", "L");
    leg->AddEntry(phase_T, "tweeter fit", "L");
    leg->AddEntry(phase_M, "mid fit", "L");
    leg->Draw();

    c->Update();
}

void amplitudeLinearFit() { // misurare frequenza di crossover dal grafico

  TCanvas *c = new TCanvas("c", "Crossover Linear Fit", 800, 600);

  TGraphErrors *g_woofer = new TGraphErrors("V_woofer.txt", "%lg %lg %lg %lg");
  TGraphErrors *g_tweeter =
      new TGraphErrors("V_tweeter.txt", "%lg %lg %lg %lg");

  // stile e colori dei marker
  g_woofer->SetMarkerStyle(1);
  g_tweeter->SetMarkerStyle(1);

  g_woofer->SetMarkerColor(kBlue);
  g_tweeter->SetMarkerColor(kGreen);

  // titoli e assi
  g_woofer->SetTitle("Filtro CrossOver;Frequenza [Hz];Ampiezza [V]");
  g_woofer->GetXaxis()->SetRangeUser(4000, 12000);
  g_woofer->Draw("APE"); // primo grafico in canvas
  g_tweeter->Draw("PE same");

  double fmin = 6900.; // ristretto range del fit
  double fmax = 7700.;

  // TF1 per le funzioni (stile linea, senza marker)
  TF1 *line_W = new TF1("line_W", "[0]+[1]*x", fmin, fmax);
  TF1 *line_T = new TF1("line_T", "[0]+[1]*x", fmin, fmax);

  // colori e stile linea
  line_W->SetLineColor(kBlue + 1);
  line_W->SetLineWidth(3);
  line_T->SetLineColor(kGreen + 1);
  line_T->SetLineWidth(3);

  // fit
  std::cout << "Woofer fit";
  g_woofer->Fit(line_W, "", "", fmin, fmax);
  std::cout << "Tweeter fit";
  g_tweeter->Fit(line_T, "", "", fmin, fmax);

  // la x d'intersezione sarà x = - (a_T - a_W)/(b_T - b_W)
  double a_W = line_W->GetParameter(0); // a intercetta, b slope
  double b_W = line_W->GetParameter(1);
  double a_T = line_T->GetParameter(0);
  double b_T = line_T->GetParameter(1);
  double da_W = line_W->GetParError(0);
  double db_W = line_W->GetParError(1);
  double da_T = line_T->GetParError(0);
  double db_T = line_T->GetParError(1);

  double denominator = b_T - b_W;
  double crossover = -(a_T - a_W) / denominator;

  double dx_da_W = -1. / denominator; // calcolo derivate
  double dx_da_T = 1. / denominator;
  double dx_db_W = -(a_T - a_W) / (denominator * denominator);
  double dx_db_T = (a_T - a_W) / (denominator * denominator);

  double delta_crossover =
      sqrt(pow(dx_da_W * da_W, 2) + pow(dx_da_T * da_T, 2) +
           pow(dx_db_W * db_W, 2) + pow(dx_db_T * db_T, 2));
  std::cout << "\n========== CROSSOVER ==========" << std::endl;
  std::cout << "Frequenza crossover: " << crossover << " Hz" << std::endl;
  std::cout << "Errore su crossover: ±" << delta_crossover << " Hz"
            << std::endl;

  // disegna le funzioni sullo stesso grafico
  line_W->Draw("same");
  line_T->Draw("same");
  TMarker *crossover_marker =
      new TMarker(crossover, line_W->Eval(crossover), 20);
  crossover_marker->SetMarkerColor(kRed);
  crossover_marker->SetMarkerSize(1.2);
  crossover_marker->Draw("same");

  // legenda
  TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->AddEntry(crossover_marker, "Crossover", "P");
  leg->AddEntry(line_W, "Woofer fit", "L");
  leg->AddEntry(line_T, "Tweeter fit", "L");
  leg->Draw();

  c->Update();
}
