#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLegend.h>
#include <cmath>
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
    double R = par[0];
    double L = par[1];

    return R / std::sqrt((R * R + W * W * L * L));
}

double V_T(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double R = par[0];
    double C = par[1];

    return R / std::sqrt(R * R + 1 / (W * W * C * C));
}

double V_M(double *f, double *par)
{
    double W = 2 * TMath::Pi() * f[0];
    double R = par[0];
    double RL = par[1];
    double L = par[2];
    double C = par[3];

    return R / std::sqrt((R + RL) * (R + RL) + (W * L - 1 / (W * C)) * (W * L - 1 / (W * C)));
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

/*double sourceTF(double *x, double *par)
{
    double f = x[0]; // qua cosa mettiamo?
    return 0.0;
}

double wooferTF(double *x, double *par)
{
    double f = x[0];
    return R_W * V_0 / std::sqrt((R_LW + R_W) * (R_LW + R_W) + (f * f * L_W * L_W));
    return 0.0;
}

double tweeterTF(double *x, double *par)
{
    double f = x[0];
    return R_T * V_0 / std::sqrt(R_T * R_T + 1 / (f * f * C_T * C_T));
    return 0.0;
}

double midTF(double *x, double *par)
{
    double f = x[0];
    return R_M * V_0 / std::sqrt((R_M + R_LM) * (R_M + R_LM) + (f * L_M - 1 / (f * C_M)) * (f * L_M - 1 / (f * C_M)));
    return 0.0;
}

*/

void plotAmplitude()
{
    // canvas
    TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

    // leggi i file di dati con errori su x e y
    TGraphErrors *g_source = new TGraphErrors("V_source.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_woofer = new TGraphErrors("V_woofer.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_tweeter = new TGraphErrors("V_tweeter.txt", "%lg %lg %lg %lg");
    TGraphErrors *g_mid = new TGraphErrors("V_mid.txt", "%lg %lg %lg %lg");

    // stile e colori dei marker
    g_source->SetMarkerStyle(20);
    g_woofer->SetMarkerStyle(20);
    g_tweeter->SetMarkerStyle(20);
    g_mid->SetMarkerStyle(20);
    g_source->SetMarkerColor(kRed);
    g_woofer->SetMarkerColor(kBlue);
    g_tweeter->SetMarkerColor(kGreen);
    g_mid->SetMarkerColor(kMagenta);

    // titoli e assi
    g_source->SetTitle("Filtro CrossOver;Frequency [Hz];#Delta V [V]");
    g_source->SetMinimum(0);
    g_source->Draw("APE"); // primo grafico in canvas
    g_woofer->Draw("PE same");
    g_tweeter->Draw("PE same");
    g_mid->Draw("PE same");

    // recupera estremo di frequenza per disegnare le funzioni
    double fmin = g_source->GetXaxis()->GetXmin();
    double fmax = g_source->GetXaxis()->GetXmax();

    // TF1 per le funzioni (stile linea, senza marker)
    TF1 *f_S = new TF1("f_S", V_S, fmin, fmax, 9);
    TF1 *f_W = new TF1("f_W", [&, f_S](double *f, double *par)
                       { return f_S->Eval(f[0]) * V_W(f, par); }, fmin, fmax, 2);
    TF1 *f_T = new TF1("f_T", [&, f_S](double *f, double *par)
                       { return f_S->Eval(f[0]) * V_T(f, par); }, fmin, fmax, 2);
    TF1 *f_M = new TF1("f_M", [&, f_S](double *f, double *par)
                       { return f_S->Eval(f[0]) * V_M(f, par); }, fmin, fmax, 4);

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
    f_S->SetParameters(V_0, L_1, C_1, R_1, L_2, C_2, R_L2, R_2, R_gen);
    f_W->SetParameters(R_1, L_1);
    f_T->SetParameters(R_1, C_1);
    f_M->SetParameters(R_2, R_L2, L_2, C_2);
    f_S->SetParNames("V_0", "L_1", "C_1", "R_1", "L_2", "C_2", "R_L2", "R_2", "R_gen");
    f_W->SetParNames("R_1", "L_1");
    f_T->SetParNames("R_1", "C_1");
    f_M->SetParNames("R_2", "R_L2", "L_2", "C_2");
    f_W->SetParLimits(0, R_1 - 5 * 0.55, R_1 + 5 * 0.55);
    f_W->SetParLimits(1, L_1 - 5 * 5.3e-4, L_1 + 5 * 5.3e-4);

    // fit
    g_source->Fit(f_S);
    g_woofer->Fit(f_W);
    g_tweeter->Fit(f_T);
    g_mid->Fit(f_M);

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

void plotPhase()
{
    // canvas
    TCanvas *c = new TCanvas("c", "Crossover Analysis", 800, 600);

    // leggi i file di dati con errori su x e y
    TGraphErrors *p_source = new TGraphErrors("P_source.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_woofer = new TGraphErrors("P_woofer.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_tweeter = new TGraphErrors("P_tweeter.txt", "%lg %lg %lg %lg");
    TGraphErrors *p_mid = new TGraphErrors("P_mid.txt", "%lg %lg %lg %lg");

    // stile e colori dei marker
    p_source->SetMarkerStyle(20);
    p_woofer->SetMarkerStyle(20);
    p_tweeter->SetMarkerStyle(20);
    p_mid->SetMarkerStyle(20);
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
    double fmin = p_source->GetXaxis()->GetXmin();
    double fmax = p_source->GetXaxis()->GetXmax();

    // TF1 per le funzioni (stile linea, senza marker)
    // TF1 *phase_S = new TF1("phase_S", "[0]", fmin, fmax, 1);
    TF1 *phase_W = new TF1("phase_W", p_W, fmin, fmax, 2);
    TF1 *phase_T = new TF1("phase_T", p_T, fmin, fmax, 2);
    TF1 *phase_M = new TF1("phase_M", p_M, fmin, fmax, 4);

    // colori e stile linea
    // phase_S->SetLineColor(kRed + 1);
    // phase_S->SetLineWidth(3);
    phase_W->SetLineColor(kBlue + 1);
    phase_W->SetLineWidth(3);
    phase_T->SetLineColor(kGreen + 1);
    phase_T->SetLineWidth(3);
    phase_M->SetLineColor(kMagenta + 1);
    phase_M->SetLineWidth(3);

    // parametri delle funzioni
    // phase_S->SetParameters(0.);
    phase_W->SetParameters(R_1, L_1);
    phase_T->SetParameters(R_1, C_1);
    phase_M->SetParameters(R_2, R_L2, L_2, C_2);
    phase_W->SetParNames("R_1", "L_1");
    phase_T->SetParNames("R_1", "C_1");
    phase_M->SetParNames("R_2", "R_L2", "L_2", "C_2");

    // fit
    // p_source->Fit(phase_S);
    p_woofer->Fit(phase_W);
    p_tweeter->Fit(phase_T);
    p_mid->Fit(phase_M);

    // disegna le funzioni sullo stesso grafico
    // phase_S->Draw("same");
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
    // leg->AddEntry(phase_S, "source TF", "L");
    leg->AddEntry(phase_W, "woofer TF", "L");
    leg->AddEntry(phase_T, "tweeter TF", "L");
    leg->AddEntry(phase_M, "mid TF", "L");
    leg->Draw();

    c->Update();
}
