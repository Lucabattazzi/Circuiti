#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>
#include <TComplex.h>
#include <TMath.h>

void plotWaveforms(double freq = 19000.)
{
    // canvas
    TCanvas *c = new TCanvas("c", "Forme d'onda", 800, 600);

    TString fname[4] = {"data/wave_source",
                        "data/wave_woofer",
                        "data/wave_tweeter",
                        "data/wave_mid"};
    for (int i{0}; i < 4; ++i)
    {
        fname[i] += freq;
        fname[i] += ".txt";
    }

    // leggi i file di dati con errori su x e y
    TGraphErrors *g_source = new TGraphErrors(fname[0], "%lg %lg %lg");
    TGraphErrors *g_woofer = new TGraphErrors(fname[1], "%lg %lg %lg");
    TGraphErrors *g_tweeter = new TGraphErrors(fname[2], "%lg %lg %lg");
    TGraphErrors *g_mid = new TGraphErrors(fname[3], "%lg %lg %lg");

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
    g_source->SetTitle(";Tempo [s];Tensione [V]");
    g_source->Draw("APE"); // primo grafico in canvas
    g_woofer->Draw("PE same");
    g_tweeter->Draw("PE same");
    g_mid->Draw("PE same");

    // recupera estremo di frequenza per disegnare le funzioni
    double fmin = g_source->GetXaxis()->GetXmin();
    double fmax = g_source->GetXaxis()->GetXmax();
    // double fmax{0.23e-3};

    // TF1 per le funzioni (stile linea, senza marker)
    TF1 *f_S = new TF1("f_S", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x+[2])", 0., fmax);
    TF1 *f_W = new TF1("f_W", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x+[2])", 0., fmax);
    TF1 *f_T = new TF1("f_T", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x+[2])", 0., fmax);
    TF1 *f_M = new TF1("f_M", "[0]*TMath::Sin(2*TMath::Pi()*[1]*x+[2])", 0., fmax);

    // colori e stile linea
    f_S->SetLineColor(kRed + 1);
    f_W->SetLineColor(kBlue + 1);
    f_T->SetLineColor(kGreen + 1);
    f_M->SetLineColor(kMagenta + 1);

    TF1 *f[4] = {f_S, f_W, f_T, f_M};
    double A[4] = {2.5, 1., 2., 1.};
    for (int i = 0; i < 4; ++i)
    {
        f[i]->SetParNames("A", "f", "phi");
        f[i]->SetParameters(A[i], freq, 0.);

        f[i]->SetParLimits(0, 0., 3.);
        f[i]->SetParLimits(1, freq - 1000., freq + 1000.);
        f[i]->SetParLimits(2, -TMath::Pi(), TMath::Pi());

        f[i]->SetLineWidth(3);
        f[i]->SetNpx(10000);
    }

    g_source->Fit(f_S);
    g_woofer->Fit(f_W);
    g_tweeter->Fit(f_T);
    g_mid->Fit(f_M);

    for (int i = 0; i < 4; ++i)
    {
        std::cout << "Function " << f[i]->GetName() << " reduced chi2: " << f[i]->GetChisquare() / f[i]->GetNDF() << '\n';
        // disegna le funzioni sullo stesso grafico
        f[i]->Draw("same");
    }

    // legenda
    TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90);
    leg->SetBorderSize(1);
    leg->AddEntry(g_source, "source", "lep");
    leg->AddEntry(g_woofer, "woofer", "lep");
    leg->AddEntry(g_tweeter, "tweeter", "lep");
    leg->AddEntry(g_mid, "mid", "lep");
    leg->AddEntry(f_S, "source fit", "L");
    leg->AddEntry(f_W, "woofer fit", "L");
    leg->AddEntry(f_T, "tweeter fit", "L");
    leg->AddEntry(f_M, "mid fit", "L");
    leg->Draw();

    c->Update();
}