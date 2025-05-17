#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLegend.h>
#include <cmath>

double R_LW = 9.99E+02;    
double R_W = 9.99E+02;    
double L_W = 5.34E-02;   
double R_T = 9.99E+02;
double C_T = 1.01E-08;
double R_LM = 2.05E+02;
double R_M = 1.50E+03;
double C_M= 1.01E-08;
double L_M= 4.64E-02;
double V_0 = 5;

double sourceTF(double *x, double *par) {
    double f = x[0]; // qua cosa mettiamo?
    return 0.0;
}

double wooferTF(double *x, double *par) {
    double f = x[0];
    return R_W*V_0/std::sqrt((R_LW+R_W)*(R_LW+R_W)+(f*f*L_W*L_W));
    return 0.0;
}

double tweeterTF(double *x, double *par) {
    double f = x[0];
    return R_T*V_0/std::sqrt(R_T*R_T+1/(f*f*C_T*C_T));
    return 0.0;
}
double midTF(double *x, double *par) {
    double f = x[0];
    return R_M*V_0/std::sqrt((R_M+R_LM)*(R_M+R_LM)+(f*L_M-1/(f*C_M))*(f*L_M-1/(f*C_M)));
    return 0.0;
}

void plotCrossOver() {
    // canvas
    TCanvas *c = new TCanvas("c","Crossover Analysis",800,600);
    
    // 3) leggi i file di dati con errori su x e y
    TGraphErrors *g_source  = new TGraphErrors("V_source.txt",  "%lg %lg");
    TGraphErrors *g_woofer  = new TGraphErrors("V_woofer.txt",  "%lg %lg");
    TGraphErrors *g_tweeter = new TGraphErrors("V_tweeter.txt", "%lg %lg");
    TGraphErrors *g_mid     = new TGraphErrors("V_mid.txt",     "%lg %lg");

    // stile e colori dei marker
    g_source->SetMarkerColor(kRed);
    g_woofer->SetMarkerColor(kBlue);
    g_tweeter->SetMarkerColor(kGreen+2);
    g_mid->SetMarkerColor(kMagenta);

    // titoli e assi
    g_source->SetTitle("Filtro CrossOver;Frequency [Hz];#Delta V");
    g_source->Draw("AP");  // primo grafico in canvas
    g_woofer->Draw("P same");
    g_tweeter->Draw("P same");
    g_mid->Draw("P same");

    // recupera estremo di frequenza per disegnare le funzionis
    double fmin = g_source->GetXaxis()->GetXmin();
    double fmax = g_source->GetXaxis()->GetXmax();

    // 5) TF1 per le funzioni (stile linea, senza marker)
    TF1 *f_src  = new TF1("f_src",  sourceTF,  fmin, fmax, 0);
    TF1 *f_woof = new TF1("f_woof", wooferTF,  fmin, fmax, 0);
    TF1 *f_twee = new TF1("f_twee", tweeterTF, fmin, fmax, 0);
    TF1 *f_midf = new TF1("f_midf", midTF,     fmin, fmax, 0);

    // colori e stile linea
    f_src ->SetLineColor(kRed);    f_src ->SetLineWidth(4); f_src ->SetLineStyle(2);
    f_woof->SetLineColor(kBlue);   f_woof->SetLineWidth(4); f_woof->SetLineStyle(2);
    f_twee->SetLineColor(kGreen+2);f_twee->SetLineWidth(4); f_twee->SetLineStyle(2);
    f_midf->SetLineColor(kMagenta);f_midf->SetLineWidth(4); f_midf->SetLineStyle(2);

    // disegna le funzioni sullo stesso grafico
    f_src ->Draw("same");
    f_woof->Draw("same");
    f_twee->Draw("same");
    f_midf->Draw("same");

    // legenda
    TLegend *leg = new TLegend(0.65,0.65,0.90,0.90);
    leg->SetBorderSize(0);
    leg->AddEntry(g_source, "source",  "P");
    leg->AddEntry(g_woofer, "woofer",  "P");
    leg->AddEntry(g_tweeter,"tweeter", "P");
    leg->AddEntry(g_mid,    "mid",     "P");
    leg->AddEntry(f_src,    "source TF","L");
    leg->AddEntry(f_woof,   "woofer TF","L");
    leg->AddEntry(f_twee,   "tweeter TF","L");
    leg->AddEntry(f_midf,   "mid TF",   "L");
    leg->Draw();

    c->Update();
}
