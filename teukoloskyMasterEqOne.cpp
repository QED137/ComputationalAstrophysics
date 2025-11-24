
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

// ===========================================================
// COMPLEX STRUCT
// ===========================================================
struct ComplexPsi {
    double re;
    double im;
};

// ===========================================================
// PHYSICAL PARAMETERS
// ===========================================================
struct Params {
    double M = 1.0;
    double a = 0.7;   // spin parameter
    double omega = 0.5;
    int s = -2;
    int l = 2, m = 2;
};

// ===========================================================
// GRID STRUCTURE
// ===========================================================
struct Grid {
    int Nr = 250;
    int Nth = 180;

    double rMin = 1.1;    // Just inside horizon (Kerr–Schild works)
    double rMax = 80.0;

    vector<double> r;         // physical r
    vector<double> rStar;     // tortoise coordinate
    vector<double> theta;
    double drStar, dth;

    vector<vector<ComplexPsi>> psi;
    vector<vector<ComplexPsi>> psiPrev;
};

// ===========================================================
// KERR FUNCTIONS
// ===========================================================
double Delta(double r, const Params& P){
    return r*r - 2*P.M*r + P.a*P.a;
}

double tortoise_dr(double r, const Params& P){
    double D = Delta(r, P);
    return (r*r + P.a*P.a) / D;
}

// Effective Teukolsky potential (s = -2)
// Simplified stable form
double Veff(double r, double th, const Params& P){
    double D = Delta(r, P);
    double s = -2;

    double lambda = P.l*(P.l+1) - s*(s+1);
    return (D / (r*r + P.a*P.a)) * lambda;
}

// ===========================================================
// INITIALIZATION
// ===========================================================
void initialize(Grid &G, const Params &P){
    // Build θ-grid
    G.theta.resize(G.Nth);
    G.dth = M_PI/(G.Nth-1);
    for(int j=0;j<G.Nth;j++){
        G.theta[j] = j*G.dth;
    }

    // Build r-grid + tortoise grid
    G.r.resize(G.Nr);
    G.rStar.resize(G.Nr);
    double dr = (G.rMax - G.rMin)/(G.Nr - 1);

    double rstar_acc = 0;
    for(int i=0;i<G.Nr;i++){
        double r = G.rMin + i*dr;
        G.r[i] = r;

        if(i>0){
            double drs = tortoise_dr(G.r[i], P) * (G.r[i]-G.r[i-1]);
            rstar_acc += drs;
        }
        G.rStar[i] = rstar_acc;
    }

    G.drStar = G.rStar[1] - G.rStar[0];

    // Allocate field arrays
    G.psi.resize(G.Nr, vector<ComplexPsi>(G.Nth));
    G.psiPrev.resize(G.Nr, vector<ComplexPsi>(G.Nth));

    // Initial Gaussian pulse
    double r0 = 20.0;
    double sigma = 3.0;

    for(int i=0;i<G.Nr;i++){
        for(int j=0;j<G.Nth;j++){
            double r = G.r[i];
            double th = G.theta[j];
            double amp = exp(-pow((r-r0)/sigma,2)) * sin(th);

            G.psi[i][j] = {amp, 0.0};
            G.psiPrev[i][j] = {amp, 0.0};
        }
    }

    cout<<"Initialization complete.\n";
}

// ===========================================================
// PDE OPERATORS
// ===========================================================

// Radial 2nd derivative in r* (tortoise)
ComplexPsi radial_term(int i,int j,const Grid& G,const Params&P){
    ComplexPsi out;
    double d2r_re = (G.psi[i+1][j].re - 2*G.psi[i][j].re + G.psi[i-1][j].re)
                      /(G.drStar*G.drStar);
    double d2r_im = (G.psi[i+1][j].im - 2*G.psi[i][j].im + G.psi[i-1][j].im)
                      /(G.drStar*G.drStar);
    out.re = d2r_re;
    out.im = d2r_im;
    return out;
}

// Angular term (full Teukolsky spin-weighted θ operator)
ComplexPsi angular_term(int i,int j,const Grid& G,const Params&P){
    double th = G.theta[j];
    double s = -2;

    ComplexPsi out;
    double d2th_re =
        (G.psi[i][j+1].re - 2*G.psi[i][j].re + G.psi[i][j-1].re) / (G.dth*G.dth);
    double d2th_im =
        (G.psi[i][j+1].im - 2*G.psi[i][j].im + G.psi[i][j-1].im) / (G.dth*G.dth);

    double cot = 1.0/tan(th);
    double corr = s*s*cot*cot - s;

    out.re = d2th_re + corr*G.psi[i][j].re;
    out.im = d2th_im + corr*G.psi[i][j].im;
    return out;
}

// Frame dragging mixing term
ComplexPsi frame_drag(int i,int j,const Grid& G,const Params&P){
    double r = G.r[i];
    double D = Delta(r,P);

    ComplexPsi out;
    out.re = -4 * P.a * (r - P.M)/D * G.psi[i][j].im;
    out.im = +4 * P.a * (r - P.M)/D * G.psi[i][j].re;
    return out;
}

// ===========================================================
// TIME EVOLUTION
// ===========================================================
void evolve(Grid& G, const Params& P,double dt){
    vector<vector<ComplexPsi>> psiNew = G.psi;

    for(int i=1;i<G.Nr-1;i++){
        for(int j=1;j<G.Nth-1;j++){
            auto R = radial_term(i,j,G,P);
            auto Th = angular_term(i,j,G,P);
            auto C = frame_drag(i,j,G,P);

            double V = Veff(G.r[i], G.theta[j], P);

            psiNew[i][j].re =
                2*G.psi[i][j].re - G.psiPrev[i][j].re
                + dt*dt*(R.re + Th.re + C.re - V*G.psi[i][j].re);

            psiNew[i][j].im =
                2*G.psi[i][j].im - G.psiPrev[i][j].im
                + dt*dt*(R.im + Th.im + C.im - V*G.psi[i][j].im);
        }
    }

    // Absorbing outer boundary
    int w = 8;
    for(int k=0;k<w;k++){
        double damp = exp(-pow((w-k)/3.0,2));

        for(int j=0;j<G.Nth;j++){
            psiNew[k][j].re *= damp;
            psiNew[k][j].im *= damp;
            psiNew[G.Nr-1-k][j].re *= damp;
            psiNew[G.Nr-1-k][j].im *= damp;
        }
    }

    // Rotate buffers
    G.psiPrev = G.psi;
    G.psi = psiNew;
}

// ===========================================================
// OUTPUT
// ===========================================================
void dump_csv(const Grid& G,int step){
    string name = "frame_" + to_string(step) + ".csv";
    ofstream f(name);
    f<<"rStar,theta,re,im\n";

    for(int i=0;i<G.Nr;i++){
        for(int j=0;j<G.Nth;j++){
            f<<G.rStar[i]<<","<<G.theta[j]<<","<<G.psi[i][j].re<<","<<G.psi[i][j].im<<"\n";
        }
    }
    f.close();
    cout<<"Saved "<<name<<"\n";
}

// ===========================================================
// MAIN
// ===========================================================
int main(){
    Params P;
    Grid G;

    initialize(G,P);

    double dt = 0.02;
    int steps = 1500;

    for(int t=0;t<steps;t++){
        evolve(G,P,dt);

        if(t%50==0){
            dump_csv(G,t);
            cout<<"t="<<t<<"\n";
        }
    }

    cout<<"Simulation complete.\n";
    return 0;
}
