#include <bits/stdc++.h>
using namespace std;
#define M_PI 3.14159265358979323846

// Parameters

struct Params {
    
    double L   = 1.0;   
    double Q   = 1.0;   
    double T   = 100.0; 
    double c0  = 10.0;  
    int    N   = 301;
    double dx  = 0.0;
    double alpha = 0.5; 
    double dt  = 0.0;
    int    Nt  = 2000;

    
    double gamma = 0.0; 
    double Fext  = 0.0;
    double Omega = 0.0;
    int    Omegak = 0;  
    int    ip   = 0;   
    double uinit = 0.01;
    double xA    = 0.5;
    double vinit = 0.0;
    double sigma = 0.05;
    int    ic_type = 0; 

    
    int    bc   = 0;

    
    int Kmodes  = 0;

    
    double rho_jump = 1.0; 

    string tag = "run";
};

static string getArg(int argc, char** argv, const string& key)
{
    for (int i = 1; i < argc; ++i) {
        string s(argv[i]);
        auto pos = s.find('=');
        if (pos != string::npos && s.substr(0, pos) == key)
            return s.substr(pos+1);
    }
    return "";
}

inline double sqr(double x) { return x*x; }



void build_rho_c(const Params& P, vector<double>& rho, vector<double>& c)
{
    int N = P.N;
    rho.assign(N, P.Q);
    c.assign(N, 0.0);

    for (int i = 0; i < N; ++i) {
        double x = P.dx * i;
        if (P.rho_jump > 1.0 && x > P.L/2.0) {
            rho[i] = P.Q * P.rho_jump; // heavy half
        }
        c[i] = sqrt(P.T / rho[i]);
    }
}




void apply_BC(const Params& P, vector<double>& u)
{
    int N = P.N;
    if (P.bc == 0) {
        
        u[0]     = 0.0;
        u[N-1]   = 0.0;
    } else {
        u[0]     = u[1];
        u[N-1]   = u[N-2];
    }
}



void compute_energies(const Params& P,
                      const vector<double>& rho,
                      const vector<double>& u,
                      const vector<double>& v,
                      double& Ekin, double& Epot)
{
    int N = P.N;
    double dx = P.dx;

    Ekin = 0.0;
    for (int i = 1; i <= N-2; ++i) {
        Ekin += dx * rho[i] * 0.5 * sqr(v[i]);
    }

    Epot = 0.0;
    for (int i = 1; i <= N-1; ++i) {
        double du = (u[i] - u[i-1]) / dx;
        Epot += dx * P.T * 0.5 * du * du;
    }
}


// spectral coefficients b_k, d_k, e_k 

void compute_modes(const Params& P,
                   const vector<double>& u,
                   const vector<double>& v,
                   vector<double>& b,
                   vector<double>& d,
                   vector<double>& e)
{
    int N = P.N;
    int Km = P.Kmodes;
    b.assign(Km+1, 0.0);
    d.assign(Km+1, 0.0);
    e.assign(Km+1, 0.0);

    if (Km == 0) return;

    double norm = 2.0 / (N - 1); // from Eq.(13),(15)

    for (int k = 1; k <= Km; ++k) {
        double bk = 0.0, dk = 0.0;
        for (int i = 1; i <= N-2; ++i) {
            double s = sin(k * M_PI * i / (N - 1));  
            bk += u[i] * s;  
            dk += v[i] * s;
        }
        bk *= norm;
        dk *= norm;
        b[k] = bk;
        d[k] = dk;

        // Eq.(19), 
        double dx = P.dx;
        double term1 = P.T * dx * sqr(sin(k * M_PI / (N - 1))) / sqr(dx);
        double term2 = P.Q * dx;
        e[k] = (N - 1) / 4.0 * (term1 * sqr(bk) + term2 * sqr(dk));
    }
}

// discrete dispersion omega_k – Eq.(20)
double omega_k(const Params& P, int k)
{
    return (P.c0 / P.dx) * sqrt(2.0 * (1.0 - cos(k * M_PI / (P.N - 1))));
}


// initialise displacement and velocity

void init_IC(const Params& P,
             vector<double>& u,
             vector<double>& v)
{
    int N = P.N;
    u.assign(N, 0.0);
    v.assign(N, 0.0);

    if (P.ic_type == 0) {
        // Gaussian packet – Eq.(21,22) 
        for (int i = 0; i < N; ++i) {
            double x = P.dx * i;
            double arg = x - P.xA; // t=0
            double u0 = P.uinit * exp(-arg*arg / (2.0 * P.sigma * P.sigma));
            u[i] = u0;
            v[i] = P.vinit * (x - P.xA) / (P.sigma * P.sigma) * u0;  
        }
    } else {
        //  (Task 5)
        for (int i = 0; i < N; ++i) {
            double x = P.dx * i;
            double u0 = sin(M_PI * x / P.L)
                      + sin(2.0 * M_PI * x / P.L)
                      + sin(3.0 * M_PI * x / P.L);
            u[i] = P.uinit * u0;
            v[i] = 0.0;
        }
    }
    apply_BC(P, u);
}




int main(int argc, char** argv)
{
    ios::sync_with_stdio(false);

    Params P;

    
    if (auto s = getArg(argc, argv, "--Nt");      !s.empty()) P.Nt      = stoi(s);
    if (auto s = getArg(argc, argv, "--gamma");   !s.empty()) P.gamma   = stod(s);
    if (auto s = getArg(argc, argv, "--Fext");    !s.empty()) P.Fext    = stod(s);
    if (auto s = getArg(argc, argv, "--Omega");   !s.empty()) P.Omega   = stod(s);
    if (auto s = getArg(argc, argv, "--Omegak");  !s.empty()) P.Omegak  = stoi(s);
    if (auto s = getArg(argc, argv, "--ip");      !s.empty()) P.ip      = stoi(s);
    if (auto s = getArg(argc, argv, "--bc");      !s.empty()) P.bc      = stoi(s); // 0=fixed,1=open
    if (auto s = getArg(argc, argv, "--uinit");   !s.empty()) P.uinit   = stod(s);
    if (auto s = getArg(argc, argv, "--xA");      !s.empty()) P.xA      = stod(s);
    if (auto s = getArg(argc, argv, "--vinit");   !s.empty()) P.vinit   = stod(s);
    if (auto s = getArg(argc, argv, "--sigma");   !s.empty()) P.sigma   = stod(s);
    if (auto s = getArg(argc, argv, "--Kmodes");  !s.empty()) P.Kmodes  = stoi(s);
    if (auto s = getArg(argc, argv, "--ic");      !s.empty()) P.ic_type = stoi(s);
    if (auto s = getArg(argc, argv, "--rho_jump");!s.empty()) P.rho_jump= stod(s);
    if (auto s = getArg(argc, argv, "--tag");     !s.empty()) P.tag     = s;

   
    P.dx  = P.L / (P.N - 1);
    P.c0  = sqrt(P.T / P.Q);
    P.dt  = P.alpha * P.dx / P.c0;

    if (P.Omegak > 0) {
        P.Omega = omega_k(P, P.Omegak);
    }

    
    vector<double> rho, c, u, v, a, vhalf;
    build_rho_c(P, rho, c);
    init_IC(P, u, v);

    int N = P.N;
    a.assign(N, 0.0);
    vhalf.assign(N, 0.0);

    // --- output files ---
    string fe_name = "energies_" + P.tag + ".csv";
    string fm_name = "map_"      + P.tag + ".csv";
    string fmode_name = "modes_" + P.tag + ".csv";

    ofstream fe(fe_name);
    fe << "t,Ekin,Epot,Etot\n";

    ofstream fm(fm_name);
    fm << "t,x,u\n";

    ofstream fmode;
    if (P.Kmodes > 0) {
        fmode.open(fmode_name);
        fmode << "t";
        for (int k = 1; k <= P.Kmodes; ++k) fmode << ",b" << k;
        for (int k = 1; k <= P.Kmodes; ++k) fmode << ",d" << k;
        for (int k = 1; k <= P.Kmodes; ++k) fmode << ",e" << k;
        fmode << ",Etot\n";
    }

    //  time evolution 
    for (int it = 1; it <= P.Nt; ++it) {
        double t = (it - 1) * P.dt;

        // acceleration at time t (Eq.6)
        for (int i = 1; i <= N-2; ++i) {
            double lap = (u[i+1] - 2.0*u[i] + u[i-1]) / (P.dx * P.dx);
            double force = 0.0;
            if (i == P.ip) force = P.Fext * sin(P.Omega * t);
            a[i] = c[i]*c[i]*lap - P.gamma * v[i] + force;
        }

        // v at half 
        for (int i = 1; i <= N-2; ++i) {
            vhalf[i] = v[i] + 0.5 * P.dt * a[i];
        }

        // new u
        for (int i = 1; i <= N-2; ++i) {
            u[i] += vhalf[i] * P.dt;
        }
        apply_BC(P, u);

        // acceleration  t+dt
        double t_next = t + P.dt;
        for (int i = 1; i <= N-2; ++i) {
            double lap = (u[i+1] - 2.0*u[i] + u[i-1]) / (P.dx * P.dx);
            double force = 0.0;
            if (i == P.ip) force = P.Fext * sin(P.Omega * t_next);
            a[i] = c[i]*c[i]*lap - P.gamma * vhalf[i] + force;
        }

        // new v
        for (int i = 1; i <= N-2; ++i) {
            v[i] = vhalf[i] + 0.5 * P.dt * a[i];
        }

        // energies
        double Ek, Ep;
        compute_energies(P, rho, u, v, Ek, Ep);
        double Et = Ek + Ep;
        fe << t_next << "," << Ek << "," << Ep << "," << Et << "\n";

        // map output
        for (int i = 0; i < N; ++i) {
            double x = P.dx * i;
            fm << t_next << "," << x << "," << u[i] << "\n";
        }
        fm << "\n";

        
        if (P.Kmodes > 0) {
            vector<double> b, d, e;
            compute_modes(P, u, v, b, d, e);
            fmode << t_next;
            for (int k = 1; k <= P.Kmodes; ++k) fmode << "," << b[k];
            for (int k = 1; k <= P.Kmodes; ++k) fmode << "," << d[k];
            double Etot_modes = 0.0;
            for (int k = 1; k <= P.Kmodes; ++k) {
                fmode << "," << e[k];
                Etot_modes += e[k];
            }
            fmode << "," << Etot_modes << "\n";
        }
    }

    cerr << "DONE " << P.tag << "\n";
    return 0;
}
