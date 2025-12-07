#include <bits/stdc++.h>
using namespace std;


struct Params {
    double L = 10.0;
    int    N = 51;
    double dt = 10.0;
    double tmax = 1e4;
    double D = 0.1;

    double h = 0.0;
    double hw = 0.0;
    double T_inf = 273.0;

    double xa = 8.0, xb = 8.8;
    double ya = 2.0, yb = 4.0;
    double Smax = 2.0;

    double ywin1 = 6.0, ywin2 = 9.0;

    double xc = 2.0, yc = 8.0;
    double Tlow = 293.0;
    double Thigh = 1e4;

    int Kmax = 30;
    double TOL = 1e-8;

    string tag = "case";
};

static string getArg(int argc, char** argv, const string& k)
{
    for (int i = 1; i < argc; i++) {
        string s(argv[i]);
        auto p = s.find('=');
        if (p != string::npos && s.substr(0, p) == k)
            return s.substr(p + 1);
    }
    return "";
}

inline int idx(int i, int j, int N) { return i + N * j; }


void build_source(const Params& P, double dx, vector<double>& S) {
    int N = P.N;
    S.assign(N * N, 0.0);
    int ia = int(round(P.xa / dx)), ib = int(round(P.xb / dx));
    int ja = int(round(P.ya / dx)), jb = int(round(P.yb / dx));
    for (int j = ja; j <= jb && j < N; j++)
        for (int i = ia; i <= ib && i < N; i++)
            S[idx(i, j, N)] = P.Smax;
}


void apply_BC(const Params& P, double dx, vector<double>& T) {
    int N = P.N;

    auto alpha_from_h = [&](double hval) {
        return hval * dx / P.D; // h * Δ / D
    };

    // WEST (i = 0)
    for (int j = 1; j <= N - 2; j++) {
        double alpha = alpha_from_h(P.h);
        int k0 = idx(0, j, N), k1 = idx(1, j, N);
        T[k0] = (alpha * P.T_inf + T[k1]) / (1.0 + alpha);
    }

    // NORTH (j = N-1)
    for (int i = 1; i <= N - 2; i++) {
        double alpha = alpha_from_h(P.h);
        int kN = idx(i, N - 1, N), k2 = idx(i, N - 2, N);
        T[kN] = (alpha * P.T_inf + T[k2]) / (1.0 + alpha);
    }

    // EAST (i = N-1), window vs wall
    for (int j = 1; j <= N - 2; j++) {
        double y = j * dx;
        bool window = (y >= P.ywin1 && y <= P.ywin2);
        double h_use = window ? P.hw : P.h;
        double alpha = alpha_from_h(h_use);

        int kR = idx(N - 1, j, N);
        int k2 = idx(N - 2, j, N);
        T[kR] = (alpha * P.T_inf + T[k2]) / (1.0 + alpha);
    }

    // SOUTH (j = 0)
    for (int i = 1; i <= N - 2; i++) {
        double alpha = alpha_from_h(P.h);
        int k0 = idx(i, 0, N), k1 = idx(i, 1, N);
        T[k0] = (alpha * P.T_inf + T[k1]) / (1.0 + alpha);
    }

    // CORNERS (simple averages of neighbours)
    T[idx(0,      0,      N)] = 0.5 * (T[idx(0,      1,      N)] + T[idx(1,      0,      N)]);
    T[idx(0,      N - 1,  N)] = 0.5 * (T[idx(0,      N - 2,  N)] + T[idx(1,      N - 1,  N)]);
    T[idx(N - 1,  N - 1,  N)] = 0.5 * (T[idx(N - 2,  N - 1,  N)] + T[idx(N - 1,  N - 2,  N)]);
    T[idx(N - 1,  0,      N)] = 0.5 * (T[idx(N - 2,  0,      N)] + T[idx(N - 1,  1,      N)]);
}


// BUILD RHS R_ij  (Eq.21)

void build_R(const Params& P, double dx,
             const vector<double>& T,
             const vector<double>& S, int w,
             vector<double>& R)
{
    int N = P.N;
    R.assign(N * N, 0.0);

    double coef = P.D * P.dt / (2.0 * dx * dx);

    for (int j = 1; j <= N - 2; j++) {
        for (int i = 1; i <= N - 2; i++) {
            int k = idx(i, j, N);
            double lap = (T[idx(i + 1, j, N)] + T[idx(i - 1, j, N)]
                        + T[idx(i, j + 1, N)] + T[idx(i, j - 1, N)]
                        - 4.0 * T[k]);
            R[k] = T[k] + coef * lap + 0.5 * P.dt * w * S[k];
        }
    }
}


// GAUSS–SEIDEL CN STEP (Eq.22–23, with source)

double gauss_seidel_CN(const Params& P, double dx,
                       vector<double>& T,
                       const vector<double>& S,
                       const vector<double>& R,
                       int w)
{
    int N = P.N;
    double coef  = P.D * P.dt / (2.0 * dx * dx);
    double denom = 1.0 + 2.0 * P.D * P.dt / (dx * dx);

    double resNorm = 0.0;

    for (int it = 0; it < P.Kmax; it++) {
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 1; i <= N - 2; i++) {
                int k = idx(i, j, N);
                double sumNb = T[idx(i + 1, j, N)] + T[idx(i - 1, j, N)]
                             + T[idx(i, j + 1, N)] + T[idx(i, j - 1, N)];

                // CN update with S^{n+1} term:
                T[k] = (coef * sumNb + 0.5 * P.dt * w * S[k] + R[k]) / denom;
            }
        }

        apply_BC(P, dx, T);

        // residual (optional)
        resNorm = 0.0;
        for (int j = 1; j <= N - 2; j++) {
            for (int i = 1; i <= N - 2; i++) {
                int k = idx(i, j, N);
                double sumNb = T[idx(i + 1, j, N)] + T[idx(i - 1, j, N)]
                             + T[idx(i, j + 1, N)] + T[idx(i, j - 1, N)];
                double lapNew = sumNb - 4.0 * T[k];
                double L = T[k] - coef * lapNew - 0.5 * P.dt * w * S[k];
                double diff = L - R[k];
                resNorm += diff * diff;
            }
        }
        resNorm = sqrt(resNorm) * dx;
        if (resNorm < P.TOL) break;
    }
    return resNorm;
}


// ENERGIES (Eq.27 and Eq.29)

double energy_supplied(const Params& P, double dx,
                       int w, int wold, const vector<double>& S)
{
    double sumS = 0.0;
    for (double s : S) sumS += s;
    return 0.5 * (w + wold) * dx * dx * P.dt * sumS;
}

double energy_window(const Params& P, double dx,
                     const vector<double>& T)
{
    int N = P.N;
    int j1 = int(round(P.ywin1 / dx));
    int j2 = int(round(P.ywin2 / dx));
    double E = 0.0;
    for (int j = j1; j <= j2 && j < N; j++) {
        int k = idx(N - 1, j, N);
        E += P.hw * (T[k] - P.T_inf) * dx * P.dt;
    }
    return E;
}


int main(int argc, char** argv) {
    Params P;

    if (auto s = getArg(argc, argv, "--h");     !s.empty()) P.h     = stod(s);
    if (auto s = getArg(argc, argv, "--hw");    !s.empty()) P.hw    = stod(s);
    if (auto s = getArg(argc, argv, "--Thigh"); !s.empty()) P.Thigh = stod(s);
    if (auto s = getArg(argc, argv, "--Tinf");  !s.empty()) P.T_inf = stod(s);
    if (auto s = getArg(argc, argv, "--tag");   !s.empty()) P.tag   = s;

    int N = P.N;
    double dx = P.L / (N - 1);

    int ic = int(round(P.xc / dx));
    int jc = int(round(P.yc / dx));

    vector<double> T(N * N, P.T_inf);
    vector<double> S, R;

    build_source(P, dx, S);
    apply_BC(P, dx, T);

    
    system("mkdir data");
    system("mkdir figs");

    string sf = "data/sensor_" + P.tag + ".csv";
    ofstream fs(sf);
    fs << "t,T_sensor,E_sup,E_win\n";

    
    vector<double> snaps;
    if (P.tag == "task1" || P.tag == "task2" || P.tag == "task3") {
        snaps = {10.0, 100.0, 1000.0, 10000.0};
    } else if (P.tag == "task4") {
        snaps = {500.0, 1600.0, 2600.0, 4500.0};
    } else {
        snaps = {}; // tasks 5 & 6: no maps
    }

    double t = 0.0;
    int steps = int(P.tmax / P.dt);

    int w = 1, wold = 1;
    bool control_by_sensor = (P.tag != "task6"); // Task6

    for (int n = 0; n <= steps; n++) {
        // (1) R uses old heater state w (at time n)
        build_R(P, dx, T, S, w, R);

        // (2) save old weight
        wold = w;

        // (3) change heater state 
        if (control_by_sensor) {
            double Tsensor = T[idx(ic, jc, N)];
            if (Tsensor < P.Tlow)       w = 1;
            else if (Tsensor > P.Thigh) w = 0;
        } else {
            w = 1; // Task 6: heater always ON
        }

        // (4) Gauss–Seidel using new w (for S^{n+1})
        double r = gauss_seidel_CN(P, dx, T, S, R, w);
        (void)r; // not used further

        // (6) update time
        t = (n + 1) * P.dt;

        // (5) energies at this step
        double Es = energy_supplied(P, dx, w, wold, S);
        double Ew = energy_window(P, dx, T);

        // (7) save time series
        fs << t << "," << T[idx(ic, jc, N)] << "," << Es << "," << Ew << "\n";

        // (8) snapshots for maps
        for (double ts : snaps) {
            if (fabs(t - ts) < 0.5 * P.dt) {
                char buf[256];
                sprintf(buf, "data/field_%s_t%.0f.csv", P.tag.c_str(), ts);
                ofstream ff(buf);
                ff << "x,y,T\n";
                for (int j = 0; j < N; j++) {
                    double y = j * dx;
                    for (int i = 0; i < N; i++) {
                        double x = i * dx;
                        ff << x << "," << y << "," << T[idx(i, j, N)] << "\n";
                    }
                    ff << "\n";
                }
            }
        }
    }

    cerr << "DONE " << P.tag << "\n";
    return 0;
}
