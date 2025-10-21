#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

const double g = 9.81;
const double T0 = 293.0;

//  task1 (no drag)  xmax ,tmax
double analytical_tmax(double v0, double theta) {
    return 2.0 * v0 * sin(theta) / g;
}
double analytical_xmax(double v0, double theta) {
    return v0 * v0 * sin(2.0 * theta) / g;
}

// Forward Euler 
double simulate(double dt, double v0, double theta, double D, double a, double m,
                double alpha, const string &filename) {
    double t = 0.0;
    double x = 0.0, y = 0.0;
    double vx = v0 * cos(theta);
    double vy = v0 * sin(theta);

    ofstream fout(filename);
    fout << "t,x,y\n";

    double xold, yold, told;
    bool RUN = true;

    while (RUN) {
        xold = x;
        yold = y;
        told = t;

        double v = sqrt(vx * vx + vy * vy);
        double density_factor = pow((1.0 - a * y / T0), alpha);
        if (density_factor < 0) density_factor = 0; //

        double Fx = -(D / m) * v * vx * density_factor;
        double Fy = -(D / m) * v * vy * density_factor - g;

        x += dt * vx;
        y += dt * vy;
        vx += dt * Fx;
        vy += dt * Fy;
        t += dt;

        if (y < 0.0) {
            double r = yold / (yold - y);
            x = xold + (x - xold) * r;
            y = 0.0;
            t = told + r * dt;
            RUN = false;
        }
        fout << t << "," << x << "," << y << "\n";
    }
    fout.close();
    return x;
}

int main() {
    double m = 1.0;
    double v0 = 100.0;
    double theta_deg = 45.0;
    double theta = theta_deg * M_PI / 180.0;
    double D = 0.0, a = 0.0, alpha = 2.5;

    //  Task 2: No drag  
    double tmax = analytical_tmax(v0, theta);
    double xmax = analytical_xmax(v0, theta);
    cout << "Analytical tmax=" << tmax << " s, xmax=" << xmax << " m\n";

    int n_values[] = {10, 20, 50, 100, 200, 500};
    ofstream ferr("errors.csv");
    ferr << "dt,Eglob\n";
    for (int n : n_values) {
        double dt = tmax / n;
        string filename = "traj_n" + to_string(n) + ".csv";
        double xnum = simulate(dt, v0, theta, D, a, m, alpha, filename);
        double Eglob = xmax - xnum;
        cout << "n=" << n << " dt=" << dt << " xnum=" << xnum
             << " error=" << Eglob << "\n";
        ferr << dt << "," << Eglob << "\n";
    }
    ferr.close();

    // Task 3: Drag force (a=0) 
    double dt_small = 0.01;
    vector<double> Dvals = {0.0, 1e-4, 2e-4, 5e-4, 1e-3};
    for (double Dval : Dvals) {
        string fname = "drag_traj_D" + to_string((int)(Dval * 1e6)) + ".csv";
        simulate(dt_small, v0, theta, Dval, 0.0, m, alpha, fname);
    }

    //  (D=0, 1e-3, 2e-3)
    vector<double> Dscan = {0.0, 1e-3, 2e-3};
    for (double Dval : Dscan) {
        ofstream fout("range_D" + to_string((int)(Dval * 1000)) + ".csv");
        fout << "theta_deg,xmax\n";
        for (int ang = 15; ang <= 65; ang++) {
            double th = ang * M_PI / 180.0;
            double xnum = simulate(dt_small, v0, th, Dval, 0.0, m, alpha, "tmp.csv");
            fout << ang << "," << xnum << "\n";
        }
        fout.close();
    }

    //  Task 4: Drag + altitude correction 
    m = 20.0; v0 = 700.0; D = 1e-3; a = 6.5e-3;
    vector<int> angles = {35, 45};
    for (int ang : angles) {
        double th = ang * M_PI / 180.0;
        string f1 = "altcorr_angle" + to_string(ang) + "_a0.csv";
        string f2 = "altcorr_angle" + to_string(ang) + "_aPos.csv";
        simulate(dt_small, v0, th, D, 0.0, m, alpha, f1);    // a=0
        simulate(dt_small, v0, th, D, a, m, alpha, f2);      // a>0
    }

    cout << "Simulation Done.\n";
    return 0;
}

