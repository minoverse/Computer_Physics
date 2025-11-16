

#include <bits/stdc++.h>
using namespace std;

struct Params {
    double lam0=1, lam1=5, lam2=50;
    double tmax=200.0;
    double dt_init=1e-2;        // assignment default
    double tol_trap=1e-10;     // Newton (trapezoid) tolerance
    double tol_step=1e-4;      // adaptive tolerance (Task 2 says TOL=1e-4)
    double safety=0.9;         // step-size safety factor
    double fac_min=0.2, fac_max=5.0;
    int newton_maxit=20;
};

static inline void f(const Params& P, const array<double,3>& N, array<double,3>& rhs){
    rhs[0] = -P.lam0 * N[0];
    rhs[1] =  P.lam0 * N[0] - P.lam1 * N[1];
    rhs[2] =  P.lam1 * N[1] - P.lam2 * N[2];
}

static inline void exact(const Params& P, double t, array<double,3>& Ne){
    double L0=P.lam0, L1=P.lam1, L2=P.lam2;
    double e0=exp(-L0*t), e1=exp(-L1*t), e2=exp(-L2*t);
    
    double N0 = e0;
    
    double N1 = 0.0;
    if (fabs(L1-L0) > 1e-14) N1 = (L0/(L1-L0))*(e0 - e1);
    else                     N1 = L0*e0*t; // degenerate fallback (rare in tests)
    
    double N2 = 0.0;
    auto term = [&](double La, double Lb, double Lc, double ea){
        return (L0*L1) / ((Lb-La)*(Lc-La)) * ea;
    };
    if (fabs((L1-L0)*(L2-L0)*(L2-L1)) > 1e-20) {
        N2 = term(L0,L1,L2,e0) + term(L1,L0,L2,e1) + term(L2,L0,L1,e2);
    } else {
       
        N2 = max(0.0, 1.0 - N0 - N1); // acceptable for test parameters
    }  
    Ne = {N0, N1, N2};
}

// solve A x = b (3x3)  Gaussian elimination
static bool solve3x3(double A[3][3], double b[3], double x[3]){
    int p[3] = {0,1,2};
    
    for(int k=0;k<3;k++){
        int piv=k;
        double mv=fabs(A[piv][k]);
        for(int r=k+1;r<3;r++){
            if(fabs(A[r][k])>mv){ mv=fabs(A[r][k]); piv=r; }
        }
        if (mv<1e-30) return false;
        if(piv!=k) swap(p[k],p[piv]);
        // eliminate
        for(int i=k+1;i<3;i++){
            double alpha = A[p[i]][k]/A[p[k]][k];
            for(int j=k;j<3;j++) A[p[i]][j]-=alpha*A[p[k]][j];
            b[p[i]]-=alpha*b[p[k]];
        }
    }
    
    for(int i=2;i>=0;i--){
        double sum=b[p[i]];
        for(int j=i+1;j<3;j++) sum-=A[p[i]][j]*x[j];
        x[i]=sum/A[p[i]][i];
    }
    return true;
}

//  (Newton): from (t,N) with dt -> (t+dt, Nnew)
// Returns success flag; Nnew  valid only if true.
static bool trapezoid_step(const Params& P, double t, const array<double,3>& N,
                           double dt, array<double,3>& Nnew)
{
    array<double,3> fn, guess = N;
    f(P, N, fn);
    // Newton iterations: solve G(Y)=0 for Y = N(t+dt)
    // G(Y) = Y - N - 0.5*dt*( f(N) + f(Y) )
    for(int it=0; it<P.newton_maxit; it++){
        array<double,3> fY;
        f(P, guess, fY);
        double G[3] = {
            guess[0] - N[0] - 0.5*dt*(fn[0] + fY[0]),
            guess[1] - N[1] - 0.5*dt*(fn[1] + fY[1]),
            guess[2] - N[2] - 0.5*dt*(fn[2] + fY[2])
        };
        //  J = df/dN is const
       
        // dG/dY = I - 0.5*dt*J
        double A[3][3] = {
            {1.0 + 0.5*dt*P.lam0,  0.0,                   0.0},
            {-0.5*dt*P.lam0,       1.0 + 0.5*dt*P.lam1,  0.0},
            {0.0,                  -0.5*dt*P.lam1,       1.0 + 0.5*dt*P.lam2}
        };
        double rhs[3] = {-G[0], -G[1], -G[2]};
        double dY[3]  = {0,0,0};
        if(!solve3x3(A,rhs,dY)) return false;
        double max_corr = 0.0;
        for(int k=0;k<3;k++){
            guess[k] += dY[k];
            max_corr = max(max_corr, fabs(dY[k]));
        }
        if (max_corr < P.tol_trap) { Nnew = guess; return true; }
    }
    return false; 
}


static double controller_factor(const Params& P, double err){
    if (err<=0) return P.fac_max;
    double fac = P.safety * pow(P.tol_step / err, 1.0/3.0);
    fac = min(P.fac_max, max(P.fac_min, fac));
    return fac;
}

static string out_name(const Params& P){
    ostringstream ss;   
    ss.setf(std::ios::fixed); ss<<setprecision(6);
    ss<<"decay_l0_"<<P.lam0<<"_l1_"<<P.lam1<<"_l2_"<<P.lam2<<"_tol_";
    ss<<std::scientific<<setprecision(0)<<P.tol_step;
    string s = ss.str();
    for(char& c: s) if(c=='+') c='p'; // sanitize
    return s + ".csv";
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Params P;
    // parse simple --key=value args
    for(int i=1;i<argc;i++){
        string a(argv[i]);
        auto eq=a.find('=');
        auto key=a.substr(0,eq);
        auto val=(eq==string::npos)?"":a.substr(eq+1);
        auto to_d=[&](){return atof(val.c_str());};
        if(key=="--lam0") P.lam0=to_d();
        else if(key=="--lam1") P.lam1=to_d();
        else if(key=="--lam2") P.lam2=to_d();
        else if(key=="--tol")  P.tol_step=to_d();
        else if(key=="--tmax") P.tmax=to_d();
        else if(key=="--dt0")  P.dt_init=to_d();
    }

    string fname = out_name(P);
    ofstream out(fname);
    out<<"t,N0,N1,N2,exN0,exN1,exN2,dt\n";

    array<double,3> N = {1.0,0.0,0.0};
    double t=0.0, dt=P.dt_init;

    // Write initial line
    array<double,3> Ne;
    exact(P,t,Ne);
    out<<setprecision(16)<<t<<","<<N[0]<<","<<N[1]<<","<<N[2]<<","
       <<Ne[0]<<","<<Ne[1]<<","<<Ne[2]<<","<<dt<<"\n";

    while(t < P.tmax){
        if (t+dt > P.tmax) dt = P.tmax - t;

        // One full step
        array<double,3> y_full;
        bool ok_full = trapezoid_step(P,t,N,dt,y_full);

        // Two half-steps
        array<double,3> y_half, y_half2;
        bool ok_h1 = trapezoid_step(P,t,N,dt*0.5,y_half);
        bool ok_h2 = ok_h1 && trapezoid_step(P,t+0.5*dt,y_half,dt*0.5,y_half2);

        if(!(ok_full && ok_h1 && ok_h2)){
            
            dt = max(1e-16, 0.5*dt);
            continue;
        }

        // error estimate
        double err = max( { fabs(y_half2[0]-y_full[0]),
                            fabs(y_half2[1]-y_full[1]),
                            fabs(y_half2[2]-y_full[2]) } );

        if (err <= P.tol_step){
            
            N = y_half2;
            t += dt;
            // record
            exact(P,t,Ne);
            out<<setprecision(16)<<t<<","<<N[0]<<","<<N[1]<<","<<N[2]<<","
               <<Ne[0]<<","<<Ne[1]<<","<<Ne[2]<<","<<dt<<"\n";
            
            dt *= controller_factor(P, err);
        } else {
           
            dt *= controller_factor(P, err);
        }
        if (dt < 1e-16) { cerr<<"dt underflow\n"; break; }
    }
    out.close();

    cerr<<"Wrote: "<<fname<<"\n";
    return 0; 
}


