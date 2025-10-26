

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

struct Params {
    double k   = 1.0;
    double m   = 1.0;
    double alpha = 0.0;
    double F0    = 0.0;
    double Omega = 0.0;
    double x0  = 1.0;
    double v0  = 0.0;
    double tmax = 50.0;
    long long N = 10000;
};

struct State { double x, v; };

// RHS u' = f(t,u)
static inline State rhs(double t, const State& u, const Params& p){
    State f;
    f.x = u.v;
    f.v = -(p.k/p.m)*u.x - (p.alpha/p.m)*u.v + (p.F0/p.m)*std::sin(p.Omega*t);
    return f;
}

// RK4 step
static inline State rk4_step(double t, double dt, const State& u, const Params& p){
    State k1 = rhs(t, u, p);
    State u2{u.x + 0.5*dt*k1.x, u.v + 0.5*dt*k1.v};
    State k2 = rhs(t + 0.5*dt, u2, p);
    State u3{u.x + 0.5*dt*k2.x, u.v + 0.5*dt*k2.v};
    State k3 = rhs(t + 0.5*dt, u3, p);
    State u4{u.x + dt*k3.x, u.v + dt*k3.v};
    State k4 = rhs(t + dt, u4, p);
    State out;
    out.x = u.x + (dt/6.0)*(k1.x + 2*k2.x + 2*k3.x + k4.x);
    out.v = u.v + (dt/6.0)*(k1.v + 2*k2.v + 2*k3.v + k4.v);
    return out;
}

static inline double Ekin(double m, double v){ return 0.5*m*v*v; }
static inline double Epot(double k, double x){ return 0.5*k*x*x; }

// Underdamped exact (F0=0, alpha>0); NaN for critical/overdamped
static inline double exact_x(double t, const Params& p){
    double omega0 = std::sqrt(p.k/p.m);
    if(p.alpha<=0) return p.x0*std::cos(omega0*t) + (p.v0/omega0)*std::sin(omega0*t);
    double tau = 2.0*p.m/p.alpha;
    double disc = omega0*omega0 - 1.0/(tau*tau);
    if(disc<=0) return std::numeric_limits<double>::quiet_NaN();
    double omega_d = std::sqrt(disc);
    double A = std::sqrt(p.x0*p.x0 + std::pow((p.x0 + p.v0*tau)/(omega_d*tau),2));
    double phi = std::atan2((p.x0 + p.v0*tau), (p.x0*omega_d*tau));
    return A * std::exp(-t/tau) * std::cos(omega_d*t + phi);
}

static void simulate_one(const std::string& tag, const Params& p){
    double dt = p.tmax / double(p.N);
    std::ofstream f("trajectory_"+tag+".csv");
    f << "t,x,v,Ekin,Epot,Etot\n";
    State u{p.x0,p.v0};
    double t=0;
    for(long long i=0;i<=p.N;i++){
        double ek=Ekin(p.m,u.v), ep=Epot(p.k,u.x);
        f<<std::setprecision(12)<<t<<","<<u.x<<","<<u.v<<","<<ek<<","<<ep<<","<<(ek+ep)<<"\n";
        u = rk4_step(t, dt, u, p);
        t += dt;
    }
}

static void simulate_compare_exact(const std::string& tag, Params p){
    double dt = p.tmax / double(p.N);
    std::ofstream f("compare_"+tag+".csv");
    f << "t,x_num,x_exact\n";
    State u{p.x0,p.v0};
    double t=0;
    for(long long i=0;i<=p.N;i++){
        double xe = exact_x(t,p);
        f<<std::setprecision(12)<<t<<","<<u.x<<","<<xe<<"\n";
        u = rk4_step(t, dt, u, p);
        t += dt;
    }
}

static void sweep_resonance(const Params& base_params, double omega_lo, double omega_hi, double domega){
    std::ofstream f("resonance.csv");
    f<<"Omega,alpha,amp\n";
    for(double Om=omega_lo; Om<=omega_hi+1e-12; Om+=domega){
        Params p=base_params; p.Omega=Om;
        double dt = p.tmax / double(p.N);
        State u{p.x0,p.v0};
        double t=0;
        double x_prev = u.x, x_curr = u.x, x_next;
        double last_peak = std::fabs(u.x);
        for(long long i=0;i<=p.N;i++){
            x_next = rk4_step(t, dt, u, p).x; // look-ahead
            if(x_curr >= x_prev && x_curr > x_next) last_peak = std::fabs(x_curr);
            u = rk4_step(t, dt, u, p);
            t += dt;
            x_prev = x_curr; x_curr = u.x;
        }
        f<<std::setprecision(12)<<Om<<","<<p.alpha<<","<<last_peak<<"\n";
    }
}

int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Params p;
    std::string task="2a";
    std::vector<std::string> args(argv+1, argv+argc);
    double omega_lo=0.1, omega_hi=2.0, domega=0.01;
    std::vector<double> alphas_task3{1e-4, 0.1, 0.5, 1.95};
    std::vector<double> alphas_task4{0.01, 0.1, 0.5, 1.0};

    for(const auto& a: args){
        auto pos=a.find('=');
        if(pos==std::string::npos){ if(a=="help"){ std::cerr<<"see source header for examples\n"; } continue; }
        std::string k=a.substr(0,pos), v=a.substr(pos+1);
        auto to_d = [](const std::string& s){ return std::stod(s); };
        if(k=="task") task=v;
        else if(k=="k") p.k=to_d(v);
        else if(k=="m") p.m=to_d(v);
        else if(k=="alpha") p.alpha=to_d(v);
        else if(k=="F0") p.F0=to_d(v);
        else if(k=="Omega") p.Omega=to_d(v);
        else if(k=="x0") p.x0=to_d(v);
        else if(k=="v0") p.v0=to_d(v);
        else if(k=="tmax") p.tmax=to_d(v);
        else if(k=="N") p.N=std::stoll(v);
        else if(k=="sweep"){
            double a_,b_,c_; char colon;
            std::stringstream ss(v);
            ss>>a_>>colon>>b_>>colon>>c_;
            omega_lo=a_; omega_hi=b_; domega=c_;
        } else if(k=="alphas3"){
            alphas_task3.clear();
            std::stringstream ss(v); std::string item;
            while(std::getline(ss,item,',')) alphas_task3.push_back(std::stod(item));
        } else if(k=="alphas4"){
            alphas_task4.clear();
            std::stringstream ss(v); std::string item;
            while(std::getline(ss,item,',')) alphas_task4.push_back(std::stod(item));
        }
    }

    if(task=="2a"){
        simulate_one("2a",p);
    } else if(task=="2b"){
        simulate_one("2b",p);
        simulate_compare_exact("2b",p);
    } else if(task=="3"){
        for(double a: alphas_task3){
            Params q=p; q.alpha=a;
            simulate_one(std::string("3_alpha_")+std::to_string(a),q);
        }
    } else if(task=="4"){
        std::ofstream fout("resonance_all.csv");
        fout<<"Omega,alpha,amp\n";
        for(double a: alphas_task4){
            Params q=p; q.alpha=a;
            sweep_resonance(q, omega_lo, omega_hi, domega);
            std::ifstream rin("resonance.csv");
            std::string line; std::getline(rin,line); // skip header
            while(std::getline(rin,line)) fout<<line<<"\n";
        }
    } else {
        std::cerr<<"Unknown task. Use task=2a|2b|3|4\n";
        return 1;
    }
    return 0;
}
