

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>

struct State { double x, y, vx, vy; };
const double PI = 3.141592653589793;
const double GMs = 4.0 * PI * PI;       // AU^3 / yr^2

// ---- Neri 4th-order symplectic constants ----
const double s  = cbrt(2.0);
const double a1 = 1.0 / (2.0 * (2.0 - s));
const double a2 = (1.0 - s) / (2.0 * (2.0 - s));
const double a3 = a2, a4 = a1;
const double b1 = 1.0 / (2.0 - s);
const double b2 = -s / (2.0 - s);
const double b3 = b1, b4 = 0.0;

// acceleration from Eq.(7)
inline void acc(double x,double y,double alpha,double& ax,double& ay){
    double r = std::sqrt(x*x + y*y);
    double f = -GMs / (r*r*r) * (1.0 + alpha/(r*r));
    ax = f * x; ay = f * y;
}

// main symplectic stepper 
void simulate(double tmax,double dt,double alpha,const std::string& outname){
    State s;
    double a = 0.387098, e = 0.206;
    double rmax = a*(1+e);
    double vmin = 2*PI*std::sqrt((1-e)/(a*(1+e)));  // aphelion
    s.x=rmax; s.y=0; s.vx=0; s.vy=vmin;

    std::ofstream f(outname);
    f<<"t,x,y,vx,vy\n";

    int N=(int)(tmax/dt);
    double t=0;
    for(int i=0;i<N;i++){
        f<<t<<","<<s.x<<","<<s.y<<","<<s.vx<<","<<s.vy<<"\n";

        double ax,ay;
        acc(s.x,s.y,alpha,ax,ay);
        s.x+=a1*s.vx*dt; s.y+=a1*s.vy*dt;
        s.vx+=b1*ax*dt;  s.vy+=b1*ay*dt;

        acc(s.x,s.y,alpha,ax,ay);
        s.x+=a2*s.vx*dt; s.y+=a2*s.vy*dt;
        s.vx+=b2*ax*dt;  s.vy+=b2*ay*dt;

        acc(s.x,s.y,alpha,ax,ay);
        s.x+=a3*s.vx*dt; s.y+=a3*s.vy*dt;
        s.vx+=b3*ax*dt;  s.vy+=b3*ay*dt;

        acc(s.x,s.y,alpha,ax,ay);
        s.x+=a4*s.vx*dt; s.y+=a4*s.vy*dt;
        // b4=0

        t+=dt;
    }
    f.close();
}

//  detect perihelia for precession 
struct Hit { double t, theta; };
std::vector<Hit> find_perihelia(const std::string& file){
    std::ifstream f(file);
    std::string line; getline(f,line);
    std::vector<double> t,x,y,r;
    while(std::getline(f,line)){
        double tt,xx,yy,vx,vy; char c;
        std::stringstream ss(line);
        ss>>tt>>c>>xx>>c>>yy>>c>>vx>>c>>vy;
        t.push_back(tt); x.push_back(xx); y.push_back(yy);
        r.push_back(std::sqrt(xx*xx+yy*yy));
    }
    std::vector<Hit> hits;
    for(size_t i=1;i+1<r.size();++i)
        if(r[i]<r[i-1] && r[i]<r[i+1])
            hits.push_back({t[i], std::atan2(y[i],x[i])});
    return hits;
}

double omega_from_two_perihelia(const std::string& file){
    auto h=find_perihelia(file);
    if(h.size()<2) return 0.0;
    double dtheta=h[1].theta-h[0].theta;
    if(dtheta> M_PI) dtheta-=2*M_PI;
    if(dtheta<-M_PI) dtheta+=2*M_PI;
    return dtheta/(h[1].t-h[0].t);
}


int main(){
    double TM=0.240846;
    double dt=1e-4;

    //  Task 2  short test (≈1 orbit)
    simulate(0.95*TM,dt,0.0,"orbit_short.csv");

    //  Task 3  long stability (100 orbits)
    simulate(100*TM,1e-5,0.0,"orbit_long.csv");

    // Task 4  precession visible (α=0.01)
    simulate(4*TM,dt,0.01,"orbit_precession.csv");

    //  Task 5 ω(α) sweep + fit data
    std::ofstream out("alpha_omega.csv");
    out<<"alpha,omega\n";
    double alphamax=1e-3, tmax=3.0, dt2=1e-5;
    for(int j=0;j<=6;j++){
        double alpha=alphamax/std::pow(2.0,j);
        std::string fn="run_alpha_"+std::to_string(j)+".csv";
        simulate(tmax,dt2,alpha,fn);
        double om=omega_from_two_perihelia(fn);
        out<<std::setprecision(12)<<alpha<<","<<om<<"\n";
    }
    out.close();

    // ---- Task 5 (d) 
    
    std::ifstream fin("alpha_omega.csv");
    std::string header; getline(fin,header);
    double a1, w1, a2, w2; fin>>a1; fin.ignore(1,','); fin>>w1;
    for(int k=0;k<6;k++){ fin>>a2; fin.ignore(1,','); fin>>w2; }
    fin.close();
    double slope = (w1-w2)/(a1-a2);
    double alpha_real=1.1e-8;
    double omega_real=slope*alpha_real;              
    double arcsec_cent = omega_real*(180*3600/PI)*100;
    std::cout<<"Estimated ω_real = "<<omega_real<<" rad/year\n";
    std::cout<<"≈ "<<arcsec_cent<<" arcsec/century\n";
    std::cout<<"Expected: 42.98 arcsec/century\n";
    std::cout<<"Simulations complete.\n";
}
