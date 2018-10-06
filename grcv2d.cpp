//
// Created by 服部篤樹 on 2018/05/16.
//
#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <cmath>

void spline(int nin, std::vector<double>& x, std::vector<double>& y, std::vector<double>& fdp);
void clst1(int n, double dx1, std::vector<double>& x);
double clst1_sub(double e, double dx1, int n1);
void clst2(int n, double sp1, double sp2, std::vector<double>& xout);
void stret(std::vector<double>& x, std::vector<double>& y, double s1, double s2, int n);
double fasin(double y);
double fasinh(double y);
double speval(int nin, std::vector<double>& x, std::vector<double>& y, std::vector<double>& fdp, double xx, double f);


void grcv2d(int nin, std::vector<double>& xin, std::vector<double>& yin, int n, double k, double dl1, double dl2,
            std::vector<double>& xout, std::vector<double>& yout){
    std::vector<double> arcs(301), cspx(301), cspy(301), sn(301);
    for(int j=1; j<nin; j++){
        arcs[j] = arcs[j-1] + std::sqrt(std::pow(xin[j]-xin[j-1], 2.0) + std::pow(yin[j] - yin[j-1],2.0));
    }

    /*    for(auto itr=arcs.begin(); itr<arcs.end(); itr++){
        std::cout << *itr << std::endl;
    }*/

    if(nin > 2){
        spline(nin, arcs, xin, cspx);
        spline(nin, arcs, yin, cspy);
    }

/*    for(auto itr=cspx.begin(); itr<cspx.end(); itr++){
        std::cout << *itr << std::endl;
    }*/

    /*
    * CLUSTERING
    */
    double sl1 = dl1 / arcs[nin - 1];
    double sl2 = dl2 / arcs[nin - 1];

    if(k != 0){
        sl1 = dl1 / static_cast<double>(nin - 1);
        sl2 = dl2 / static_cast<double>(nin - 1);
    }

    if(sl2 <= 0){
        clst1(n, sl1, sn);
    }else{
        clst2(n, sl1, sl2, sn);
    }

    /*
     * ARCS ---> XOUT, YOUT
     */


    if(nin == 2){
        for(int i = 0; i<n; i++){
            xout[i] = (1 - sn[i])*xin[0] + sn[i] * xin[1];
            yout[i] = (1 - sn[i])*yin[0] + sn[i] * yin[1];


        }
    }else{
        for(int i = 0; i<n; i++){
            double snd = sn[i] * arcs[nin];
            xout[i] = speval(nin, arcs, xin, cspx, snd, xout[i]);
            yout[i] = speval(nin, arcs, yin, cspx, snd, yout[i]);


        }
    }

    return;
}

void spline(int nin, std::vector<double>& x, std::vector<double>& y, std::vector<double>& fdp){
    const int n1 = nin - 1;
    std::vector<double> a(301), b(301), c(301), r(301);

    c[0] = x[1] - x[0];
    for(int i = 1; i < n1; i++){
        c[i] = x[i+1] - x[i];
        a[i] = c[i-1];
        b[i] = 2.0 * (a[i] + c[i]);
        r[i] = 6.0 * ((y[i+1] - y[i])/c[i] - (y[i] - y[i-1])/c[i-1]);
    }
    b[1] = b[1] + c[0];
    b[n1 - 1] = b[n1 - 1] + c[n1 - 1];

    for (int i = 2; i < n1; i++){
        double t = a[i]/b[i-1];
        b[i] = b[i] - t*c[i-1];
        r[i] = r[i] - t*r[i-1];
    }
    fdp[n1 - 1] = r[n1 - 1]/b[n1 - 1];

    for (int i = 1; i < nin - 2; i++){
        int nmi = nin - i;
        fdp[nmi - 1] = (r[nmi - 1] - c[nmi - 1] * fdp[nmi])/b[nmi - 1];
    }
    fdp[0] = fdp[1];
    fdp[nin-1] = fdp[n1 - 1];
    return;
}

void clst1(int n, double dx1, std::vector<double>& x){
    int n1 = n - 1;
    double dx = 1.0 / static_cast<double>(n1);
    double e;

    if (std::fabs(dx - dx1) < 1.0e-5){
        e = 1;
    }
    else{
        e = clst1_sub(1.2, dx1, n1);
        while(std::abs(e - 1.0) < 1.0e-5){
            double e1 = std::pow(e, 2.0);
            e = clst1_sub(e1, dx1, n1);
        }
        //例外処理
        if(e < 0){
            return;
        }
    }
    x[0] = 0;
    dx = dx1;
    for (int i = 1; i < n; i++){
        x[i] = x[i-1] + dx;
        dx = dx * e;
    }
    return;
}

double clst1_sub(double e, double dx1, int n1){
    int iter = 0;
    double ep;
    while(1){
        ep = e;
        e = ep - (std::pow(dx1 * ep,n1) - ep + 1.0 - dx1)/(std::pow(dx1 * static_cast<double>(n1) * e, n1 - 1.0));
        if(std::abs(e - ep) < 1.0e-5){
            return e;
        }
        //例外処理（ひとまずの妥協案）
        if (iter > 100){
            std::cout << "Iteration failed in clst1_sub()!" << std::endl;
            return -1;
        }
    }
}

void clst2(int n, double sp1, double sp2, std::vector<double>& xout){
    std::vector<double> xin(1000);
    if (sp1 < 0 || sp2 < 0){
        std::cout << "Wrong Values in sp1, sp2 :" << sp1 << ", " << sp2 << std::endl;
        std::cout << "sp1, sp2 must be 0 < sp1,sp2 < 0.3" << std::endl;
    }
    double smal = std::min(sp1,sp2);
    int n1 = n - 1;
    for (int i = 0; i < n; i++){
        xin[i] = static_cast<double>(i)/ static_cast<double>(n1);
    }
    double s11 = 1.0 / (sp1 * static_cast<double>(n1));
    double s21 = 1.0 / (sp2 * static_cast<double>(n1));

    stret(xin, xout, s11, s21, n);

    double dx1 = xout[1];
    double dx2 = xout[n-1] - xout[n1 - 1];

    //Improve accuracy of end-spacing control
    if(fabs(sp1 - dx1) <= smal && abs(sp2 - dx2) <= smal){
        //continue
    }else{
        for(int i = 0; i<20; i++){
            double ds1 = s11 * 0.1;
            double ds2 = s21 * 0.1;
            double s12 = s11 + ds1;
            double s22 = s21 + ds2;
            stret(xin, xout, s12, s21, n);

            double dx1ds1 = (xout[1] - dx1)/ds1;
            double dx2ds1 = (xout[n-1] - xout[n1 - 1] - dx2)/ds1;
            stret(xin, xout, s11, s22, n);

            double dx1ds2 = (xout[1] - dx1)/ds2;
            double dx2ds2 = (xout[n-1] - xout[n1 - 1] - dx2)/ds2;

            double df1 = sp1 - dx1;
            double df2 = sp2 - dx2;
            double d = dx1ds1*dx2ds2 - dx1ds2*dx2ds1;

            ds1 = (df1*dx2ds2 - df2*dx1ds2)/d;
            ds2 = (df2*dx1ds1 - df1*dx2ds1)/d;
            s12 = s11 + ds1;
            s22 = s21 + ds2;
            stret(xin, xout, s12, s22, n);

            dx1 = xout[1];
            dx2 = xout[n - 1] - xout[n1];
            if(fabs(sp1 - dx1) <= smal && fabs(sp2 - dx2) <= smal){
                break;
            }
            s11 = s12;
            s21 = s22;
        }
    }
    xout[0] = 0;
    xout[n - 1] = 1;
    return;
}

void stret(std::vector<double>& x, std::vector<double>& y, double s1, double s2, int n){
    double b = std::sqrt(s1 * s2);
    double a = b/s2;

    if(b < 0.999){
        double dz = fasin(b);
        double tanx;
        for (int j = 1; j<n; j++){
            tanx = std::tan(dz*x[j]);
            y[j] = tanx/(a*sin(dz) + (1 - a*cos(dz))*tanx);
        }
    } else if(b > 1.001){
        double dz = fasinh(b);
        double tanhx;
        for (int j = 0; j<n; j++){
            tanhx = tanh(dz*x[j]);
            y[j] = tanhx/(a*sinh(dz) + (1 - a*cosh(dz))*tanhx);

        }
    }else{
        for(int j = 0; j<n ; j++){
            double u = x[j]*(1 + 2*(b-1)*(x[j] - 0.5)*(1 - x[j]));
            y[j] = u/(a + (1 - a)*u);
        }
    }
    return;
}

//多分おk
// Inversion of y = sin(x)/x
double fasin(double y){
    const double a3 = -2.6449341;
    const double a4 = 6.794732;
    const double a5 = 13.205501;
    const double a6 = 11.726095;
    const double b2 = 0.057321429;
    const double b3 = 0.048774238;
    const double b4 = -0.053337753;
    const double b5 = 0.075845134;

    if (y <= 0.26938972){
        return (M_PI*((((((a6*y + a5)*y + a4)*y + a3)*y + 1)*y - 1)*y + 1));
    }
    else {
        double yb = 1 - y;
        return (std::sqrt(6*yb)*(((((b5*yb + b4)*yb + b3)*yb + b2)*yb + 0.15)*yb + 1));
    }
}

//多分おk
//Inversion of y = sinh(x)
double fasinh(double y){
    const double a2 = 0.057321429;
    const double a3 = -0.024907295;
    const double a4 = 0.0077424461;
    const double a5 = -0.0010794123;
    const double b0 = -0.02041793;
    const double b1 = 0.24902722;
    const double b2 = 1.9496443;
    const double b3 = -2.6294547;
    const double b4 = 8.56795911;

    if (y <= 2.7829681){
        double yb = y - 1;
        return (std::sqrt(6 * yb)*(((((a5*yb + a4)*yb + a3)*yb + a2)*yb -0.15)* yb + 1.0));
    }
    else{
        double v = std::log(y);
        double w = 1/y - 0.028527431;
        return (v + std::log(2*v)*(1+1/v)+(((b4*w + b3)*w + b2)*w + b1)*w + b0);
    }
}

double speval(int nin, std::vector<double>& x, std::vector<double>& y, std::vector<double>& fdp, double xx, double f){
    for (int i = 0; i < nin - 2; i++){
        if(xx <= x[i+1]){
            double dxm = xx - x[i];
            double dxp = x[i + 1] - xx;
            double del = x[i + 1] - x[i];
            return (fdp[i]*dxp*(dxp*dxp/del - del)/6 + y[i]*dxp/del + fdp[i+1]*dxm*(dxm*dxm/del - del)/6 + y[i+1]*dxm/del);
        }else{
            //continue
        }
    }
    return -1;
}