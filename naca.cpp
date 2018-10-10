//
// Created by 服部篤樹 on 2018/10/05.
//
#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "grcv2d.h"

double ynaca(double t, double xx);

void tinput(int jdim, int kdim, int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y){
    const int ld = 301;

    //一時的な出力
    std::vector<double> x1(ld), y1(ld), x2(ld), y2(ld);

    //jmaxは奇数出ないといけない
    const int jwake = 14;
    const double ds1 = 0.03;
    const double ds2 = 0.8;
    const double rout = 5.0;
    const double xdw = 5.0;
    const int mtout = 11;

    const int jlead = (jmax + 1)/2;
    const int jtail1 = jwake + 1;
    const int jtail2 = jmax + 1 - jtail1;
    const int nwall = jmax  - (jtail1 - 1)*2;

    //<<NACA0012>>
    const double tnaca = 0.12;
    const int jleadw = 301;

    /*
     * Airfoil Lower Surface
     */
    for (int i = 0; i<jleadw; i++){
        x1[i] = 1.0 - static_cast<double>(i - 1)/ static_cast<double>(jleadw - 1);
        y1[i] = -ynaca(tnaca, x1[i]);
    }
    const int nwallh = (nwall + 1)/2;
    grcv2d(jleadw, x1, y1, nwallh, 1, 0.5, 0.5, x2, y2);
    //std::cout << jtail1 << " " << jlead << " " << std::endl;

/*    //TEST
    //ファイルチェック
    //csvファイルを出力
    std::ofstream ofscsv("naca0012_1.csv");    //, std::ios::app
    if (!ofscsv) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for (int i = 0; i<x1.size() ; i++){
        ofscsv << x1[i] << "," << y1[i] << std::endl;
    }

    std::ofstream ofscsv2("naca0012_2.csv");    //, std::ios::app
    if (!ofscsv2) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for (int i = 0; i<x2.size() ; i++){
        ofscsv2 << x2[i] << "," << y2[i] << std::endl;
    }*/

    for(int j = jtail1 - 1; j<jlead; j++){
        int jj = j - jtail1 + 1;
        x[j][0] = x2[jj];
        y[j][0] = y2[jj];
        //std::cout << x[j][0] << std::endl;
    }

    /*
     * Airfoil Upper Surface
     */
    for(int i = 0; i<jleadw; i++){
        x2[i] = x1[jleadw - i - 1];
        y2[i] = ynaca(tnaca, x2[i]);
    }
    grcv2d(jleadw, x2, y2, nwallh, 1, 0.5, 0.5, x1, y1);
    for(int j = jlead; j < jtail2; j++){
        int jj = j - jlead;             //正しい？
        x[j][0] = x1[jj];
        y[j][0] = y1[jj];
        //std::cout << x[j][0] << std::endl;
    }

    /*
     * Grid on Inner Boundary : WAKE
     */
    if(jwake >= 1){
        double dsx = sqrt(pow(x[jtail1 - 1][0] - x[jtail1][0], 2.0) + pow(y[jtail1 - 1][0] - y[jtail1][0], 2.0));
        x1[1] = xdw;
        y1[1] = 0;
        x1[0] = 1;
        y1[0] = 0;
        grcv2d(2, x1, y1, jtail1, 0, dsx, 0, x2, y2);
        for(int j = 0; j<jtail1-1; j++){
            int jj = jtail1 - 1 - j;
            //std::cout << x2[jj] << std::endl;
            x[j][0] = x2[jj];
            x[jmax-j-1][0] = x[j][0];
            y[j][0] = 0;
            y[jmax-j-1][0] = y[j][0];
        }
    }

    /*
     * Grid on Outer Boundary
     */
    const int nwk2 = 50;
    if(jwake >= 1){
        x1[0] = xdw;
        y1[0] = -rout;
        x1[1] = 1;
        y1[1] = -rout;
        for(int j = 2; j<nwk2-2; j++){
            double th = 1.5 * M_PI - M_PI/ static_cast<double>(nwk2 - 2 - 1)* static_cast<double>(j - 1);
            x1[j] = 1 + rout * cos(th);
            y1[j] = rout * sin(th);
        }
        x1[nwk2 - 2] = 1;
        y1[nwk2 - 2] = rout;
        x1[nwk2 - 1] = xdw;
        y1[nwk2 - 1] = rout;
    }else{
        for(int j = 0; j<nwk2; j++){
            double th = -2 * M_PI / static_cast<double>(nwk2 - 1)/ static_cast<double>(j);
            x1[j] = rout * cos(th) + 0.5;
            y1[j] = rout * sin(th);
        }
    }
    grcv2d(nwk2, x1, y1, jmax, 1, 1, 1, x2, y2);
    for(int j = 0; j<jmax; j++){
        x[j][kmax - 1] = x2[j];
        y[j][kmax - 1] = y2[j];
        //std::cout << y[j][kmax] << std::endl;
    }


    std::ofstream ofscsv2("naca0012_2.csv");    //, std::ios::app
    if (!ofscsv2) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for(int j = 0; j<jmax; j++){
        ofscsv2 << x2[j] << "," << y2[j] << std::endl;
    }

    /*
     * 下流境界
     */
    for(int j = 0; j < jmax; j += jmax - 1){
        x1[0] = x[j][0];
        y1[0] = y[j][0];
        x1[1] = x[j][kmax - 1];
        y1[1] = y[j][kmax - 1];
        //std::cout << ds1 << " " << ds2 << std::endl;
        grcv2d(2, x1, y1, kmax, 0, ds1, ds2, x2, y2);
        for(int k = 1; k<kmax - 1; k++){
            x[j][k] = x2[k];
            y[j][k] = y2[k];
            //std::cout << y[j][k] << std::endl;
        }
    }
}

double ynaca(double t, double xx){
    double x = xx*1.00893;
    return (t/ 0.2 * (0.2696 * sqrt(x) - 0.126*x - 0.3516*pow(x, 2.0) + 0.2843*pow(x, 3.0) - 0.1015*pow(x, 4.0)));
}

void tcore(int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y,
           std::vector<std::vector<double>>& x1, std::vector<std::vector<double>>& y1){

    std::vector<double> s(300);
    s[0] = 0;
    for(int k = 1; k<kmax; k++) {
        s[k] = s[k - 1] + sqrt(pow(x[0][k] - x[0][k - 1], 2.0) + pow(y[0][k] - y[0][k - 1], 2.0));
    }
    for (int j = 0; j < jmax; j++) {
        for (int k = 0; k < kmax; k++) {
            double beta = 1 - s[k] / s[kmax - 1];
            x1[j][k] = beta * x[j][0] + (1 - beta) * x[j][kmax - 1];
            y1[j][k] = beta * y[j][0] + (1 - beta) * y[j][kmax - 1];
        }
    }

    s[0] = 0;
    for (int j = 1; j < jmax; j++) {
        s[j] = s[j - 1] + sqrt(pow(x[j][0] - x[j - 1][0], 2.0) + pow(y[j][0] - y[j - 1][0], 2.0));
    }
    for (int j = 0; j < jmax; j++) {
        for (int k = 0; k < kmax; k++) {
            double alpha = 1 - s[j] / s[jmax - 1];
            x[j][k] = x1[j][k] + alpha * (x[0][k] - x1[0][k]) + (1 - alpha) * (x[jmax - 1][k] - x1[jmax - 1][k]);
            y[j][k] = y1[j][k] + alpha * (y[0][k] - y1[0][k]) + (1 - alpha) * (y[jmax - 1][k] - y1[jmax - 1][k]);
            //std::cout << x1[j][k] << std::endl;
        }
    }
}

void gio(int id, int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y){
    //gnuplot用にファイル出力
    std::ofstream ofs("mygrid.dat");
    if (!ofs) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }
    //TecPLot用に書式を整える
    ofs << "VARIABLES=\"X\",\"Y\"" << std::endl;
    ofs << "ZONE T=\"Mygrid\",I=" << kmax << ",J=" << jmax << ",F=POINT" << std::endl;
    ofs << "DT=(DOUBLE DOUBLE)" << std::endl;
    for (int j = 0; j<jmax ; j++){
        for(int k = 0; k<kmax; k++){
            ofs << x[j][k] << " " << y[j][k] << std::endl;
        }
    }
}