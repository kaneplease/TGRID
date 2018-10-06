//
// Created by 服部篤樹 on 2018/10/05.
//

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "grcv2d.h"

double ynaca(double t, double xx);

double tinput(){
    const int jdim = 301;
    const int kdim = 100;
    const int ld = 301;

    //最終的な出力
    std::vector<std::vector<double>> x(jdim, std::vector<double>(kdim));
    std::vector<std::vector<double>> y(jdim, std::vector<double>(kdim));

    //一時的な出力
    std::vector<double> x1(ld), y1(ld), x2(ld), y2(ld);

    //jmaxは奇数出ないといけない
    const int jmax = 81;
    const int kmax = 20;
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

    //TEST
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
    }

}

double ynaca(double t, double xx){
    double x = xx*1.00893;
    return (t/ 0.2 * (0.2696 * sqrt(x) - 0.126*x - 0.3516*pow(x, 2.0) + 0.2843*pow(x, 3.0) - 0.1015*pow(x, 4.0)));
}