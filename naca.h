//
// Created by 服部篤樹 on 2018/10/06.
//

#ifndef MAKING_LATTICE_TEST_NACA_H
#define MAKING_LATTICE_TEST_NACA_H

void tinput(int jdim, int kdim, int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y);
void tcore(int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y,
           std::vector<std::vector<double>>& x1, std::vector<std::vector<double>>& y1);
void gio(int id, int jmax, int kmax, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y);
#endif //MAKING_LATTICE_TEST_NACA_H
