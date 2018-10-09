#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "grcv2d.h"
#include "naca.h"

std::vector<std::string> split(std::string& input, char delimiter);

int main() {
    /*
     * ファイル読み込み
     */
    std::ifstream ifs("n63412-il.csv");
    std::string line;
    std::vector<double> x_read;
    std::vector<double> y_read;
    while (std::getline(ifs, line)) {
        std::vector<std::string> strvec = split(line, ',');
        for (int i=0; i<strvec.size();i++){
            if(i%2 == 0){
                x_read.push_back(std::stod(strvec[i]));
            }
            else{
                y_read.push_back(std::stod(strvec[i]));
            }
        }
    }
/*    for(auto itr=y_read.begin(); itr<y_read.end(); itr++){
        std::cout << *itr << std::endl;
    }*/

    ////////////////////
    //ここから本番
    ///////////////////
    const int jdim = 301;
    const int kdim = 100;
    const int jmax = 81;
    const int kmax = 20;
    //最終的な出力
    std::vector<std::vector<double>> x(jdim, std::vector<double>(kdim));
    std::vector<std::vector<double>> y(jdim, std::vector<double>(kdim));
    std::vector<std::vector<double>> xtr(jdim, std::vector<double>(kdim));
    std::vector<std::vector<double>> ytr(jdim, std::vector<double>(kdim));

    tinput(jdim, kdim, jmax, kmax, x, y);
    tcore(jmax, kmax, x, y, xtr, ytr);
    gio(1, jmax, kmax, x, y);
    return 0;
}

//CSVファイル読み込み用
std::vector<std::string> split(std::string& input, char delimiter)
{
    std::istringstream stream(input);
    std::string field;
    std::vector<std::string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}
