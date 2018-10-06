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
    std::ifstream ifs("n63412-il.csv");
    std::string line;
    std::vector<double> x;
    std::vector<double> y;


    while (std::getline(ifs, line)) {

        std::vector<std::string> strvec = split(line, ',');

        for (int i=0; i<strvec.size();i++){
            if(i%2 == 0){
                x.push_back(std::stod(strvec[i]));
            }
            else{
                y.push_back(std::stod(strvec[i]));
            }
        }
    }
/*    for(auto itr=y.begin(); itr<y.end(); itr++){
        std::cout << *itr << std::endl;
    }*/
    auto nin = static_cast<int> (x.size());

    //////////////////////////////////////////
    int n, k;
    double dl1, dl2;
    std::vector<double> xout, yout;
    //////////////////////////////////////////


    //grcv2d(nin, x, y, n, k, dl1, dl2, xout, yout);
    tinput();

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
