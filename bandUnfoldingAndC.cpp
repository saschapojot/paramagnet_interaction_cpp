#include "s2Eigen.hpp"


int main(int argc, char **argv) {
    std::string name0="T0.111616";
    std::string name3="T2.175757575757576";

    auto rd=reader(1,name3);
    auto model=dbExchangeModel(rd.T);
    rd.searchFiles();
    rd.sortFiles();
    rd.parseCHiFile();
    rd.parse_sAllDir();
    rd.parse_EAllDir();
    rd.parse_muAllDir();
    rd.fillZeWeights(model);


}
