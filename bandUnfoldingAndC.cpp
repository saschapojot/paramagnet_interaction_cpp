#include "s2Eigen.hpp"


int main(int argc, char **argv) {
    std::string name0="T0.111616";
    std::string name3="T2.175757575757576";
    std::string name4="T3.100000";
    auto rd=reader(1,name0);
    auto model=dbExchangeModel(rd.T);


    rd.searchFiles();
    rd.sortFiles();
    rd.parseCHiFile();
    rd.parse_sAllDir();
    rd.parse_EAllDir();
//    rd.parse_muAllDir();
    rd.fillZeWeights(model);

    rd.computeAllWEAllWE2();
    std::cout<<rd.dbeta_epsilon(10)<<std::endl;


}
