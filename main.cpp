
#include "dbExchangeModel.hpp"
#include <fstream>



int main(int argc, char **argv) {

//   auto db=dbExchangeModel(2.1);
//   db.executionMC();
    std::vector<double> EAll{0.1,-10,9.1,4.3};
    std::ofstream outF("EAll",std::ios::out| std::ios::binary);

    msgpack::pack(outF,EAll);
    outF.seekp(0);
    outF.close();


    return 0;
}
