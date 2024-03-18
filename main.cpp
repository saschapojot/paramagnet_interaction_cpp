#include "dbExchangeModel.hpp"


int main(int argc, char **argv) {
    if (argc != 4) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }

    double T = std::stod(argv[1]);
    int groupNum=std::stoi(argv[2]);
    int rowNum=std::stoi(argv[3]);
    auto model = dbExchangeModel(T,groupNum,rowNum);
    model.parseCSV(groupNum,rowNum);
    bool ferro = false;
    int lag = -1;
    int totalLoopEq = -1;

    model.reachEqMC(ferro, lag, totalLoopEq);
    if (!ferro and lag > 0) {
        model.executionMC(lag, totalLoopEq);
    }



    // Display the loaded data

    return 0;
}
