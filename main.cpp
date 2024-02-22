#include "dbExchangeModel.hpp"


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }

    double T = std::stod(argv[1]);
    auto model = dbExchangeModel(T);
    bool ferro = false;
    int lag = -1;
    int totalLoopEq = -1;

    model.reachEqMC(ferro, lag, totalLoopEq);
    if (!ferro and lag > 0) {
        model.executionMC(lag, totalLoopEq);
    }
//    model.data2File(model.record);


    // Display the loaded data

    return 0;
}
