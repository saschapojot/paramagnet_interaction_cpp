#include "dbExchangeModel.hpp"

int main(int argc, char **argv) {

    double T=3.1;
    auto model= dbExchangeModel(T);
    std::vector<double>sCurr{1,1,1,1,1,-1,1,1,-1,-1};
    model.oneStepEig(sCurr);


}