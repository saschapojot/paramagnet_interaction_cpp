#include <iostream>
#include "dbExchangeModel.hpp"




int main(int argc, char **argv) {


    Eigen::MatrixXd X{{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9,0}};
    auto Y=X*1.2;



    return 0;
}
