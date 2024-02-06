//
// Created by polya on 2/6/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
#define PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>


using mat10c= Eigen::Matrix< std::complex<double>, 10, 10 >;
using mat2c=Eigen::Matrix< std::complex<double>, 2, 2 >;
using mat20c=Eigen::Matrix< std::complex<double>, 20, 20 >;
class dbExchangeModel{
public:
    dbExchangeModel(double temperature){
        this->T=temperature;

    }


public:
    int part=0; // a group of computations
    int L=0;// length of a supercell
    int M=20;// number of supercells

    double Ne=static_cast<double>(M);// electron number

    double t=0.4;// hopping coefficient
    double J=-1;// exchange interaction

    double g=0;// coupling coefficient

    double T=0;// temperature

//    using mat10c= Eigen::Matrix< std::complex<double>, 10, 10 >;
//    using mat2c=Eigen::Matrix< std::complex<double>, 2, 2 >;
    mat2c I2{{1,0},{0,1}};
//    using mat20c=Eigen::Matrix< std::complex<double>, 20, 20 >;

    mat20c constructhPart();// construct one part of Hamltonian
    mat20c  hEig(const int &j, std::vector<float>);














};





#endif //PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
