//
// Created by polya on 2/6/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
#define PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <numbers>
#include<complex>
using namespace std::complex_literals;
using mat10c = Eigen::Matrix<std::complex<double>, 10, 10>;
using mat2c = Eigen::Matrix<std::complex<double>, 2, 2>;
using mat20c = Eigen::Matrix<std::complex<double>, 20, 20>;
const auto PI=std::numbers::pi;

class dbExchangeModel {
public:
    dbExchangeModel(double temperature) {
        this->T = temperature;

        //construct SBZ values
        for (int j = 0; j < this->M; j++) {
            this->KSupValsAll.push_back(2 * PI * double(j) / (double(this->L * this->M)));
        }

        this->constructhPart();

        this->I10upup= this->kron(mat10c::Identity(),upup)*this->g;
        this->I10downdown=this->kron(mat10c::Identity(),downdown)*this->g;


    }


public:
    int part = 0; // a group of computations
    int L = 10;// length of a supercell
    int M = 20;// number of supercells

    double Ne = static_cast<double>(M);// electron number

    double t = 0.4;// hopping coefficient
    double J = -1;// exchange interaction

    double g = 0;// coupling coefficient

    double T = 1;// temperature
    std::vector<double> KSupValsAll;//all the values in SBZ
    std::vector<mat20c> preComputedHamiltonianPart;//part of the SBZ Hamiltonian that can be precomputed
    mat20c I10upup; // electron spin interaction part
    mat20c I10downdown;// electron spin interaction part




//    using mat10c= Eigen::Matrix< std::complex<double>, 10, 10 >;
//    using mat2c=Eigen::Matrix< std::complex<double>, 2, 2 >;
    mat2c I2{{1, 0},
             {0, 1}};
    mat2c upup{{1,0},{0,0}};
    mat2c downdown{{0,0},{0,1}};
    mat20c I20=mat20c::Identity()*J;
//    using mat20c=Eigen::Matrix< std::complex<double>, 20, 20 >;
public:
    mat20c kron(const mat10c &, const mat2c &);// perform Kronecker product
    void constructhPart();// construct the precomputable part of Hamiltonian
    mat20c hEig(std::vector<double> &s, const int &j); //using s value from each MC step to construct SBZ Hamiltonian, and compute its eigenvalue problem


};


#endif //PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
