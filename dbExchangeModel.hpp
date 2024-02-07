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
#include <tuple>
#include<thread>
#include <vector>
#include <thread>
#include <future>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>

using namespace std::complex_literals;
using mat10c = Eigen::Matrix<std::complex<double>, 10, 10>;
using mat2c = Eigen::Matrix<std::complex<double>, 2, 2>;
using mat20c = Eigen::Matrix<std::complex<double>, 20, 20>;
const auto PI=std::numbers::pi;
using eigVal20=Eigen::SelfAdjointEigenSolver<mat20c>::RealVectorType;
using vecVal20=Eigen::SelfAdjointEigenSolver<mat20c>::EigenvectorsType;

class dataholder{
public:
    std::vector<double>EAvgAll;//to be stored
    std::vector<std::vector<double>>sAll;//to be stored
    std::vector<double>chemPotAll;//to be stored
    std::vector<std::tuple<int,eigVal20 ,vecVal20>> eigRstAll;//to be stored


};


class dbExchangeModel {
public:
    dbExchangeModel(double temperature) {
        this->T = temperature;

        //construct SBZ values
        for (int j = 0; j < this->M; j++) {
            this->KSupValsAll.push_back(2 * PI * double(j) / (double(this->L * this->M)));
            this->KSupIndsAll.push_back(j);
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

    double g = 0.01;// coupling coefficient

    double T = 1;// temperature
    double beta=1/T;
    std::vector<double> KSupValsAll;//all the values in SBZ
    std::vector<int> KSupIndsAll;//[0,1,...,M-1]
    std::vector<mat20c> preComputedHamiltonianPart;//part of the SBZ Hamiltonian that can be precomputed
    mat20c I10upup; // electron spin interaction part
    mat20c I10downdown;// electron spin interaction part
    int burnInEst=90000;
    Eigen::SelfAdjointEigenSolver<mat20c> eigSolution;// solver for hermitian matrices
    std::vector<double>EAvgAll;//to be stored
    std::vector<std::vector<double>>sAll;//to be stored
    std::vector<double>chemPotAll;//to be stored
    std::vector<std::tuple<int,eigVal20 ,vecVal20>> eigRstAll;//to be stored




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


    ///
/// @param s spin values for a MC step
/// @param j index of one SBZ value
/// @return j, eigenvals, eigenvects
   std::tuple<int,eigVal20 ,vecVal20>  hEig(const std::vector<double> &s, const int &j); //using s value from each MC step to construct SBZ Hamiltonian, and compute its eigenvalue problem



    ///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
   std::vector<std::tuple<int,eigVal20 ,vecVal20>> s2Eig(const std::vector<double> &s);//(parallel version) using s value from each MC step to construct SBZ Hamiltonian,
   // and compute its eigenvalue problem for each K in SBZ
   ///
   /// @param s spin values in a MC step
   /// @return eigenvalues and eigenvectors for all values in SBZ
   std::vector<std::tuple<int,eigVal20 ,vecVal20>> s2EigSerial(const std::vector<double> &s);//serial version of s2Eig()

   ///
   /// @param eigResults eigenvalue results from s2EigSerial()
   /// @return all eigenvalues placed in a vector
   std::vector<double> combineFromEig(const std::vector<std::tuple<int,eigVal20 ,vecVal20>>& eigResults);

   ///
   /// @param func a monotonically increasing function
   /// @return find root func=0
   static double bisection_method(std::function<double(double)> func);

   ///
   /// @param EVec a vector of all eigenvalues given one s
   /// @return chemical potential
   double chemicalPotential(const std::vector<double>& EVec);

   ///
   /// @param EVec a vector of all eigenvalues given one s
   /// @return
   double avgEnergy(const std::vector<double>& EVec);

   void executionMC();



};


#endif //PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
