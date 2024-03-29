//
// Created by polya on 2/6/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
#define PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
#include <iostream>
#include <Eigen/Dense>
//#include <unsupported/Eigen/KroneckerProduct>

#include <Eigen/Eigenvalues>
#include <numbers>
#include<complex>

#include <tuple>

#include <vector>

#include <memory>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>
#include <random>
#include <chrono>
#include <msgpack.hpp>

#include <fstream>
#include <boost/filesystem.hpp>

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <cstdlib>
#include <regex>
#include <sstream>

using namespace std::complex_literals;
//using mat10c = Eigen::Matrix<std::complex<double>, 10, 10>;
//using mat2c = Eigen::Matrix<std::complex<double>, 2, 2>;
//using mat20c = Eigen::Matrix<std::complex<double>, 20, 20>;
using vec20c=Eigen::Vector<std::complex<double>,20>;
const auto PI=std::numbers::pi;
//using eigVal20=Eigen::SelfAdjointEigenSolver<mat20c>::RealVectorType;
//using vecVal20=Eigen::SelfAdjointEigenSolver<mat20c>::EigenvectorsType;
using matc=Eigen::MatrixXcd;
//using vecc=Eigen::VectorXcd;
using eigVal=Eigen::SelfAdjointEigenSolver<matc>::RealVectorType;
using vecVal=Eigen::SelfAdjointEigenSolver<matc>::EigenvectorsType;

namespace fs = boost::filesystem;
class oneEigSolution{
public:
    oneEigSolution(const int&ind,const std::vector<double>& vals, const std::vector<std::complex<double>> &vecs){
        this->j=ind;
        this->eigVals=std::vector<double>(vals);
        this->eigVecs=std::vector<std::complex<double>>(vecs);
    }
private:
    friend class boost::serialization::access; // Allow serialization to access private members

//    template<class Archive>
//    void serialize(Archive& ar, const unsigned int version) {
//        ar & BOOST_SERIALIZATION_NVP(j);
//        ar & BOOST_SERIALIZATION_NVP(eigVals);
//        ar & BOOST_SERIALIZATION_NVP(eigVecs);
//    }

public:
    //the data structure holding eigenvalue solution for one mc step
    int j;
    std::vector<double > eigVals;
    std::vector<std::complex<double>> eigVecs;


};


class dataholder{
public:

    std::vector<std::vector<double>>sAll;//to be stored
    std::vector<double>EAll;//to be stored
    std::vector<double>muAll;//to be stored

//    std::vector<std::vector<std::tuple<int,eigVal20 ,vecVal20>>> eigRstAll;//to be stored after converting to flattenedEigSolution
//    std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>>flattenedEigSolution;
//    std::vector<oneEigSolution> multipleSolutions;

public:
    ///
/// @return flattened value, Eigen datatypes to std for serialization
    void flattenEigData();

    ///
    /// @param filename xml file
    /// name of eigen-solutions
//     void saveEigToXML(const std::string &filename);


    ///
    /// @param filename xml file name of vec
    ///@param vec vector to be saved
   static  void saveVecToXML(const std::string &filename,const std::vector<double> &vec);

    ///
    /// @param filename  xml file name of vecvec
    /// @param vecvec vector<vector> to be saved
    static void saveVecVecToXML(const std::string &filename,const std::vector<std::vector<double>> &vecvec);







};


class dbExchangeModel {
public:
    dbExchangeModel(double temperature,const int &groupNum,const int& rowNum) {
        this->T = temperature;
        this->beta=1/temperature;
        this->group=groupNum;
        this->row=rowNum;
        this->parseCSV(groupNum,rowNum);
        this->Ne=static_cast<double >(this->M)*2*L*f;
        this->I2L=matc::Identity(2*L,2*L);


        //construct SBZ values
        for (int j = 0; j < this->M; j++) {
            this->KSupValsAll.push_back(2 * PI * double(j) / (double(this->L * this->M)));
            this->KSupIndsAll.push_back(j);
        }

        this->constructhPart();


        this->I10upup= this->kron(matc::Identity(L,L),upup)*this->g;
        this->I10downdown=this->kron(matc::Identity(L,L),downdown)*this->g;
//        std::cout<<"L="<<L<<", M="<<M<<", Ne="<<Ne<<", t="<<t
//        <<", J="<<J<<", g="<<g<<", T="<<T<<", f="<<f<<", Ne="<<Ne<<std::endl;



    }


public:
//    dataholder record;
    int group = 0; // a group of computations
    int row=0;//select a row in csv
    int L = 0;// length of a supercell
    int M = 0;// number of supercells
//    std::string outDir;

    double Ne =0;// static_cast<double>(M);// electron number

    double t = 0;// hopping coefficient
    double J = 0;// exchange interaction

    double g = 0;// coupling coefficient
    double f=0;//filling ratio
    double T;// temperature
    int lastFileNum=0;
    double beta;
    std::vector<double> KSupValsAll;//all the values in SBZ
    std::vector<int> KSupIndsAll;//[0,1,...,M-1]
    std::vector<matc> preComputedHamiltonianPart;//part of the SBZ Hamiltonian that can be precomputed
    matc I10upup; // electron spin interaction part
    matc I10downdown;// electron spin interaction part


    std::vector<double>sRange{-1,1};
    int sweepNumInOneFlush=3000;// flush the results to python every sweepNumInOneFlush*L iterations
    int flushMaxNum=700;
    int dataNumTotal=8000;
    Eigen::SelfAdjointEigenSolver<matc> eigSolution;// solver for hermitian matrices
//    std::vector<double>EAvgAll;//to be stored
//    std::vector<std::vector<double>>sAll;//to be stored
//    std::vector<double>chemPotAll;//to be stored
//    std::vector<std::tuple<int,eigVal20 ,vecVal20>> eigRstAll;//to be stored




//    using mat10c= Eigen::Matrix< std::complex<double>, 10, 10 >;
//    using mat2c=Eigen::Matrix< std::complex<double>, 2, 2 >;
    matc I2{{1, 0},
             {0, 1}};

    matc upup{{1,0},{0,0}};
    matc downdown{{0,0},{0,1}};
    matc I2L;

//    using mat20c=Eigen::Matrix< std::complex<double>, 20, 20 >;
public:
    matc kron(const matc &, const matc &);// perform Kronecker product, RHS is 2 by 2
    void constructhPart();// construct the precomputable part of Hamiltonian//executed in constructor

    /// initialize parameters using  values from  csv
    /// @param groupNum group number of csv
    /// @param rowNum row number in csv
    void parseCSV(const int groupNum,const int& rowNum);
    ///
/// @param s spin values for a MC step
/// @param j index of one SBZ value
/// @return j, eigenvals, eigenvects

   std::tuple<int,eigVal ,vecVal>  hEig(const std::vector<double> &s, const int &j); //using s value from each MC step to construct SBZ Hamiltonian, and compute its eigenvalue problem



    ///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
   std::vector<std::tuple<int,eigVal ,vecVal>> s2Eig(const std::vector<double> &s);//(parallel version) using s value from each MC step to construct SBZ Hamiltonian,
   // and compute its eigenvalue problem for each K in SBZ
   ///
   /// @param s spin values in a MC step
   /// @return eigenvalues and eigenvectors for all values in SBZ
public:
   std::vector<std::tuple<int,eigVal ,vecVal>> s2EigSerial(const std::vector<double> &s);//serial version of s2Eig()

   ///
   /// @param eigResults eigenvalue results from s2EigSerial()
   /// @return all eigenvalues placed in a vector
   std::vector<double> combineFromEig(const std::vector<std::tuple<int,eigVal ,vecVal>>& eigResults);

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
   std::vector<double> avgEnergy(const std::vector<double>& EVec);

   ///
   /// @param cmd python execution string
   /// @return signal from the python
   static std::string execPython(const char* cmd);

   ///
   /// @param lag decorrelation length
   /// @param loopEq total loop numbers in reaching equilibrium
   void executionMC(const int& lag,const int & loopEq);// mc simulation without inquiring equilibrium after reaching equilibrium

   ///
   /// @param ferro is ferromagnetic
   /// @param lag decorrelation length
   /// @param loopTotal total mc steps
   void reachEqMC(bool& ferro, int &lag, int&loopTotal);// mc simulation while communicating with python to inquire equilibrium

   ///
   /// @param record holding data
//   void data2File(const dataholder& record);//data serialization

//    template<typename T>void serializationViaFStream(const T& values, const std::string & fileName);
//
//    void serializationViaFStream(const std::vector<std::vector<double>>& vecvec,const std::string & fileName);
//
//    void serializationViaFStream(const std::vector<double>& vec,const std::string & fileName);
//
    void serializationViaFStream(const std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>>& vecvectuple,const std::string & fileName);

    template<class T>
    static void printVec(const std::vector<T>& vec){
        for(int i=0;i<vec.size()-1;i++){
            std::cout<<vec[i]<<",";
        }
        std::cout<<vec[vec.size()-1]<<std::endl;
    }

    /// eigenvalues and mu for a s vector
    /// @param sCurr s vector
    void oneStepEig(const std::vector<double>& sCurr);


};


#endif //PARAMAGNET_INTERACTION_CPP_DBEXCHANGEMODEL_HPP
