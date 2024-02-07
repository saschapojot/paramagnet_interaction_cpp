//
// Created by polya on 2/6/24.
//
#include "dbExchangeModel.hpp"


///
/// @param LHS electron part of the matrix
/// @param RHS spin part of the matrix
/// @return Kronecker product of the 2 matrices
mat20c dbExchangeModel::kron(const mat10c &LHS, const mat2c &RHS) {

    mat20c retMat = mat20c::Zero();
    for (int i = 0; i < LHS.rows(); i++) {
        for (int j = 0; j < LHS.cols(); j++) {
            retMat.block<2, 2>(i * RHS.rows(), j * RHS.rows()) = LHS(i, j) * RHS;
        }
    }

    return retMat;

}


void dbExchangeModel::constructhPart() {

    for (int j = 0; j < this->M; j++) {
        double K = this->KSupValsAll[j];
        mat20c retMat = mat20c::Zero();

        mat10c tmp0 = mat10c::Zero();
        tmp0(this->L - 1, 0) = -t * std::exp(L * K * 1i);

        retMat += this->kron(tmp0, I2);

        mat10c tmp1 = mat10c::Zero();
        tmp1(0, this->L - 1) = -t * std::exp(-L * K * 1i);

        retMat += this->kron(tmp1, I2);

        for (int l = 1; l <= this->L - 1; l++) {
            mat10c tmp = mat10c::Zero();
            tmp(l - 1, l) = -t;
            tmp(l, l - 1) = -t;
            retMat += this->kron(tmp, I2);
        }

        this->preComputedHamiltonianPart.push_back(retMat);


    }


}


///
/// @param s spin values for a MC step
/// @param j index of one SBZ value
/// @return j, eigenvals, eigenvects
std::tuple<int, eigVal20, vecVal20> dbExchangeModel::hEig(const std::vector<double> &s, const int &j) {
    mat20c part = this->preComputedHamiltonianPart[j];
    double sSum = 0;
    for (int l = 0; l < L; l++) {
        sSum += s[l] * s[(l + 1) % L];

    }

    part += this->I20 * sSum;


    mat20c I10upupCopy = this->I10upup;
    for (int j = 0; j < L; j++) {
        I10upupCopy(2 * j, 2 * j) *= s[j];
        I10upupCopy(2 * j + 1, 2 * j + 1) *= s[j];
    }

    mat20c I10downdownCopy = this->I10downdown;
    for (int j = 0; j < L; j++) {
        I10downdownCopy(2 * j, 2 * j) *= s[j];
        I10downdownCopy(2 * j + 1, 2 * j + 1) *= s[j];
    }


    mat20c wholeh = part + I10upupCopy - I10downdownCopy;// h(K,s)
    this->eigSolution.compute(wholeh);
    eigVal20 vals = eigSolution.eigenvalues();
    vecVal20 vecs = eigSolution.eigenvectors();

    return std::make_tuple(j, vals, vecs);


}


///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
std::vector<std::tuple<int, eigVal20, vecVal20>> dbExchangeModel::s2Eig(const std::vector<double> &s) {
//TODO: there are mistakes in this function
    std::vector<std::tuple<int, eigVal20, vecVal20>> retVec;

    std::vector<std::future<std::tuple<int, eigVal20, vecVal20>>> futAll(this->M);

    for (auto j: this->KSupIndsAll) {
        futAll[j] = std::async(std::launch::async, [this, &s, j]() {
            return this->hEig(s, j);
        });
    }

    for (auto i = 0; i < futAll.size(); i++) {
        std::tuple<int, eigVal20, vecVal20> rst = futAll[i].get();
        retVec.push_back(rst);
    }

    return retVec;
}


///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
std::vector<std::tuple<int, eigVal20, vecVal20>> dbExchangeModel::s2EigSerial(const std::vector<double> &s) {

    std::vector<std::tuple<int, eigVal20, vecVal20>> retVec;

    for (auto j: this->KSupIndsAll) {
        retVec.push_back(this->hEig(s, j));
    }

    return retVec;


}


///
/// @param eigResults eigenvalue results from s2EigSerial()
/// @return all eigenvalues placed in a vector
std::vector<double> dbExchangeModel::combineFromEig(
        const std::vector<std::tuple<int, eigVal20, vecVal20>> &eigResults) {

    std::vector<double> retEVec;

    for (auto elem: eigResults) {
        eigVal20 oneVecTmp = std::get<1>(elem);
        for (auto e: oneVecTmp) {
            retEVec.push_back(e);
        }
    }

    return retEVec;


}


///
/// @param func a monotonically increasing function
/// @return find root func=0
double dbExchangeModel::bisection_method(std::function<double(double)> func) {

    double a = -1; //left searching point
    double b = 1;// right searching point
    double leftSearchLenth = 1;
    double rightSearchLength = 1;
    //search for left end
    while (func(a) > 0) {
        a -= leftSearchLenth;
        leftSearchLenth *= 2;
    }
    //search for right end
    while (func(b) < 0) {
        b += rightSearchLength;
        rightSearchLength *= 2;
    }

    int maxiter = 10000;
    double tol = 1e-8;


    int counter = 0;

    //===================================
    //bisection search
    while (counter <= maxiter) {
        double midpoint = (a + b) / 2.0;
        double midVal = func(midpoint);
        if (std::abs(midVal) < 1e-11 or (b - a) / 2 < tol) {

            return midpoint;
        }

        if (func(a) * midVal < 0) {
            b = midpoint;
        } else {
            a = midpoint;
        }
        counter++;


    }
    //==========================================
    return (a + b) / 2.0;


}














///
/// @param EVec a vector of all eigenvalues given one s
/// @return chemical potential
double dbExchangeModel::chemicalPotential(const std::vector<double>& EVec) {

    // local function
   auto muf=[&EVec,this](double mu)-> double {

       std::vector<double> occ;
       for(auto e:EVec){
           double tmp=1/(std::exp(this->beta*(e-mu))+1);
           occ.push_back(tmp);
       }
       double  sum_of_elems = std::accumulate(occ.begin(), occ.end(),
                                           decltype(occ)::value_type(0.0));
       return sum_of_elems;

   };

   double muVal=this->bisection_method(muf);

   return muVal;


}


double dbExchangeModel::avgEnergy(const std::vector<double> &EVec) {
    double muVal= this->chemicalPotential(EVec);

    double weightedE;
    for (auto e: EVec){
        double tmp=1/(std::exp(this->beta*(e-muVal))+1)*e;
        weightedE+=tmp;

    }
    return weightedE;


}


