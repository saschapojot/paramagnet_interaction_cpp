//
// Created by polya on 2/6/24.
//
#include "dbExchangeModel.hpp"


///
/// @param LHS electron part of the matrix
/// @param RHS spin part of the matrix
/// @return Kronecker product of the 2 matrices
mat20c dbExchangeModel::kron(const mat10c & LHS, const mat2c & RHS) {

    mat20c retMat=mat20c::Zero();
    for (int i=0;i<LHS.rows();i++){
        for (int j=0;j<LHS.cols();j++){
            retMat.block<2,2>(i*RHS.rows(),j*RHS.rows())=LHS(i,j)*RHS;
        }
    }

    return retMat;

}


void dbExchangeModel::constructhPart(){

    for(int j=0;j<this->M;j++){
        double K=this->KSupValsAll[j];
        mat20c retMat=mat20c::Zero();

        mat10c tmp0=mat10c::Zero();
        tmp0(this->L-1,0)=-t*std::exp(L*K*1i);

        retMat+=this->kron(tmp0,I2);

        mat10c tmp1=mat10c::Zero();
        tmp1(0,this->L-1)=-t*std::exp(-L*K*1i);

        retMat+=this->kron(tmp1,I2);

        for(int l=1;l<=this->L-1;l++){
            mat10c tmp=mat10c::Zero();
            tmp(l-1,l)=-t;
            tmp(l,l-1)=-t;
            retMat+=this->kron(tmp,I2);
        }

        this->preComputedHamiltonianPart.push_back(retMat);


    }





}




mat20c dbExchangeModel::hEig(std::vector<double> &s, const int &j) {
    mat20c part=this->preComputedHamiltonianPart[j];
    double sSum=0;
    for(int l=0;l<L;l++){
        sSum+=s[l]*s[(l+1)%L];

    }

    part+=this->I20*sSum;






}









































