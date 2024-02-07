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



///
/// @param s spin values for a MC step
/// @param j index of one SBZ value
/// @return j, eigenvals, eigenvects
std::tuple<int,eigVal20 ,vecVal20 > dbExchangeModel::hEig(const std::vector<double> &s, const int &j) {
    mat20c part=this->preComputedHamiltonianPart[j];
    double sSum=0;
    for(int l=0;l<L;l++){
        sSum+=s[l]*s[(l+1)%L];

    }

    part+=this->I20*sSum;


    mat20c I10upupCopy=this->I10upup;
    for(int j=0;j<L;j++){
        I10upupCopy(2*j,2*j)*=s[j];
        I10upupCopy(2*j+1,2*j+1)*=s[j];
    }

    mat20c I10downdownCopy=this->I10downdown;
    for (int j=0;j<L;j++){
        I10downdownCopy(2*j,2*j)*=s[j];
        I10downdownCopy(2*j+1,2*j+1)*=s[j];
    }


   mat20c wholeh=part+I10upupCopy-I10downdownCopy;// h(K,s)
   this->eigSolution.compute(wholeh);
   eigVal20 vals=eigSolution.eigenvalues();
   vecVal20 vecs=eigSolution.eigenvectors();

   return std::make_tuple(j,vals,vecs);









}




///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
std::vector<std::tuple<int,eigVal20 ,vecVal20>> dbExchangeModel::s2Eig(const std::vector<double> &s) {
//TODO: there are mistakes in this function
    std::vector<std::tuple<int,eigVal20 ,vecVal20>> retVec;

    std::vector<std::future<std::tuple<int,eigVal20 ,vecVal20>>> futAll(this->M);

    for (auto j: this->KSupIndsAll){
        futAll[j]=std::async(std::launch::async,[this,&s,j](){
            return this->hEig(s, j);
        });
    }

   for (auto i=0;i<futAll.size();i++){
       std::tuple<int,eigVal20 ,vecVal20>  rst=futAll[i].get();
       retVec.push_back(rst);
   }

    return retVec;
}



///
/// @param s spin values in a MC step
/// @return eigenvalues and eigenvectors for all values in SBZ
std::vector<std::tuple<int,eigVal20 ,vecVal20>> dbExchangeModel::s2EigSerial(const std::vector<double> &s) {

    std::vector<std::tuple<int,eigVal20 ,vecVal20>> retVec;

    for(auto j: this->KSupIndsAll){
        retVec.push_back(this->hEig(s,j));
    }

    return retVec;


}



































