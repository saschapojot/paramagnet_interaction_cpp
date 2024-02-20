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
//std::vector<std::tuple<int, eigVal20, vecVal20>> dbExchangeModel::s2Eig(const std::vector<double> &s) {
////TODO: there are mistakes in this function
//    std::vector<std::tuple<int, eigVal20, vecVal20>> retVec;
//
//    std::vector<std::future<std::tuple<int, eigVal20, vecVal20>>> futAll(this->M);
//
//    for (auto j: this->KSupIndsAll) {
//        futAll[j] = std::async(std::launch::async, [this, &s, j]() {
//            return this->hEig(s, j);
//        });
//    }
//
//    for (auto i = 0; i < futAll.size(); i++) {
//        std::tuple<int, eigVal20, vecVal20> rst = futAll[i].get();
//        retVec.push_back(rst);
//    }
//
//    return retVec;
//}


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
        if (std::abs(midVal) < 1e-9 or (b - a) / 2 < tol) {

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
double dbExchangeModel::chemicalPotential(const std::vector<double> &EVec) {

    // local function
    auto muf = [&EVec, this](double mu) -> double {

        std::vector<double> occ;
        for (auto e: EVec) {
            double tmp = 1 / (std::exp(this->beta * (e - mu)) + 1);
            occ.push_back(tmp);
        }
        double sum_of_elems = std::accumulate(occ.begin(), occ.end(),
                                              decltype(occ)::value_type(0.0));
        return sum_of_elems - this->Ne;

    };

    double muVal = bisection_method(muf);

    return muVal;


}


std::vector<double> dbExchangeModel::avgEnergy(const std::vector<double> &EVec) {
    double muVal = this->chemicalPotential(EVec);

    double weightedE;
    for (auto e: EVec) {
        double tmp = 1 / (std::exp(this->beta * (e - muVal)) + 1) * e;
        weightedE += tmp;

    }
    std::vector<double> retVec{weightedE, muVal};
    return retVec;


}


void dbExchangeModel::executionMC() {
    //init
//    dataholder record=dataholder();//records all data

    std::random_device rd;
    std::uniform_int_distribution<int> indsAll(0, 1);
    std::uniform_int_distribution<int> flipInds(0, L - 1);


    std::vector<double> sCurr;//init s
    for (int i = 0; i < this->L; i++) {
        sCurr.push_back(this->sRange[indsAll(rd)]);
    }
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);


    auto tripleCurr = this->s2EigSerial(sCurr);//init eig result
    std::vector<double> EVec = this->combineFromEig(tripleCurr);//init EVec
    auto EAndMuCurr = this->avgEnergy(EVec);// init E and mu
    double EAvgCurr = EAndMuCurr[0];
    double muCurr = EAndMuCurr[1];
    int flipNum = 0;
    int noFlipNum = 0;


    //start of MC
//    int sweep=20000;
//    int counter=sweep*this->L;


    const auto tMCStart{std::chrono::steady_clock::now()};

    int maxSweepIter = this->sweepNumInOneFlush * this->flushMaxNum;
    for (int swp = 0; swp < maxSweepIter; swp++) {
        std::unique_ptr<dataholder> record_ptr;
        int loopStart = swp * L;

        for (int i = 0; i < this->sweepNumInOneFlush * this->L; i++) {
            //perform a flip
            auto sNext = std::vector<double>(sCurr);
            int flipTmpInd = flipInds(rd);
            sNext[flipTmpInd] *= -1;
            auto tripleNext = this->s2EigSerial(sNext);
            auto EVecNext = this->combineFromEig(tripleNext);
            auto EAndMuNext = this->avgEnergy(EVecNext);
            double EAvgNext = EAndMuNext[0];
            double muNext = EAndMuCurr[1];
            double DeltaE = (EAvgNext - EAvgCurr) / this->M;
            //decide if flip is accepted
            if (DeltaE <= 0) {
                sCurr = std::vector<double>(sNext);
                tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                EAvgCurr = EAvgNext;
                muCurr = muNext;
                flipNum++;

            } else {
                double r = distUnif01(e2);
                if (r < std::exp(-this->beta * DeltaE)) {
                    sCurr = std::vector<double>(sNext);
                    tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                    EAvgCurr = EAvgNext;
                    muCurr = muNext;
                    flipNum++;

                } else {
                    noFlipNum++;
                }


            }

            record_ptr->sAll.push_back(sCurr);
            record_ptr->EAll.push_back(EAvgCurr);
            record_ptr->muAll.push_back(muCurr);
            record_ptr->eigRstAll.push_back(tripleCurr);


        }
        int loopEnd = loopStart + this->sweepNumInOneFlush * this->L - 1;

        record_ptr->flattenEigData();


    }
//    for (int i=0;i<counter;i++) {
//
//        //perform a flip
//
//        auto sNext =std::vector<double>(sCurr);
//
//        int flipTmpInd = flipInds(rd);
//        sNext[flipTmpInd] *= -1;
//        auto tripleNext = this->s2EigSerial(sNext);
//        auto EVecNext = this->combineFromEig(tripleNext);
//
//        auto EAndMuNext = this->avgEnergy(EVecNext);
//        double EAvgNext = EAndMuNext[0];
//        double muNext = EAndMuCurr[1];
//        double DeltaE = (EAvgNext - EAvgCurr) / this->M;
//        //decide if flip is accepted
//        if (DeltaE <= 0) {
//            sCurr = std::vector<double>(sNext);
//            tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
//            EAvgCurr = EAvgNext;
//            muCurr = muNext;
//            flipNum++;
//
//        } else {
//            double r = distUnif01(e2);
//            if (r < std::exp(-this->beta * DeltaE)) {
//                sCurr = std::vector<double>(sNext);
//                tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
//                EAvgCurr = EAvgNext;
//                muCurr = muNext;
//                flipNum++;
//
//            } else {
//                noFlipNum++;
//            }
//
//
//        }
//        record_ptr->sAll.push_back(sCurr);
//        record_ptr->EAll.push_back(EAvgCurr);
//        record_ptr->muAll.push_back(muCurr);
//        record_ptr->eigRstAll.push_back(tripleCurr);
//
//        if (i%5000==0){
//            std::cout<<"flip "<<i<<std::endl;
//            const auto tMC5000{std::chrono::steady_clock::now()};
//            const std::chrono::duration<double> elapsed_seconds5000{tMC5000 - tMCStart};
//            std::cout<<elapsed_seconds5000.count()/3600.0<<" h"<<std::endl;
//        }
//
//
//    }

    //end of mc
//    record.flattenEigData();


    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    std::cout << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    std::cout << "flip number: " << flipNum << std::endl;
    std::cout << "no flip number: " << noFlipNum << std::endl;


}

///
/// @return flattened value, Eigen datatypes to std for serialization
void dataholder::flattenEigData() {


    for (auto const &vecAllForOneS: this->eigRstAll) {
//        std::vector<std::tuple<int, std::vector<double>, std::vector<std::complex<double>> >> flattenedVecAllForOneS;// eigensolutions for all j=0,1,...,M-1


        for (auto const &tp: vecAllForOneS) {
            int j = std::get<0>(tp);
            eigVal20 eigValsTmp = std::get<1>(tp);
            vecVal20 eigVecsTmp = std::get<2>(tp);//stored in column major format

            std::vector<double> stdEigValsTmp;
            std::vector<std::complex<double>> stdEigVecsTmp;

            for (auto const &x: eigValsTmp) {
                stdEigValsTmp.emplace_back(x);
            }
            for (auto const &x: eigVecsTmp.reshaped()) {
                stdEigVecsTmp.emplace_back(x);
            }
//            std::tuple<int, std::vector<double>, std::vector<std::complex<double>>> flattenedTuple = std::make_tuple(j,
//            stdEigValsTmp,
//                    stdEigVecsTmp);
            this->multipleSolutions.push_back(oneEigSolution(j, stdEigValsTmp, stdEigVecsTmp));


        }


    }


}

//template<typename T>
//void dbExchangeModel::serializationViaFStream(const T &values, const std::string & fileName) {
//    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
//
//    msgpack::pack(outFTmp, values);
//    outFTmp.seekp(0);
//    outFTmp.close();
//
//
//}

void
dbExchangeModel::serializationViaFStream(const std::vector<std::vector<double>> &vecvec, const std::string &fileName) {
    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
    msgpack::pack(outFTmp, vecvec);
    outFTmp.seekp(0);
    outFTmp.close();


}


void dbExchangeModel::serializationViaFStream(const std::vector<double> &vec, const std::string &fileName) {
    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
    msgpack::pack(outFTmp, vec);
    outFTmp.seekp(0);
    outFTmp.close();


}

void dbExchangeModel::serializationViaFStream(
        const std::vector<std::vector<std::tuple<int, std::vector<double>, std::vector<std::complex<double>> >>> &vecvectuple,
        const std::string &fileName) {
    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
    msgpack::pack(outFTmp, vecvectuple);
    outFTmp.seekp(0);
    outFTmp.close();


}


void dbExchangeModel::data2File(const dataholder &record) {
    namespace fs = std::filesystem;

    //output folder
//    std::string outDir="./part"+std::to_string(this->part);
//    if(! fs::is_directory(outDir)|| !fs::exists(outDir)){
//        fs::create_directory(outDir);
//    }

    std::string outSubDir = "./part" + std::to_string(this->part) + "/";

    if (!fs::is_directory(outSubDir) || !fs::exists(outSubDir)) {
        fs::create_directory(outSubDir);
    }


    std::string prefix =
            "T" + std::to_string(this->T) + "t" + std::to_string(this->t) + "J" + std::to_string(this->J) + "g" +
            std::to_string(this->g) + "part" + std::to_string(this->part);
    //output sAll
    std::string outsAllName = outSubDir + prefix + ".sAll";
    serializationViaFStream(record.sAll, outsAllName);

    //output EAll
    std::string outEAllName = outSubDir + prefix + ".EAll";
    serializationViaFStream(record.EAll, outEAllName);

    //output muAll
    std::string outMuAllName = outSubDir + prefix + ".muAll";
    serializationViaFStream(record.muAll, outMuAllName);

    //output flattened solution
//    std::string outFlEigSolName=outSubDir+prefix+".flattenedEigSolution";
//    serializationViaFStream(record.flattenedEigSolution,outFlEigSolName);


}

//template void dbExchangeModel::serializationViaFStream<std::vector<std::vector<double>>>(const std::vector<std::vector<double>> &values, const std::string &fileName);
//
//template void dbExchangeModel::serializationViaFStream<std::vector<double>>(const std::vector<double> &values, const std::string &fileName);
//
//template void dbExchangeModel::serializationViaFStream<std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>>>(const std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>> &values, const std::string &fileName);