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
//
//    int nm=11;
//    double val=vals[nm];
//    auto vec=vecs.col(nm);
////    auto diff=wholeh*vec-val*vec;
//    std::cout<<"norm="<<vec.norm()<<std::endl;

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
        if (std::abs(midVal) < 1e-8 or (b - a) / 2 < tol) {

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
//    std::cout<<"Ne="<<Ne<<std::endl;
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
/// initialize parameters using  values from  csv
/// @param groupNum group number of csv
/// @param rowNum row number in csv
void dbExchangeModel::parseCSV(const int groupNum,const int& rowNum) {
    std::string commandToReadCSV = "python3 parseCSV.py " + std::to_string(groupNum) + " " + std::to_string(rowNum);

    std::string result = this->execPython(commandToReadCSV.c_str());

    std::regex patternL("L(\\d+)M");
    std::smatch matchL;
    if (std::regex_search(result, matchL, patternL)) {
        this->L = std::stoi(matchL[1].str());
    }

    std::regex patternM("M(\\d+)t");
    std::smatch matchM;
    if (std::regex_search(result, matchM, patternM)) {
        this->M = std::stoi(matchM[1].str());
    }

    std::regex pattern_t("t([+-]?\\d+(\\.\\d+)?)J");
    std::smatch match_t;
    if (std::regex_search(result, match_t, pattern_t)) {
        this->t = std::stod(match_t[1].str());
    }

    std::regex patternJ("J([+-]?\\d+(\\.\\d+)?)g");
    std::smatch  matchJ;
    if(std::regex_search(result,matchJ,patternJ)){
        this->J=std::stod(matchJ[1].str());
    }

    std::regex pattern_g("g([+-]?\\d+(\\.\\d+)?)");
    std::smatch match_g;
    if(std::regex_search(result,match_g,pattern_g)){
        this->g=std::stod(match_g[1].str());
    }

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

///
/// @param ferro is ferromagnetic
/// @param lag decorrelation length
/// @param loopTotal total mc steps
void dbExchangeModel::reachEqMC(bool &ferro, int &lag, int &loopTotal) {
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
//    namespace fs = boost::filesystem;

    //output directory
    std::ostringstream sObjT;
    sObjT << std::fixed;
    sObjT << std::setprecision(10);
    sObjT << T;
    std::string TStr = sObjT.str();
    std::string outDir =
            "./group" + std::to_string(this->group) + "data" + "/row" + std::to_string(this->row) + "/T" + TStr + "/";
    std::string outEAllSubDir = outDir + "EAll/";
    std::string outMuAllSubDir = outDir + "muAll/";
    std::string outSAllSubDir = outDir + "sAll/";
    std::string outEigAllSubDir = outDir + "eigAll/";


    if (!fs::is_directory(outEAllSubDir) || !fs::exists(outEAllSubDir)) {
        fs::create_directories(outEAllSubDir);
    }

    if (!fs::is_directory(outMuAllSubDir) || !fs::exists(outMuAllSubDir)) {
        fs::create_directories(outMuAllSubDir);
    }
    if (!fs::is_directory(outSAllSubDir) || !fs::exists(outSAllSubDir)) {
        fs::create_directories(outSAllSubDir);
    }

    if (!fs::is_directory(outEigAllSubDir) || !fs::exists(outEigAllSubDir)) {
        fs::create_directories(outEigAllSubDir);
    }



    //start of MC
//    int sweep=20000;
//    int counter=sweep*this->L;



    const auto tMCStart{std::chrono::steady_clock::now()};
    int counter = 0;
    int fls = 0;
    bool active = true;
//    std::regex continueRegex("continue");
    std::regex stopRegex("stop");
    std::regex wrongRegex("wrong");
    std::regex ErrRegex("Err");
    std::regex lagRegex("lag=\\s*(\\d+)");
    std::regex fileNumRegex("fileNum=\\s*(\\d+)");
    std::regex ferroRegex("ferro");
    std::regex eqRegex("equilibrium");

//    std::smatch matchContinue;
    std::smatch matchEAvgStop;
    std::smatch matchEAvgWrong;
    std::smatch matchEAvgErr;
    std::smatch matchEAvgLag;
    std::smatch matchFileNum;
    std::smatch matchEAvgFerro;
    std::smatch matchEAvgEq;

    std::smatch matchSAvgStop;
    std::smatch matchSAvgWrong;
    std::smatch matchSAvgErr;
    std::smatch matchSAvgFileNum;
    std::smatch matchSAvgLag;
    std::smatch matchSAvgFerro;
    std::smatch matchSAvgEq;



//    std::cout<<"reachEq part"<<std::endl;
    while (fls < this->flushMaxNum and active == true) {
        std::unique_ptr<dataholder> record_ptr = std::make_unique<dataholder>();
        int loopStart = fls * this->sweepNumInOneFlush * this->L;

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
            double DeltaE = (EAvgNext - EAvgCurr) / static_cast<double>(this->M);
            //decide if flip is accepted
            if (DeltaE <= 0) {
                sCurr = std::vector<double>(sNext);
                tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                EAvgCurr = EAvgNext;
                muCurr = muNext;
                flipNum++;

            } else {
                double r = distUnif01(e2);
                if (r <= std::exp(-this->beta * DeltaE)) {
                    sCurr = std::vector<double>(sNext);
                    tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                    EAvgCurr = EAvgNext;
                    muCurr = muNext;
                    flipNum++;

                } else {
                    noFlipNum++;
                }


            }
//            std::cout<<"========================"<<std::endl;
//            std::cout<<"step "<<counter<<std::endl;
//            std::cout<<"EAvg Curr="<<EAvgCurr<<std::endl;
//            std::cout<<"mu Curr="<<muCurr<<std::endl;
//            std::cout<<"sCurr=";
//            printVec(sCurr);
//            std::cout<<"============"<<std::endl;
            record_ptr->sAll.push_back(sCurr);
            record_ptr->EAll.push_back(EAvgCurr);
            record_ptr->muAll.push_back(muCurr);

//            record_ptr->eigRstAll.push_back(tripleCurr);
            counter += 1;


        }
        int loopEnd = loopStart + this->sweepNumInOneFlush * this->L - 1;

//        record_ptr->flattenEigData();
//        std::ostringstream  sObjT;
//        sObjT<<std::fixed;
//        sObjT<<std::setprecision(10);
//        sObjT<<T;
//        std::string TStr=sObjT.str();
        std::string filenameMiddle = "loopStart" + std::to_string(loopStart) +
                                     "loopEnd" + std::to_string(loopEnd) + "T" + TStr;

        std::string outEFileTmp = outEAllSubDir + filenameMiddle + ".EAll.xml";

        record_ptr->saveVecToXML(outEFileTmp, record_ptr->EAll);

        std::string outMuFileTmp = outMuAllSubDir + filenameMiddle + ".muAll.xml";

        record_ptr->saveVecToXML(outMuFileTmp, record_ptr->muAll);

        std::string outSFileTmp = outSAllSubDir + filenameMiddle + ".sAll.xml";

        record_ptr->saveVecVecToXML(outSFileTmp, record_ptr->sAll);


//        std::string outEigFileTmp = outEigAllSubDir + filenameMiddle + ".eigAll.bin";

//        this->serializationViaFStream(record_ptr->flattenedEigSolution, outEigFileTmp);

//        record_ptr->saveEigToXML(outEigFileTmp);

        const auto tflushEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{tflushEnd - tMCStart};
        std::cout << "flush " << fls << std::endl;
        std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;

        //communicate with python to inquire equilibrium

        //inquire equilibrium of EAvg
        std::string commandEAvg = "python3 checkVec.py " + outEAllSubDir;
        //inquire equilibrium of s
        std::string commandSAvg = "python3 checksVecVec.py " + outSAllSubDir;
        std::string resultEAvg;
        std::string resultSAvg;
        if (fls % 6 == 5) {
            try {
                resultEAvg = this->execPython(commandEAvg.c_str());
                resultSAvg = this->execPython(commandSAvg.c_str());

                std::cout << "EAvg message from python: " << resultEAvg << std::endl;
                std::cout << "sAvg message from python: " << resultSAvg << std::endl;

            } catch (const std::exception &e) {
                std::cerr << "Error: " << e.what() << std::endl;
                std::exit(10);
            }
            catch (...) {
                // Handle any other exceptions
                std::cerr << "Error" << std::endl;
                std::exit(11);
            }

            // parse result
            if (std::regex_search(resultEAvg, matchEAvgErr, ErrRegex)) {
                std::cout << "error encountered" << std::endl;
                std::exit(12);
            }

            if (std::regex_search(resultEAvg, matchEAvgWrong, wrongRegex)) {
                std::exit(13);
            }

//        bool ferro= false;

            if (std::regex_search(resultEAvg, matchEAvgStop, stopRegex) and
                std::regex_search(resultSAvg, matchSAvgStop, stopRegex)) {
                if (std::regex_search(resultEAvg, matchEAvgFerro, ferroRegex) and
                    std::regex_search(resultSAvg, matchSAvgFerro, ferroRegex)) {
                    active = false;
                    ferro = true;
                    std::regex_search(resultEAvg, matchFileNum, fileNumRegex);
                    std::string fileNumStr = matchFileNum.str(1);
//                    std::cout<<"fileNumStr: "<<fileNumStr<<std::endl;
                    this->lastFileNum = std::stoi(fileNumStr);
                }


            }


            if (std::regex_search(resultEAvg, matchEAvgEq, eqRegex) and
                std::regex_search(resultSAvg, matchSAvgEq, eqRegex)) {
                if (std::regex_search(resultEAvg, matchEAvgLag, lagRegex) and
                    std::regex_search(resultSAvg, matchSAvgLag, lagRegex)) {
                    std::string lagStrEAvg = matchEAvgLag.str(1);
                    std::string lagStr_s = matchSAvgLag.str(1);
                    int lag_s = std::stoi(lagStr_s);

                    int lagEAvg = std::stoi(lagStrEAvg);


                    lag = lag_s > lagEAvg ? lag_s : lagEAvg;
                    std::cout << "lag_s=" << lag_s << std::endl;
                    std::cout << "lagEAvg=" << lagEAvg << std::endl;
                    std::cout << "lag=" << lag << std::endl;
                    std::regex_search(resultEAvg, matchFileNum, fileNumRegex);
                    std::string fileNumStr = matchFileNum.str(1);
//                    std::cout<<"fileNumStr: "<<fileNumStr<<std::endl;
                    this->lastFileNum = std::stoi(fileNumStr);
                    active = false;
                }

            }

            if (std::regex_search(resultEAvg, matchEAvgStop, stopRegex) and
                std::regex_search(resultSAvg, matchSAvgEq, eqRegex)) {
                std::regex_search(resultSAvg, matchSAvgLag, lagRegex);
                std::string lagStr_s = matchSAvgLag.str(1);
                int lag_s = std::stoi(lagStr_s);
                lag = lag_s;
                std::regex_search(resultSAvg, matchSAvgFileNum, fileNumRegex);
                std::string fileNumStr = matchSAvgFileNum.str(1);
                this->lastFileNum = std::stoi(fileNumStr);
                active = false;

            }

        }//end if
        fls++;
    }//end of while loop




    std::ofstream outSummary(outDir + "summary.txt");
    loopTotal = flipNum + noFlipNum;


    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    outSummary << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    outSummary << "total sweep number: " << static_cast<int>(loopTotal / this->L) << std::endl;
    outSummary << "total loop number: " << loopTotal << std::endl;

    outSummary << "flip number: " << flipNum << std::endl;
    outSummary << "no flip number: " << noFlipNum << std::endl;
    outSummary << "lastFileNum=" << lastFileNum << std::endl;
    outSummary << "equilibrium reached: " << !active << std::endl;
    outSummary << "ferro: " << ferro << std::endl;
    outSummary << "lag=" << lag << std::endl;
    outSummary.close();
}//end of function






///
/// @param lag decorrelation length
/// @param loopEq total loop numbers in reaching equilibrium
void dbExchangeModel::executionMC(const int &lag, const int &loopEq) {

    int counter=0;
    int remainingDataNum = this->dataNumTotal - this->lastFileNum*sweepNumInOneFlush*L;

    int remainingLoopNum = remainingDataNum * lag;

    if (remainingLoopNum <= 0) {
        return;
    }
    double remainingLoopNumDB = static_cast<double>(remainingLoopNum);
    double LDB = static_cast<double>(this->L);

    int remainingSweepNum = std::ceil(remainingLoopNumDB / LDB);
    double remainingSweepNumDB = static_cast<double >(remainingSweepNum);

    double sweepNumInOneFlushDB = static_cast<double >(sweepNumInOneFlush);

    double remainingFlushNumDB = std::ceil(remainingSweepNumDB / sweepNumInOneFlushDB);
    int remainingFlushNum = static_cast<int>(remainingFlushNumDB);

    //init
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
    namespace fs = boost::filesystem;

    //output directory
    std::ostringstream  sObjT;
    sObjT<<std::fixed;
    sObjT<<std::setprecision(10);
    sObjT<<T;
    std::string TStr=sObjT.str();
    std::string outDir = "./group" + std::to_string(this->group)+"data"+"/row"+std::to_string(this->row) + "/T" + TStr + "/";
    std::string outEAllSubDir = outDir + "EAll/";
    std::string outMuAllSubDir = outDir + "muAll/";
    std::string outSAllSubDir = outDir + "sAll/";
    std::string outEigAllSubDir = outDir + "eigAll/";
    const auto tMCStart{std::chrono::steady_clock::now()};

    std::cout<<"remaining flush number: "<<remainingFlushNum<<std::endl;

    for (int fls = 0; fls < remainingFlushNum; fls++) {
        std::unique_ptr<dataholder> record_ptr = std::make_unique<dataholder>();
        int loopStart =loopEq+ fls * this->sweepNumInOneFlush * this->L;
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
            double DeltaE = (EAvgNext - EAvgCurr) / static_cast<double>(this->M);
            //decide if flip is accepted
            if (DeltaE <= 0) {
                sCurr = std::vector<double>(sNext);
                tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                EAvgCurr = EAvgNext;
                muCurr = muNext;
                flipNum++;

            } else {
                double r = distUnif01(e2);
                if (r <= std::exp(-this->beta * DeltaE)) {
                    sCurr = std::vector<double>(sNext);
                    tripleCurr = std::vector<std::tuple<int, eigVal20, vecVal20>>(tripleNext);
                    EAvgCurr = EAvgNext;
                    muCurr = muNext;
                    flipNum++;

                } else {
                    noFlipNum++;
                }


            }
//            std::cout<<"========================"<<std::endl;
//            std::cout<<"step "<<counter<<std::endl;
//            std::cout<<"EAvg Curr="<<EAvgCurr<<std::endl;
//            std::cout<<"mu Curr="<<muCurr<<std::endl;
//            std::cout<<"sCurr=";
//            printVec(sCurr);
//            std::cout<<"============"<<std::endl;


            record_ptr->sAll.push_back(sCurr);
            record_ptr->EAll.push_back(EAvgCurr);
            record_ptr->muAll.push_back(muCurr);
//            record_ptr->eigRstAll.push_back(tripleCurr);
            counter += 1;


        }//end of sweeps in 1 flush


        int loopEnd = loopStart + this->sweepNumInOneFlush * this->L - 1;
//        record_ptr->flattenEigData();
        std::string filenameMiddle = "loopStart" + std::to_string(loopStart) +
                                     "loopEnd" + std::to_string(loopEnd) + "T" + std::to_string(this->T)  + "AfterEq";

        std::string outEFileTmp = outEAllSubDir + filenameMiddle + ".EAll.xml";

        record_ptr->saveVecToXML(outEFileTmp, record_ptr->EAll);

        std::string outMuFileTmp = outMuAllSubDir + filenameMiddle + ".muAll.xml";

        record_ptr->saveVecToXML(outMuFileTmp, record_ptr->muAll);

        std::string outSFileTmp = outSAllSubDir + filenameMiddle + ".sAll.xml";

        record_ptr->saveVecVecToXML(outSFileTmp, record_ptr->sAll);


//        std::string outEigFileTmp = outEigAllSubDir + filenameMiddle + ".eigAll.bin";
//
//        this->serializationViaFStream(record_ptr->flattenedEigSolution, outEigFileTmp);

//        record_ptr->saveEigToXML(outEigFileTmp);

        const auto tflushEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{tflushEnd - tMCStart};
        std::cout << "flush " << fls << std::endl;
        std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;


    }//end of flush loop
    std::ofstream outSummary(outDir + "summaryAfterEq.txt");
    int loopTotal = flipNum + noFlipNum;


    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    outSummary << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    outSummary << "total sweep number: " << static_cast<int>(loopTotal / this->L) << std::endl;
    outSummary << "total loop number: " << loopTotal << std::endl;

    outSummary << "flip number: " << flipNum << std::endl;
    outSummary << "no flip number: " << noFlipNum << std::endl;


    outSummary.close();


}//end of function executionMC()
















std::string dbExchangeModel::execPython(const char *cmd) {
    std::array<char, 4096> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output



}

///
/// @return flattened value, Eigen datatypes to std for serialization
//void dataholder::flattenEigData() {
//
//
//    for (auto const &vecAllForOneS: this->eigRstAll) {
//        std::vector<std::tuple<int, std::vector<double>, std::vector<std::complex<double>> >> eigFor1s;
//
//        for (auto const &tp: vecAllForOneS) {
//            int j = std::get<0>(tp);
//            eigVal20 eigValsTmp = std::get<1>(tp);
//            vecVal20 eigVecsTmp = std::get<2>(tp);//stored in column major format
//
//            std::vector<double> stdEigValsTmp;
//            std::vector<std::complex<double>> stdEigVecsTmp;
//
//            for (auto const &x: eigValsTmp) {
//                stdEigValsTmp.emplace_back(x);
//            }
//            for (auto const &x: eigVecsTmp.reshaped()) {
//                stdEigVecsTmp.emplace_back(x);
//            }
//            eigFor1s.push_back(std::make_tuple(j, stdEigValsTmp, stdEigVecsTmp));
//
////            this->multipleSolutions.push_back(oneEigSolution(j, stdEigValsTmp, stdEigVecsTmp));
////              this->flattenedEigSolution.push_back(std::make_tuple(j,))
//
//        }
//
//        this->flattenedEigSolution.push_back(eigFor1s);
//
//    }
//
//
//}

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

//void
//dbExchangeModel::serializationViaFStream(const std::vector<std::vector<double>> &vecvec, const std::string &fileName) {
//    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
//    msgpack::pack(outFTmp, vecvec);
//    outFTmp.seekp(0);
//    outFTmp.close();
//
//
//}


//void dbExchangeModel::serializationViaFStream(const std::vector<double> &vec, const std::string &fileName) {
//    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
//    msgpack::pack(outFTmp, vec);
//    outFTmp.seekp(0);
//    outFTmp.close();
//
//
//}

void dbExchangeModel::serializationViaFStream(
        const std::vector<std::vector<std::tuple<int, std::vector<double>, std::vector<std::complex<double>> >>> &vecvectuple,
        const std::string &fileName) {
    std::ofstream outFTmp(fileName, std::ios::out | std::ios::binary);
    msgpack::pack(outFTmp, vecvectuple);
    outFTmp.seekp(0);
    outFTmp.close();


}


//void dbExchangeModel::data2File(const dataholder &record) {
//    namespace fs = std::filesystem;
//
//    //output folder
////    std::string outDir="./part"+std::to_string(this->part);
////    if(! fs::is_directory(outDir)|| !fs::exists(outDir)){
////        fs::create_directory(outDir);
////    }
//
//    std::string outSubDir = "./part" + std::to_string(this->part) + "/";
//
//    if (!fs::is_directory(outSubDir) || !fs::exists(outSubDir)) {
//        fs::create_directory(outSubDir);
//    }
//
//
//    std::string prefix =
//            "T" + std::to_string(this->T) + "t" + std::to_string(this->t) + "J" + std::to_string(this->J) + "g" +
//            std::to_string(this->g) + "part" + std::to_string(this->part);
//    //output sAll
//    std::string outsAllName = outSubDir + prefix + ".sAll";
//    serializationViaFStream(record.sAll, outsAllName);
//
//    //output EAll
//    std::string outEAllName = outSubDir + prefix + ".EAll";
//    serializationViaFStream(record.EAll, outEAllName);
//
//    //output muAll
//    std::string outMuAllName = outSubDir + prefix + ".muAll";
//    serializationViaFStream(record.muAll, outMuAllName);
//
//    //output flattened solution
////    std::string outFlEigSolName=outSubDir+prefix+".flattenedEigSolution";
////    serializationViaFStream(record.flattenedEigSolution,outFlEigSolName);
//
//
//}
//
/////
///// @param filename xml file name of eigen-solutions
//void dataholder::saveEigToXML(const std::string &filename) {
//    std::ofstream ofs(filename);
//    boost::archive::xml_oarchive oa(ofs);
////    oa & BOOST_SERIALIZATION_NVP(this->multipleSolutions);
//      for (size_t j=0;j<this->multipleSolutions.size();j++){
//          auto eig=multipleSolutions[j];
//          oa &BOOST_SERIALIZATION_NVP(eig);
//      }
////    oa.put("<\\boost_serialization>\r\n");
//}
//
//
///
/// @param filename xml file name of vec
///@param vec vector to be saved
void dataholder::saveVecToXML(const std::string &filename, const std::vector<double> &vec) {

    std::ofstream ofs(filename);
    boost::archive::xml_oarchive oa(ofs);
    oa & BOOST_SERIALIZATION_NVP(vec);
//    oa.put("<\\boost_serialization>\r\n");


}


///
/// @param filename  xml file name of vecvec
/// @param vecvec vector<vector> to be saved
void dataholder::saveVecVecToXML(const std::string &filename, const std::vector<std::vector<double>> &vecvec) {


    std::ofstream ofs(filename);
    boost::archive::xml_oarchive oa(ofs);
    oa & BOOST_SERIALIZATION_NVP(vecvec);
//    oa.put("<\\boost_serialization>\r\n");



}



void dbExchangeModel::oneStepEig(const std::vector<double> &sCurr) {

auto tripleCurr=this->s2EigSerial(sCurr);
    std::vector<double> EVec = this->combineFromEig(tripleCurr);
    auto EAndMuCurr = this->avgEnergy(EVec);
    double EAvgCurr = EAndMuCurr[0];
    double muCurr = EAndMuCurr[1];
    std::cout<<"sCurr=";
    printVec(sCurr);

    std::cout<<"EAvgCurr="<<EAvgCurr<<std::endl;
    std::cout<<"muCurr="<<muCurr<<std::endl;


}