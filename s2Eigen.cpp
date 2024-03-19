//
// Created by polya on 2/29/24.
//

#include "s2Eigen.hpp"

///search Ell, muAll, sAll folder's files
void reader::searchFiles() {
    this->EPath = this->TDir + "/EAll/";
    this->muPath = this->TDir + "/muAll/";
    this->sAllPath = this->TDir + "/sAll/";

    for (const auto &entry: fs::directory_iterator(EPath)) {
        this->EFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(muPath)) {
        this->muFilesAll.push_back(entry.path().string());
    }
    for (const auto &entry: fs::directory_iterator(sAllPath)) {
        this->sAllFilesAll.push_back(entry.path().string());
    }
//    std::cout<<"finishing search files"<<std::endl;


}


///
/// @param path the path containing xml files
/// @return sorted xml files by starting loop
std::vector <std::string> reader::sortOneDir(const std::vector <std::string> &allFiles) {
    std::vector<int> startingLoopsAll;
    for (const std::string &name: allFiles) {
        std::regex startPattern("loopStart(\\d+)loopEnd");
        std::smatch matchPattern;
        if (std::regex_search(name, matchPattern, startPattern)) {
            startingLoopsAll.push_back(std::stoi(matchPattern.str(1)));
        }

    }

    std::vector <size_t> inds = this->argsort<int>(startingLoopsAll);;

    std::vector <std::string> sortedFiles;
    for (const auto &i: inds) {
        sortedFiles.push_back(allFiles[i]);
    }

    return sortedFiles;

}


///sort files by starting loop
void reader::sortFiles() {

    this->sortedEFilesAll = this->sortOneDir(this->EFilesAll);
    this->sortedmuFilesAll = this->sortOneDir(this->muFilesAll);
    this->sortedsAllFilesAll = this->sortOneDir(this->sAllFilesAll);


}


///parse Txxx.xxxchi.txt to extract lag, lastElemNum, lastFilesNum
void reader::parseCHiFile() {
    std::string chiFileName;
    //search chi file
    std::regex chiFilePattern("T[+-]?\\d*(\\.\\d+)?chi");
    std::smatch matchChi;
    for (const auto &entry: fs::directory_iterator(this->TDir)) {
        std::string fileName = entry.path().string();
        if (std::regex_search(fileName, matchChi, chiFilePattern)) {
            chiFileName = fileName;
            break;

        }


    }

    std::regex lagPattern("lag=([+-]?\\d+)");
    std::regex lastFilesNumPattern("lastFilesNum=(\\d+)");
    std::regex lastElemNumPattern("lastElemNum=(\\d+)");

    std::smatch matchLag;
    std::smatch matchFileNum;
    std::smatch matchElemNum;

    std::ifstream chiIn(chiFileName);
    for (std::string line; std::getline(chiIn, line);) {

        //extract lag value
        if (std::regex_search(line, matchLag, lagPattern)) {
            this->lag = std::stoi(matchLag.str(1));
//            std::cout<<"lag="<<lag<<std::endl;
        }
        //extract lastFilesNum
        if (std::regex_search(line, matchFileNum, lastFilesNumPattern)) {

            this->lastFilesNum = std::stoi(matchFileNum.str(1));

//            std::cout<<"lastFilesNum="<<lastFilesNum<<std::endl;
        }

        //extract lastElemNum
        if (std::regex_search(line, matchElemNum, lastElemNumPattern)) {
            this->lastElemNum = std::stoi(matchElemNum.str(1));

        }


    }

//    std::cout<<"lag="<<lag<<std::endl;

    if (this->lag == -1) {
        this->ferro = true;
//        std::cout<<"c'est ferro"<<std::endl;

    }


}


/// parse sVec values in sAll directory
void reader::parse_sAllDir() {
//    std::cout<<"enterinng parse s"<<std::endl;
//std::cout<<"T="<<T<<std::endl;
//std::cout<<"ferro="<<ferro<<std::endl;
    if (this->ferro == true) {// ferro case, read last file
        std::vector <std::vector<double>> sAll;
//        std::cout<<"entering ferro"<<std::endl;
//      std::cout<<"file num="<<sortedsAllFilesAll.size()<<std::endl;
//      std::cout<<"file name: "<<this->sortedsAllFilesAll[sortedsAllFilesAll.size()-1]<<std::endl;
        std::ifstream ifs(this->sortedsAllFilesAll[sortedsAllFilesAll.size() - 1]);
        if (!ifs.is_open()) {
            std::cerr << "cannot open" << std::endl;
            return;
        }
        boost::archive::xml_iarchive ia(ifs);

        ia >> BOOST_SERIALIZATION_NVP(sAll);

        for (int i = sAll.size() - lastElemNum; i < sAll.size(); i++) {
            this->sSelected.push_back(sAll[i]);
        }
//        std::cout<<sSelected.size()<<std::endl;





    }// end of ferro case
    else {//paramagnetic case, read last lastFilesNum files
//        std::cout<<"entering para"<<std::endl;
        std::vector <std::string> selectedFiles_sAll;
        std::vector <std::vector<double>> sVecsCombined;
        for (int i = this->sortedsAllFilesAll.size() - lastFilesNum; i < this->sortedsAllFilesAll.size(); i++) {
            selectedFiles_sAll.push_back(this->sortedsAllFilesAll[i]);

        }
//        std::cout<<selectedFiles_sAll[selectedFiles_sAll.size()-2]<<std::endl;
//        std::cout<<"sorted file num="<<sortedsAllFilesAll.size()<<std::endl;
//        std::cout<<"last files num="<<lastFilesNum<<std::endl;
//        printVec(selectedFiles_sAll);
        for (const auto &oneName: selectedFiles_sAll) {//read  xml files in sVec
            std::vector <std::vector<double>> sAll;
            std::ifstream ifs(oneName);
            if (!ifs.is_open()) {
                std::cerr << "cannot open" << std::endl;
                return;
            }
            boost::archive::xml_iarchive ia(ifs);
            ia >> BOOST_SERIALIZATION_NVP(sAll);
            for (const auto &vec: sAll) {
                sVecsCombined.push_back(vec);

            }


        }//end of reading  sVec xml files

        int onePartLengh = static_cast<int>(sVecsCombined.size() / 2);
//        std::cout<<"one part length="<<onePartLengh<<std::endl;
//        std::cout<<"lag="<<lag<<std::endl;
        for (int i = 0; i < onePartLengh; i += lag) {
            this->sSelected.push_back(sVecsCombined[i]);
        }
        int i = onePartLengh;
        for (; i < 2 * onePartLengh; i += lag) {
            this->sSelected.push_back(sVecsCombined[i]);
        }
//            std::cout<<"s last ind="<<i-lag<<std::endl;
//std::cout<<"s num="<<sSelected.size()<<std::endl;

//            int i1=1000;
//            int i2=int(sSelected.size()*0.75);
//            std::cout<<i1<<" th vec = ";
//        printVec(sSelected[i1]);
//        std::cout<<i2<<" th vec = ";
//        printVec(sSelected[i2]);
    }//end of paramagnetic case

//std::cout<<"data point num="<<sSelected.size()<<std::endl;

}

///parse EAvg values in EAll directory
void reader::parse_EAllDir() {
    if (this->ferro == true) {//ferro case, read last file
        std::vector<double> EInLastFile;
        std::ifstream ifs(this->sortedEFilesAll[sortedEFilesAll.size() - 1]);
        if (!ifs.is_open()) {
            std::cerr << "cannot open" << std::endl;
            return;
        }
        boost::archive::xml_iarchive ia(ifs);
        ia >> BOOST_SERIALIZATION_NVP(EInLastFile);
        for (int i = EInLastFile.size() - lastElemNum; i < EInLastFile.size(); i++) {
            this->ESelected.push_back(EInLastFile[i]);

        }
//        printVec(ESelected);

    }// end of ferro case
    else {//paramagnetic case, read last lastFilesNum files
        std::vector <std::string> selectedFilesEAll;
        std::vector<double> EAllCombined;
        for (int i = this->sortedEFilesAll.size() - lastFilesNum; i < this->sortedEFilesAll.size(); i++) {
            selectedFilesEAll.push_back(this->sortedEFilesAll[i]);


        }
//        std::cout<<selectedFilesEAll[selectedFilesEAll.size()-2]<<std::endl;
//        std::cout<<"selected file num="<<selectedFilesEAll.size()<<std::endl;

        for (const auto &oneName: selectedFilesEAll) {//read  xml files in EAll
            std::vector<double> EPerFile;
            std::ifstream ifs(oneName);
            if (!ifs.is_open()) {
                std::cerr << "cannot open" << std::endl;
                return;
            }
            boost::archive::xml_iarchive ia(ifs);
            ia >> BOOST_SERIALIZATION_NVP(EPerFile);
            for (const auto &val: EPerFile) {
                EAllCombined.push_back(val);
            }


        }//end of reading  xml files in EAll

        int onePartLength = static_cast<int>(EAllCombined.size() / 2);
        for (int i = 0; i < onePartLength; i += lag) {
            this->ESelected.push_back(EAllCombined[i]);
        }
        int i = onePartLength;
        for (; i < 2 * onePartLength; i += lag) {
            this->ESelected.push_back(EAllCombined[i]);
        }
//        std::cout << "E last ind=" << i - lag << std::endl;

    }//end of paramagnetic case

//    int i1=1000;
//int i2=int(ESelected.size()*0.75);
//std::cout<<i1<<" th value="<<ESelected[i1]<<std::endl;
//std::cout<<i2<<" th value="<<ESelected[i2]<<std::endl;
//std::cout<<"E num="<<ESelected.size()<<std::endl;

}


///parse mu values in muAll directory
void reader::parse_muAllDir() {
    if (this->ferro == true) {//ferro case, read last file
        std::vector<double> muInLastFile;
        std::ifstream ifs(this->sortedmuFilesAll[sortedmuFilesAll.size() - 1]);
        if (!ifs.is_open()) {
            std::cerr << "cannot open" << std::endl;
            return;
        }
        boost::archive::xml_iarchive ia(ifs);
        ia >> BOOST_SERIALIZATION_NVP(muInLastFile);
        for (int i = muInLastFile.size() - lastElemNum; i < muInLastFile.size(); i++) {
            this->muSelected.push_back(muInLastFile[i]);
        }

    }// end of ferro case
    else {//paramagnetic case, read last lastFilesNum files
        std::vector <std::string> selectedMuFilesAll;
        std::vector<double> muAllCombined;
        for (int i = this->sortedmuFilesAll.size() - lastFilesNum; i < this->sortedmuFilesAll.size(); i++) {
            selectedMuFilesAll.push_back(this->sortedmuFilesAll[i]);
        }
//        std::cout<<selectedMuFilesAll[selectedMuFilesAll.size()-2]<<std::endl;

        for (const auto &oneName: selectedMuFilesAll) {//read xml files in muAll
            std::vector<double> muPerFile;
            std::ifstream ifs(oneName);
            if (!ifs.is_open()) {
                std::cerr << "cannot open" << std::endl;
                return;
            }
            boost::archive::xml_iarchive ia(ifs);
            ia >> BOOST_SERIALIZATION_NVP(muPerFile);
            for (const auto &val: muPerFile) {
                muAllCombined.push_back(val);
            }

        }//end of reading xml files in muAll
        int onePartLength = static_cast<int>(muAllCombined.size() / 2);
        for (int i = 0; i < onePartLength; i += lag) {
            this->muSelected.push_back(muAllCombined[i]);
        }
        int i = onePartLength;
        for (; i < 2 * onePartLength; i += lag) {
            this->muSelected.push_back(muAllCombined[i]);
        }
//        std::cout<<"mu last ind="<<i-lag<<std::endl;


    }//end of paramagnetic case
//    std::cout<<"mu num="<<this->muSelected.size()<<std::endl;

}


///compute Ze weights and projections
void reader::fillZeWeights(dbExchangeModel &model) {

//check consistency from reading and new computation
//std::cout<<"E num="<<this->ESelected.size()<<std::endl;
//std::cout<<"s num="<<this->sSelected.size()<<std::endl;
//std::cout<<"mu num="<<this->muSelected.size()<<std::endl;
//        int i=int(ESelected.size()*0.75);
//        double muTmp=this->muSelected[i];
//        double ETmp=this->ESelected[i];
//        std::vector<double> one_s=this->sSelected[i];
//
//        auto tripleTmp=model.s2EigSerial(one_s);
//        std::vector<double> EVecTmp=model.combineFromEig(tripleTmp);
//        auto EAndMuTmp=model.avgEnergy(EVecTmp);
//        double EAvgTmp=EAndMuTmp[0];
//        double mu2Tmp=EAndMuTmp[1];
//        std::cout<<"T="<<model.T<<std::endl;
//        std::cout<<"T reader="<<this->T<<std::endl;
//
//        std::cout<<"s from reading:";
//    printVec(one_s);
//
//        std::cout<<"E from reading="<<ETmp<<std::endl;
//        std::cout<<"E from computing="<<EAvgTmp<<std::endl;
//
//        std::cout<<"mu from reading="<<muTmp<<std::endl;
//        std::cout<<"mu from computing="<<mu2Tmp<<std::endl;
    int P = this->sSelected.size();
//    std::cout<<"P="<<std::endl;
    for (int p = 0; p < P; p++) {
        auto s = this->sSelected[p];
        auto tripleTmp = model.s2EigSerial(s);
        std::vector<double> EInOne_s = model.combineFromEig(tripleTmp);
        auto EAnd_muTmp = model.avgEnergy(EInOne_s);
        //compute one A mat
//        std::cout<<"len="<<tripleTmp.size()<<std::endl;
        //projection

        for (int m = 0; m < M; m++) {
            for (int a = 0; a < L; a++) {
                std::vector <std::vector<double>> oneATmp;
                int rowCount = 0;
                for (int j = 0; j < 2; j++) {
                    for (int b = 0; b < 2 * L; b++) {

                        vec20c yEigen(this->yVecsAll[m][a][j].data());
                        double ETmp = std::get<1>(tripleTmp[m])[b];
                        vec20c zb = std::get<2>(tripleTmp[m]).col(b);
                        double cbja = std::abs(zb.dot(yEigen));
//                        this->AMatsAll[p][m][a][rowCount]=std::vector<double>({ETmp,cbja});
//                        std::cout<<"ETmp="<<ETmp<<", cbja="<<cbja<<std::endl;
                        oneATmp.push_back({ETmp, cbja});
                        rowCount++;
                    }
                }

                int sortBy = 0;
                std::sort(oneATmp.begin(), oneATmp.end(),
                          [sortBy](const std::vector<double> &a, const std::vector<double> &b) {
                              return a[sortBy] < b[sortBy];
                          });

                this->AMatsAll[p][m][a] = oneATmp;
//                std::cout<<"oneATmp:"<<oneATmp[0][1]<<std::endl;
            }//a end
        }//m end





        double EAvgTmp = EAnd_muTmp[0];
        double muTmp = EAnd_muTmp[1];
//        std::cout<<"EAvgTmp="<<EAvgTmp<<", muTmp="<<muTmp<<std::endl;
//        if(ferro==true){
//            std::cout<<"E="<<EAvgTmp<<std::endl;
//            std::cout<<"mu="<<muTmp<<std::endl;
//        }
        this->muRecomputed.push_back(muTmp);
        this->epsilonSelected.push_back(EAvgTmp / static_cast<double >(M));
        std::vector <std::vector<double>> weightForOne_s;
        std::vector <std::vector<double>> eigsForOne_s;
        double sumForOne_s = 0;
        for (const auto &tp: tripleTmp) {//for all K
            std::vector<double> weightForOne_sOneK;
            std::vector<double> eigsForOne_sOneK;
//            int j=std::get<0>(tp);
//            double KTmp=model.KSupValsAll[j];
            eigVal20 eigValsTmp = std::get<1>(tp);
            for (const double &e: eigValsTmp) {//for all eigenvalues in one K
                double expVal = std::exp(beta * (e - muTmp));
                double wtTmp = std::pow(1 / (expVal + 1), 2) * expVal;
//                if (ferro==true){
//
//                    std::cout<<"wtTmp="<<wtTmp<<std::endl;
//                }
                weightForOne_sOneK.push_back(wtTmp);
                eigsForOne_sOneK.push_back(e);
            }
            sumForOne_s += std::accumulate(weightForOne_sOneK.begin(), weightForOne_sOneK.end(), 0.0);
            weightForOne_s.push_back(weightForOne_sOneK);
            eigsForOne_s.push_back(eigsForOne_sOneK);

        }//end of for all K
        this->ZeVals.push_back(sumForOne_s);
//        std::cout<<"sumForOne_s"<<sumForOne_s<<std::endl;

        for (auto &vec: weightForOne_s) {
            for (auto &val: vec) {
                val /= sumForOne_s;
            }
        }
        this->eigValsRecomputed.push_back(eigsForOne_s);
        this->ZeWeightsAll.push_back(weightForOne_s);
//        if (ferro==true){
//            printVec(ZeWeightsAll[0][0]);
//        }

//        double tmp=0;
//        for(const auto &vec:weightForOne_s){
//            tmp+=std::accumulate(vec.begin(),vec.end(),0.0);
//        }
//        std::cout<<"total weight: "<<tmp<<std::endl;





    }//end of for all s (p)


}//end of fillZeWeights


///parse last file under EAll,muAll, sAll
void reader::parseLast(dbExchangeModel &model) {

    int ind = 0;
    std::string lastEFile = "/home/polya/Documents/cppCode/paramagnet_interaction_cpp/part1/T3.100000/EAll/loopStart660000loopEnd689999T3.100000t0.400000J-1.000000g0.010000part1L10M20AfterEq.EAll.xml";

    std::vector<double> EInLastFile;
    std::ifstream ifsE(lastEFile);
    boost::archive::xml_iarchive iaE(ifsE);
    iaE >> BOOST_SERIALIZATION_NVP(EInLastFile);
    ind = EInLastFile.size() - 2;
//    std::cout << "E=" << EInLastFile[ind] << std::endl;

    std::string lastMuFile = "/home/polya/Documents/cppCode/paramagnet_interaction_cpp/part1/T3.100000/muAll/loopStart660000loopEnd689999T3.100000t0.400000J-1.000000g0.010000part1L10M20AfterEq.muAll.xml";
    std::vector<double> muInLastFile;
    std::ifstream ifsMu(lastMuFile);
    boost::archive::xml_iarchive iaMu(ifsMu);
    iaMu >> BOOST_SERIALIZATION_NVP(muInLastFile);
//    std::cout << "mu=" << muInLastFile[ind] << std::endl;


    std::string lastSFile = "/home/polya/Documents/cppCode/paramagnet_interaction_cpp/part1/T3.100000/sAll/loopStart660000loopEnd689999T3.100000t0.400000J-1.000000g0.010000part1L10M20AfterEq.sAll.xml";
    std::vector <std::vector<double>> sInLastFile;
    std::ifstream ifsS(lastSFile);
    boost::archive::xml_iarchive iaS(ifsS);
    iaS >> BOOST_SERIALIZATION_NVP(sInLastFile);
//    std::cout << "s=";
//    printVec(sInLastFile[ind]);
    std::vector<double> lastS = sInLastFile[ind];
    auto tripleTmp = model.s2EigSerial(lastS);
    auto EVecLast = model.combineFromEig(tripleTmp);
    auto EAndMuLast = model.avgEnergy(EVecLast);
    double EAvgLast = EAndMuLast[0];
//    double muLast = EAndMuLast[1];
//    std::cout << "computation E=" << EAvgLast << std::endl;
//    std::cout << "computation mu=" << muLast << std::endl;

//    std::cout << "model t=" << model.t << ", J=" << model.J << ", g=" << model.g << ", T=" << model.T << std::endl;

}


/// @param s_ind index of s vector
/// @return WE for each s
double reader::computeWE(const int &s_ind) {

    std::vector <std::vector<double>> &weightsTmp = this->ZeWeightsAll[s_ind];
    std::vector <std::vector<double>> &eigsTmp = this->eigValsRecomputed[s_ind];
    double sumTmp = 0;
    for (int i = 0; i < weightsTmp.size(); i++) {
        auto &vec_wt = weightsTmp[i];
        auto &vec_eigs = eigsTmp[i];
        sumTmp += std::inner_product(vec_wt.begin(), vec_wt.end(), vec_eigs.begin(), 0.0);
    }

    return sumTmp;

}


/// @param s_ind index of s vector
/// @return WE2 for each s
double reader::computeWE2(const int &s_ind) {
    std::vector <std::vector<double>> &weightsTmp = this->ZeWeightsAll[s_ind];
    std::vector <std::vector<double>> &eigsTmp = this->eigValsRecomputed[s_ind];
    double sumTmp = 0;
    for (int i = 0; i < weightsTmp.size(); i++) {
        auto &vec_wt = weightsTmp[i];
        auto &vec_eigs = eigsTmp[i];
        std::vector<double> vec_eigs2;
        for (const auto &e: vec_eigs) {
            vec_eigs2.push_back(std::pow(e, 2));
        }
        sumTmp += std::inner_product(vec_wt.begin(), vec_wt.end(), vec_eigs2.begin(), 0.0);


    }
    return sumTmp;

}

///compute all WE and WE2
void reader::computeAllWEAllWE2() {
    for (int i = 0; i < this->sSelected.size(); i++) {
        double tmp = this->computeWE(i);
        double tmp2 = this->computeWE2(i);
//        std::cout<<"tmp="<<tmp<<std::endl;
//        std::cout<<"tmp2="<<tmp2<<std::endl;
        this->WE.push_back(tmp);
        this->WE2.push_back(tmp2);
    }


}


///
/// @param i index of s, to be deleted in computation
/// @return pseudovalue of \partial_{\beta}\braket{\epsilon} with si deleted
double reader::dbeta_epsilon(const int &i) {
    std::vector<double> epsilonDeleted;
    std::vector<double> ZeDeleted;
    std::vector<double> WEDeleted;
    std::vector<double> WE2Deleted;
//    if(T==1.066667){
//        std::cout<<this->epsilonSelected.size()<<std::endl;
//std::cout<<this->ZeVals.size()<<std::endl;
//std::cout<<this->WE.size()<<std::endl;
//std::cout<<this->WE2.size()<<std::endl;
//    }
//std::cout<<this->epsilonSelected.size()<<std::endl;
//std::cout<<this->ZeVals.size()<<std::endl;
//std::cout<<this->WE.size()<<std::endl;
//std::cout<<this->WE2.size()<<std::endl;

    for (int j = 0; j < this->epsilonSelected.size(); j++) {
        if (j == i) {
            continue;
        } else {
            epsilonDeleted.push_back(this->epsilonSelected[j]);
            ZeDeleted.push_back(this->ZeVals[j]);
            WEDeleted.push_back(this->WE[j]);
            WE2Deleted.push_back(this->WE2[j]);
        }
    }//end for

    //compute part1
    std::vector<double> oneMinus_betaEpsilonVec;
    for (const auto &val: epsilonDeleted) {
        oneMinus_betaEpsilonVec.push_back(1 - this->beta * val);
    }

    std::vector<double> W2EMinusWE2Vec;
    for (int l = 0; l < WEDeleted.size(); l++) {
        double tmp = std::pow(WEDeleted[l], 2) - WE2Deleted[l];
        W2EMinusWE2Vec.push_back(tmp);

    }
    std::vector<double> part1_vec;
    for (int l = 0; l < ZeDeleted.size(); l++) {
        double tmp = oneMinus_betaEpsilonVec[l] * ZeDeleted[l] * W2EMinusWE2Vec[l] / static_cast<double>(M);
        part1_vec.push_back(tmp);
    }
    double part1 = std::accumulate(part1_vec.begin(), part1_vec.end(), 0.0);
    part1 /= static_cast<double >(part1_vec.size());



//compute part2
    double part2 = std::accumulate(epsilonDeleted.begin(), epsilonDeleted.end(), 0.0);
    part2 /= static_cast<double >(epsilonDeleted.size());


//compute part3
    std::vector<double> part3_vec;
    for (int l = 0; l < ZeDeleted.size(); l++) {
        double tmp = ZeDeleted[l] * W2EMinusWE2Vec[l];
        part3_vec.push_back(tmp);
    }

    double part3 = std::accumulate(part3_vec.begin(), part3_vec.end(), 0.0);
    part3 /= static_cast<double>(part3_vec.size());
    part3 *= beta;
    part3 /= static_cast<double>(M);

//compute part4

    std::vector<double> part4_vec;
    for (const auto &val: epsilonDeleted) {
        part4_vec.push_back(std::pow(val, 2));
    }

    double part4 = std::accumulate(part4_vec.begin(), part4_vec.end(), 0.0);
    part4 /= static_cast<double>(part4_vec.size());
//
//    if(T==1.066667){
//        std::cout<<"part1="<<part1<<", part2="<<part2<<", part3="<<part3<<", part4="<<part4<<std::endl;
//    }
    double rst = part1 + part2 * part3 + std::pow(part2, 2) - part4;
    return rst;


}// end dbeta_epsilon



///
/// @param ps pseudovalue of dbeta_epsilon
/// @param sd standard deviation of dbeta_epsilon
void reader::pseudoValueOfC(double &ps, double &sd) {

    std::vector<double> C_All;
    for (int i = 0; i < this->sSelected.size(); i++) {
        C_All.push_back(-dbeta_epsilon(i) / std::pow(this->T, 2));

    }
//    std::cout<<"sSelected len="<<sSelected.size()<<std::endl;
//    std::cout<<"C_All len="<<C_All.size()<<std::endl;
//    std::cout<<C_All[10]<<std::endl;
//    if(T==1.066667){
//        printVec(C_All);
//    }
    ps = std::accumulate(C_All.begin(), C_All.end(), 0.0);
    ps /= static_cast<double >(C_All.size());

    std::vector<double> varVec;
    for (int i = 0; i < C_All.size(); i++) {
        double tmp = std::pow(ps - C_All[i], 2) / static_cast<double >(C_All.size() - 1);
        varVec.push_back(tmp);
    }
//    if(T==1.066667){
//        printVec(varVec);
//    }

    double var = std::accumulate(varVec.begin(), varVec.end(), 0.0);

    sd = std::sqrt(var / static_cast<double >(C_All.size()));
//    std::cout<<"C="<<ps<<std::endl;
//    std::cout<<"sd="<<sd<<std::endl;






}


///compute C and write to file
void reader::C2File() {
    double C_ps = 0, sd = 0;
    this->pseudoValueOfC(C_ps, sd);
//    if(T==1.066667){
//        std::cout<<"C="<<C_ps<<std::endl;
//        std::cout<<"sd="<<sd<<std::endl;
//    }
//std::cout<<"==================="<<std::endl;
//    std::cout<<"T="<<T<<std::endl;
//    std::cout<<"C="<<C_ps<<std::endl;
//    std::cout<<"sd="<<sd<<std::endl;
//    std::cout<<"==================================="<<std::endl;
    std::string outCDir = "./group" + std::to_string(this->groupNum)+"data/row"+std::to_string(rowNum) +"/"+ "CAll/";
    namespace fs = boost::filesystem;
    if (!fs::is_directory(outCDir) || !fs::exists(outCDir)) {
        fs::create_directories(outCDir);
    }
    std::string outCFileName = outCDir + "T" + std::to_string(this->T) + "C.txt";
    std::ofstream ofs(outCFileName);
    ofs << "C=" << C_ps << std::endl;
    ofs << "halfLength=" << 1.960 * sd << std::endl;
    ofs.close();


}
/////
///// @tparam T
///// @param vec
///// @return norm of vector
//template<class valType>
//double reader::vectorNorm(const std::vector<valType>&vec){
//
//    double sum = 0.0;
//    for (const auto& val : vec) {
//        sum += std::pow(std::abs(val), 2);
//    }
//    return std::sqrt(sum);
//
//}


/// construct all y vectors
/// @param model dbExchangeModel object
void reader::construct_yAll(const dbExchangeModel &model) {
    //y vectors are indexed by m,a,j
    for (int m = 0; m < M; m++) {
        double Km = model.KSupValsAll[m];
        std::vector < std::vector < std::vector < std::complex < double >> >> yVecForOne_m;

        for (int a = 0; a < L; a++) {
            std::vector < std::vector < std::complex < double>>> yVecForOne_a;
            for (int j = 0; j < xVecsAll.size(); j++) {
                std::vector <std::complex<double>> one_yVec;
                for (int l = 0; l < L; l++) {
                    one_yVec.emplace_back(xVecsAll[j][0]);
                    one_yVec.emplace_back(xVecsAll[j][1]);
                }
                for (int l = 0; l < L; l++) {
                    std::complex<double> coefTmp = std::exp(1i * static_cast<double>(l) * Km) * std::exp(
                            1i * 2.0 * PI * static_cast<double>(l) / static_cast<double>(L) * static_cast<double >(a));
                    one_yVec[2 * l] *= coefTmp;
                    one_yVec[2 * l + 1] *= coefTmp;

                }//finishing constructing one y
                double nm = vectorNorm(one_yVec);//normaliztion
                for (auto &val: one_yVec) {
                    val /= nm;
                }
                yVecForOne_a.push_back(one_yVec);

            }//end of j loop
            yVecForOne_m.push_back(yVecForOne_a);
        }//end of a loop
        this->yVecsAll.push_back(yVecForOne_m);
    }//end of m loop

//    printVec(yVecsAll[2][2][1]);


}


///initialize all A matrices, AMatsAll indexed by s,m,a,j
void reader::initAMatsAll() {
    std::vector <std::vector<double>> oneAMat;
    oneAMat.reserve(2 * L);
    for (int b = 0; b < 2 * L; b++) {
        for (int j = 0; j < 2; j++) {
            oneAMat.push_back({0.0, 0.0});

        }
    }

    //a
    std::vector < std::vector < std::vector < double>>> AMat_a;
    AMat_a.reserve(L);
    for (int a = 0; a < L; a++) {
        AMat_a.push_back(oneAMat);
    }

    //m
    std::vector < std::vector < std::vector < std::vector < double >> >> AMat_m;
    AMat_m.reserve(M);
    for (int m = 0; m < M; m++) {
        AMat_m.push_back(AMat_a);
    }

    int P = this->sSelected.size();
    //s
    this->AMatsAll.reserve(P);
    for (int p = 0; p < P; p++) {
        this->AMatsAll.push_back(AMat_m);
    }
//std::cout<<"A size"<<sizeof(AMatsAll)<<std::endl;

}


///compute the average and error of eigenvalue
void reader::computeMeanE() {
    for (int m = 0; m < M; m++) {
        std::vector <std::vector<double>> meanFor_one_m;
        std::vector <std::vector<double>> hfLengthForOne_m;
        for (int a = 0; a < L; a++) {

            std::vector <std::vector<double>> EForAll_s_One_m_One_a;


//            std::cout<<"P="<<AMatsAll.size()<<std::endl;
            for (const auto &all_m_a: this->AMatsAll) {//for each s
                const auto &ATmp = all_m_a[m][a];
                std::vector<double> EVecTmp;
                for (const auto &row: ATmp) {
                    EVecTmp.push_back(row[0]);//elements in EVecTmp are in ascending order
//                    std::cout<<row[0]<<std::endl;
                }
                EForAll_s_One_m_One_a.push_back(EVecTmp);
            }//end for each s

            //collapse EForAll_s_One_m_One_a by row sum
            std::vector<double> rowMean(EForAll_s_One_m_One_a[0].size(), 0.0);

            for (const auto &vec: EForAll_s_One_m_One_a) {
                addRightVecToLeft(rowMean, vec);
            }

            //mean
            int P = EForAll_s_One_m_One_a.size();
            for (auto &val: rowMean) {
                val /= static_cast<double >(P);
            }

            meanFor_one_m.push_back(rowMean);

            //half length
            std::vector <std::vector<double>> hfLengthForAll_s_One_m_One_a;
            for (const auto &all_m_a: this->AMatsAll) {//for each s
                const auto &ATmp = all_m_a[m][a];
                std::vector<double> EVecTmp;
                for (const auto &row: ATmp) {
                    EVecTmp.push_back(row[0]);//elements in EVecTmp are in ascending order
//
                }
                hfLengthForAll_s_One_m_One_a.push_back(EVecTmp);
            }// end for each s

            //substract rowMean from each row of hfLengthForAll_s_One_m_One_a, take square
            std::vector <std::vector<double>> diffSquared;
            for (const auto &vec: hfLengthForAll_s_One_m_One_a) {
                auto rowTmp = leftMinusRightVec(vec, rowMean);
                auto rowTmp2 = vectorSquared(rowTmp);
                diffSquared.push_back(rowTmp2);

            }

            std::vector<double> var(diffSquared[0].size(), 0.0);
            for (const auto &vec: diffSquared) {
                addRightVecToLeft(var, vec);
            }

            for (auto &val: var) {
                val /= static_cast<double >(P - 1);
            }

            std::vector<double> hfLengthForOne_a;
            for (const auto &val: var) {
                hfLengthForOne_a.push_back(1.96 * std::sqrt(val / static_cast<double >(P)));
            }
            hfLengthForOne_m.push_back(hfLengthForOne_a);

        }//a end
        this->meanE.push_back(meanFor_one_m);
        this->hfLengthE.push_back(hfLengthForOne_m);

    }// m end



}


/// add the value of RHS to LHS
/// @param LHS vector
/// @param RHS vector
void reader::addRightVecToLeft(std::vector<double> &LHS, const std::vector<double> &RHS) {
    int length = RHS.size();
    //LHS and RHS must have the same size

    for (int i = 0; i < length; i++) {
        LHS[i] += RHS[i];
    }


}


///
/// @param LHS vector
/// @param RHS vector
/// @return LHS-RHS
std::vector<double> reader::leftMinusRightVec(const std::vector<double> &LHS, const std::vector<double> &RHS) {
    //LHS size must == RHS size
    int length = LHS.size();
    std::vector<double> rst(length, 0.0);
    for (int i = 0; i < length; i++) {
        rst[i] = LHS[i] - RHS[i];
    }

    return rst;


}


///
/// @param vec vector
/// @return elementwise squared vector
std::vector<double> reader::vectorSquared(const std::vector<double> &vec) {

    int length = vec.size();
    std::vector<double> retVec(length, 0.0);
    for (int i = 0; i < length; i++) {
        retVec[i] = std::pow(vec[i], 2);
    }

    return retVec;
}


///compute the marker size of unfolded band
void reader::computeMarkerSize() {
    for (int m = 0; m < M; m++) {
        std::vector <std::vector<double>> sizeForOne_m;
        for (int a = 0; a < L; a++) {
            std::vector <std::vector<double>> sizeForAll_s_One_m_One_a;
            int P = AMatsAll.size();
            for (const auto &all_m_a: this->AMatsAll) {//for each s
                const auto &ATmp = all_m_a[m][a];
                std::vector<double> sizeVecTmp;
                for (const auto &row: ATmp) {
                    sizeVecTmp.push_back(std::pow(row[1], 2));
                }

                sizeForAll_s_One_m_One_a.push_back(sizeVecTmp);

            }//end for each s

            //mean of marker size
            std::vector<double> rowMean(sizeForAll_s_One_m_One_a[0].size(), 0.0);
            for (const auto &vec: sizeForAll_s_One_m_One_a) {

                addRightVecToLeft(rowMean, vec);
            }
            for (auto &val: rowMean) {
                val /= static_cast<double >(P);
            }

            sizeForOne_m.push_back(rowMean);
        }//end a
        this->markerSize.push_back(sizeForOne_m);
    }//end m


}


///
/// @return flattened data
std::vector<double> reader::flattenData(const std::vector<std::vector<std::vector<double>>

> &data) {
    std::vector<double> outVec;
    for (
        const auto &vecvec
            : data) {
        for (
            const auto &vec
                : vecvec) {
            for (
                const auto &val
                    : vec) {
                outVec.push_back(val);
            }
        }
    }
    return outVec;
}


///write meanE, hfLengthE, markerSize to file
void reader::bandToFile() {
    std::string outDir = this->TDir + "/";
    std::string suffix = std::to_string(M) + "L" + std::to_string(L) + "vecLength" + std::to_string(4 * L) + ".xml";
    std::string EMeanFile = outDir + "EMeanM" + suffix;
    std::string EHfLength = outDir + "EHfLengthM" + suffix;
    std::string mkSize = outDir + "markerSizeM" + suffix;

    auto flatttenedEMean = flattenData(meanE);
    std::ofstream ofs0(EMeanFile);
    boost::archive::xml_oarchive oa0(ofs0);
    oa0 & BOOST_SERIALIZATION_NVP(flatttenedEMean);


    auto flattenedEHfLength = flattenData(hfLengthE);
    std::ofstream ofs1(EHfLength);
    boost::archive::xml_oarchive oa1(ofs1);
    oa1 & BOOST_SERIALIZATION_NVP(flattenedEHfLength);


    auto flattenedSize = flattenData(markerSize);
    std::ofstream ofs2(mkSize);
    boost::archive::xml_oarchive oa2(ofs2);
    oa2 & BOOST_SERIALIZATION_NVP(flattenedSize);


}
