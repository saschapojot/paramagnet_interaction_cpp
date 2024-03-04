//
// Created by polya on 2/29/24.
//

#include "s2Eigen.hpp"

///search Ell, muAll, sAll folder's files
void reader::searchFiles() {
    this-> EPath = this->TDir + "/EAll/";
    this-> muPath = this->TDir + "/muAll/";
    this-> sAllPath = this->TDir + "/sAll/";

    for (const auto &entry: fs::directory_iterator(EPath)) {
        this->EFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(muPath)) {
        this->muFilesAll.push_back(entry.path().string());
    }
    for (const auto &entry: fs::directory_iterator(sAllPath)) {
        this->sAllFilesAll.push_back(entry.path().string());
    }


}


///
/// @param path the path containing xml files
/// @return sorted xml files by starting loop
std::vector<std::string> reader::sortOneDir(const std::vector<std::string> & allFiles) {
    std::vector<int> startingLoopsAll;
    for (const std::string &name: allFiles) {
        std::regex startPattern("loopStart(\\d+)loopEnd");
        std::smatch matchPattern;
        if (std::regex_search(name, matchPattern, startPattern)) {
            startingLoopsAll.push_back(std::stoi(matchPattern.str(1)));
        }

    }

    std::vector<size_t> inds=this->argsort<int>(startingLoopsAll);;

    std::vector<std::string> sortedFiles;
    for(const auto &i: inds){
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
        }
        //extract lastFilesNum
        if (std::regex_search(line, matchFileNum, lastFilesNumPattern)) {
            this->lastFilesNum =1;// std::stoi(matchFileNum.str(1));
        }

        //extract lastElemNum
        if (std::regex_search(line, matchElemNum, lastElemNumPattern)) {
            this->lastElemNum = std::stoi(matchElemNum.str(1));
        }


    }

//    std::cout<<"lag="<<lag<<std::endl;
//    std::cout<<"lastFilesNum="<<lastFilesNum<<std::endl;
//    std::cout<<"lastElemNum="<<lastElemNum<<std::endl;
    if (this->lag == -1) {
        this->ferro = true;
//        std::cout<<"ferro: "<<ferro<<std::endl;
    }


}


/// parse sVec values in sAll directory
void reader::parse_sAllDir() {
//    std::cout<<"ferro="<<ferro<<std::endl;
    if (this->ferro == true) {// ferro case, read last file
        std::vector<std::vector<double>> sAll;
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
        std::vector<std::string> selectedFiles_sAll;
        std::vector<std::vector<double>> sVecsCombined;
        for (int i = this->sortedsAllFilesAll.size() - lastFilesNum; i < this->sortedsAllFilesAll.size(); i++) {
            selectedFiles_sAll.push_back(this->sortedsAllFilesAll[i]);

        }
//        std::cout<<selectedFiles_sAll[selectedFiles_sAll.size()-2]<<std::endl;
//        std::cout<<"sorted file num="<<sortedsAllFilesAll.size()<<std::endl;
//        std::cout<<"last files num="<<lastFilesNum<<std::endl;
        for (const auto &oneName: selectedFiles_sAll) {//read  xml files in sVec
            std::vector<std::vector<double>> sAll;
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

        int onePartLengh=static_cast<int>(sVecsCombined.size()/2);
//        std::cout<<"one part length="<<onePartLengh<<std::endl;
//        std::cout<<"lag="<<lag<<std::endl;
            for(int i=0;i<onePartLengh;i+=lag){
                this->sSelected.push_back(sVecsCombined[i]);
            }
            int i=onePartLengh;
            for(;i<2*onePartLengh;i+=lag){
                this->sSelected.push_back(sVecsCombined[i]);
            }
            std::cout<<"s last ind="<<i-lag<<std::endl;

//            int i1=1000;
//            int i2=int(sSelected.size()*0.75);
//            std::cout<<i1<<" th vec = ";
//        printVec(sSelected[i1]);
//        std::cout<<i2<<" th vec = ";
//        printVec(sSelected[i2]);
    }//end of paramagnetic case



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
        std::vector<std::string> selectedFilesEAll;
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
        int i=onePartLength;
        for (; i < 2 * onePartLength; i += lag) {
            this->ESelected.push_back(EAllCombined[i]);
        }std::cout<<"E last ind="<<i-lag<<std::endl;

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
        std::vector<std::string> selectedMuFilesAll;
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
        int i=onePartLength;
        for (; i < 2 * onePartLength; i += lag) {
            this->muSelected.push_back(muAllCombined[i]);
        }
        std::cout<<"mu last ind="<<i-lag<<std::endl;


    }//end of paramagnetic case
//    std::cout<<"mu num="<<this->muSelected.size()<<std::endl;

}


///compute Ze weights
void reader::fillZeWeights( dbExchangeModel& model){
//    for(const auto &s: this->sSelected){
//        auto tripleTmp=model.s2EigSerial(s);
//        std::vector<double>EInOne_s=model.combineFromEig(tripleTmp);
//        auto EAnd_muTmp=model.avgEnergy(EInOne_s);
//
//
//
//    }
        int i=this->ESelected.size()-1;
        double muTmp=this->muSelected[i];
        double ETmp=this->ESelected[i];
        std::vector<double> one_s=this->sSelected[i];

        auto tripleTmp=model.s2EigSerial(one_s);
        std::vector<double> EVecTmp=model.combineFromEig(tripleTmp);
        auto EAndMuTmp=model.avgEnergy(EVecTmp);
        double EAvgTmp=EAndMuTmp[0];
        double mu2Tmp=EAndMuTmp[1];
        std::cout<<"T="<<model.T<<std::endl;
        std::cout<<"T reader="<<this->T<<std::endl;

        std::cout<<"s from reading:";
    printVec(one_s);

        std::cout<<"E from reading="<<ETmp<<std::endl;
        std::cout<<"E from computing="<<EAvgTmp<<std::endl;

        std::cout<<"mu from reading="<<muTmp<<std::endl;
        std::cout<<"mu from computing="<<mu2Tmp<<std::endl;

}