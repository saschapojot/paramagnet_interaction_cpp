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
            this->lastFilesNum = std::stoi(matchFileNum.str(1));
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
//        std::cout<<"sorted file num="<<sortedsAllFilesAll.size()<<std::endl;
//        std::cout<<"last files num="<<lastFilesNum<<std::endl;
        for (const auto &oneName: selectedFiles_sAll) {//read all xml files in sVec
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


        }//end of reading all sVec xml files

        int onePartLengh=static_cast<int>(sVecsCombined.size()/2);
//        std::cout<<"one part length="<<onePartLengh<<std::endl;
//        std::cout<<"lag="<<lag<<std::endl;
            for(int i=0;i<onePartLengh;i+=lag){
                this->sSelected.push_back(sVecsCombined[i]);
            }
            for(int i=onePartLengh;i<2*onePartLengh;i+=lag){
                this->sSelected.push_back(sVecsCombined[i]);
            }

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
//        std::cout<<"selected file num="<<selectedFilesEAll.size()<<std::endl;

        for (const auto &oneName: selectedFilesEAll) {//read all xml files in EAll
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


        }//end of reading all xml files in EAll

        int onePartLength = static_cast<int>(EAllCombined.size() / 2);
        for (int i = 0; i < onePartLength; i += lag) {
            this->ESelected.push_back(EAllCombined[i]);
        }
        for (int i = onePartLength; i < 2 * onePartLength; i += lag) {
            this->ESelected.push_back(EAllCombined[i]);
        }

    }//end of paramagnetic case

//    int i1=1000;
//int i2=int(ESelected.size()*0.75);
//std::cout<<i1<<" th value="<<ESelected[i1]<<std::endl;
//std::cout<<i2<<" th value="<<ESelected[i2]<<std::endl;

}