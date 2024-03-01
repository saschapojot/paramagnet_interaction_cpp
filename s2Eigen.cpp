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
void reader::parse_sAllDir(){
//    std::cout<<"ferro="<<ferro<<std::endl;
  if(this->ferro==true){// ferro case, read last file
      std::vector<std::vector<double>> sAll;
//      std::cout<<"file num="<<sortedsAllFilesAll.size()<<std::endl;
//      std::cout<<"file name: "<<this->sortedsAllFilesAll[sortedsAllFilesAll.size()-1]<<std::endl;
      std::ifstream ifs(this->sortedsAllFilesAll[sortedsAllFilesAll.size()-1]);
      if(!ifs.is_open()){
          std::cerr<<"cannot open"<<std::endl;
          return;
      }
      boost::archive::xml_iarchive ia(ifs);
//      vecvec sAllVecvec;
      ia>> BOOST_SERIALIZATION_NVP(sAll);
//      for(const auto& vec: sAllVecvec.data){
//          sAll.push_back(vec);
//      }
        for (int i=sAll.size()-lastElemNum;i<sAll.size();i++){
            this->sSelected.push_back(sAll[i]);
        }
//        std::cout<<sSelected.size()<<std::endl;





  }// end of ferro case
  else{//paramagnetic case, read last lastFilesNum files
      std::vector<std::string> selectedFiles_sAll;
      for (int i=this->sortedsAllFilesAll.size()-lastFilesNum;i<this->sortedsAllFilesAll.size();i++){
      selectedFiles_sAll.push_back(this->sortedsAllFilesAll[i]);

      }
      std::vector<std::vector<double>> sAll;





  }//end of paramagnetic case



}