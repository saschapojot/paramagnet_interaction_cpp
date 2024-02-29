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