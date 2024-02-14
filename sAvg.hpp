//
// Created by polya on 2/13/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_SAVG_HPP
#define PARAMAGNET_INTERACTION_CPP_SAVG_HPP
#include <string>
#include <iostream>
#include <filesystem>
#include <vector>
#include <regex>
#include <msgpack.hpp>
#include <numeric>
#include <algorithm>
namespace fs = std::filesystem;
class loaderAndComputer{
    //load data and compute physical observables
public:
    std::vector<std::string> EFilesAll;
    std::vector<std::string>muFilesAll;
    std::vector<std::string>sAllFilesAll;

    std::vector<std::string> sortedEFilesAll;
    std::vector<std::string>sortedmuFilesAll;
    std::vector<std::string>sortedsAllFilesAll;
    int part=0;

    std::string searchPath="./part"+std::to_string(part)+"/";

public:
    void searchFiles();//search .sAll, .EAll, .muAll files

    ///
    /// @param fileVector a vector of filenames
    /// @return filenames sorted by T
    std::vector<std::string> sortFiles(const std::vector<std::string> & fileVector);//sort .sAll, .EAll, .muAll files by T

    void filesAllSorted();

    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v);

};

#endif //PARAMAGNET_INTERACTION_CPP_SAVG_HPP
