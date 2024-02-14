//
// Created by polya on 2/13/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_SAVG_HPP
#define PARAMAGNET_INTERACTION_CPP_SAVG_HPP
#include "dbExchangeModel.hpp"
#include <string>
#include <regex>
#include <numeric>

#include <sstream>




namespace fs = std::filesystem;
class dataStorage{
public:
    std::vector<std::vector<double>>sAll;
    std::vector<double>EAll;
    std::vector<double>muAll;
    std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>>flattenedEigSolution;


};



class loaderAndComputer{
    //load data and compute physical observables
public:
    std::vector<std::string> EFilesAll;
    std::vector<std::string>muFilesAll;
    std::vector<std::string>sAllFilesAll;
    std::vector<std::string> eigFilesAll;

    std::vector<std::string> sortedEFilesAll;
    std::vector<std::string>sortedmuFilesAll;
    std::vector<std::string>sortedsAllFilesAll;
    std::vector<std::string> sortedEigFilesAll;
    int part=0;
    std::vector<dataStorage> storages;





    std::string searchPath="./part"+std::to_string(part)+"/";

public:
    void searchFiles();//search .sAll, .EAll, .muAll, .flattenedEigSolution files

    ///
    /// @param fileVector a vector of filenames
    /// @return filenames sorted by T
    std::vector<std::string> sortFiles(const std::vector<std::string> & fileVector);//sort .sAll, .EAll, .muAll, .flattenedEigSolution files by T

    void filesAllSorted();//get the sorted file names

    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v);//sort and return index


    ///
    /// @param fileName a file that contains serialized vector<double>
    std::vector<double> readVecDouble(const std::string &fileName);

    ///
    /// @param fileName a file that contains serialized vector<vector<double>>
    std::vector<std::vector<double>> readVecVecDouble(const std::string &fileName);

    ///
    /// @param fileName a file that contains serialized vector<vector<tuple>>
    std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>> readVecVecTuple(const std::string &fileName);

    void fillIntodataStorage();// fill  deserialized values into  storages

};


#endif //PARAMAGNET_INTERACTION_CPP_SAVG_HPP
