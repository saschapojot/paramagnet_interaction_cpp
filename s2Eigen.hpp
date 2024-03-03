//
// Created by polya on 2/29/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
#define PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
#include "dbExchangeModel.hpp"



#include <string>


namespace fs = std::filesystem;
//Read from the xml files and solve eigenvalue problem, then map the folded bands to unfolded bands
//Also, compute specific heat
//needs to read E, mu, and sAvg

class reader {

public:
    reader(const int &p, const std::string &TDir) {
        this->part = p;
        this->TDir = "./part" + std::to_string(p) + "/" + TDir;

        std::regex TPattern("T([+-]?\\d*(\\.\\d+)?)");
        std::smatch T_match;
        if (std::regex_search(TDir, T_match, TPattern)) {

            this->T = std::stod(T_match.str(1));
        }
//        std::cout<<TDir<<std::endl;
//        std::cout<<T<<std::endl;

    }

    ///search Ell, muAll, sAll folder's files
    void searchFiles();

    ///
    /// @param path the path containing xml files
    /// @return sorted xml files by starting loop
    std::vector<std::string> sortOneDir(const std::vector<std::string> &allFiles);


    ///sort files by starting loop
    void sortFiles();



    template<class T>
    std::vector<size_t> argsort(const std::vector<T> &v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] <= v[i2]; });
        return idx;
    }

    ///parse Txxx.xxxchi.txt to extract lag, lastElemNum, lastFilesNum
    void parseCHiFile();

    /// parse sVec values in sAll directory
    void parse_sAllDir();

    template<class T>
            static void printVec(const std::vector<T>& vec){
                for(int i=0;i<vec.size()-1;i++){
                    std::cout<<vec[i]<<",";
                }
                std::cout<<vec[vec.size()-1]<<std::endl;
            }


    ///parse EAvg values in EAll directory
    void parse_EAllDir();


public:
    std::vector<std::string> EFilesAll;
    std::vector<std::string> muFilesAll;
    std::vector<std::string> sAllFilesAll;
    std::vector<std::string> sortedEFilesAll;
    std::vector<std::string> sortedmuFilesAll;
    std::vector<std::string> sortedsAllFilesAll;
    int part = 1;

    double T = 0;
    std::string TDir;
    std::string EPath;
    std::string muPath;
    std::string sAllPath;
    double C = 0;
    int L = 10;
    int M = 20;

    bool ferro;
    int lag=0;
    int lastFilesNum=0;//used for paramagnetic case
    int lastElemNum=0;//used for ferromagnetic case

    std::vector<std::vector<double>>sSelected;// selected s vectors for computation of unfolded bands
    std::vector<double>ESelected;//selected E values for computation of unfolded bands and C


};





#endif //PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
