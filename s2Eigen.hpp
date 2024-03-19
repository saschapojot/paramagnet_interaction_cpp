//
// Created by polya on 2/29/24.
//

#ifndef PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
#define PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
#include "dbExchangeModel.hpp"



#include <string>



//Read from the xml files and solve eigenvalue problem, then map the folded bands to unfolded bands
//Also, compute specific heat
//needs to read E, mu, and sAvg

class reader {

public:
    reader(const int &groupNum,const int& rowNum, const std::string &TDir) {
        this->groupNum = groupNum;
        this->rowNum=rowNum;
        this->TDir ="./group"+std::to_string(groupNum)+"data/row"+std::to_string(rowNum)+"/"+ TDir;

        std::regex TPattern("T([+-]?\\d*(\\.\\d+)?)");
        std::smatch T_match;
        if (std::regex_search(TDir, T_match, TPattern)) {

            this->T = std::stod(T_match.str(1));
//            std::cout<<T<<std::endl;
            this->beta = 1.0 / T;
        }
//        std::cout<<TDir<<std::endl;
//        std::cout<<T<<std::endl;
        std::vector<double> x0{1, 0};
        std::vector<double> x1{0, 1};
        this->xVecsAll.push_back(x0);
        this->xVecsAll.push_back(x1);
//        std::cout<<"finishing constructor"<<std::endl;

    }

    ///search Ell, muAll, sAll folder's files
    void searchFiles();

    ///
    /// @param path the path containing xml files
    /// @return sorted xml files by starting loop
    std::vector <std::string> sortOneDir(const std::vector <std::string> &allFiles);


    ///sort files by starting loop
    void sortFiles();


    template<class T>
    std::vector <size_t> argsort(const std::vector <T> &v) {
        std::vector <size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] <= v[i2]; });
        return idx;
    }

    ///parse Txxx.xxxchi.txt to extract lag, lastElemNum, lastFilesNum
    void parseCHiFile();

    /// parse sVec values in sAll directory
    void parse_sAllDir();

    template<class T>
    static void printVec(const std::vector <T> &vec) {
        for (int i = 0; i < vec.size() - 1; i++) {
            std::cout << vec[i] << ",";
        }
        std::cout << vec[vec.size() - 1] << std::endl;
    }


    ///parse EAvg values in EAll directory
    void parse_EAllDir();

    ///parse mu values in muAll directory
    void parse_muAllDir();

    ///compute Ze weights and projections
    void fillZeWeights(dbExchangeModel &model);

    ///parse last file under EAll,muAll, sAll
    void parseLast(dbExchangeModel &model);


    ///
    /// @param s_ind index of s vector
    /// @return WE for each s
    double computeWE(const int &s_ind);

    //compute WE2 for each s


    /// @param s_ind index of s vector
    /// @return WE2 for each s
    double computeWE2(const int &s_ind);


    ///compute all WE and WE2
    void computeAllWEAllWE2();

    ///
    /// @param i index of s, to be deleted in computation
    /// @return pseudovalue of \partial_{\beta}\braket{\epsilon} with si deleted
    double dbeta_epsilon(const int &i);

    ///
    /// @param ps pseudovalue of C
    /// @param sd standard deviation of C
    void pseudoValueOfC(double &ps, double &sd);

    ///compute C and write to file
    void C2File();

    /// construct all y vectors
    /// @param model dbExchangeModel object
    void construct_yAll(const dbExchangeModel &model);

    ///
    /// @tparam valType
    /// @param vec
    /// @return norm of vector
    template<class valType>
   static double vectorNorm(const std::vector <valType> &vec) {

        double sum = 0.0;
        for (const auto &val: vec) {
            sum += std::pow(std::abs(val), 2);
        }
        return std::sqrt(sum);
    }
    ///initialize all A matrices, AMatsAll indexed by s,m,a,j
    void initAMatsAll();

    ///compute the average and error of eigenvalue
    void computeMeanE();


    ///compute the marker size of unfolded band
    void computeMarkerSize();

    /// add the value of RHS to LHS
    /// @param LHS vector
    /// @param RHS vector
    static void addRightVecToLeft(std::vector<double>& LHS, const std::vector<double> & RHS);

    ///
    /// @param LHS vector
    /// @param RHS vector
    /// @return LHS-RHS
    static std::vector<double> leftMinusRightVec(const std::vector<double> &LHS,const std::vector<double> &RHS);


    ///
    /// @param vec vector
    /// @return elementwise squared vector
    static std::vector<double> vectorSquared(const std::vector<double> &vec);

    ///
    /// @return flattened data
    static std::vector<double> flattenData(const std::vector<std::vector<std::vector<double>>>& data);

    ///write meanE, hfLengthE, markerSize to file
    void bandToFile();

public:
    std::vector <std::string> EFilesAll;
    std::vector <std::string> muFilesAll;
    std::vector <std::string> sAllFilesAll;
    std::vector <std::string> sortedEFilesAll;
    std::vector <std::string> sortedmuFilesAll;
    std::vector <std::string> sortedsAllFilesAll;
//    int part = 1;
    int groupNum;
    int rowNum;

    double T = 0;
    double beta = 0;
    std::string TDir;
    std::string EPath;
    std::string muPath;
    std::string sAllPath;
//    double C = 0;
    int L = 10;
    int M = 20;

    bool ferro=false;
    int lag = 0;
    int lastFilesNum = 0;//used for paramagnetic case
    int lastElemNum = 0;//used for ferromagnetic case

    std::vector <std::vector<double>> sSelected;// selected s vectors for computation of unfolded bands
    std::vector<double> ESelected;//selected E values for computation of unfolded bands and C
    std::vector<double> muSelected;//selected mu values for computation of unfolded bands and C
//    std::vector<double> KSupValsAll;//all the values in SBZ
//    std::vector<int> KSupIndsAll;//[0,1,...,M-1]

    std::vector<double> muRecomputed;// the reading from muAll is not accurate, we recompute mu.
    std::vector<double> epsilonSelected;// energy per supercell
    std::vector <std::vector<std::vector < double>>>
    ZeWeightsAll;// weights for all s, for all values in SBZ, for all eigenvalues given one value in SBZ
    std::vector<double> ZeVals;// Ze values for each s
    std::vector <std::vector<std::vector < double>>>
    eigValsRecomputed;// recomputed values for eigenvalue,  for all s, for all values in SBZ, for all eigenvalues
    std::vector<double> WE;// W(E(\beta,s)) for each s
    std::vector<double> WE2;// W(E^{2}(\beta,s)) for each s

    std::vector <std::vector<double>> xVecsAll;
    std::vector <std::vector<std::vector < std::vector < std::complex < double>>>>>
    yVecsAll;//m,a,j,vec

  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> AMatsAll;//s,m,a,
    std::vector<std::vector<std::vector<double>>> meanE;//m,a
    std::vector<std::vector<std::vector<double>>>hfLengthE;//m,a

    std::vector<std::vector<std::vector<double>>> markerSize;//m,a
};





#endif //PARAMAGNET_INTERACTION_CPP_S2EIGEN_HPP
