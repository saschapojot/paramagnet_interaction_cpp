#include "s2Eigen.hpp"


std::vector<std::string> scanFiles(const int &part){
    std::string searchPath="./part"+std::to_string(part)+"/";
    std::vector<std::string> TDirs;
    if(fs::exists(searchPath) && fs::is_directory(searchPath)){
        for(const auto &entry:fs::directory_iterator(searchPath)){
            TDirs.push_back(entry.path().filename());
        }
    }
//    for(const auto&s:TDirs)
//    {
//        std::cout<<s<<std::endl;
//    }

    return TDirs;

}


int main(int argc, char **argv) {

    int part=1;
    std::vector<std::string> TDirs= scanFiles(part);
//    for (const auto& s:TDirs) {
//        const auto tCStart{std::chrono::steady_clock::now()};
//        auto rd = reader(part, s);
//
//        auto model = dbExchangeModel(rd.T);
//
//
//        rd.searchFiles();
//        rd.sortFiles();
//        rd.parseCHiFile();
//        rd.parse_sAllDir();
//        rd.parse_EAllDir();
//
//        rd.fillZeWeights(model);
//
//        rd.computeAllWEAllWE2();
//
//        rd.C2File();
//        const auto tCEnd{std::chrono::steady_clock::now()};
//        const std::chrono::duration<double> elapsed_seconds{tCEnd - tCStart};
//        std::cout << "time elapsed: " << elapsed_seconds.count() << " s" << std::endl;
//    }
//
////////////////////////one file

    std::string name0="T0.111616";
    auto rd=reader(1,name0);
    auto model=dbExchangeModel(rd.T);

    rd.searchFiles();
    rd.sortFiles();
    rd.parseCHiFile();
    rd.parse_sAllDir();
    rd.parse_EAllDir();

    rd.construct_yAll(model);
    rd.initAMatsAll();

    rd.fillZeWeights(model);
    int p=10;
    int m=7;
    int a=2;
    rd.printVec(rd.AMatsAll[p][m][a][0]);
    rd.printVec(rd.AMatsAll[p][m][a][1]);
    rd.printVec(rd.AMatsAll[p][m][a][2]);
    rd.printVec(rd.AMatsAll[p][m][a][3]);
    rd.printVec(rd.AMatsAll[p][m][a][4]);
    rd.printVec(rd.AMatsAll[p][m][a][5]);


    rd.computeAllWEAllWE2();

    rd.C2File();







}
