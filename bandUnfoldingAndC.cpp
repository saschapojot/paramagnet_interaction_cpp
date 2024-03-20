#include "s2Eigen.hpp"


std::vector<std::string> scanFiles(const int &groupNum, const int& rowNum){
    std::string searchPath="./group"+std::to_string(groupNum)+"data/row"+std::to_string(rowNum)+"/";
    std::vector<std::string> TDirs;
    if(fs::exists(searchPath) && fs::is_directory(searchPath)){
        for(const auto &entry:fs::directory_iterator(searchPath)){
            if(entry.path().filename().string()[0]=='T'){
                TDirs.push_back(entry.path().filename().string());
            }

        }
    }
//    for(const auto&s:TDirs)
//    {
//        std::cout<<s<<std::endl;
//    }

    return TDirs;

}


int main(int argc, char **argv) {
    if(argc!=3){
        std::cerr<<"wrong number of arguments"<<std::endl;
        exit(1);
    }
    int groupNum=std::stoi(argv[1]);
    int rowNum=std::stoi(argv[2]);
    std::vector<std::string> TDirs = scanFiles(groupNum,rowNum);
    reader::printVec(TDirs);
    reader::printVec(TDirs);
    for (const auto& s:TDirs) {
        std::cout<<"file is "<<s<<std::endl;
        const auto tCStart{std::chrono::steady_clock::now()};
        auto rd = reader(groupNum,rowNum, s);

        auto model = dbExchangeModel(rd.T,groupNum,rowNum);


        rd.searchFiles();
        rd.sortFiles();
        rd.parseCHiFile();
        rd.parse_sAllDir();
        rd.parse_EAllDir();
        rd.construct_yAll(model);
    rd.initAMatsAll();

        rd.fillZeWeights(model);

        rd.computeAllWEAllWE2();

        rd.C2File();
        rd.computeMeanE();

        rd.computeMarkerSize();
        rd.bandToFile();
        const auto tCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{tCEnd - tCStart};
        std::cout << "time elapsed: " << elapsed_seconds.count() << " s" << std::endl;
    }

////////////////////////one file
//    const auto tCStart{std::chrono::steady_clock::now()};
//    std::string name0 = "T0.0000100000";
//    auto rd = reader(groupNum,rowNum, name0);
//    auto model = dbExchangeModel(rd.T,groupNum,rowNum);
//
//    rd.searchFiles();
//    rd.sortFiles();
//    rd.parseCHiFile();
//    rd.parse_sAllDir();
//    rd.parse_EAllDir();
//
//    rd.construct_yAll(model);
//    rd.initAMatsAll();
//
//    rd.fillZeWeights(model);
//
//
//    rd.computeAllWEAllWE2();
//
//    rd.C2File();
//
//
//    rd.computeMeanE();
//
//    rd.computeMarkerSize();
//    rd.bandToFile();
//
//    const auto tCEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_seconds{tCEnd - tCStart};
//    std::cout << "time elapsed: " << elapsed_seconds.count() << " s" << std::endl;


}
