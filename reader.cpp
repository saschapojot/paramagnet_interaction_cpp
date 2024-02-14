#include "sAvg.hpp"


int main(int argc, char **argv){
    auto loader=loaderAndComputer();
    loader.searchFiles();
    loader.filesAllSorted();

    for(const auto &str:loader.sortedsAllFilesAll){
        std::cout<<str<<std::endl;
    }





}