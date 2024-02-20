#include "sAvg.hpp"



int main(int argc, char **argv){
    const auto tStart{std::chrono::steady_clock::now()};
    auto loader=loaderAndComputer();
    loader.searchFiles();
    loader.filesAllSorted();
    loader.fillIntodataStorage();
    loader.data2PhysicalQuantities();

    loader.physicalQuantities2csv();

    loader.diagostics(loader.storages.size()-1);

    const auto tEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};

    std::cout<<"total time: "<<elapsed_secondsAll.count()/3600.0<<" h"<<std::endl;







}




