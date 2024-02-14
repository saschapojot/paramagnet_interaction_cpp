#include "sAvg.hpp"


int main(int argc, char **argv){
    auto loader=loaderAndComputer();
    loader.searchFiles();
    loader.filesAllSorted();
    loader.fillIntodataStorage();

    double mu=loader.storages[1].muAll[7];
    double E=loader.storages[2].EAll[1];

    std::vector<double> s=loader.storages[0].sAll[5];
    std::tuple<int,std::vector<double>,std::vector<std::complex<double>> > onetp=loader.storages[1].flattenedEigSolution[3][7];

    std::cout<<"mu="<<mu<<std::endl;
    std::cout<<"E="<<E<<std::endl;

    std::cout<<"==================="<<std::endl;
    std::cout<<"s: ";
    for(auto val: s){
        std::cout<<val<<", ";
    }
    std::cout<<std::endl;

    std::cout<<"eig="<<std::get<1>(onetp)[9]<<std::endl;


std::cout<<"end"<<std::endl;




}