#include "sAvg.hpp"

//#include <cstdio>
//#include <iostream>
//#include <memory>
//#include <array>
//
//std::string exec(const char* cmd) {
//    std::array<char, 128> buffer; // Buffer to store command output
//    std::string result; // String to accumulate output
//
//    // Open a pipe to read the output of the executed command
//    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
//    if (!pipe) {
//        throw std::runtime_error("popen() failed!");
//    }
//
//    // Read the output a chunk at a time and append it to the result string
//    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
//        result += buffer.data();
//    }
//
//    return result; // Return the accumulated output
//}
//
//int main() {
//    std::string command = "python3 your_script.py";
//    try {
//        std::string result = exec(command.c_str());
//        std::cout << "Output from Python script: " << result;
//    } catch (const std::runtime_error& e) {
//        std::cerr << "Error: " << e.what() << std::endl;
//    }
//    return 0;
//}
//


int main(int argc, char **argv){
//    const auto tStart{std::chrono::steady_clock::now()};
//    auto loader=loaderAndComputer();
//    loader.searchFiles();
//    loader.filesAllSorted();
//    loader.fillIntodataStorage();
//    loader.data2PhysicalQuantities();
//
//    loader.physicalQuantities2csv();
//
////    loader.diagostics(loader.storages.size()-1);
//
//    const auto tEnd{std::chrono::steady_clock::now()};
//    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};
//
//    std::cout<<"total time: "<<elapsed_secondsAll.count()/3600.0<<" h"<<std::endl;
//



    namespace fs = boost::filesystem;
    int part=1;
    //output directory
    double T=1;
    std::string outEAllSubDir = "./part" + std::to_string(part) + "/T"+std::to_string(T)+"/EAll/";
    std::string outMuAllSubDir = "./part" + std::to_string(part) + "/T"+std::to_string(T)+"/muAll/";
    std::string outSAllSubDir = "./part" + std::to_string(part) + "/T"+std::to_string(T)+"/sAll/";
    std::string outEigAllSubDir = "./part" + std::to_string(part) + "/T"+std::to_string(T)+"/eigAll/";


    if(!fs::is_directory(outEAllSubDir) || !fs::exists(outEAllSubDir)){
        fs::create_directories(outEAllSubDir);
    }




}




