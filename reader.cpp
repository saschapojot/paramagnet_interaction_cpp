#include "sAvg.hpp"

#include <cstdio>
#include <iostream>
#include <memory>
#include <array>
#include <exception>
std::string exec(const char* cmd) {
    std::array<char, 512> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output
}

int main() {
    std::string command = "python3 checkVec.py ./part1/T2.100000/EAll/";
    std::string result="";
    try {
       result = exec(command.c_str());
        std::cout << "Output from Python script: " << result;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    catch (...) {
        // Handle any other exceptions
        std::cerr << "Error" << std::endl;
    }
    return 0;
}



//int main(int argc, char **argv){
////    const auto tStart{std::chrono::steady_clock::now()};
////    auto loader=loaderAndComputer();
////    loader.searchFiles();
////    loader.filesAllSorted();
////    loader.fillIntodataStorage();
////    loader.data2PhysicalQuantities();
////
////    loader.physicalQuantities2csv();
////
//////    loader.diagostics(loader.storages.size()-1);
////
////    const auto tEnd{std::chrono::steady_clock::now()};
////    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};
////
////    std::cout<<"total time: "<<elapsed_secondsAll.count()/3600.0<<" h"<<std::endl;
////
//
//
//
//
//
//
//
//}
//
//
//
//


//#include <iostream>
//#include <regex>
//#include <string>
//
//int main() {
//    std::string text = "Example text with some keywords: apple, banana, orange";
//    std::regex keyword_regex("(apple|banana|orange)"); // Define your keywords here
//
//    std::smatch match;
//    while (std::regex_search(text, match, keyword_regex)) {
//        std::cout << "Found keyword: " << match.str(0) << std::endl;
//        text = match.suffix().str(); // Proceed to the rest of the string
//    }
//
//    return 0;
//}
