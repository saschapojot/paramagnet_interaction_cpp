#include "sAvg.hpp"



void loaderAndComputer::searchFiles() {

    for (const auto &entry: fs::directory_iterator(this->searchPath)) {
        std::regex ptEAll("EAll");
        std::regex ptsAll("sAll");
        std::regex ptsmuAll("muAll");

        std::string fileTmp = entry.path().string();
        if (std::regex_search(fileTmp, ptEAll)) {
            this->EFilesAll.push_back(fileTmp);
        }

        if (std::regex_search(fileTmp, ptsAll)) {
            this->sAllFilesAll.push_back(fileTmp);
        }
        if (std::regex_search(fileTmp, ptsmuAll)) {
            this->muFilesAll.push_back(fileTmp);
        }


    }





}

template <typename T>
std::vector<size_t> loaderAndComputer::sort_indexes(const std::vector<T> &v){

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;


}



std::vector<std::string>  loaderAndComputer::sortFiles(const std::vector<std::string> & fileVector) {
    std::vector<double> temperatures;

    std::regex patternTt("T(-?\\d+\\.\\d+)t");
    std::smatch match;
    for (const auto &oneFile: fileVector) {

        if (std::regex_search(oneFile, match, patternTt)) {
            std::string Tstr = match[1];
            temperatures.push_back(std::stod(Tstr));
        }


    }


    std::vector<size_t> idxTemperature = this->sort_indexes(temperatures);

    std::vector<std::string> sortedFileVector;
    for (const auto &id: idxTemperature) {
        sortedFileVector.push_back(fileVector[id]);
    }

    return sortedFileVector;



}



void loaderAndComputer::filesAllSorted(){

    this->sortedsAllFilesAll=this->sortFiles(sAllFilesAll);
    this->sortedEFilesAll=this->sortFiles(EFilesAll);
    this->sortedmuFilesAll=this->sortFiles(muFilesAll);

}


