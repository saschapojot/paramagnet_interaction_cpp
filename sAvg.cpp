#include "sAvg.hpp"



void loaderAndComputer::searchFiles() {

    for (const auto &entry: fs::directory_iterator(this->searchPath)) {
        std::regex ptEAll("EAll");
        std::regex ptsAll("sAll");
        std::regex ptsmuAll("muAll");
        std::regex ptEigAll("flattenedEigSolution");

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

        if (std::regex_search(fileTmp,ptEigAll)){
            this->eigFilesAll.push_back(fileTmp);
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
    this->sortedEigFilesAll=this->sortFiles(eigFilesAll);

}

///
/// @param fileName a file that contains serialized vector<double>
std::vector<double> loaderAndComputer::readVecDouble(const std::string &fileName) {
    std::ifstream inF(fileName, std::ios::in | std::ios::binary);
    std::stringstream ssTmp;
    ssTmp << inF.rdbuf();
    inF.close();

    std::string contentAll = ssTmp.str();
    msgpack::object_handle oh = msgpack::unpack(contentAll.data(), contentAll.size());
    msgpack::object deserialized = oh.get();
    auto vec = deserialized.as<std::vector<double>>();

    return vec;


}


///
/// @param fileName a file that contains serialized vector<vector<double>>
std::vector<std::vector<double>> loaderAndComputer::readVecVecDouble(const std::string &fileName){
    std::ifstream inF(fileName,std::ios::in|std::ios::binary);
    std::stringstream ssTmp;
    ssTmp << inF.rdbuf();
    inF.close();

    std::string contentAll = ssTmp.str();
    msgpack::object_handle oh = msgpack::unpack(contentAll.data(), contentAll.size());
    msgpack::object deserialized = oh.get();

    auto vecvec=deserialized.as<std::vector<std::vector<double>>>();

    return vecvec;




}


///
/// @param fileName a file that contains serialized vector<vector<tuple>>
std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>> loaderAndComputer::readVecVecTuple(const std::string &fileName){
    std::ifstream inF(fileName,std::ios::in|std::ios::binary);
    std::stringstream ssTmp;
    ssTmp << inF.rdbuf();
    inF.close();

    std::string contentAll = ssTmp.str();
    msgpack::object_handle oh = msgpack::unpack(contentAll.data(), contentAll.size());
    msgpack::object deserialized = oh.get();
    auto vecvectp=deserialized.as<std::vector<std::vector<std::tuple<int,std::vector<double>,std::vector<std::complex<double>> >>>>();

    return vecvectp;


}

// fill  deserialized values into  storages
void loaderAndComputer::fillIntodataStorage(){
    size_t numEFiles=this->sortedEFilesAll.size();
    size_t numMuFiles=this->sortedmuFilesAll.size();
    size_t numsFiles=this->sortedsAllFilesAll.size();
    size_t numEigFiles=this->sortedEigFilesAll.size();

    int diff0=numEFiles-numMuFiles;
    int diff1=numEFiles-numsFiles;
    int diff2=numEFiles-numEigFiles;

    int n2=diff0*diff0+diff1*diff1+diff2*diff2;
    if (n2!=0){
        std::cout<<"file numberd unmatch"<<std::endl;
        std::exit(3);
    }

    for (size_t i=0;i<numEFiles;i++){
        dataStorage stgTmp;
        stgTmp.EAll=this->readVecDouble(this->sortedEFilesAll[i]);
        stgTmp.muAll=this->readVecDouble(this->sortedmuFilesAll[i]);
        stgTmp.sAll=this->readVecVecDouble(this->sortedsAllFilesAll[i]);
        stgTmp.flattenedEigSolution=this->readVecVecTuple(this->sortedEigFilesAll[i]);
        this->storages.push_back(stgTmp);


    }

};
