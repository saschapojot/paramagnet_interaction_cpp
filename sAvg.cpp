#include "sAvg.hpp"



void loaderAndComputer::searchFiles() {

    for (const auto &entry: fs::directory_iterator(this->searchPath)) {
        std::regex ptEAll("EAll");
        std::regex ptsAll("sAll");
        std::regex ptsmuAll("muAll");
//        std::regex ptEigAll("flattenedEigSolution");

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

//        if (std::regex_search(fileTmp,ptEigAll)){
//            this->eigFilesAll.push_back(fileTmp);
//        }


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
    std::vector<double> sortedTemperatures;
    for (const auto &id: idxTemperature) {
        sortedFileVector.push_back(fileVector[id]);
        sortedTemperatures.push_back(temperatures[id]);
    }

    this->TAll=std::vector<double>(sortedTemperatures);//TODO: performed multiple times




    return sortedFileVector;



}



void loaderAndComputer::filesAllSorted(){

    this->sortedsAllFilesAll=this->sortFiles(sAllFilesAll);
    this->sortedEFilesAll=this->sortFiles(EFilesAll);
    this->sortedmuFilesAll=this->sortFiles(muFilesAll);
//    this->sortedEigFilesAll=this->sortFiles(eigFilesAll);

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
//    size_t numEigFiles=this->sortedEigFilesAll.size();

    int diff0=numEFiles-numMuFiles;
    int diff1=numEFiles-numsFiles;
//    int diff2=numEFiles-numEigFiles;

    int n2=diff0*diff0+diff1*diff1;
            //+diff2*diff2;
    if (n2!=0){
        std::cout<<"file numbers unmatch"<<std::endl;
        std::exit(3);
    }

    for (size_t i=0;i<numEFiles;i++){
        dataStorage stgTmp;
//        std::cout<<"reading file "<<this->sortedEFilesAll[i]<<std::endl;
        stgTmp.EAll=this->readVecDouble(this->sortedEFilesAll[i]);

//        std::cout<<"reading file "<<this->sortedmuFilesAll[i]<<std::endl;

        stgTmp.muAll=this->readVecDouble(this->sortedmuFilesAll[i]);
//        std::cout<<"reading file "<<this->sortedsAllFilesAll[i]<<std::endl;
        stgTmp.sAll=this->readVecVecDouble(this->sortedsAllFilesAll[i]);
//        std::cout<<"reading file "<<this->sortedEigFilesAll[i]<<std::endl;
//        stgTmp.flattenedEigSolution=this->readVecVecTuple(this->sortedEigFilesAll[i]);
        this->storages.push_back(stgTmp);


    }


};

///
/// @param oneData one object of dataStorage
/// @param T temperature of oneData
/// @return mean s, chi, C
void loaderAndComputer::computePhysicalQuantities(const dataStorage& oneData,const double&T) {
    int startingPosition = 10000;
    int sep = 20;

    //compute s
    ////////////////////////////////////////////
    std::vector<std::vector<double>> sSelected;

    if (startingPosition < oneData.sAll.size()) {
        auto iter = oneData.sAll.begin() + startingPosition;

        while (iter < oneData.sAll.end()) {
            std::vector<double> oneS = *iter;
            sSelected.push_back(oneS);
            std::advance(iter, sep);
        }
    }

    std::vector<double> sAvg;//average value of spins for one mc step
    for (const auto &vec: sSelected) {
        sAvg.push_back(this->average(vec));
    }
    std::vector<double> sAvgAbs;//absolute value of average value of spins for one mc step
    for (const auto &val: sAvg) {
        sAvgAbs.push_back(std::abs(val));
    }

    double absAvgS = average(sAvgAbs);
    this->sAvgAbsAll.push_back(absAvgS);

    ////////////////////////////////////////////////

    //compute chi
    ///////////////////////////////////////////////
    double meanAllS = this->average(sAvg);
    std::vector<double> sSquared;
    for (const auto &val: sAvg) {
        sSquared.push_back(val * val);
    }

    double meanS2 = this->average(sSquared);
    double chiTmp = (meanS2 - meanAllS * meanAllS) / T;

    this->chiAll.push_back(chiTmp);

    /////////////////////////////////////////////////

    //compute C
    std::vector<double> ESelected;

    if (startingPosition < oneData.EAll.size()) {
        auto iter = oneData.EAll.begin() + startingPosition;
        while (iter < oneData.EAll.end()) {
            double ETmp = *iter;
            ESelected.push_back(ETmp);
            std::advance(iter, sep);

        }


    }

    double meanE=this->average(ESelected);
    std::vector<double> ESelected2;
    for (auto const &val:ESelected){
        ESelected2.push_back(val*val);
    }
    double meanE2=this->average(ESelected2);

    double CTmp=(meanE2-meanE*meanE)/(T*T);

    this->spHeatAll.push_back(CTmp);


}


///
/// @param vec a vector of double
/// @return average of values in vec
double loaderAndComputer::average(const std::vector<double> &vec){
    if(vec.empty()){
        return 0;
    }
    double const count = static_cast<double>(vec.size());
    return std::reduce(vec.begin(), vec.end()) / count;


}



///perform computePhysicalQuantities on each element of storages
void loaderAndComputer::data2PhysicalQuantities(){
    size_t num=this->storages.size();
    for (size_t i=0;i<num;i++){
        this->computePhysicalQuantities(this->storages[i],this->TAll[i]);

    }


}


///store the physical quantities in a csv file
void loaderAndComputer::physicalQuantities2csv(){
    std::string outCsvName=searchPath+"reduced.csv";
    std::ofstream csvOut(outCsvName);

    csvOut<<"T,s,chi,C"<<std::endl;
    size_t  vecLength=this->TAll.size();
    for (size_t i=0;i<vecLength;i++){
        std::vector<double> oneRow{this->TAll[i],this->sAvgAbsAll[i],this->chiAll[i],this->spHeatAll[i]};
        this->write1Line2Csv(oneRow,csvOut);
    }

    csvOut.close();



}

///
/// @param oneRow a row of numbers to output
/// @param outF
void loaderAndComputer::write1Line2Csv(const std::vector<double> &oneRow, std::ofstream& outF){
    size_t num=oneRow.size();
    for(size_t i=0;i<num-1;i++){
        outF<<oneRow[i]<<",";
    }
    outF<<oneRow[num-1]<<std::endl;



}