

#include "dbExchangeModel.hpp"
#include <boost/math/distributions/kolmogorov_smirnov.hpp>

#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <Eigen/Dense>

// Include Eigen headers
#include <Eigen/Dense>

// Define a sample data structure with an Eigen matrix
struct MyData {
    Eigen::Matrix2f matrix;
    std::string name;

    // Serialization function for Eigen matrix
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & (const int  &)matrix.rows();
        ar & (const int  &)matrix.cols();
        ar & boost::serialization::make_array(matrix.data(), matrix.size());
        ar & name;
    }
};

int main() {
    // Create an instance of MyData with an Eigen matrix
    MyData data;
    data.matrix = Eigen::Matrix2f::Random(3, 3);
    data.name = "Eigen Matrix Serialization";

    // Serialize the data to a file
    {
        std::ofstream ofs("data.txt");
        boost::archive::text_oarchive ar(ofs);
        ar << data;
    }

//     Deserialize the data from the file
    MyData loadedData;
    {
        std::ifstream ifs("data.txt");
        boost::archive::text_iarchive ar(ifs);
        ar >> loadedData;
    }

    // Display the loaded data
    std::cout << "Loaded data:\n" << loadedData.matrix << "\nName: " << loadedData.name << std::endl;

    return 0;
}
