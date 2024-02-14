
#include "dbExchangeModel.hpp"






int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }

    double T = std::stod(argv[1]);
    auto model = dbExchangeModel(T);
    model.executionMC();
    model.data2File(model.record);


    // Display the loaded data
    std::cout << "Loaded data:\n" << loadedData.matrix << "\nName: " << loadedData.name << std::endl;

    return 0;
}
