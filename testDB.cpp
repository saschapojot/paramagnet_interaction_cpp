#include "dbExchangeModel.hpp"

int main(int argc, char **argv) {

    int i = 2;
    for (int j = 0; j < 10; j++) {
        if (j == i) {
            continue;
        } else { std::cout << j << std::endl; }
    }


}