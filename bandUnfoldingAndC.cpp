#include "s2Eigen.hpp"


int main(int argc, char **argv) {

    auto rd=reader(1,"T2.175758");
    rd.searchFiles();
    rd.sortFiles();

}
