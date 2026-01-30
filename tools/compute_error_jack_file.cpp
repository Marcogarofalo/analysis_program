#define CONTROL

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>

#include "global.hpp"
#include "resampling_new.hpp"

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " jack_file" << std::endl;
        return 1;
    }
    std::ifstream file(argv[1]);
    std::string riga;
    int line_count = 0;

    if (file.is_open()) {
        while (std::getline(file, riga)) {
            line_count++;
        }
        file.close();
        std::cout << "lines: " << line_count << "  jack(+1 of the mean): " << line_count - 1 << std::endl;
    }
    else {
        std::cerr << "Impossibile aprire il file." << std::endl;
    }


    resampling_f* myj = new resampling_jack(line_count - 2);


    double* tmp = (double*)malloc(sizeof(double) * myj->Njack);
    myj->read_jack_from_file(tmp, argv[1]);
    printf("mean = %.12f  +-  %.12f\n", myj->mean(tmp), myj->comp_error(tmp));

}