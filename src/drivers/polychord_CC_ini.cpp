#include "interfaces.hpp"
#include "CC_ini_likelihood.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc > 1) {
        std::string input_file = argv[1];
        set_ini(input_file);
        run_polychord(loglikelihood,setup_loglikelihood,input_file) ;
        return 0;
    }
    else{
        std::cerr << "PolyChord should be called with at most one argument, the input file" << std::endl;
        return 1;
    }
}

