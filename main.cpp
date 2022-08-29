#include "src/Driver.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"

#include <iostream>
#include <string>
#include <chrono>

int main(int argc, char** argv)
{
    CommandLineArguments cmd(argc, argv);

    // initialize the driver 
    InputParser ip;
    std::string fname = argv[1];

    ParameterPack pack;
    ip.ParseFile(fname,pack);

    Driver d(pack, cmd);
    d.run();
}