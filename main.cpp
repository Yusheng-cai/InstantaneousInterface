#include "src/Driver.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    CommandLineArguments cmd(argc, argv);
    InputParser ip;
    std::string fname = argv[1];

    ParameterPack pack;
    ip.ParseFile(fname,pack);

    Driver d(pack, cmd);

    d.update();
    d.calculate();
}