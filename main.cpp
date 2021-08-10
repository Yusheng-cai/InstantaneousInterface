#include "src/Driver.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"

#include <iostream>
#include <string>
#include <chrono>

int main(int argc, char** argv)
{
    CommandLineArguments cmd(argc, argv);
    InputParser ip;
    std::string fname = argv[1];

    ParameterPack pack;
    ip.ParseFile(fname,pack);

    Driver d(pack, cmd);

    auto start = std::chrono::high_resolution_clock::now();
    d.update();
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "Time it took for update is " << diff.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    d.calculate();
    end = std::chrono::high_resolution_clock::now();
    auto diff2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "Time it took for calculate is " << diff2.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    d.finishCalculate();
    end = std::chrono::high_resolution_clock::now();
    auto diff3 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "Time it took for finish calc is " << diff3.count() << std::endl;

    d.printOutputfile();
}