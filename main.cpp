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

    for (int i=0;i<d.getNumFrames();i++)
    {
        if (d.CheckValidStep(i))
        {
            std::cout << "Frame = " << i << std::endl;
            d.readFrameXdr(i);

            d.update();

            auto start = std::chrono::high_resolution_clock::now();

            d.calculate();

            d.printOutputfileIfOnStep();
            auto end = std::chrono::high_resolution_clock::now();
            auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
            std::cout << "Time it took for calculate is " << diff.count() << " milliseconds." << std::endl;
        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    d.finishCalculate();
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "Time it took for finish calculate is " << diff.count() << std::endl;

    d.printFinalOutput();
}