#include "wavobser_class_version.h"
#include "iostream"


int main(int argc, char* argv[]) {
    WaveSimulation sim;
    printf("Starting wave simulation with %d iterations...\n", argc > 1 && argc > 0 ? std::atoi(argv[1]) : 10);
    sim.run(argc > 1 ? std::atoi(argv[1]) : 10); // Run with a limit if provided, default to 10
    return 0;
}