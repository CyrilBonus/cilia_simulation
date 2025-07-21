#include "reading.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <set>
#include <iostream>

bool Reader::read_csv(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) return false;
    std::string line;
    std::set<std::pair<int, int>> seen_coords; // To track duplicates
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y;
        char comma;
        if (!(iss >> x >> comma >> y)) continue; // skip invalid lines
        int xi = static_cast<int>(std::round(x));
        std::cout << "Read coordinate: (" << xi << ", " << y << ")\n";
        int yi = static_cast<int>(std::round(y));
        if (seen_coords.count({xi, yi})) {
            std::cout << "Duplicate coordinate found: (" << xi << ", " << yi << ")\n";
        }
        seen_coords.emplace(xi, yi);
        coords.emplace_back(xi, yi);
    }
    return true;
}

// tester s'il y a deux coordonnées qui se superposent

// clustering meanshift
