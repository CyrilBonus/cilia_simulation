#pragma once
#include <vector>
#include <string>
#include <utility>

class Reader {
public:
    std::vector<std::pair<int, int>> coords;
    bool read_csv(const std::string& filename);
};