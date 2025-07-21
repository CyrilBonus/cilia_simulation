#pragma once
#include <vector>
#include <cmath>
#include <cstdio> 

struct GreenFunctions {
    static double G11(double, double, double, double, double, double, double, double);
    static double G12(double, double, double, double, double, double, double, double);
    static double G13(double, double, double, double, double, double, double, double);
    static double G21(double, double, double, double, double, double, double, double);
    static double G22(double, double, double, double, double, double, double, double);
    static double G23(double, double, double, double, double, double, double, double);
    static double G31(double, double, double, double, double, double, double, double);
    static double G32(double, double, double, double, double, double, double, double);
    static double G33(double, double, double, double, double, double, double, double);
};

class WaveSimulation {
public:
    double omega, omega1, omega2, OMEGA, A1, B1, B2, A2, qx, qy, a;
    double eta, friction, theta, kay, R0;
    double g2[3], h2[3];
    int we;

    WaveSimulation();

    double force(double p);

    void run(int limit = 10);

private:
    void initialize_arrays(int num_x, int num_y);
    void fill_active(const std::vector<std::pair<int, int>>& coords, int min_x, int min_y);
    void write_output(FILE* fp, double t, const std::vector<std::vector<double>>& phi, const std::vector<std::vector<double>>& p2, const std::vector<std::vector<int>>& active);
};