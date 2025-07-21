//quel type de représentation ? Notre simulation contrairement à celles des documents possède des trous
//quelle dimension ? automatique ou pré-définie ?

// réunion jeudi 14h 03/07

#include "wavobser_class_version.h"
#include "reading.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <climits>
#include <ctime>
#include <cstdio>

#define Pi 3.14159265
#define tnum 2000000
#define lx 1
#define ly 1
#define q1 2.0*Pi/nu2
#define o2 10000
#define h 1

WaveSimulation::WaveSimulation()
    :       omega(0), omega1(0), omega2(0), OMEGA(1), A1(0.0), B1(0.0), B2(0.5), A2(0.5), qx(0), qy(0), a(lx / 20.0),
      eta(1.0), friction(0), theta(0.0), kay(0.0), R0(0.2), we(1)
{
    friction = 4.0 * Pi * eta * a;
}

double WaveSimulation::force(double p) {
    return (1.0 + A1 * cos(p) + A2 * cos(2.0 * p) + B1 * sin(p) + B2 * sin(2.0 * p));
}

double GreenFunctions::G11(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = 1.0 / d1 + deltax * deltax / pow(d1, 3) - 1.0 / d2 - deltax * deltax / pow(d2, 3) +
        2.0 * z1 * z2 * (3.0 * deltax * deltax / pow(d2, 5) - 1.0 / pow(d2, 3));
    return (G);
}

double GreenFunctions::G12(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltax*deltay / pow(d1, 3) - deltax*deltay / pow(d2, 3) + 2.0*z1*z2*3.0*deltax*deltay / pow(d2, 5);

    return (G);
}

double GreenFunctions::G13(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltax*deltaz_1 / pow(d1, 3) - deltax*deltaz_2 / pow(d2, 3) +
        2.0*z2*deltax / pow(d2, 3) - 2.0*3.0*z1*z2 * (deltax*deltaz_2 / pow(d2, 5));

    return (G);
}

double GreenFunctions::G21(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltax*deltay / pow(d1, 3) - deltax*deltay / pow(d2, 3) + 2.0*z1*z2*3.0*deltax*deltay / pow(d2, 5);

    return (G);
}

double GreenFunctions::G22(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = 1.0 / d1 + deltay*deltay / pow(d1, 3) - 1.0 / d2 - deltay*deltay / pow(d2, 3) +
        2.0*z1*z2*(3.0*deltay*deltay / pow(d2, 5) - 1.0 / pow(d2, 3));

    return (G);
}

double GreenFunctions::G23(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltay*deltaz_1 / pow(d1, 3) - deltay*deltaz_2 / pow(d2, 3) +
        2.0*z2*deltay / pow(d2, 3) - 2.0*3.0*z1*z2 * (deltay*deltaz_2 / pow(d2, 5));

    return (G);
}

double GreenFunctions::G31(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltax*deltaz_1 / pow(d1, 3) - deltax*deltaz_2 / pow(d2, 3) +
        2.0*z2*deltax / pow(d2, 3) + 2.0*3.0*z1*z2 * (deltax*deltaz_2 / pow(d2, 5));

    return (G);
}

double GreenFunctions::G32(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = deltay*deltaz_1 / pow(d1, 3) - deltay*deltaz_2 / pow(d2, 3) +
        2.0*z2*deltay / pow(d2, 3) + 2.0*3.0*z1*z2 * (deltay*deltaz_2 / pow(d2, 5));

    return (G);
}

double GreenFunctions::G33(double deltax, double deltay, double deltaz_1, double deltaz_2, double d1, double d2, double z1, double z2) {
    double G = 1.0 / d1 + deltaz_1 * deltaz_1 / pow(d1, 3) - 1.0 / d2 - deltaz_2 * deltaz_2 / pow(d2, 3) +
        2.0 * z1 * z2*(-3.0*deltaz_2* deltaz_2 / pow(d2, 5) + 1.0 / pow(d2, 3));

    return (G);
}

// Main simulation logic
void WaveSimulation::run(int limit) {
    Reader reader;
    if (!reader.read_csv("selected_points.csv")) {
        std::cerr << "Failed to read file.\n";
        return;
    }

    // Find bounds
    int min_x = INT_MAX, max_x = INT_MIN, min_y = INT_MAX, max_y = INT_MIN;
    for (const auto& [x, y] : reader.coords) {
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    std::cout << "un max x de " << max_x << " et un max y de " << max_y << std::endl;

    int num_x = max_x - min_x + 1;
    int num_y = max_y - min_y + 1;

    std::cout << "Grid size: " << num_x << " x " << num_y << std::endl;

    // Allocate simulation arrays dynamically
    std::vector<std::vector<double>> phi(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> ph2(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<int>> active(num_x, std::vector<int>(num_y, 0));
    std::vector<std::vector<double>> t1x(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> t1y(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> t1z(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> t2x(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> t2y(num_x, std::vector<double>(num_y, 0.0));
    std::vector<std::vector<double>> t2z(num_x, std::vector<double>(num_y, 0.0));

    // Fill active matrix using coordinates
    for (const auto& [x, y] : reader.coords) {
        int xi = x - min_x;
        int yi = y - min_y;
        active[xi][yi] = 1;
    }

    // Count active elements
    size_t n_active = 0;
    for (int i = 0; i < num_x; ++i)
        for (int j = 0; j < num_y; ++j)
            if (active[i][j]) ++n_active;

    std::cout << "We found " << n_active << " active elements." << std::endl;

    int i, j, l, m, o;

    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time); // Record start time

    g2[0] = -cos(theta); g2[1] = -sin(theta); g2[2] = cos(kay);
    
    h2[0] = sin(kay)*sin(theta); h2[1] = -sin(kay)*cos(theta);

    double deltay, deltaz_1, deltaz_2, d1, d2;

    double t, dt = 0.002, x1, y1, z1, x2, y2, z2, deltaphi; char filename[256] = { 0 };

    for (i = 0; i < num_x; i++)
    {
        for (j = 0; j < num_y; j++)
        {
            if (!active[i][j]) continue; // CUSTOM
            phi[i][j] = 2.0 * Pi*rand() / (RAND_MAX + 1.0); 

            //phi[i][j] = 0.0;
        }
    }

    FILE *fp = fopen("o.dat", "w");

    for (l = 1; l <= tnum; l++)
    {
        t = l*dt;

        if (l==1)
        {
            sprintf(filename, "%d.dat", l-1);

            fp = fopen(filename, "w");

        //custom
        fprintf(fp, "# A1 = %lf\n", A1);
        fprintf(fp, "# B1 = %lf\n", B1);
        fprintf(fp, "# A2 = %lf\n", A2);
        fprintf(fp, "# B2 = %lf\n\n", B2);
        //custom

            
            for (i = 0; i < num_x; i++)
            {
                int wrote = 0;
                for (j = 0; j < num_y; j++)
                {
                    if (!active[i][j]) continue; // CUSTOM
                    fprintf(fp, "%lf %d %d %lf %lf\n", t, i, j, phi[i][j], phi[i][j]-A2/2.0*sin(2.0*phi[i][j])+B2/2.0*cos(2.0*phi[i][j]));
                    wrote = 1;
                }

                if (wrote) fprintf(fp, "\n");
            }

            fclose(fp);
        }

        for (o = 1; o <= o2; o++)
        {
            if (l == tnum / o2 * o)
            {
                sprintf(filename, "%d.dat", o);

                fp = fopen(filename, "w");
                // CUSTOM: Print elapsed time when limit th file is written
                if (o == limit +1 ) {
                    clock_gettime(CLOCK_MONOTONIC, &end_time);
                    double elapsed = (end_time.tv_sec - start_time.tv_sec) +
                                     (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
                    printf("Time to generate %d files: %.3f seconds\n",limit, elapsed);
                    return; // Stop after limit files
                }
                // END CUSTOM
            }
        }
        
        for (i = 0; i < num_x; i++)
        {
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                OMEGA = this->force(phi[i][j]) / friction / R0;

                t1x[i][j] = g2[0] * sin(phi[i][j]) + h2[0] * cos(phi[i][j]);

                t1y[i][j] = g2[1] * sin(phi[i][j]) + h2[1] * cos(phi[i][j]);

                t1z[i][j] = g2[2] * cos(phi[i][j]);

                for (m = 0; m < num_x; m++)
                {
                    double deltax2[num_y], deltay2[num_y], deltaz_12[num_y], deltaz_23[num_y], conti[num_y], d12[num_y], d23[num_y];

#pragma omp parallel for

                    for (int m2 = 0; m2 < num_y; m2++)
                    {
                        if (!active[m][m2]) continue; // CUSTOM
                        if(i == m && j ==m2)
                        {
                            conti[m2] = 0.0;

                            continue;
                        }

                        t2x[m][m2] = g2[0] * sin(phi[m][m2]) + h2[0] * cos(phi[m][m2]);

                        t2y[m][m2] = g2[1] * sin(phi[m][m2]) + h2[1] * cos(phi[m][m2]);

                        t2z[m][m2] = g2[2] * cos(phi[m][m2]);

                        deltay2[m2] = ly*j - ly*m2 
                                            + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                                                    R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));

                        deltax2[m2] = lx * i - lx * m
                                + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay)) -
                                     R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));

                        if (i - m > (num_x - 1) / 2.0 + 0.1)
                        {
                            deltax2[m2] = lx * i - lx * (num_x + m) 
                                            + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay))-
                                                R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));
                        }

                        else if (i - m < -(num_x - 1) / 2.0 - 0.1)
                        {
                            deltax2[m2] = lx * i - lx * (-num_x + m)
                                    + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay)) -
                                        R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));
                        }

                        if (j - m2 > (num_y - 1) / 2.0 + 0.1)
                        {
                            deltay2[m2] = ly*j - ly*(num_y+m2) 
                                            + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                                                    R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));
                        }

                        else if (j - m2 < - (num_y - 1) / 2.0 - 0.1)
                        {
                            deltay2[m2] = ly*j - ly*(-num_y+m2)
                                            + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                                                    R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));
                	    }


                        deltaz_12[m2] = R0 * (cos(kay) * sin(phi[i][j]) - cos(kay) *sin(phi[m][m2]));

                        deltaz_23[m2] = R0 * (cos(kay) *sin(phi[i][j]) + cos(kay) *sin(phi[m][m2])) + 2.0 * h;

                        d12[m2] = sqrt(deltax2[m2] * deltax2[m2] + deltay2[m2] * deltay2[m2] + deltaz_12[m2] * deltaz_12[m2]);

                        d23[m2] = sqrt(deltax2[m2] * deltax2[m2] + deltay2[m2] * deltay2[m2] + deltaz_23[m2] * deltaz_23[m2]);

                        conti[m2] = GreenFunctions::G11(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G12(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G13(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2z[m][m2];

                        conti[m2] += GreenFunctions::G21(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G22(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G23(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2z[m][m2];

                        conti[m2] += GreenFunctions::G31(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G32(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G33(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2z[m][m2];
                    }
                    
                    for (int m2 = 0; m2 < num_y; m2++)
                    {
                        if (!active[m][m2]) continue; // CUSTOM
                        OMEGA += 1.0 / R0*force(phi[m][m2])*conti[m2] / 8.0 / eta /Pi;
                    }
                }

                deltaphi = OMEGA*dt;

                ph2[i][j] = phi[i][j] + deltaphi;
            }
        }
        
        double p2[num_x][num_y]; 

        for (i = 0; i < num_x; i++)
                {
                        for (j = 0; j < num_y; j++)
                        {
                                if (!active[i][j]) continue; // CUSTOM
                                p2[i][j] = ph2[i][j] -A2/2.0*sin(2.0*phi[i][j])+B2/2.0*cos(2.0*phi[i][j]);

                                while (p2[i][j] >= 2.0 * Pi)
                                {
                                        p2[i][j] -= 2.0 * Pi;
                                }

                                while (p2[i][j] < 0)
                                {
                                        p2[i][j] += 2.0 * Pi;
                                }
                        }
                }

        for (i = 0; i < num_x; i++)
        {
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                phi[i][j] = ph2[i][j]; 

                while (phi[i][j] >= 2.0 * Pi)
                {
                    phi[i][j] -= 2.0 * Pi; 
                }

                while (phi[i][j] < 0)
                {
                    phi[i][j] += 2.0 * Pi;
                }
            }
        }


//		if (l > tnum * 0.90)
        {
            for (o = 1; o <= o2; o++)
            {
                if (l == tnum / o2 * o)
                {
                        for (i = 0; i < num_x; i++)
                        {
                            int wrote = 0;
                            for (j = 0; j < num_y; j++)
                            {
                                if (!active[i][j]) continue; // CUSTOM
                                //printf("%lf %d %d %lf\n", t, i, j, phi[i][j]); 

                                fprintf(fp, "%lf %d %d %lf %lf\n", t, i, j, phi[i][j], p2[i][j]);
                                wrote = 1;
                            }
                            
                            if (wrote) fprintf(fp, "\n");

                        }

                        fclose(fp); 
                }
            }
        }
    }

    for (l = 1; l <= o2; l++)
    {
        t = (tnum+l) * dt;
        sprintf(filename, "%d.dat", l + o2);
    
        fp = fopen(filename, "w");

        for (i = 0; i < num_x; i++)
        {
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                OMEGA = this->force(phi[i][j]) / friction / R0;

                t1x[i][j] = g2[0] * sin(phi[i][j]) + h2[0] * cos(phi[i][j]);

                t1y[i][j] = g2[1] * sin(phi[i][j]) + h2[1] * cos(phi[i][j]);

                t1z[i][j] = g2[2] * cos(phi[i][j]);

                for (m = 0; m < num_x; m++)
                {
                    double deltax2[num_y], deltay2[num_y], deltaz_12[num_y], deltaz_23[num_y], conti[num_y], d12[num_y], d23[num_y];

#pragma omp parallel for

                    for (int m2 = 0; m2 < num_y; m2++)
                    {
                        if (!active[m][m2]) continue; // CUSTOM

                        if (i == m && j == m2)
                        {
                            conti[m2] = 0.0;

                            continue;
                        }

                        t2x[m][m2] = g2[0] * sin(phi[m][m2]) + h2[0] * cos(phi[m][m2]);

                        t2y[m][m2] = g2[1] * sin(phi[m][m2]) + h2[1] * cos(phi[m][m2]);

                        t2z[m][m2] = g2[2] * cos(phi[m][m2]);

                        deltay2[m2] = ly * j - ly * m2
                            + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                            R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));

                        deltax2[m2] = lx * i - lx * m
                            + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay)) -
                            R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));

                        if (i - m > (num_x - 1) / 2.0 + 0.1)
                        {
                            deltax2[m2] = lx * i - lx * (num_x + m)
                                + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay)) -
                                R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));
                        }

                        else if (i - m < -(num_x - 1) / 2.0 - 0.1)
                        {
                            deltax2[m2] = lx * i - lx * (-num_x + m)
                                + R0 * (cos(phi[i][j]) * cos(theta) + sin(phi[i][j]) * sin(theta) * sin(kay)) -
                                R0 * (cos(phi[m][m2]) * cos(theta) + sin(phi[m][m2]) * sin(theta) * sin(kay));
                        }

                        if (j - m2 > (num_y - 1) / 2.0 + 0.1)
                        {
                            deltay2[m2] = ly * j - ly * (num_y + m2)
                                + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                                R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));
                        }

                        else if (j - m2 < -(num_y - 1) / 2.0 - 0.1)
                        {
                            deltay2[m2] = ly * j - ly * (-num_y + m2)
                                + R0 * (cos(phi[i][j]) * sin(theta) - sin(phi[i][j]) * cos(theta) * sin(kay)) -
                                R0 * (cos(phi[m][m2]) * sin(theta) - sin(phi[m][m2]) * cos(theta) * sin(kay));
                        }


                        deltaz_12[m2] = R0 * (cos(kay) * sin(phi[i][j]) - cos(kay) *sin(phi[m][m2]));

                        deltaz_23[m2] = R0 * (cos(kay) *sin(phi[i][j]) + cos(kay) *sin(phi[m][m2])) + 2.0 * h;

                        d12[m2] = sqrt(deltax2[m2] * deltax2[m2] + deltay2[m2] * deltay2[m2] + deltaz_12[m2] * deltaz_12[m2]);

                        d23[m2] = sqrt(deltax2[m2] * deltax2[m2] + deltay2[m2] * deltay2[m2] + deltaz_23[m2] * deltaz_23[m2]);

                        conti[m2] = GreenFunctions::G11(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G12(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G13(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1x[i][j] * t2z[m][m2];

                        conti[m2] += GreenFunctions::G21(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G22(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G23(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1y[i][j] * t2z[m][m2];

                        conti[m2] += GreenFunctions::G31(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2x[m][m2];

                        conti[m2] += GreenFunctions::G32(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2y[m][m2];

                        conti[m2] += GreenFunctions::G33(deltax2[m2], deltay2[m2], deltaz_12[m2], deltaz_23[m2], d12[m2], d23[m2], R0*sin(phi[i][j]) + h, R0*sin(phi[m][m2]) + h)*t1z[i][j] * t2z[m][m2];
                    }

                    for (int m2 = 0; m2 < num_y; m2++)
                    {
                        if (!active[m][m2]) continue; // CUSTOM
                        OMEGA += 1.0 / R0 * force(phi[m][m2])*conti[m2] / 8.0 / eta / Pi;
                    }
                }

                deltaphi = OMEGA * dt;

                ph2[i][j] = phi[i][j] + deltaphi;
            }
        }

        double p2[num_x][num_y];

        for (i = 0; i < num_x; i++)
        {
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                p2[i][j] = ph2[i][j] - A2 / 2.0*sin(2.0*phi[i][j]) + B2 / 2.0*cos(2.0*phi[i][j]);

                while (p2[i][j] >= 2.0 * Pi)
                {
                    p2[i][j] -= 2.0 * Pi;
                }

                while (p2[i][j] < 0)
                {
                    p2[i][j] += 2.0 * Pi;
                }
            }
        }

        for (i = 0; i < num_x; i++)
        {
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                phi[i][j] = ph2[i][j];

                while (phi[i][j] >= 2.0 * Pi)
                {
                    phi[i][j] -= 2.0 * Pi;
                }

                while (phi[i][j] < 0)
                {
                    phi[i][j] += 2.0 * Pi;
                }
            }
        }

        for (i = 0; i < num_x; i++)
        {
            int wrote  = 0;
            for (j = 0; j < num_y; j++)
            {
                if (!active[i][j]) continue; // CUSTOM
                //printf("%lf %d %d %lf\n", t, i, j, phi[i][j]); 

                fprintf(fp, "%lf %d %d %lf %lf\n", t, i, j, phi[i][j], p2[i][j]);
                wrote = 1;
            }

             if (wrote) fprintf(fp, "\n");
        }

        fclose(fp);
    }
}

void WaveSimulation::initialize_arrays(int num_x, int num_y) {
   
}

void WaveSimulation::fill_active(const std::vector<std::pair<int, int>>& coords, int min_x, int min_y) {
}

void WaveSimulation::write_output(FILE* fp, double t, const std::vector<std::vector<double>>& phi, const std::vector<std::vector<double>>& p2, const std::vector<std::vector<int>>& active) {
    for (int i = 0; i < phi.size(); ++i) {
        int wrote = 0;
        for (int j = 0; j < phi[0].size(); ++j) {
            if (!active[i][j]) continue;
            fprintf(fp, "%lf %d %d %lf %lf\n", t, i, j, phi[i][j], p2[i][j]);
            wrote = 1;
        }
        if (wrote) fprintf(fp, "\n");
    }
}

// we need to stop generating files when the 10th file is written
// and print the elapsed time for generating these files but we already do that


