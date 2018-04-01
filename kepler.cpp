#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <time.h>
using namespace std;

#include "vector.hpp"
#include "odeint.hpp"
using namespace cpl;

const double pi = 4 * atan(1.0);
const double GmPlusM = 4 * pi * pi;
const double MEarth = 5.972; //*e24 kg
const double MSun = 1.989e6;
const double MMercury = 0.3285;
const double G = 6.674e-11;
const double me =0.165e-6 ;
const double mm = 3e-6;
bool switch_t_with_y = false;    //  to interpolate to y = 0

//  Derivative vector for Newton's law of gravitation
Vector f(const Vector& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4], r_x1 = x[5], r_y1 = x[6], v_x1 = x[7], v_y1 = x[8];
    double rSquared = r_x*r_x + r_y*r_y;
    double rSquared1 = r_x1*r_x1 + r_y1*r_y1;
    double rEarthMercury = pow((r_x-r_x1)*(r_x-r_x1)+(r_y-r_y1)*(r_y-r_y1), 1.5);
    double rCubed = rSquared * sqrt(rSquared);
    double rCubed1 = rSquared1 * sqrt(rSquared1);
    
    Vector f(9);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[5] = v_x1;
    f[6] = v_y1;
    f[3] = - GmPlusM * r_x/rCubed - GmPlusM/MSun*MMercury*(r_x-r_x1)/rEarthMercury ;
    f[4] = - GmPlusM * r_y / rCubed - GmPlusM/MSun*MMercury*(r_y-r_y1)/rEarthMercury;
    f[7] =- GmPlusM * r_x1/rCubed1 - GmPlusM/MSun*MEarth*(r_x1-r_x)/rEarthMercury ;
    f[8] = - GmPlusM * r_y1 / rCubed1 - GmPlusM/MSun*MEarth*(r_y1-r_y)/rEarthMercury;

    if (switch_t_with_y) {
        //  use y as independent variable
        for (int i = 0; i < 5; i++)
            f[i] /= v_y;
    }
    return f;
}

//  Change independent variable from t to y and step back to y = 0
void interpolate_crossing(Vector x, int& crossing) {
    ++crossing;
    switch_t_with_y = true;
    RK4Step(x, -x[2], f);
    cout << " crossing " << crossing << "\t t = " << x[0]
         << "\t x = " << x[1] << endl;
    switch_t_with_y = false;
}

int main() {
    cout << " Kepler orbit comparing fixed and adaptive Runge-Kutta\n"
         << " -----------------------------------------------------\n"
         << " Enter aphelion distance in AU, and eccentricity: ";
    double r_ap, eccentricity, a, T, v0, v01, b, eccentricity1, r_ap1;
    r_ap1 = 0.467;
    eccentricity1 = 0.205;
    cin >> r_ap >> eccentricity;
    a = r_ap / (1 + eccentricity);
    b = r_ap1 /(1 + eccentricity1);
    T = pow(a, 1.5);
    v0 = sqrt(GmPlusM * (2 / r_ap - 1 / a));
    v01 =sqrt(GmPlusM * (2 / r_ap1 - 1 / b));
    cout << " Enter number of periods, step size, and adaptive accuracy: ";
    double periods, dt, accuracy;
    cin >> periods >> dt >> accuracy;
    Vector x0(8);
    x0[0] = 0;  x0[1] = r_ap;  x0[2] = 0;  x0[3] = 0;  x0[4] = v0; x0[5] = r_ap1; x0[6] = 0; x0[7] = 0; x0[8] = v01;

    ofstream dataFile("fixed.data");
    Vector x = x0;
    clock_t t1;
    int steps = 0, crossing = 0;
    t1 = clock();
    cout << "\n Integrating with fixed step size" << endl;
    do {
        for (int i = 0; i < 8; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        double y = x[2];
        double y1 = x[6];
        RK4Step(x, dt, f);
        ++steps;
    if (y * x[2] < 0 )
           interpolate_crossing(x, crossing);
    }
    while (x[0] < periods * T);
    cout << " number of fixed size steps = " << steps << endl;
    cout << " data in file fixed.data" << endl;
    dataFile.close();
    t1 = clock()-t1;
    cout << "fixed time: "<< t1<< endl;
    dataFile.open("adaptive.data");
    x = x0;
    steps = crossing = 0;
    t1 = clock();
    double dt_max = 0, dt_min = 100;
    cout << "\n Integrating with adaptive step size" << endl;
    do {
        for (int i = 0; i < 8; i++)
            dataFile << x[i] << '\t';
        dataFile  <<'\n';
        double t_save = x[0];
        double y = x[2];
        adaptiveRK4Step(x, dt, accuracy, f);
        double step_size = x[0] - t_save;
        ++steps;
        if (step_size < dt_min) dt_min = step_size;
        if (step_size > dt_max) dt_max = step_size;
        if (y * x[2] < 0)
            interpolate_crossing(x, crossing);
    } while (x[0] < periods * T);
    t1 = clock()-t1;
    cout << " number of adaptive steps = " << steps << endl;
    cout << " step size: min = " << dt_min << "  max = " << dt_max << endl;
    cout << " data in file adaptive.data" << endl;
    cout << t1<< endl;
    dataFile.close();
}

