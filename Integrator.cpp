#include "Integrator.h"
#include <cmath>
#include <cstdio>

StateVector rk4Step(double t, const StateVector& state, double dt, DerivFunc derivs) {
    StateVector k1 = derivs(t, state);

    StateVector temp1 = state + k1 * (dt / 2.0);
    StateVector k2 = derivs(t + dt / 2.0, temp1);

    StateVector temp2 = state + k2 * (dt / 2.0);
    StateVector k3 = derivs(t + dt / 2.0, temp2);

    StateVector temp3 = state + k3 * dt;
    StateVector k4 = derivs(t + dt, temp3);

    return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

void integrate(double t0, double tend, double dt, StateVector& state, DerivFunc derivs, const char* outputFile) {
    FILE* file = nullptr;
    if (outputFile != nullptr) {
        file = fopen(outputFile, "w");
        if (file != nullptr) {
            // Заголовок CSV
            fprintf(file, "JD,X,Y,Z,VX,VY,VZ\n");
        }
    }

    double t = t0;
    //double step = fmin(t0 - t, dt);
    while (t < tend) {
        // Запись в файл
        if (file != nullptr) {
            fprintf(file, "%.16le,%.16le,%.16le,%.16le,%.16le,%.16le,%.16le\n",
                t, state.r.x, state.r.y, state.r.z,
                state.v.x, state.v.y, state.v.z);
        }

        //step = fmin(t0 - t, dt);
        state = rk4Step(t, state, dt, derivs);
        t += dt;
    }

    // Запись последней точки
    if (file != nullptr) {
        double r = state.r.norm();
        fprintf(file, "%.16le,%.16le,%.16le,%.16le,%.16le,%.16le,%.16le\n",
            t, state.r.x, state.r.y, state.r.z,
            state.v.x, state.v.y, state.v.z);
        fclose(file);
    }
}