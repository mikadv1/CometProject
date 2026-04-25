#include <cstdio>
#include <cmath>
#include <clocale>
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"
#include "FileIO.h"
#include "GaussNewton.h"

int main() {
    //std::setlocale(LC_ALL, "Russian");

    if (!initCometModel("C:\\diploma\\data\\epm2021.bsp")) return 1;

    // Состояние в момент времени t0
    StateVector state0;
    state0.r = { 1.057725183264539E+00, -5.962453673239701E+00, -2.179617339224094E+00 };
    state0.v = { -1.377811929112319E-02, 2.988753120767556E-02, 1.139417199619259E-02 };
    //NGVector ng0 = { 5.320206165314E-8, 1.148166060448E-8, -6.854491829872E-9 };
    NGVector ng0 = { 0.0, 0.0, 0.0 };

    OrbitParams init_params(state0.r, state0.v, ng0);

    double t0 = 2460800.5, tend = t0 + 350.0, dt = 0.0125;

    //Trajectory traj;
    //integrate(t0, tend, dt * 10, dt, state0, ng0, traj);
    //processObservations("C:\\diploma\\data\\obs.csv", "C:\\diploma\\data\\residuals.csv", traj);

    std::vector<Observation> obs_vector;
    obs_vector.reserve(8500);
    loadObservations("C:\\diploma\\data\\obs.csv", obs_vector);

    std::array<bool, 9> selected = { 1, 1, 1, 1, 1, 1, 0, 0, 0 };
    OrbitParams estimated = fitOrbit(init_params, obs_vector, t0, tend, dt, selected, 1.0e-10);

    //integrate(t0, tend, dt * 10, dt, { estimated.r0, estimated.v0 }, estimated.ng0, traj);
    //processObservations("C:\\diploma\\data\\obs.csv", "C:\\diploma\\data\\residuals_after.csv", traj);
   
    cleanupCometModel();
    return 0;
}