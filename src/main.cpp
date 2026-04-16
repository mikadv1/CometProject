#include <cstdio>
#include <cstring>
#include <cmath>
#include <windows.h>
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"
#include "FileIO.h"


int main() {
    SetConsoleOutputCP(65001);

    if (!initCometModel("C:\\diploma\\data\\epm2021.bsp")) return 1;

    // Состояние в момент времени t0
    StateVector state0;
    state0.r = { 1.057725183264539E+00, -5.962453673239701E+00, -2.179617339224094E+00 };
    state0.v = { -1.377811929112319E-02, 2.988753120767556E-02, 1.139417199619259E-02 };

    setNGParameters(5.320206165314E-8, 1.148166060448E-8, -6.854491829872E-9);

    double t0 = 2460800.5, tend = t0 + 350.0, grid_dt = 0.125, internal_dt = 0.0125;

    // Интегрирование
    Trajectory traj;
    integrate(t0, tend, grid_dt, internal_dt, state0, cometDerivatives, traj);

    // Обработка наблюдений
    processObservations(
        "C:\\diploma\\data\\obs.csv",
        "C:\\diploma\\data\\residuals.csv",
        traj
    );

    traj.clear();
    cleanupCometModel();
    return 0;
}