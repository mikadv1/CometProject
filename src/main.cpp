#include <cstdio>
#include <windows.h>
#include "Constants.h"
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"
#include "Reduction.h"

int main() {
    SetConsoleOutputCP(65001);

    // 1. Инициализация модели
    if (!initCometModel("C:\\diploma\\data\\epm2021.bsp")) {
        printf("[ERROR] Failed to initialize comet model!\n");
        return 1;
    }

    // 2. Начальные условия (из JPL Horizons)
    StateVector state0;
    state0.r = { 2.704338288122531E-01, -4.501623714812790E+00, 2.893195749808039E-01 };  // а.е.
    state0.v = { -1.384389231863844E-02, 3.253008104441191E-02, -1.470083979858731E-03 };  // а.е./день

    double t0 = 2460857.5;  // JD
    double grid_dt = 1.0e-1;
    double internal_dt = 1.0e-2;  // шаг интегрирования 
    double tend = t0 + 300.0; 

    Trajectory traj;
    integrate(t0, tend, grid_dt, internal_dt, state0, cometDerivatives, traj);

    Vector3 r_station(174.768, -4926.470, 4033.989);
    ReductionResult result;

    reduceObservation(
        2460858.728550,
        976840.200 / ARCSEC_PER_RAD,
        -67259.500 / ARCSEC_PER_RAD,
        r_station,
        traj,
        result
    );

    printf("JD UT:  %.8f d\n", result.jd_utc);
    printf("JD TDB: %.8f d\n", result.jd_tdb);
    printf("dRA:  %.3f arcsec\n", result.dRA * ARCSEC_PER_RAD);
    printf("dDec: %.3f arcsec\n", result.dDec * ARCSEC_PER_RAD);
    printf("RA obs:    %.3f arcsec\n", result.ra_obs * ARCSEC_PER_RAD);
    printf("RA model:  %.3f arcsec\n", result.ra_model * ARCSEC_PER_RAD);
    printf("Dec obs:   %.3f arcsec\n", result.dec_obs * ARCSEC_PER_RAD);
    printf("Dec model: %.3f arcsec\n", result.dec_model * ARCSEC_PER_RAD);


    // 6. Очистка
    traj.clear();

    cleanupCometModel();

    return 0;
}