#include <cstdio>
#include <windows.h>
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"

int main() {
    SetConsoleOutputCP(65001);

    // 1. Инициализация модели
    if (!initCometModel("C:\\diploma\\ephs\\epm2021.bsp")) {
        printf("[ERROR] Failed to initialize comet model!\n");
        return 1;
    }

    // 2. Начальные условия (из JPL Horizons)
    StateVector state0;
    state0.r = { 2.704338288122531E-01, -4.501623714812790E+00, 2.893195749808039E-01 };  // а.е.
    state0.v = { -1.384389231863844E-02, 3.253008104441191E-02, -1.470083979858731E-03 };  // а.е./день

    double t0 = 2460857.5;  // JD
    double dt = 1.0e-2;  // шаг интегрирования 
    double tend = t0 + 300.0; 

    // 3. Интегрирование
    StateVector state = state0;
    integrate(t0, tend, dt, state, cometDerivatives, "orbit_output.csv");

    // 4. Вывод результата
    printf("Final position: %.16le %.16le %.16le AU\n", state.r.x, state.r.y, state.r.z);

    cleanupCometModel();

    return 0;
}