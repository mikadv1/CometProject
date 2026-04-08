#include <cstdio>
#include <windows.h>
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"

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

    // 3. Интегрирование
    Trajectory traj;
    integrate(t0, tend, grid_dt, internal_dt, state0, cometDerivatives, traj);

    printf("Integration is done. Points: %d\n", traj.nPoints);

    // 4. Тест интерполяции: запрос положения в произвольный момент
    double query_t = 2460876.5;  // произвольное время
    if (isTimeInRange(traj, query_t)) {
        StateVector interpolated = interpolateLinear(traj, query_t);
        printf("Interpolated at JD %.8f:\n", query_t);
        printf("  Position: (%.16le, %.16le, %.16le) AU\n",
            interpolated.r.x, interpolated.r.y, interpolated.r.z);
    }

    // 5. (Опционально) Вывод траектории в файл — теперь это делает main, а не integrator
    FILE* f = fopen("data\\orbit_output.csv", "w");
    if (f) {
        fprintf(f, "JD,X,Y,Z,VX,VY,VZ\n");
        for (int i = 0; i < traj.nPoints; i++) {
            fprintf(f, "%.8f,%.16le,%.16le,%.16le,%.16le,%.16le,%.16le\n",
                traj.t[i],
                traj.state[i].r.x, traj.state[i].r.y, traj.state[i].r.z,
                traj.state[i].v.x, traj.state[i].v.y, traj.state[i].v.z);
        }
        fclose(f);
    }

    // 6. Очистка
    traj.clear();

    cleanupCometModel();

    return 0;
}