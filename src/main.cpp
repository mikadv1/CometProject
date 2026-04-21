#include <cstdio>
#include <cstring>
#include <cmath>
#include <clocale>
#include <windows.h>
#include "Types.h"
#include "Integrator.h"
#include "CometModel.h"
#include "FileIO.h"
#include "GaussNewton.h"

void testSolveCholesky() {
    // 1. Создаём симметричную положительно определённую матрицу A (9x9)
    Matrix9x9 A{};
    for (int i = 0; i < 9; i++) {
        A[i][i] = 10.0 + i; 
        for (int j = 0; j < 9; j++) {
            if (i != j) {
                A[i][j] = 0.5 * (i + j + 1);  
            }
        }
    }

    // 2. Задаём известный вектор решения x_true
    ParamVector x_true = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

    std::array<bool, 9> selected = { true, true, true, false, false, false, true, true, true };

    ParamVector b_partial{ 113.5, 138., 162.5, 248., 113., 130., 260.5, 285., 309.5 };

    ParamVector x_partial = solveCholesky(A, b_partial, selected);

    printf("\nVector b (A * x_true):\n");
    for (int i = 0; i < 9; i++) {
        printf("  b[%d] = %.4f\n", i, b_partial[i]);
    }
    printf("\n");

    for (int i = 0; i < 9; i++) {
        printf("  x[%d] = %.6f (expected: %.6f)\n", i, x_partial[i], x_true[i]);
    }
}

int main() {
    std::setlocale(LC_ALL, "Russian");

    if (!initCometModel("C:\\diploma\\data\\epm2021.bsp")) return 1;

    // Состояние в момент времени t0
    StateVector state0;
    state0.r = { 1.057725183264539E+00, -5.962453673239701E+00, -2.179617339224094E+00 };
    state0.v = { -1.377811929112319E-02, 2.988753120767556E-02, 1.139417199619259E-02 };

    NGVector ng0 = { 5.320206165314E-8, 1.148166060448E-8, -6.854491829872E-9 };

    double t0 = 2460800.5, tend = t0 + 350.0, grid_dt = 0.125, internal_dt = 0.0125;

    
    testSolveCholesky();
   
    cleanupCometModel();
    return 0;
}