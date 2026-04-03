#pragma once
#include "Types.h"

struct Trajectory {
    double* t;          // Массив времени (JD), размер nPoints
    StateVector* state; // Массив состояний, размер nPoints
    int nPoints;        // Количество точек
    double dt;          // Шаг сетки (дни)
    double t0;          // Начальное время
    double tend;        // Конечное время

    Trajectory();
    ~Trajectory();

    // Выделение памяти под траекторию
    void allocate(int n, double start_t, double step);

    // Освобождение памяти
    void clear();
};

// Тип функции для правой части ОДУ
using DerivFunc = StateVector(*)(double t, const StateVector& state);

// Заполняет траекторию на сетке [t0, tend] с шагом grid_dt
// internal_dt - внутренний шаг интегрирования
void integrate(
    double t0,
    double tend,
    double grid_dt,      
    double internal_dt,  
    const StateVector& state0,
    DerivFunc derivs,
    Trajectory& traj
);

// Линейная интерполяция состояния на произвольный момент времени
StateVector interpolateLinear(const Trajectory& traj, double t);

bool isTimeInRange(const Trajectory& traj, double t);

// Шаг метода Рунге-Кутты 4 порядка
static StateVector rk4Step(double t, const StateVector& state, double dt, DerivFunc derivs);