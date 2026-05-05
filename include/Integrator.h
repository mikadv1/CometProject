#pragma once
#include "Types.h"
#include "CometModel.h"

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

// Заполняет траекторию на сетке [t0, tend] с шагом grid_dt
// internal_dt - внутренний шаг интегрирования
void integrate(
    double t0,
    double tend,
    double grid_dt,      
    double internal_dt,  
    const StateVector& state0,
    const NGVector& ng,
    Trajectory& traj
);

// Линейная интерполяция состояния на произвольный момент времени
StateVector interpolateLinear(const Trajectory& traj, double t);

bool isTimeInRange(const Trajectory& traj, double t);

// Шаг метода Рунге-Кутты 4 порядка
template<typename StateType, typename derivFunc>
StateType rk4Step(double t, const StateType& state, double dt, const NGVector& ng, derivFunc derivs) {
    StateType k1 = derivs(t, state, ng);
    StateType k2 = derivs(t + dt / 2.0, state + k1 * (dt / 2.0), ng);
    StateType k3 = derivs(t + dt / 2.0, state + k2 * (dt / 2.0), ng);
    StateType k4 = derivs(t + dt, state + k3 * dt, ng);

    return state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}