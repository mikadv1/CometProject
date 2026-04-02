#pragma once
#include "Types.h"

// Тип функции для правой части ОДУ
using DerivFunc = StateVector(*)(double t, const StateVector& state);

// Шаг метода Рунге-Кутты 4 порядка
StateVector rk4Step(double t, const StateVector& state, double dt, DerivFunc derivs);

// Интегрирование на интервале [t0, tend]
void integrate(double t0, double tend, double dt, StateVector& state, DerivFunc derivs, const char* outputFile = nullptr);