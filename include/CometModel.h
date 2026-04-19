#pragma once
#include "Types.h"

double marsden_g(double r_au);

Vector3 computeNGAcceleration(const StateVector& state, double A1, double A2, double A3);

// Правая часть системы ОДУ для кометы
StateVector cometDerivatives(double t, const StateVector& state, const NGVector& ng);

// Инициализация модели (загрузка эфемерид планет)
bool initCometModel(const char* ephFile);

// Очистка ресурсов
void cleanupCometModel();