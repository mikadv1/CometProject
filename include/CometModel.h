#pragma once
#include "Types.h"

double marsden_g(double r_au);

Vector3 computeNGAcceleration(const StateVector& state);

// Правая часть системы ОДУ для кометы
StateVector cometDerivatives(double t, const StateVector& state);

// Инициализация модели (загрузка эфемерид планет)
bool initCometModel(const char* ephFile);

// Очистка ресурсов
void cleanupCometModel();

void setNGParameters(double a1, double a2, double a3);

void getNGParameters(double& a1, double& a2, double& a3);