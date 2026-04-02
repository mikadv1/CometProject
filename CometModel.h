#pragma once
#include "Types.h"

// Правая часть системы ОДУ для кометы
StateVector cometDerivatives(double t, const StateVector& state);

// Инициализация модели (загрузка эфемерид планет)
bool initCometModel(const char* ephFile);

// Очистка ресурсов
void cleanupCometModel();