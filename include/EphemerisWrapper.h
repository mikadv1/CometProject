#pragma once
#include "Types.h"

// Инициализация эфемерид
bool initEphemeris(const char* filePath);

// Получить положение планеты на момент JD
Vector3 getPlanetPosition(int bodyCode, double jd);

// Очистка
void cleanupEphemeris();