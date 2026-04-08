#pragma once
#include "Types.h"

// Инициализация эфемерид
bool initEphemeris(const char* filePath);

// Получить положение и скорость тела на момент JD
StateVector getBodyState(int bodyCode, double jd);

// Разница шкал TT и TDB
double TTmTDB(double jd_tdb);

// Очистка
void cleanupEphemeris();