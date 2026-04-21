#pragma once
#include <vector>
#include "Types.h"
#include "Integrator.h"
#include "Reduction.h"

// Загрузка траектории из CSV (формат JPL Horizons)
bool loadTrajectoryFromCSV(const char* filename, Trajectory& traj, double t0, double tend, double dt);

// Обработка наблюдений и запись невязок
bool processObservations(const char* obsFile, const char* outFile, const Trajectory& traj);

bool loadObservations(const char* obsFile, std::vector<Observation>& obs_vector);