#pragma once
#include "Types.h"
#include "Integrator.h"

struct ReductionResult {
    double jd_utc;
    double jd_tdb;
    double dRA;
    double dDec;
};

struct Observation {
    double jd_utc;
    double ra;
    double dec;
    Vector3 r_station_itrs;
    char type;
};

void reduceObservation(
    const Observation& obs,
    const Trajectory& traj,
    ReductionResult& result
);

double utc2tdb(double jd_utc);

Vector3 stationITRS2GCRS(double jd_utc, double jd_tt, const Vector3& r_itrs);

double solveLightTime(
    double jd_obs_tdb,
    const Vector3& r_obs,
    const Trajectory& traj
);

Vector3 Deflection(const Vector3& rho, const Vector3& r_obs);

void cartesianToSphericalArcsec(const Vector3& r, double& ra, double& dec);