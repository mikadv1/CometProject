#pragma once
#include "Types.h"
#include "Integrator.h"

struct ReductionResult {
    double jd_utc;
    double jd_tdb;
    double ra_obs;
    double dec_obs;
    double ra_model;
    double dec_model;
    double dRA;
    double dDec;
    int obsCode;
};

void reduceObservation(
    double jd_utc,
    double ra_obs,
    double dec_obs,
    const Vector3& r_station_itrs,
    const Trajectory& traj,
    ReductionResult& result
);

double utc2tdb(double jd_utc);

Vector3 stationITRS2GCRS(double jd_utc, const Vector3& r_itrs);

double solveLightTime(
    double jd_obs_tdb,
    const Vector3& r_obs,
    const Trajectory& traj
);

Vector3 applyGravDeflection(const Vector3& rho, const Vector3& r_obs, const Vector3& r_comet);

Vector3 applyAberration(const Vector3& rho, const Vector3& v_earth);

void cartesianToSpherical(const Vector3& r, double& ra, double& dec);