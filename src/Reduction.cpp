#include "Reduction.h"
#include "EphemerisWrapper.h"
#include "Constants.h"
#include "ephaccess.h"
#include <cmath>

extern "C" {
#include "sofa.h"
#include "sofam.h"
}

double utc2tdb(double jd_utc) {
    double tai1, tai2, tt1, tt2, tdb1, tdb2;
    double dtr;

    iauUtctai(jd_utc, 0.0, &tai1, &tai2);
    iauTaitt(tai1, tai2, &tt1, &tt2);
    dtr = TTmTDB(tt1 + tt2);
    iauTttdb(tt1, tt2, dtr, &tdb1, &tdb2);

    return tdb1 + tdb2;
}

Vector3 stationITRF2GCRS(double jd_utc, double jd_tt, const Vector3& r_itrs) {
    double xp = 0.0;
    double yp = 0.0;
    double c2t[3][3];
    double r_in[3], r_out[3];

    iauC2t06a(jd_tt, 0.0, jd_utc, 0.0, xp, yp, c2t);

    r_in[0] = r_itrs.x;
    r_in[1] = r_itrs.y;
    r_in[2] = r_itrs.z;

    iauTrxp(c2t, r_in, r_out);

    return Vector3(r_out[0], r_out[1], r_out[2]);
}

double solveLightTime(
    double jd_obs_tdb,
    const Vector3& r_obs,
    const Trajectory& traj
) {
    double delta = 0.0;

    for (int iter = 0; iter < 20; iter++) {
        double jd_emit = jd_obs_tdb - delta;

        StateVector comet = interpolateLinear(traj, jd_emit);
        double dist = (comet.r - r_obs).norm();
        double new_delta = dist / C_LIGHT_AU_DAY;

        if (fabs(new_delta - delta) < 1e-12) {
            break;
        }

        delta = new_delta;
    }    
    return delta;
}

Vector3 Deflection(const Vector3& rho, const Vector3& r_obs) {
    double p[3], q[3], e[3], pout[3];

    double rho_norm = rho.norm();
    double r_sun = r_obs.norm();

    p[0] = rho.x / rho_norm;
    p[1] = rho.y / rho_norm;
    p[2] = rho.z / rho_norm;

    Vector3 r_comet = rho + r_obs;
    double q_norm = r_comet.norm();
    q[0] = r_comet.x / q_norm;
    q[1] = r_comet.y / q_norm;
    q[2] = r_comet.z / q_norm;

    e[0] = r_obs.x / r_sun;
    e[1] = r_obs.y / r_sun;
    e[2] = r_obs.z / r_sun;

    double bm = 1.0;
    double em = r_sun;
    double dlim = 1e-3;     // ~2.6 degrees

    iauLd(bm, p, q, e, em, dlim, pout);

    return Vector3(pout[0] * rho_norm, pout[1] * rho_norm, pout[2] * rho_norm);
}

void cartesianToSphericalArcsec(const Vector3& r, double& ra, double& dec) {
    double p[3], theta, phi;

    p[0] = r.x;
    p[1] = r.y;
    p[2] = r.z;

    iauC2s(p, &theta, &phi);

    ra = theta * ARCSEC_PER_RAD;
    dec = phi * ARCSEC_PER_RAD;
}

void reduceObservation(
    const Observation& obs,
    const Trajectory& traj,
    ReductionResult& result
) {
    result.jd_utc = obs.jd_utc;

    double jd_tdb = utc2tdb(obs.jd_utc);
    result.jd_tdb = jd_tdb;

    if (obs.type != 'S') {
        Vector3 r_station_gcrs = stationITRF2GCRS(obs.jd_utc, jd_tdb, obs.r_station_itrf);
    }
    Vector3 r_station_gcrs = obs.r_station_itrf;
    r_station_gcrs = r_station_gcrs / AU_KM;

    StateVector earth_state = getBodyState(EPH_EARTH, jd_tdb);

    Vector3 r_obs = earth_state.r + r_station_gcrs;

    double light_time = solveLightTime(jd_tdb, r_obs, traj);
    double jd_emit = jd_tdb - light_time;

    StateVector comet = interpolateLinear(traj, jd_emit);

    Vector3 rho = comet.r - r_obs;

    rho = Deflection(rho, r_obs);

    double ra_model, dec_model;
    cartesianToSphericalArcsec(rho, ra_model, dec_model);

    double dRA = obs.ra - ra_model;
    while (dRA > PI * ARCSEC_PER_RAD) dRA -= PI * 2 * ARCSEC_PER_RAD;
    while (dRA < -PI * ARCSEC_PER_RAD) dRA += PI * 2 * ARCSEC_PER_RAD;

    double dDec = obs.dec - dec_model;

    result.dRA = dRA;
    result.dDec = dDec;

    return;
}