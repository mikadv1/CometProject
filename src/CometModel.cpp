#include "CometModel.h"
#include "EphemerisWrapper.h"
#include "Constants.h"
#include <cmath>
#include <cstdio>

StateVector cometDerivatives(double t, const StateVector& state) {
    StateVector result;
    result.r = state.v;

    // 1. Гравитация Солнца
    double r = state.r.norm();
    double a_grav = -GM_SUN / (r * r * r);
    result.v = state.r * a_grav; 

    // 2. Учет квадрупольного момента Солнца
    Vector3 p(SUN_POLE_X, SUN_POLE_Y, SUN_POLE_Z);
    Vector3 r_hat = state.r * (1 / r);
    double r_dot_p = r_hat * p;
    double coeff = 3 * J2_SUN * GM_SUN * R_SUN * R_SUN / (r * r * r * r);
    Vector3 a_obl = (r_hat * (2.5 * r_dot_p * r_dot_p - 0.5) - p * r_dot_p) * coeff;
    result.v = result.v + a_obl;

    // 3. Возмущения от планет  
    for (int i = 0; i < NUM_OBJECTS; i++) {
        StateVector body_state = getBodyState(OBJECTS[i].code, t);
        Vector3 r_rel = state.r - body_state.r;
        double d = r_rel.norm();

        double a_perturb = -OBJECTS[i].gm / (d * d * d);
        result.v = result.v + r_rel * a_perturb;
    }

    // 4. Негравитационные эффекты (будет добавлено позже)
    // ...

    return result;
}

bool initCometModel(const char* ephFile) {
    if (!initEphemeris(ephFile)) {
        return false;
    }
    return true;
}

void cleanupCometModel() {
    cleanupEphemeris();
}