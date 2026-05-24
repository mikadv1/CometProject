#include "CometModel.h"
#include "EphemerisWrapper.h"
#include "Constants.h"
#include <cmath>

// Функция активности g(r) для CO2
double marsden_g(double r_au) {
    return pow(r_au, -2.0);
}

// Вычисление негравитационного ускорения в барицентрической системе
Vector3 computeNGAcceleration(const StateVector& state, double A1, double A2, double A3) {
    double r = state.r.norm();

    Vector3 e_R = state.r / r;
    Vector3 h_vec = state.r.cross(state.v);
    double h = h_vec.norm();

    Vector3 e_N = h_vec / h;
    Vector3 e_T = e_N.cross(e_R);
    double g_val = marsden_g(r);

    return (e_R * A1 + e_T * A2 + e_N * A3) * g_val;
}

StateVector cometDerivatives(double t, const StateVector& state, const NGVector& ng) {
    StateVector result;
    result.r = state.v;

    // 1. Гравитация Солнца
    StateVector sun = getBodyState(10, t);  // Положение Солнца относительно барицентра
    Vector3 r_rel = state.r - sun.r;
    double r = r_rel.norm();
    double a_grav = -GM_SUN / (r * r * r);
    result.v = r_rel * a_grav;

    // 2. Учет квадрупольного момента Солнца
    Vector3 p(SUN_POLE_X, SUN_POLE_Y, SUN_POLE_Z);
    Vector3 r_hat = r_rel / r;
    double r_dot_p = r_hat * p;
    double coeff = 3.0 * J2_SUN * GM_SUN * R_SUN * R_SUN / (r * r * r * r);
    Vector3 a_obl = (r_hat * (2.5 * r_dot_p * r_dot_p - 0.5) - p * r_dot_p) * coeff;
    //result.v = result.v + a_obl;

    // 3. Возмущения от планет  
    for (int i = 0; i < NUM_OBJECTS; i++) {
        StateVector body_state = getBodyState(OBJECTS[i].code, t);
        Vector3 r_rel = state.r - body_state.r;
        double d = r_rel.norm();

        double a_perturb = -OBJECTS[i].gm / (d * d * d);
        result.v = result.v + r_rel * a_perturb;
    }

    // 4. Негравитационные эффекты 
    Vector3 a_ng = computeNGAcceleration(state - sun, ng.A1, ng.A2, ng.A3);
    result.v = result.v + a_ng;

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