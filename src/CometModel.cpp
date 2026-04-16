#include "CometModel.h"
#include "EphemerisWrapper.h"
#include "Constants.h"
#include <cmath>

// Параметры модели Марсдена
double A1_NG = 0.0;  // радиальная компонента (а.е./день²)
double A2_NG = 0.0;  // поперечная компонента
double A3_NG = 0.0;  // нормальная компонента

// Параметры функции активности для водяного льда
const double NG_R0 = 2.808;
const double NG_M = 2.15;
const double NG_N = 5.093;
const double NG_K = 4.6142;
const double NG_ALPHA = 0.1112624;

// Функция активности g(r) для H2O
//double marsden_g(double r_au) {
//    double x = r_au / NG_R0;
//    return NG_ALPHA * pow(x, -NG_M) * pow(1.0 + pow(x, NG_N), -NG_K);
//}

// Функция активности g(r) для CO2
double marsden_g(double r_au) {
    return pow(r_au, -4.5);
}

// Вычисление негравитационного ускорения в барицентрической системе
Vector3 computeNGAcceleration(const StateVector& state) {
    double r = state.r.norm();

    // Единичные векторы орбитальной системы (RTN)
    Vector3 e_R = state.r / r;

    Vector3 h_vec = state.r.cross(state.v);
    double h = h_vec.norm();
    Vector3 e_N = h_vec / h;

    Vector3 e_T = e_N.cross(e_R);

    // Значение функции активности
    double g_val = marsden_g(r);

    // Суммарное ускорение
    return (e_R * A1_NG + e_T * A2_NG + e_N * A3_NG) * g_val;
}

StateVector cometDerivatives(double t, const StateVector& state) {
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
    Vector3 r_hat = state.r / r;
    double r_dot_p = r_hat * p;
    double coeff = 3.0 * J2_SUN * GM_SUN * R_SUN * R_SUN / (r * r * r * r);
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

    // 4. Негравитационные эффекты 
    Vector3 a_ng = computeNGAcceleration(state);
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

// Функции для установки параметров извне
void setNGParameters(double a1, double a2, double a3) {
    A1_NG = a1;
    A2_NG = a2;
    A3_NG = a3;
}

void getNGParameters(double& a1, double& a2, double& a3) {
    a1 = A1_NG;
    a2 = A2_NG;
    a3 = A3_NG;
}