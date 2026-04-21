#include "GaussNewton.h"
#include "CometModel.h"
#include "EphemerisWrapper.h"
#include "Constants.h"
#include "ephaccess.h"
#include <cmath>
#include <algorithm>

OrbitParams::OrbitParams() : r0(), v0(), ng0() {}

OrbitParams::OrbitParams(const Vector3& r, const Vector3& v, const NGVector ng)
    : r0(r), v0(v), ng0(ng) {}

OrbitParams OrbitParams::operator+(const OrbitParams& other) const {
    return { r0 + other.r0, v0 + other.v0, ng0 + other.ng0 };
}

OrbitParams OrbitParams::operator-(const OrbitParams& other) const {
    return { r0 - other.r0, v0 - other.v0, ng0 - other.ng0 };
}

OrbitParams OrbitParams::operator*(double scalar) const {
    return { r0 * scalar, v0 * scalar, ng0 * scalar };
}

ExtendedState::ExtendedState() : state(), dxdP{} {}

void ExtendedState::initIdentity() {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 9; j++) {
            dxdP[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

ExtendedState ExtendedState::operator+(const ExtendedState& other) const {
    ExtendedState result;
    result.state = state + other.state;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 9; j++) {
            result.dxdP[i][j] = dxdP[i][j] + other.dxdP[i][j];
        }
    }
    return result;
}

ExtendedState ExtendedState::operator*(double scalar) const {
    ExtendedState result;
    result.state = state * scalar;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 9; j++) {
            result.dxdP[i][j] = dxdP[i][j] * scalar;
        }
    }
    return result;
}

ExtendedState extendedDerivs(double t, const ExtendedState& aug, const NGVector& ng) {
    ExtendedState result;

    // 1. Базовая динамика
    StateVector deriv_state = cometDerivatives(t, aug.state, ng);
    result.state = deriv_state;

    // 2. Матрица dF/dx (6x6)
    std::array<std::array<double, 6>, 6> dFdx{};

    // Блок скоростей: d(dr/dt)/dv = I_3
    for (int i = 0; i < 3; i++) {
        dFdx[i][i + 3] = 1.0;
    }

    // Гравитация Солнца
    StateVector sun = getBodyState(10, t);
    Vector3 r_rel = aug.state.r - sun.r;
    double R = r_rel.norm();
    Vector3 r_hat = r_rel / R;
    double GM_R3 = GM_SUN / (R * R * R);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double kronecker = (i == j) ? 1.0 : 0.0;
            dFdx[i + 3][j] = -GM_R3 * (kronecker - 3.0 * r_hat[i] * r_hat[j]);
        }
    }

    // Возмущения от планет и астероидов
    for (int obj = 0; obj < NUM_OBJECTS; obj++) {
        StateVector body = getBodyState(OBJECTS[obj].code, t);
        Vector3 r_body_rel = aug.state.r - body.r;
        double d = r_body_rel.norm();
        Vector3 d_hat = r_body_rel / d;
        double GM_d3 = OBJECTS[obj].gm / (d * d * d);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double kronecker = (i == j) ? 1.0 : 0.0;
                dFdx[i + 3][j] -= GM_d3 * (kronecker - 3.0 * d_hat[i] * d_hat[j]);
            }
        }
    }

    // 3. Матрица dF/dP (6x9) - только для параметров Марсдена
    std::array<std::array<double, 6>, 9> dFdP{};

    double r_helio = (aug.state.r - sun.r).norm();
    double g_val = marsden_g(r_helio);

    Vector3 h_vec = aug.state.r.cross(aug.state.v);
    double h = h_vec.norm();
    Vector3 e_R = aug.state.r / r_helio;
    Vector3 e_N = h_vec / h;
    Vector3 e_T = e_N.cross(e_R);

    for (int i = 0; i < 3; i++) {
        dFdP[6][i + 3] = g_val * e_R[i];
        dFdP[7][i + 3] = g_val * e_T[i];
        dFdP[8][i + 3] = g_val * e_N[i];
    }

    // 4. Полная производная d(dx/dP)/dt = dF/dx * dx/dP + dF/dP
    for (int row = 0; row < 6; row++) {
        for (int col = 0; col < 9; col++) {
            double sum = 0.0;
            for (int k = 0; k < 6; k++) {
                sum += dFdx[row][k] * aug.dxdP[k][col];
            }
            result.dxdP[row][col] = sum + dFdP[col][row];
        }
    }

    return result;
}

ExtendedState integrateExtended(double t0, double t_end, double internal_dt,
    const ExtendedState& init, const NGVector& ng) {
    ExtendedState aug = init;
    double t = t0;

    while (t < t_end) {
        double dt = std::min(internal_dt, t_end - t);
        aug = rk4Step(t, aug, dt, ng, extendedDerivs);
        t += dt;
    }

    return aug;
}

Matrix2x9 computeObsJacobian(const ExtendedState& aug, const Vector3& r_observer) {
    Matrix2x9 J{};

    Vector3 rho = aug.state.r - r_observer;
    double x = rho.x, y = rho.y, z = rho.z;
    double R_xy = sqrt(x * x + y * y);
    double R = rho.norm();

    // d(Obs)/d(r) матрица 2x3
    std::array<std::array<double, 3>, 2> dObs_dr{};

    dObs_dr[0][0] = -y / (x * x + y * y);
    dObs_dr[0][1] = x / (x * x + y * y);
    dObs_dr[0][2] = 0.0;

    dObs_dr[1][0] = -x * z / (R * R * R_xy);
    dObs_dr[1][1] = -y * z / (R * R * R_xy);
    dObs_dr[1][2] = R_xy / (R * R);

    // Перевод в угловые секунды
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            dObs_dr[i][j] *= ARCSEC_PER_RAD;
        }
    }

    // Цепное правило: d(Obs)/dP = d(Obs)/d(state) * d(state)/dP
    for (int obsIdx = 0; obsIdx < 2; obsIdx++) {
        for (int pIdx = 0; pIdx < 9; pIdx++) {
            double sum = 0.0;
            for (int stateIdx = 0; stateIdx < 3; stateIdx++) {
                sum += dObs_dr[obsIdx][stateIdx] * aug.dxdP[stateIdx][pIdx];
            }
            J[obsIdx][pIdx] = sum;
        }
    }

    return J;
}

OrbitParams gaussNewtonStep(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt) {

    int nObs = static_cast<int>(obs.size());
    int nParams = 9;

    std::vector<double> residuals(2 * nObs);
    std::vector<Matrix2x9> J_mats(nObs);

    double t0 = traj.t0;

    for (int i = 0; i < nObs; i++) {
        const Observation& observation = obs[i];

        ExtendedState aug;
        aug.state.r = params.r0;
        aug.state.v = params.v0;
        aug.initIdentity();

        aug = integrateExtended(t0, observation.jd_utc, internal_dt, aug, params.ng0);

        // Вычислить модельное наблюдение через reduceObservation
        ReductionResult res;
        reduceObservation(observation, traj, res);

        // Невязки
        residuals[2 * i] = res.dRA;
        residuals[2 * i + 1] = res.dDec;

        // Якобиан
        StateVector earth = getBodyState(EPH_EARTH, observation.jd_utc);
        double jd_tdb = utc2tdb(observation.jd_utc);
        Vector3 r_station_gcrs = stationITRS2GCRS(observation.jd_utc, jd_tdb, observation.r_station_itrs);
        Vector3 r_observer = earth.r + r_station_gcrs / AU_KM;

        J_mats[i] = computeObsJacobian(aug, r_observer);
    }

    // Формируем нормальные уравнения: (J^T J) * dP = J^T r
    Matrix9x9 JTJ{};
    ParamVector JTr{};

    for (int i = 0; i < nObs; i++) {
        for (int row = 0; row < 2; row++) {
            int idx = 2 * i + row;
            for (int col = 0; col < nParams; col++) {
                JTr[col] += J_mats[i][row][col] * residuals[idx];
                for (int col2 = 0; col2 < nParams; col2++) {
                    JTJ[col][col2] += J_mats[i][row][col] * J_mats[i][row][col2];
                }
            }
        }
    }

    // Решаем систему
    std::array<bool, 9> selected = {1, 1, 1, 1, 1, 1, 0, 0, 0};
    // Решаем СЛАУ с селекцией параметров
    ParamVector step = solveCholesky(JTJ, JTr, selected);

    // Обновляем параметры
    OrbitParams newParams = params;
    newParams.r0.x -= step[0];
    newParams.r0.y -= step[1];
    newParams.r0.z -= step[2];
    newParams.v0.x -= step[3];
    newParams.v0.y -= step[4];
    newParams.v0.z -= step[5];
    newParams.ng0.A1 -= step[6];
    newParams.ng0.A2 -= step[7];
    newParams.ng0.A3 -= step[8];

    return newParams;
}

double computeResidualNorm(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt) {
    double sum = 0.0;

    for (const auto& observation : obs) {
        ReductionResult res;
        reduceObservation(observation, traj, res);

        double dRA = res.dRA;
        double dDec = res.dDec;
        sum += dRA * dRA + dDec * dDec;
    }

    return sqrt(sum / static_cast<double>(obs.size()));
}

double params_diff(const OrbitParams& p1, const OrbitParams& p2) {
    double rel_err[9];

    // Позиция
    rel_err[0] = (p1.r0.x - p2.r0.x) / p1.r0.x;
    rel_err[1] = (p1.r0.y - p2.r0.y) / p1.r0.y;
    rel_err[2] = (p1.r0.z - p2.r0.z) / p1.r0.z;

    // Скорость
    rel_err[3] = (p1.v0.x - p2.v0.x) / p1.v0.x;
    rel_err[4] = (p1.v0.y - p2.v0.y) / p1.v0.y;
    rel_err[5] = (p1.v0.z - p2.v0.z) / p1.v0.z;

    // Параметры Марсдена
    rel_err[6] = (p1.ng0.A1 - p2.ng0.A1) / p1.ng0.A1;
    rel_err[7] = (p1.ng0.A2 - p2.ng0.A2) / p1.ng0.A2;
    rel_err[8] = (p1.ng0.A3 - p2.ng0.A3) / p1.ng0.A3;

    // Норма вектора относительных ошибок
    double sum = 0.0;
    for (int i = 0; i < 9; i++) {
        sum += rel_err[i] * rel_err[i];
    }

    return sqrt(sum);
}

void printParams(const OrbitParams& p, int iteration) {
    printf("\n=== Iteration %d Parameters ===\n", iteration);
    printf("r0 = (%.16le, %.16le, %.16le) AU\n", p.r0.x, p.r0.y, p.r0.z);
    printf("v0 = (%.16le, %.16le, %.16le) AU/day\n", p.v0.x, p.v0.y, p.v0.z);
    printf("NG = (%.16le, %.16le, %.16le) AU/day^2\n", p.ng0.A1, p.ng0.A2, p.ng0.A3);
}

ParamVector solveCholesky(const Matrix9x9& A, const ParamVector& b, const std::array<bool, 9>& selected) {

    // 1. Построение карты выбранных параметров
    int n_selected = 0;
    int index_map[9];  // index_map[i] = оригинальный индекс параметра

    for (int i = 0; i < 9; i++) {
        if (selected[i]) {
            index_map[n_selected] = i;
            n_selected++;
        }
    }

    // 2. Формирование сокращённой системы
    std::array<std::array<double, 9>, 9> A_red{};  // макс. размер 9x9
    std::array<double, 9> b_red{};

    for (int i = 0; i < n_selected; i++) {
        b_red[i] = b[index_map[i]];
        for (int j = 0; j < n_selected; j++) {
            A_red[i][j] = A[index_map[i]][index_map[j]];
        }
    }

    // 3. Разложение Холецкого: A = L * L^T
    std::array<std::array<double, 9>, 9> L{};

    for (int i = 0; i < n_selected; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            if (i == j) {
                // Диагональный элемент
                for (int k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                double diag = A_red[j][j] - sum;
                L[j][j] = sqrt(diag);
            }
            else {
                // Недиагональный элемент
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A_red[i][j] - sum) / L[j][j];
            }
        }
    }

    // 4. Прямая подстановка: L * y = b_red
    std::array<double, 9> y{};
    for (int i = 0; i < n_selected; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b_red[i] - sum) / L[i][i];
    }

    // 5. Обратная подстановка: L^T * x_red = y
    std::array<double, 9> x_red{};
    for (int i = n_selected - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n_selected; j++) {
            sum += L[j][i] * x_red[j];
        }
        x_red[i] = (y[i] - sum) / L[i][i];
    }

    // 6. Отображение решения обратно в полный вектор
    ParamVector x{};
    for (int i = 0; i < n_selected; i++) {
        x[index_map[i]] = x_red[i];
    }

    return x;
}

OrbitParams fitOrbit(OrbitParams init_params,
    const std::vector<Observation>& obs,
    const Trajectory& planet_traj,
    double internal_dt,
    double tolerance,
    int max_iterations) {

    OrbitParams params = init_params;

    // Вывод начальных параметров
    printParams(params, 0);

    double initial_rms = computeResidualNorm(params, obs, planet_traj, internal_dt);
    printf("Initial RMS: %.3f arcsec\n", initial_rms);

    for (int iter = 0; iter < max_iterations; iter++) {
        OrbitParams prev = params;
        params = gaussNewtonStep(params, obs, planet_traj, internal_dt);

        // Вывод параметров после итерации
        printParams(params, iter + 1);

        double rms = computeResidualNorm(params, obs, planet_traj, internal_dt);
        double delta = params_diff(params, prev);

        printf("RMS = %.3f arcsec, delta = %.3e\n", rms, delta);

        if (delta < tolerance) {
            printf("Converged after %d iterations\n", iter + 1);
            break;
        }
    }

    return params;
}