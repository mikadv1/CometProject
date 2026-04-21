#include "Integrator.h"
#include <cmath>
#include <cstdio>

Trajectory::Trajectory() : t(nullptr), state(nullptr), nPoints(0), dt(0), t0(0), tend(0) {}

Trajectory::~Trajectory() {
	clear();
}

void Trajectory::allocate(int n, double start_t, double step) {
	clear();
	nPoints = n;
	dt = step;
	t0 = start_t;
	tend = start_t + (n - 1) * step;

	t = new double[nPoints];
	state = new StateVector[nPoints];

	for (int i = 0; i < nPoints; i++) {
		t[i] = t0 + i * dt;
	}
}

void Trajectory::clear() {
	if (t) delete[] t;
	if (state) delete[] state;
	t = nullptr;
	state = nullptr;
	nPoints = 0;
	dt = 0;
	t0 = 0;
	tend = 0;
}

void integrate(
    double t0,
    double tend,
    double grid_dt,      // шаг выходной сетки
    double internal_dt,  // внутренний шаг РК4
    const StateVector& state0,
    const NGVector& ng,
    Trajectory& traj
) {
    int nPoints = static_cast<int>((tend - t0) / grid_dt) + 1;
    traj.allocate(nPoints, t0, grid_dt);

    // Начальное состояние
    StateVector state = state0;
    double t = t0;
    int gridIdx = 0;

    // Записываем начальную точку
    traj.state[0] = state0;

    // Интегрирование
    while (t < tend) {
        // Определяем шаг: min(internal_dt, next_grid_t - t)
        double next_grid_t = t0 + (gridIdx + 1) * grid_dt;
        double step = (next_grid_t - t < internal_dt) ? (next_grid_t - t) : internal_dt;

        // Шаг РК4
        state = rk4Step(t, state, step, ng, cometDerivatives);
        t += step;

        // Если достигли следующей точки выходной сетки — записываем
        if (gridIdx + 1 < nPoints && t >= traj.t[gridIdx + 1] - 1e-12) {
            gridIdx++;
            traj.state[gridIdx] = state;
        }
    }
}

StateVector interpolateLinear(const Trajectory& traj, double t) {
    // Проверка границ
    if (t <= traj.t[0] || t >= traj.t[traj.nPoints - 1]) {
        printf("Time %.16le is out of bounds of trajectory", t);
        return StateVector();
    }

    int left = (int)(t - traj.t0) / traj.dt;
    int right = left + 1;

    // Коэффициент интерполяции
    double alpha = (t - traj.t[left]) / (traj.t[right] - traj.t[left]);

    // Линейная интерполяция позиции и скорости
    return StateVector(
        traj.state[left].r + (traj.state[right].r - traj.state[left].r) * alpha,
        traj.state[left].v + (traj.state[right].v - traj.state[left].v) * alpha
    );
}

bool isTimeInRange(const Trajectory& traj, double t) {
    return (t >= traj.t[0] && t <= traj.t[traj.nPoints - 1]);
}

