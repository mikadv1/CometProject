#pragma once
#include "Types.h"
#include "Integrator.h"
#include "Reduction.h"
#include "CometModel.h"
#include <array>
#include <vector>
#include <cstdio>

using ParamVector = std::array<double, 9>;
using Matrix6x9 = std::array<std::array<double, 9>, 6>;
using Matrix2x9 = std::array<std::array<double, 9>, 2>;
using Matrix9x9 = std::array<std::array<double, 9>, 9>;

struct OrbitParams {
    Vector3 r0;
    Vector3 v0;
    NGVector ng0;

    OrbitParams();
    OrbitParams(const Vector3& r, const Vector3& v, const NGVector ng);

    OrbitParams operator+(const OrbitParams& other) const;
    OrbitParams operator-(const OrbitParams& other) const;
    OrbitParams operator*(double scalar) const;
};

struct ExtendedState {
    StateVector state;
    Matrix6x9 dxdP;

    ExtendedState();
    void initIdentity();

    ExtendedState operator+(const ExtendedState& other) const;
    ExtendedState operator*(double scalar) const;
};

// Правая часть расширенной системы
ExtendedState extendedDerivs(double t, const ExtendedState& aug, const NGVector& ng);

using extendedDerivFunc = ExtendedState(*)(double t, const ExtendedState& state, const NGVector&);

// Интегрирование с чувствительностью
ExtendedState integrateExtended(double t0, double t_end, double internal_dt, const ExtendedState& init, const NGVector& ng);

// Один шаг Гаусса-Ньютона
OrbitParams gaussNewtonStep(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt,
    const std::array<bool, 9> selected,
    ParamVector& errors);

// Полный цикл подгонки
OrbitParams fitOrbit(
    OrbitParams init_params,
    const std::vector<Observation>& obs,
    double t0,
    double tend,
    double internal_dt,
    const std::array<bool, 9> selected,
    double tolerance = 1.0e-15,
    int max_iterations = 50);

// Вспомогательные функции
double computeRMS(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt);

double computeWRMS(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt);

double computeS(const OrbitParams& params,
    const std::vector<Observation>& obs,
    const Trajectory& traj,
    double internal_dt);

void printParams(const OrbitParams& p, const ParamVector& errors, int iteration);

ParamVector solveCholesky(const Matrix9x9& A, const ParamVector& b,
    const std::array<bool, 9>& selected);