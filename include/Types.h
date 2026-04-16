#pragma once

struct Vector3 {
    double x, y, z;

    Vector3();
    Vector3(double x_, double y_, double z_);
    Vector3 operator+(const Vector3& other) const;
    Vector3 operator-(const Vector3& other) const;
    Vector3 operator*(double scalar) const;
    Vector3 operator/(double scalar) const;
    double operator*(const Vector3& other) const;
    Vector3 cross(const Vector3& other) const;
    double norm() const;
};

struct StateVector {
    Vector3 r;  // Позиция (а.е.)
    Vector3 v;  // Скорость (а.е./день)

    StateVector operator+(const StateVector& other) const;
    StateVector operator-(const StateVector& other) const;
    StateVector operator*(double scalar) const;
};