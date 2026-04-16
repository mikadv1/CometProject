#include "Types.h"
#include <cmath>

Vector3::Vector3() : x(0.0), y(0.0), z(0.0) {}

Vector3::Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vector3 Vector3::operator+(const Vector3& other) const {
    return { x + other.x, y + other.y, z + other.z };
}

Vector3 Vector3::operator-(const Vector3& other) const {
    return { x - other.x, y - other.y, z - other.z };
}

Vector3 Vector3::operator*(double scalar) const {
    return { x * scalar, y * scalar, z * scalar };
}

Vector3 Vector3::operator/(double scalar) const {
    return {x / scalar, y / scalar, z / scalar};
}

double Vector3::operator*(const Vector3& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

Vector3 Vector3::cross(const Vector3& other) const {
    return {
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    };
}

double Vector3::norm() const {
    return sqrt(x * x + y * y + z * z);
}

StateVector StateVector::operator+(const StateVector& other) const {
    return StateVector(r + other.r, v + other.v);
}

StateVector StateVector::operator-(const StateVector& other) const {
    return StateVector(r - other.r, v - other.v);
}

StateVector StateVector::operator*(double scalar) const {
    return StateVector(r * scalar, v * scalar);
}