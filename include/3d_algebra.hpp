#pragma once
#include <math.h>
#include <cstdlib>
#include <array>
#include <iostream>

namespace Algebra {


/// Determines the computing accuracy
const double EPS = 1e-4;

class Vec3;

static Vec3 operator+(const Vec3& l, const Vec3& r);
static Vec3 operator-(const Vec3& l, const Vec3& r);
static Vec3 operator*(const Vec3& l, double scal);
static Vec3 operator*(double scal, const Vec3& r);

class Vec3 : public std::array<double, 3> {
public:
    Vec3() {}
    Vec3(double v0, double v1, double v2) {
        x = v0;
        y = v1;
        z = v2;
    }
    Vec3(double* v) {
        std::copy(v, v + 3, begin());
    }

    Vec3& operator =(const Vec3& src) {
        std::copy(src.begin(), src.end(), begin());
        return *this;
    }

    Vec3& operator +=(const Vec3& src) {
        (*this)[0] += src[0];
        (*this)[1] += src[1];
        (*this)[2] += src[2];
        return *this;
    }

    Vec3& operator -=(const Vec3& src) {
        (*this)[0] -= src[0];
        (*this)[1] -= src[1];
        (*this)[2] -= src[2];
        return *this;
    }

    Vec3& operator *=(double scal) {
        (*this)[0] *= scal;
        (*this)[1] *= scal;
        (*this)[2] *= scal;
        return *this;
    }

    Vec3 operator -() const {
        return Vec3(-(*this)[0], -(*this)[1], -(*this)[2]);
    }

    Vec3 crossCanon(const Vec3& r) const {
        return Vec3(y * r.z - z * r.y, z * r.x - x * r.z, x * r.y - y * r.x);
    }

    double dotCanon(const Vec3& r) const {
        return x * r.x + y * r.y + z * r.z;
    }

    double normL2() const {
        return sqrt(sqrNormL2());
    }

    double sqrNormL2() const {
        return x * x + y * y + z * z;
    }

    Vec3 normalizedL2() const {
        return (*this) * (1 / normL2());
    }

    double& x = data()[0];
    double& y = data()[1];
    double& z = data()[2];
};

class Plane : public std::array<double, 4> {
public:
    Plane() {}
    Plane(double a, double b, double c, double d) {
        (*this)[0] = a;
        (*this)[1] = b;
        (*this)[2] = c;
        (*this)[3] = d;
    }

    Vec3 n() const {
        return Vec3(a, b, c);
    }

    double& a = data()[0];
    double& b = data()[1];
    double& c = data()[2];
    double& d = data()[3];
};


static bool isZero(double v) {
    return v > -EPS && v < EPS;
}
static bool eq(double v0, double v1) {
    return isZero(v1 - v0);
}

static Vec3 operator+(const Vec3& l, const Vec3& r) {
    return Vec3(l.x + r.x, l.y + r.y, l.z + r.z);
}
static Vec3 operator-(const Vec3& l, const Vec3& r) {
    return Vec3(l.x - r.x, l.y - r.y, l.z - r.z);
}
static Vec3 operator*(const Vec3& l, double scal) {
    return Vec3(l.x * scal, l.y * scal, l.z * scal);
}
static Vec3 operator*(double scal, const Vec3& r) {
    return r * scal;
}
static bool operator==(const Vec3& l, const Vec3& r) {
    return eq(l.x, r.x) && eq(l.y, r.y) && eq(l.z, r.z);
}
static bool operator!=(const Vec3& l, const Vec3& r) {
    return !(l == r);
}

static bool isZero(Vec3 v) {
    return isZero(v.x) && isZero(v.y) && isZero(v.z);
}


/// "Canon" means "Canonical". It can be used in a Euclidean space with an orthonormal basis
static Vec3 crossCanon(const Vec3& l, const Vec3& r) {
    return l.crossCanon(r);
}

/// "Canon" means "Canonical". It can be used in a Euclidean space with an orthonormal basis
static double dotCanon(const Vec3& l, const Vec3& r) {
    return l.dotCanon(r);
}

/// "Canon" means "Canonical". It can be used in a Euclidean space with an orthonormal basis
static Vec3 projCanon(const Vec3& to, const Vec3& vec) {
    return (dotCanon(vec, to) / to.sqrNormL2()) * to;
}

/// "Canon" means "Canonical". It can be used in a Euclidean space with an orthonormal basis
static Vec3 projCanon(const Plane& to, const Vec3& vec) {
    return vec - ((dotCanon(vec, to.n()) + to.d) / to.n().sqrNormL2()) * to.n();
}

}
