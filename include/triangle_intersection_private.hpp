#pragma once
#include "3d_algebra.hpp"

using Algebra::Vec3;
using Algebra::Plane;

namespace TriangleIntersection {

bool trianglePointIntersected(const Vec3 tr[3], Vec3 ver);

bool triangleLineSegmentIntersected(const Vec3 tr[3], const Vec3& lineR0, const Vec3& lineDir, const Plane& plane);

bool triangleTriangleIntersected(const Vec3 tr0[3], const Plane& plane0, const Vec3 tr1[3], const Plane& plane1);

bool lineSegmentPointIntersected(const Vec3& lineR0, const Vec3& lineDir, const Vec3& point);

bool lineSegmentLineSegmentIntersected(const Vec3& lineR0_a, const Vec3& lineDir_a, const Vec3& lineR0_b, const Vec3& lineDir_b);

}
