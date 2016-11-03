#include "triangle_intersection_private.hpp"

using namespace Algebra;

namespace TriangleIntersection {

////////////////////////////////////////
//// trianglePointIntersected
//// This's the core of other triangle methods
////////////////////////////////////////

static bool isParallel(Vec3 a, Vec3 b) {
    return isZero(crossCanon(a, b));
}


static double getFactorOfParallelVectors(Vec3 a, Vec3 b) {
    double t;
    if (!isZero(b.x))
        t = a.x / b.x;
    else if (!isZero(b.y))
        t = a.y / b.y;
    else
        t = a.z / b.z;

    return t;
}

class PointInTriangle {
public:
    static bool isPointInTriangle(const Vec3 tr[3], const Vec3& point)
    {
        const Vec3& tr0 = tr[0];
        const Vec3& tr1 = tr[1];
        const Vec3& tr2 = tr[2];
        if (sameSide(point, tr0, tr1, tr2) && sameSide(point, tr1, tr0, tr2) && sameSide(point, tr2, tr0, tr1))
        {
            const Vec3 trNormal = crossCanon(tr0 - tr1, tr0 - tr2);
            if (fabs(dotCanon(tr0 - point, trNormal)) < EPS) //the point on the plane checking
                return true;
        }

        return false;
    }

private:
static bool sameSide(Vec3 p0, Vec3 tr0, Vec3 tr1, Vec3 tr2) {
        Vec3 v0 = crossCanon(tr2 - tr1, p0 - tr1);
        Vec3 v1 = crossCanon(tr2 - tr1, tr0 - tr1);
        if (isZero(v0) || isZero(v1))
            return true;
        if (!isParallel(v0, v1))
            return false;
        if (!isParallel(v0, v1))
            return false;
        return getFactorOfParallelVectors(v0, v1) >= -EPS;
    }

};

bool trianglePointIntersected(const Vec3 tr[3], Vec3 point) {
    return PointInTriangle::isPointInTriangle(tr, point);
}


////////////////////////////////////////
//// triangleLineSegmentIntersected
////////////////////////////////////////

static bool getPointInPlaneIntersectedWithLine(const Plane& plane, const Vec3& lineR0, const Vec3& lineDir, Vec3& ret) {
    const double weight = dotCanon(plane.n(), lineDir);
    if (isZero(weight)) //lineDir is parallel to plane
        return false;

    ret = lineR0 - ((dotCanon(lineR0, plane.n()) + plane.d) / weight) * lineDir;
    return true;
}

static bool isLineTimeWithinSegment(double time0, double time1) {
    if (time1 > 0.) {
        if (time0 < -EPS || time0 > time1 + EPS)
            return false;
    } else
        if (time0 > EPS || time0 < time1 - EPS)
            return false;
    return true;
}

bool triangleLineSegmentIntersected(const Vec3 tr[3], const Vec3& lineR0, const Vec3& lineDir, const Plane& plane) {
    Vec3 intersectedPoint;

    if (!getPointInPlaneIntersectedWithLine(plane, lineR0, lineDir, intersectedPoint)) { //lineDir is parallel to plane
        if (projCanon(plane, lineR0) == lineR0)
            return lineSegmentLineSegmentIntersected(lineR0, lineDir, tr[0], tr[1] - tr[0]) || //intersecting line segment with triangle line segments
                   lineSegmentLineSegmentIntersected(lineR0, lineDir, tr[0], tr[2] - tr[0]) ||
                   lineSegmentLineSegmentIntersected(lineR0, lineDir, tr[1], tr[2] - tr[1]);
        return false;
    }

    const Vec3 lineToTriange = intersectedPoint - lineR0;

    if (!isZero(lineToTriange.x)) {
        if (!isLineTimeWithinSegment(lineToTriange.x, lineDir.x))
            return false;
    } else if (!isZero(lineToTriange.y)) {
        if (!isLineTimeWithinSegment(lineToTriange.y, lineDir.y))
            return false;
    } else if (!isZero(lineToTriange.z)) {
        if (!isLineTimeWithinSegment(lineToTriange.z, lineDir.z))
            return false;
    }

    return trianglePointIntersected(tr, intersectedPoint);
}



////////////////////////////////////////////////////////////////////////////////
//// ltriangleTriangleIntersected, ineSegmentPointIntersected, lineSegmentLineSegmentIntersected
////////////////////////////////////////////////////////////////////////////////

bool triangleTriangleIntersected(const Vec3 tr0[3], const Plane& plane0, const Vec3 tr1[3], const Plane& plane1) {
    return triangleLineSegmentIntersected(tr0, tr1[0], tr1[1] - tr1[0], plane0) ||
           triangleLineSegmentIntersected(tr0, tr1[0], tr1[2] - tr1[0], plane0) ||
           triangleLineSegmentIntersected(tr1, tr0[0], tr0[1] - tr0[0], plane1) ||
           triangleLineSegmentIntersected(tr1, tr0[0], tr0[2] - tr0[0], plane1);
}

bool lineSegmentPointIntersected(const Vec3& lineR0, const Vec3& lineDir, const Vec3& point) {
    const Vec3 pointDir = point - lineR0;

    if (isZero(pointDir))
        return true;

    if (!isParallel(lineDir, pointDir))
        return false;

    const double t = getFactorOfParallelVectors(pointDir, lineDir);

    return t > -EPS && t < 1. + EPS;
}


bool lineSegmentLineSegmentIntersected(const Vec3& lineR0_a, const Vec3& lineDir_a, const Vec3& lineR0_b, const Vec3& lineDir_b) {
    const Vec3 a_b_crossed = crossCanon(lineDir_a, lineDir_b);
    if (isZero(a_b_crossed)) //lines are parallel
        return lineSegmentPointIntersected(lineR0_a, lineDir_a, lineR0_b) ||
               lineSegmentPointIntersected(lineR0_a, lineDir_a, lineR0_b + lineDir_b) ||
               lineSegmentPointIntersected(lineR0_b, lineDir_b, lineR0_a) ||
               lineSegmentPointIntersected(lineR0_b, lineDir_b, lineR0_a + lineDir_a);

    const Vec3 shift_b_crossed = crossCanon(lineR0_b - lineR0_a, lineDir_b);

    if (!isZero(shift_b_crossed) && !isParallel(a_b_crossed, shift_b_crossed))
        return false;

    const double aTime = getFactorOfParallelVectors(shift_b_crossed, a_b_crossed);

    if (aTime < -EPS || aTime > 1. + EPS)
        return false;

    const Vec3 interPoint = lineR0_a + aTime * lineDir_a - lineR0_b;

    const double bTime = getFactorOfParallelVectors(interPoint, lineDir_b);

    return bTime > -EPS && bTime < 1. + EPS;
}

}
