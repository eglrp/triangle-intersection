#include "triangle_intersection.hpp"
#include "triangle_intersection_private.hpp"

#include <algorithm>

#include <iostream>

namespace TriangleIntersection {

enum ThreePointsShape { POINT = 0, TRIANGE = 1, LINE_SEGMENT = 2 };

class Classification {
public:
    static ThreePointsShape threePoints(const Vec3 tr[3]) {
        if (isTriangle(tr))
            return TRIANGE;
        else if (isPoint(tr))
            return POINT;
        else
            return LINE_SEGMENT;
    }

private:
    static bool isPoint(const Vec3 tr[3]) {
        return tr[0] == tr[1] && tr[0] == tr[2];
    }
    static bool isTriangle(const Vec3 tr[3]) {
        return tr[0] != tr[1] && tr[0] != tr[2] && tr[1] != tr[2];
    }
    static bool isLineSegment(const Vec3 tr[3]) {
        return !isPoint(tr) && !isTriangle(tr);
    }
};

namespace Reshaping {

/// It's computed by solving triple production equation (known method, you
/// can find more information on https://en.wikipedia.org/wiki/Plane_(geometry))
static Plane getPlane(const Vec3 tr[3]) {
    Plane plane;
    plane.a = tr[0].y * tr[1].z -
              tr[0].y * tr[2].z -
              tr[1].y * tr[0].z +
              tr[1].y * tr[2].z +
              tr[2].y * tr[0].z -
              tr[2].y * tr[1].z;
    plane.b = -tr[0].x * tr[1].z +
              tr[0].x * tr[2].z +
              tr[1].x * tr[0].z -
              tr[1].x * tr[2].z -
              tr[2].x * tr[0].z +
              tr[2].x * tr[1].z;
    plane.c = tr[0].x * tr[1].y -
              tr[0].x * tr[2].y -
              tr[1].x * tr[0].y +
              tr[1].x * tr[2].y +
              tr[2].x * tr[0].y -
              tr[2].x * tr[1].y;
    plane.d = -tr[0].x * tr[1].y * tr[2].z +
              tr[0].x * tr[2].y * tr[1].z +
              tr[1].x * tr[0].y * tr[2].z -
              tr[1].x * tr[2].y * tr[0].z -
              tr[2].x * tr[0].y * tr[1].z +
              tr[2].x * tr[1].y * tr[0].z;

    return plane;
}

static void getLineFromThreePoins(const Vec3 tr[3], Vec3& lineR0, Vec3& lineDir) {
    if (tr[0] != tr[2]) { // != isn't binary
        lineR0 = tr[0];
        lineDir = tr[2] - tr[0];
    } else {
        lineR0 = tr[0];
        lineDir = tr[1] - tr[0];
    }
}

}

/// Used to improve a computing accuracy
static void linearReWeighing(double& d0, double& d1, double& d2, double& d3, double& d4, double& d5) {
    using std::max;

    const double max_ = max(max(max(max(max(fabs(d0), fabs(d1)), fabs(d2)), fabs(d3)), fabs(d4)), fabs(d5));
    //double max_ = vmax(fabs(d0), fabs(d1), fabs(d2), fabs(d3), fabs(d4), fabs(d5));

    const double weight = max_ < Algebra::EPS ? 1. : std::max(150. / max_, 1e-5);

    d0 *= weight;
    d1 *= weight;
    d2 *= weight;
    d3 *= weight;
    d4 *= weight;
    d5 *= weight;
}
static void linearReWeighing(Vec3 tr0[3], Vec3 tr1[3]) {
    linearReWeighing(tr0[0].x, tr0[1].x, tr0[2].x, tr1[0].x, tr1[1].x, tr1[2].x);
    linearReWeighing(tr0[0].y, tr0[1].y, tr0[2].y, tr1[0].y, tr1[1].y, tr1[2].y);
    linearReWeighing(tr0[0].z, tr0[1].z, tr0[2].z, tr1[0].z, tr1[1].z, tr1[2].z);
}

static bool isfinite(double data[9]) {
    for (int i = 0; i < (sizeof(data) / sizeof(*data)); i++)
        if (!std::isfinite(data[i]))
            return false;
    return true;
}

bool intersected(double tr0_p[9], double tr1_p[9]) {
    if (!isfinite(tr0_p) || !isfinite(tr1_p))
        return false;

    Vec3 tr0[3] = {Vec3(tr0_p), Vec3(tr0_p + 3), Vec3(tr0_p + 6)};
    Vec3 tr1[3] = {Vec3(tr1_p), Vec3(tr1_p + 3), Vec3(tr1_p + 6)};

    ThreePointsShape shape0 = Classification::threePoints(tr0);
    ThreePointsShape shape1 = Classification::threePoints(tr1);

    linearReWeighing(tr0, tr1);

    if (shape0 == shape1) {
        if (shape0 == TRIANGE)
            return triangleTriangleIntersected(tr0, Reshaping::getPlane(tr0), tr1, Reshaping::getPlane(tr1));
        else if (shape0 == LINE_SEGMENT) {
            Vec3 lineDir_a;
            Vec3 lineR0_a;
            Reshaping::getLineFromThreePoins(tr0, lineR0_a, lineDir_a);
            Vec3 lineDir_b;
            Vec3 lineR0_b;
            Reshaping::getLineFromThreePoins(tr1, lineR0_b, lineDir_b);
            return lineSegmentLineSegmentIntersected(lineR0_a, lineDir_a, lineR0_b, lineDir_b);
        }
        else
            return tr0[0] == tr1[0];
    } else {
        Vec3* tr_p_ordered[3] = {nullptr, nullptr, nullptr};
        tr_p_ordered[shape0] = tr0;
        tr_p_ordered[shape1] = tr1;

        if (tr_p_ordered[TRIANGE] && tr_p_ordered[LINE_SEGMENT]) {
            Vec3 lineDir;
            Vec3 lineR0;
            Reshaping::getLineFromThreePoins(tr_p_ordered[LINE_SEGMENT], lineR0, lineDir);
            return triangleLineSegmentIntersected(tr_p_ordered[TRIANGE], lineR0, lineDir, Reshaping::getPlane(tr_p_ordered[TRIANGE]));
        }
        if (tr_p_ordered[TRIANGE] && tr_p_ordered[POINT])
            return trianglePointIntersected(tr_p_ordered[TRIANGE], tr_p_ordered[POINT][0]);
        /*if (tr_p_ordered[LINE_SEGMENT] && tr_p_ordered[POINT])*/ {
            Vec3 lineDir;
            Vec3 lineR0;
            Reshaping::getLineFromThreePoins(tr_p_ordered[LINE_SEGMENT], lineR0, lineDir);
            return lineSegmentPointIntersected(lineR0, lineDir, tr_p_ordered[POINT][0]);
        }
    }
}

}
