## Synopsis

The library implements triangles (lines, points) intersection check in three-dimensional Euclidean space with a orthonormal right handed basis.

## Assembling

Can be built with a UNIX compiler (change compiler flags to build in Windows)

## Interface

*   include/triangle_intersection.hpp:
```
namespace TriangleIntersection {
	bool intersected(double tr0[9], double tr1[9]);
}
```