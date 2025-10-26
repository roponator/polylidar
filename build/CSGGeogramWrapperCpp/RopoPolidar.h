#pragma once

#include "CommonStuff.h"



extern "C" __declspec(dllexport) NGonedFaceInteropMesh* __Cpp_ConvertTriangulatedMeshToNGonedFacedPolygon(
	InteropMesh* Mesh,
	double PlaneComp_DotEps,
	double PlaneComp_PosEps,
	double MaxTriangleEdgeLengthInSegment,
	double MinimumNumberOfTrianglesInAPlanarTriangleSegment,
	double MinimumNumberOfVerticesForAnInteriorHoleInAPolygon,
	double MaxPointToPlaneDistDuringRegionGrowing, // A value of 0.0 ignores this constraint.
	double NormThresh,
	double NormThreshMin
);