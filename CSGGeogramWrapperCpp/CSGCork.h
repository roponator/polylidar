#pragma once

#include "CommonStuff.h"


struct ConfigCork
{
public:

	float Tolerance_RelationTest = 0.0001f;
	float Tolerance_TryMakeVertex = 0.0001f;

	float PlaneComp_DotEps = 0.0001f;
	float PlaneComp_PosEps = 0.0001f;

	bool bMeshValidation = true;

	bool simplify_coplanar_facets_enabled;
	double simplify_coplanar_facets_angle_tolerance;

	bool use_delaunay;

	bool detect_intersecting_neighbors;

	bool fast_union;

	bool bMeshRepair_Enabled;
	double MeshRepair_VertWeldRelativeToLongestSeg;

	bool bForceSinglethreadedAndDeterministic;
};


extern "C" __declspec(dllexport) void CORK__Cpp_Bool_ResetAndStart(
	float PlaneComp_DotEps,
	float PlaneComp_PosEps,
	float Tolerance_RelationTest,
	float Tolerance_TryMakeVertex,
	bool bMeshValidation,
	bool simplify_coplanar_facets_enabled,
	double simplify_coplanar_facets_angle_tolerance,
	bool use_delaunay,
	bool detect_intersecting_neighbors,
	bool fast_union,
	bool bMeshRepair_Enabled,
	double MeshRepair_VertWeldRelativeToLongestSeg,
	bool bForceSinglethreadedAndDeterministic
);

//extern "C" __declspec(dllexport) void __Cpp_Terminate();
extern "C" __declspec(dllexport) void CORK__Cpp_Bool_AddBrush(InteropMesh* Mesh, BrushType BrushMode);
extern "C" __declspec(dllexport) InteropMesh* CORK__Cpp_Bool_GetMesh_AllocsMesh();
extern "C" __declspec(dllexport) void CORK__Cpp_Bool_Cleanup();
extern "C" __declspec(dllexport) void CORK__Cpp_Bool_FreeMesh(InteropMesh* Mesh);

