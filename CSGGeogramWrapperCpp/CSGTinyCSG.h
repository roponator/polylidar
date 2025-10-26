#pragma once

#include "CommonStuff.h"



struct ConfigBoolean
{
public:

	float Tolerance_RelationTest = 0.0001f;
	float Tolerance_TryMakeVertex = 0.0001f;

	float PlaneComp_DotEps = 0.0001f;
	float PlaneComp_PosEps = 0.0001f;
};


extern "C" __declspec(dllexport) void __Cpp_Bool_ResetAndStart(
	float PlaneComp_DotEps,
	float PlaneComp_PosEps,
	float Tolerance_RelationTest,
	float Tolerance_TryMakeVertex
);

//extern "C" __declspec(dllexport) void __Cpp_Terminate();
extern "C" __declspec(dllexport) void __Cpp_Bool_AddBrush(InteropMesh* Mesh, BrushType BrushMode);
extern "C" __declspec(dllexport) InteropMesh* __Cpp_Bool_GetMesh_AllocsMesh();
extern "C" __declspec(dllexport) void __Cpp_Bool_Cleanup();
extern "C" __declspec(dllexport) void __Cpp_Bool_FreeMesh(InteropMesh* Mesh);

