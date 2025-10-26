
// TODO: KEEO ALL UNIONS IN SAME SCOPE UNTIL A SUBTRACTION COMES? MIGHT IMPROVE HOLES!

//#define ROPO_LOG 1

#include "CSGCork.h"

#include <Windows.h>
#include <string>
#include <filesystem>
#include <iostream>
#include <thread>


#include <cork.h>


using namespace std;


class CSGCork
{
public:


	static CorkTriMesh AccumMesh;

	static int GlobalMeshCounter;

	static ConfigCork CurrentCfg;

	static void __Cpp_ResetAndStart(
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
	)
	{
		CurrentCfg.Tolerance_RelationTest = Tolerance_RelationTest;
		CurrentCfg.Tolerance_TryMakeVertex = Tolerance_TryMakeVertex;
		CurrentCfg.PlaneComp_DotEps = PlaneComp_DotEps;
		CurrentCfg.PlaneComp_PosEps = PlaneComp_PosEps;
		CurrentCfg.bMeshValidation = bMeshValidation;
		CurrentCfg.simplify_coplanar_facets_enabled = simplify_coplanar_facets_enabled;
		CurrentCfg.simplify_coplanar_facets_angle_tolerance = simplify_coplanar_facets_angle_tolerance;
		CurrentCfg.use_delaunay = use_delaunay;
		CurrentCfg.detect_intersecting_neighbors = detect_intersecting_neighbors;
		CurrentCfg.fast_union = fast_union;
		CurrentCfg.bMeshRepair_Enabled = bMeshRepair_Enabled;
		CurrentCfg.MeshRepair_VertWeldRelativeToLongestSeg = MeshRepair_VertWeldRelativeToLongestSeg;
		CurrentCfg.bForceSinglethreadedAndDeterministic = bForceSinglethreadedAndDeterministic;


		__Cpp_Cleanup();



#if ROPO_LOG
		GEO::Logger::out("") << "__Cpp_ResetAndStart 1" << std::endl;
#endif




#if ROPO_LOG
		GEO::Logger::out("") << "__Cpp_ResetAndStart End" << std::endl;
#endif

	}

	static CorkTriMesh InteropMeshToCorkMesh(InteropMesh* M)
	{
		// TODO CLEAN MESH, SEE WORK MESH CLEAN FUNC!
		CorkTriMesh CM;
		CM.n_triangles = M->NumIndices / 3;
		CM.n_vertices = M->NumVertices;
	

		InteropVector3Float* Verts = (InteropVector3Float*)M->Float3DVertices;

		if (M->NumIndices > 0 && M->NumVertices > 0)
		{
			CM.triangles = new uint[M->NumIndices];
			CM.vertices = new float[M->NumVertices*3];

			for (int i = 0; i < M->NumVertices; ++i)
			{
				CM.vertices[i * 3 + 0] = Verts[i].x;
				CM.vertices[i * 3 + 1] = Verts[i].y;
				CM.vertices[i * 3 + 2] = Verts[i].z;

			}

			for (int i = 0; i < M->NumIndices; ++i)
			{
				CM.triangles[i] = M->TriangleIndicesInt32[i];
			}
		}

		return CM;
	}

	static InteropMesh* CorkMeshToInteropMesh(CorkTriMesh CM)
	{
		// TODO CLEAN MESH, SEE WORK MESH CLEAN FUNC!
			// todo cleanup!
		InteropMesh* M = new InteropMesh();

		auto NumVerts = CM.n_vertices;
		auto NumTris = CM.n_triangles;

		M->NumVertices = NumVerts;
		M->NumIndices = NumTris * 3;

		M->Float3DVertices = nullptr;
		M->TriangleIndicesInt32 = nullptr;

		if (NumVerts > 0)
		{
			M->Float3DVertices = new float[NumVerts * 3];
			M->TriangleIndicesInt32 = new int[NumTris * 3];

			InteropVector3Float* Verts = (InteropVector3Float*)M->Float3DVertices;

			for (int i = 0; i < NumTris; ++i)
			{
				M->TriangleIndicesInt32[i * 3 + 0] = CM.triangles[i * 3 + 0];
				M->TriangleIndicesInt32[i * 3 + 1] = CM.triangles[i * 3 + 1];
				M->TriangleIndicesInt32[i * 3 + 2] = CM.triangles[i * 3 + 2];
			}

			for (int i = 0; i < NumVerts; ++i)
			{
				Verts[i].x = CM.vertices[i * 3 + 0];
				Verts[i].y = CM.vertices[i * 3 + 1];
				Verts[i].z = CM.vertices[i * 3 + 2];

			}
		}

		return M;
	}

	static InteropMesh* __Cpp_GetMesh_AllocsMesh()
	{
		if (AccumMesh.n_vertices == 0 || AccumMesh.n_triangles == 0)
		{
			return nullptr;
		}

		return CorkMeshToInteropMesh(AccumMesh);

	}

	static void __Cpp_Cleanup()
	{
#if ROPO_LOG
		if (bInit)
		{
			GEO::Logger::out("") << "__Cpp_Cleanup 1" << std::endl;
		}
#endif

		// TODO: PROPERLY FREE BASED ON WHO OWNS IT!
		AccumMesh = CorkTriMesh();
		AccumMesh.n_triangles = 0;
		AccumMesh.n_vertices = 0;
		AccumMesh.triangles = nullptr;
		AccumMesh.vertices = nullptr;

		GlobalMeshCounter = 0;



#if ROPO_LOG
		if (bInit)
		{
			GEO::Logger::out("") << "__Cpp_Cleanup 2" << std::endl;
		}
#endif
	}


	static void __Cpp_AddBrush(InteropMesh* Mesh, BrushType BrushMode)
	{
		GlobalMeshCounter += 1;

#if ROPO_LOG
		GEO::Logger::out("") << "__Cpp_AddBrush 1. Num brushes: " << GlobalMeshCounter << std::endl;
#endif

		CorkTriMesh NewInBrush = InteropMeshToCorkMesh(Mesh);

		bool bIsFirstBrush = (AccumMesh.n_triangles == 0);

		CorkTriMesh TempRes;

		if (bIsFirstBrush)
		{
			// TODO SPECIAL CASE OF OWNERSHIP FOR MEMORY SINCE IT'S MY NOT CORKS!
			TempRes = NewInBrush;
		}
		else
		{

			if (BrushMode == BrushType::Positive)
			{
				computeUnion(AccumMesh, NewInBrush, &TempRes);
			}
			else if (BrushMode == BrushType::Subtractive)
			{
				computeDifference(AccumMesh, NewInBrush, &TempRes);
			}
			else if (BrushMode == BrushType::Intersective)
			{
				Helpers::ShowErrorMsg("Intersection not yet supported, though should be trivial imo, just too lazy now");
			}
			else
			{
				Helpers::ShowErrorMsg("unknown csg mode");
			}

		}

		// todo cleanup meshes based on who owns it!
		AccumMesh = TempRes;
		
	}


	static void __Cpp_FreeMesh(InteropMesh* Mesh)
	{
#if ROPO_LOG
		GEO::Logger::out("") << "__Cpp_FreeMesh 1" << std::endl;
#endif

		Mesh->Cleanup();

#if ROPO_LOG
		GEO::Logger::out("") << "__Cpp_FreeMesh 2" << std::endl;
#endif

	}
	

};

CorkTriMesh CSGCork::AccumMesh = CorkTriMesh();
int CSGCork::GlobalMeshCounter = 0;
ConfigCork CSGCork::CurrentCfg = ConfigCork();


static void CORK__Cpp_Bool_ResetAndStart(
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
)
{
	CSGCork::__Cpp_ResetAndStart(
		PlaneComp_DotEps,
		PlaneComp_PosEps,
		Tolerance_RelationTest,
		Tolerance_TryMakeVertex,
		bMeshValidation,
		simplify_coplanar_facets_enabled,
		simplify_coplanar_facets_angle_tolerance,
		use_delaunay,
		detect_intersecting_neighbors,
		fast_union,
		bMeshRepair_Enabled,
		MeshRepair_VertWeldRelativeToLongestSeg,
		bForceSinglethreadedAndDeterministic);
}


InteropMesh* CORK__Cpp_Bool_GetMesh_AllocsMesh()
{
	return CSGCork::__Cpp_GetMesh_AllocsMesh();
}

void CORK__Cpp_Bool_Cleanup()
{
	CSGCork::__Cpp_Cleanup();
}


void CORK__Cpp_Bool_AddBrush(InteropMesh* Mesh, BrushType BrushMode)
{
	return CSGCork::__Cpp_AddBrush(Mesh, BrushMode);

}


void CORK__Cpp_Bool_FreeMesh(InteropMesh* Mesh)
{
	return CSGCork::__Cpp_FreeMesh(Mesh);
}