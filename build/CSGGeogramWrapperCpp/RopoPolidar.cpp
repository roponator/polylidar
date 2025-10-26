
// TODO: KEEO ALL UNIONS IN SAME SCOPE UNTIL A SUBTRACTION COMES? MIGHT IMPROVE HOLES!

//#define ROPO_LOG 1

#include "RopoPolidar.h"

#include <Windows.h>
#include <string>
#include <filesystem>
#include <iostream>
#include <thread>
#include <string>

#include "Polylidar/Polylidar.hpp"

using namespace Polylidar;
using namespace std;


class CSGGeogramBoolean
{
public:

	struct vec3
	{
		double x, y, z;

		vec3(double x, double y, double z)
		{
			this->x = x;
			this->y = y;
			this->z = z;

		};

		vec3()
		{
			this->x = 0;
			this->y = 0;
			this->z = 0;

		};

		vec3 operator*(const vec3& other)const
		{
			return vec3(
				this->x * other.x,
				this->y * other.y,
				this->z * other.z
			);
		}

		vec3 operator/(const vec3& other)const
		{
			return vec3(
				this->x / other.x,
				this->y / other.y,
				this->z / other.z
			);
		}

		vec3 operator+(const vec3& other)const
		{
			return vec3(
				this->x + other.x,
				this->y + other.y,
				this->z + other.z
			);
		}

		vec3 operator-(const vec3& other)const
		{
			return vec3(
				this->x - other.x,
				this->y - other.y,
				this->z - other.z
			);
		}

		vec3 operator*(double other)const
		{
			return vec3(
				this->x * other,
				this->y * other,
				this->z * other
			);
		}

		vec3 operator/(double other)const
		{
			return vec3(
				this->x / other,
				this->y / other,
				this->z / other
			);
		}

		vec3 operator+(double other)const
		{
			return vec3(
				this->x + other,
				this->y + other,
				this->z + other
			);
		}

		vec3 operator-(double other)const
		{
			return vec3(
				this->x - other,
				this->y - other,
				this->z - other
			);
		}
	};


	static double dot(const vec3& a, const vec3& b)
	{
		return
			a.x * b.x +
			a.y * b.y +
			a.z * b.z;
	}

	static vec3 cross(const vec3& x, const vec3& y)
	{
		return vec3(
			x.y * y.z - y.y * x.z,
			x.z * y.x - y.z * x.x,
			x.x * y.y - y.x * x.y
		);
	}

	static double length(const vec3& a)
	{
		return sqrt(dot(a, a));
	}

	static vec3 normalize(const vec3& a)
	{
		double L = length(a);
		return a * (1.0 / L);
	}

	struct plane_t
	{
		vec3 normal;
		double     offset;
	};

	static plane_t make_plane(
		const vec3& point,
		const vec3& normal)
	{
		return plane_t{ normal, -dot(point, normal) };
	}

	static vec3 VecToGlm(const InteropVector3Float& V)
	{
		vec3 o;
		o.x = V.x;
		o.y = V.y;
		o.z = V.z;
		return o;
	}

	static InteropVector3Float GlmToVec(const vec3& V)
	{
		return { (float)V.x, (float)V.y, (float)V.z };
	}


	static bool ArePlanesEqual(const plane_t& p1, const plane_t& p2, float DotEps = 0.0001f, float EpsPos = 0.0001f)
	{
		// First normalize both (important so offset comparison is meaningful)
		vec3 n1 = normalize(p1.normal);
		vec3 n2 = normalize(p2.normal);
		float d1 = p1.offset / length(p1.normal);
		float d2 = p2.offset / length(p2.normal);

		if (dot(n1, n2) < (1.0f - DotEps))
		{
			return false;
		}

		return fabs(d1 - d2) < EpsPos;
	}

	static plane_t TriangleToPlane(const vec3& a, const vec3& b, const vec3& c)
	{
		plane_t p;
		p.normal = normalize(cross(b - a, c - a));
		p.offset = -dot(p.normal, a);  // todo: replace a with center?
		return p;
	}

	static std::vector<plane_t> ExtractUniquePlanes(const std::vector<plane_t>& Planes, double PlaneComp_DotEps, double PlaneComp_PosEps)
	{
		vector<plane_t> UniquePlanes;
		UniquePlanes.reserve(Planes.size());

		for (int i = 0; i < Planes.size(); ++i)
		{
			const auto& A = Planes[i];
			bool bUnique = true;

			// Loop against already added ones, not original array!
			for (int j = 0; j < UniquePlanes.size(); ++j)
			{
				const auto& B = UniquePlanes[j];

				if (ArePlanesEqual(A, B, PlaneComp_DotEps, PlaneComp_PosEps))
				{
					bUnique = false;
				}
			}

			if (bUnique)
			{
				UniquePlanes.push_back(A);
			}
		}

		return UniquePlanes;
	}

	static std::vector<plane_t> InteropMeshToPlanes(InteropMesh* M)
	{
		int NumFaces = M->NumIndices / 3;

		if ((M->NumIndices % 3) != 0)
		{
			Helpers::ShowErrorMsg("Indices must be div 3");
		}

		std::vector<plane_t> Result;
		Result.reserve(NumFaces);

		InteropVector3Float* Verts = (InteropVector3Float*)M->Float3DVertices;

		for (int i = 0; i < NumFaces; ++i)
		{
			if (true)
			{
				auto Plane = TriangleToPlane(
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 0]]),
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 1]]),
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 2]])
				);

				Result.push_back(Plane);
			}
			else
			{
				auto Plane = TriangleToPlane(
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 2]]),
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 1]]),
					VecToGlm(Verts[M->TriangleIndicesInt32[i * 3 + 0]])
				);

				Result.push_back(Plane);
			}


		}

		return Result;
	}

	static Polylidar::MeshHelper::HalfEdgeTriangulation InteropMeshToPolylidarMesh(const InteropMesh* Mesh)
	{

		size_t NumVerts = Mesh->NumVertices;
		size_t NumTris = Mesh->NumIndices / 3;
		size_t NumIndices = Mesh->NumIndices;

		if ((Mesh->NumIndices % 3) != 0)
		{
			Helpers::ShowErrorMsg("Indices must be div 3");
		}

		// TODO CLEANUP
		std::vector<double> RawVerts;
		RawVerts.reserve(NumVerts * 3);

		std::vector <size_t> RawIndices;
		RawIndices.reserve(NumIndices);

		for (int i = 0; i < NumVerts; ++i)
		{
			const float* MeshVert = &Mesh->Float3DVertices[i * 3];

			RawVerts.push_back(MeshVert[0]);
			RawVerts.push_back(MeshVert[1]);
			RawVerts.push_back(MeshVert[2]);
		}

		for (int i = 0; i < NumIndices; ++i)
		{
			RawIndices.push_back((size_t)Mesh->TriangleIndicesInt32[i]);
		}

		
		/** @brief Vertices in the mesh, N X 2 or N X 3 */
		Polylidar::Matrix<double> Vertices(std::move(RawVerts), NumVerts, 3);

		/** @brief Triangles in the mesh, K X 3 */
		Polylidar::Matrix<size_t> Triangles(std::move(RawIndices), NumTris, 3);
	
		auto Res = Polylidar::MeshHelper::HalfEdgeTriangulation(std::move(Vertices), std::move(Triangles));
		
		Res.ComputeTriangleNormals(); // MUST BE DONE OR IT CRASHES!

		return Res;
	}

	static NGonedFaceInteropMesh* ConvertTriangulatedMeshToNGonedFacedPolygon(InteropMesh* Mesh,
		double PlaneComp_DotEps,
		double PlaneComp_PosEps,
		double MaxTriangleEdgeLengthInSegment,
		double MinimumNumberOfTrianglesInAPlanarTriangleSegment,
		double MinimumNumberOfVerticesForAnInteriorHoleInAPolygon,
		double MaxPointToPlaneDistDuringRegionGrowing, // A value of 0.0 ignores this constraint.
		double NormThresh,
		double NormThreshMin //    Minimum value of the dot product between a triangle and surface normal being extracted
	)
	{

		double __alpha_NOT_USED_2D_CLOUD_ONLY = PL_DEFAULT_ALPHA;
		int _task_threads = 32;

		/*
			 * @brief Construct a new Polylidar 3D object.
	 *
	 * @param _alpha              Maximum circumcircle radius (only applicable for 2D point sets). Leave at 0.0 for 3D
	 * Data.
	 * @param _lmax               Maximum triangle edge length in a planar triangle segment
	 * @param _min_triangles      Minimum number of triangles in a planar triangle segment
	 * @param _min_hole_vertices  Minimum number of vertices for an interior hole in a polygon
	 * @param _z_thresh           Maximum point to plane distance during region growing (3D only). A value of 0.0
	 * ignores this constraint.
	 * @param _norm_thresh        IGNORE - will be deprecated or repurposed (3D only)
	 * @param _norm_thresh_min    Minimum value of the dot product between a triangle and surface normal being extracted
	 * (3D only)
	 * @param _task_threads       Number of threads for task-based parallelism (Marl library)
	 */

	 // We must find unique planes which are used by polylidar.
	 // It's ok if a concave mesh has duplicate planes for topologicaly disconnected components
	 // since polylidar takes topology into account.
		auto RawMeshPlanes = InteropMeshToPlanes(Mesh);
		auto UniquePlanes = ExtractUniquePlanes(RawMeshPlanes, PlaneComp_DotEps, PlaneComp_PosEps);

		// todo cleanup

		int NumPlanes = UniquePlanes.size();
		double* RawPlanes = new double[NumPlanes * 3];

		for (int i = 0; i < NumPlanes; ++i)
		{
			double* RawPlaneNormalVec = &RawPlanes[i * 3];
			const plane_t* Plane = &UniquePlanes[i];

			RawPlaneNormalVec[0] = Plane->normal.x;
			RawPlaneNormalVec[1] = Plane->normal.y;
			RawPlaneNormalVec[2] = Plane->normal.z;
		}

		Polylidar::Matrix<double> LibPlanes(RawPlanes, UniquePlanes.size(), 3);

		auto PolyMesh = InteropMeshToPolylidarMesh(Mesh);

		Polylidar::Polylidar3D pl(
			__alpha_NOT_USED_2D_CLOUD_ONLY,
			MaxTriangleEdgeLengthInSegment,
			MinimumNumberOfTrianglesInAPlanarTriangleSegment,
			MinimumNumberOfVerticesForAnInteriorHoleInAPolygon,
			MaxPointToPlaneDistDuringRegionGrowing,
			NormThresh,
			NormThreshMin,
			_task_threads = 32);

		Polylidar::PlanesGroup planes_group;
		Polylidar::PolygonsGroup polygons_group;
		std::tie(planes_group, polygons_group) = pl.ExtractPlanesAndPolygonsOptimized(PolyMesh, LibPlanes);

		// Copy result to my interop ngon
		{
			NGonedFaceInteropMesh* ResultNGon = new NGonedFaceInteropMesh();

			ResultNGon->NumVertices = Mesh->NumVertices;

			ResultNGon->Float3DVertices = new float[Mesh->NumVertices * 3];

			for (int i = 0; i < ResultNGon->NumVertices; ++i)
			{
				ResultNGon->Float3DVertices[i * 3 + 0] = Mesh->Float3DVertices[i * 3 + 0];
				ResultNGon->Float3DVertices[i * 3 + 1] = Mesh->Float3DVertices[i * 3 + 1];
				ResultNGon->Float3DVertices[i * 3 + 2] = Mesh->Float3DVertices[i * 3 + 2];
			}


			// Pass 1: count num polys to create
			int NumPolysToCreate = 0;
			for (const auto& IterPolysForThisDominantPlane : polygons_group)
			{
				NumPolysToCreate += IterPolysForThisDominantPlane.size();
			}

			// Pass 2: create 

			ResultNGon->NumPolygons = NumPolysToCreate;
			ResultNGon->Polygons = new InteropNGonPoly*[NumPolysToCreate];

			int IdxPoly = 0;
			for (const auto& IterPolysForThisDominantPlane : polygons_group)
			{
				for (const auto& IterPoly : IterPolysForThisDominantPlane)
				{
					if (IterPoly.holes.size() > 0)
					{
						Helpers::ShowErrorMsg("Error: a polygon has holes, not sure if legit, bug? Comment out any try again if needed");
					}

					ResultNGon->Polygons[IdxPoly] = new InteropNGonPoly; // todo optimize away after I ensure struct interop is good, remove double ptr
					InteropNGonPoly* NewPoly = ResultNGon->Polygons[IdxPoly];

				
					const auto& PolyOutline = IterPoly.shell;
					
					int NumPolyIndices = PolyOutline.size();

					NewPoly->NumIndices = NumPolyIndices;
					NewPoly->PolyIndicesInt32 = new int[NumPolyIndices];

				
					for (int i = 0; i < NumPolyIndices; ++i)
					{
						NewPoly->PolyIndicesInt32[i] = PolyOutline[i];
					}

					++IdxPoly;
				}
			}

			return ResultNGon;

		} // End copy block


	}


};
/*

static InteropMesh __Cpp_Bool_ResetAndStart(
	float PlaneComp_DotEps,
	float PlaneComp_PosEps,
	float Tolerance_RelationTest,
	float Tolerance_TryMakeVertex
)
{
	CSGGeogramBoolean::__Cpp_ResetAndStart(
		PlaneComp_DotEps,
		PlaneComp_PosEps,
		Tolerance_RelationTest,
		Tolerance_TryMakeVertex);
}

*/

NGonedFaceInteropMesh* __Cpp_ConvertTriangulatedMeshToNGonedFacedPolygon(
	InteropMesh* Mesh,
	double PlaneComp_DotEps,
	double PlaneComp_PosEps,
	double MaxTriangleEdgeLengthInSegment,
	double MinimumNumberOfTrianglesInAPlanarTriangleSegment,
	double MinimumNumberOfVerticesForAnInteriorHoleInAPolygon,
	double MaxPointToPlaneDistDuringRegionGrowing, // A value of 0.0 ignores this constraint.
	double NormThresh,
	double NormThreshMin 
)
{
	return CSGGeogramBoolean::ConvertTriangulatedMeshToNGonedFacedPolygon(
		Mesh,
		PlaneComp_DotEps,
		PlaneComp_PosEps,
		MaxTriangleEdgeLengthInSegment,
		MinimumNumberOfTrianglesInAPlanarTriangleSegment,
		MinimumNumberOfVerticesForAnInteriorHoleInAPolygon,
		MaxPointToPlaneDistDuringRegionGrowing, // A value of 0.0 ignores this constraint.
		NormThresh,
		NormThreshMin
	);
}