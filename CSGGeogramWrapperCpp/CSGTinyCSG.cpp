
// TODO: KEEO ALL UNIONS IN SAME SCOPE UNTIL A SUBTRACTION COMES? MIGHT IMPROVE HOLES!

//#define ROPO_LOG 1

#include "CSGTinyCSG.h"

#include <Windows.h>
#include <string>
#include <filesystem>
#include <iostream>
#include <thread>
#include "Exporter.h"
#include <csg.hpp>
#include <string>

using namespace csg;
using namespace std;
using namespace glm;

class CSGGeogramBoolean
{
public:

	static constexpr volume_t AIR = csg::AIR;
	static constexpr volume_t SOLID = csg::SOLID;

	static world_t* world;

	static int GlobalMeshCounter;
	static bool bIntersectionMode;


	static ConfigBoolean CurrentCfg;

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
		return { V.x, V.y, V.z };
	}

	static bool ArePlanesEqual(const plane_t& p1, const plane_t& p2, float DotEps = 0.0001f, float EpsPos = 0.0001f)
	{
		// First normalize both (important so offset comparison is meaningful)
		glm::vec3 n1 = glm::normalize(p1.normal);
		glm::vec3 n2 = glm::normalize(p2.normal);
		float d1 = p1.offset / glm::length(p1.normal);
		float d2 = p2.offset / glm::length(p2.normal);

		if (glm::dot(n1, n2) < (1.0f - DotEps))
		{
			return false;
		}

		return fabs(d1 - d2) < EpsPos;
	}

	static plane_t TriangleToPlane(const vec3& a, const vec3& b, const vec3& c)
	{
		plane_t p;
		p.normal = glm::normalize(glm::cross(b - a, c - a));
		p.offset = -glm::dot(p.normal, a);  // todo: replace a with center?
		return p;
	}

	static void __Cpp_ResetAndStart(
		float PlaneComp_DotEps,
		float PlaneComp_PosEps,
		float Tolerance_RelationTest,
		float Tolerance_TryMakeVertex
	)
	{
		CurrentCfg = ConfigBoolean();

		CurrentCfg.Tolerance_RelationTest = Tolerance_RelationTest;
		CurrentCfg.Tolerance_TryMakeVertex = Tolerance_TryMakeVertex;
		CurrentCfg.PlaneComp_DotEps = PlaneComp_DotEps;
		CurrentCfg.PlaneComp_PosEps = PlaneComp_PosEps;

		__Cpp_Cleanup();

		bIntersectionMode = false;

		world = new csg::world_t();

		world_t::Tolerance_RelationTest = CurrentCfg.Tolerance_RelationTest;
		world_t::Tolerance_TryMakeVertex = CurrentCfg.Tolerance_TryMakeVertex;

		// get/set the void volume
		volume_t void_volume = world->get_void_volume();
		world->set_void_volume(AIR);



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

	static void ExportBrushesToObj()
	{


		Exporter Ex;
		
		int FaceIdx = 0;
		int BrushIdx = 0;
		int NumAddedMeshes = 0;

		brush_t* brush = world->first();
		while (brush)
		{
			auto faces = brush->get_faces();
			for (face_t& face : faces)
			{
				int FragIdx = 0;
				for (fragment_t& fragment : face.fragments)
				{
					const vector<vertex_t>& vertices = fragment.vertices;
					vector<triangle_t> tris = triangulate(fragment);

					int NumVerts = vertices.size();
					int NumFaces = tris.size();


					// todo leaks mem
					Vector3f* Verts = new Vector3f[vertices.size()];
					Vector3f* Normals = new Vector3f[vertices.size()];
					Vector2f* UVs = new Vector2f[vertices.size()];
					int* Indices = new int[NumFaces * 3];

					for (int Vi = 0; Vi < vertices.size(); ++Vi)
					{
						const auto& V = vertices[Vi];
						Verts[Vi] = Vector3f(V.position.x, V.position.y, V.position.z);
						Normals[Vi] = Vector3f(1, 0, 0);
						UVs[Vi] = Vector2f(0, 0);
					}


					Ex.add_mesh(Vector3f(0, 0, 0), Vector3f(0, 0, 0), Vector3f(1, 1, 1), Verts, Normals, UVs, NumVerts, NumVerts, NumVerts);

					{
						int FaceCounter = 0;
						for (triangle_t& tri : tris)
						{
							int I0 = tri.i;
							int I1 = tri.k;
							int I2 = tri.j;

							Indices[FaceCounter * 3 + 0] = I0;
							Indices[FaceCounter * 3 + 1] = I1;
							Indices[FaceCounter * 3 + 2] = I2;

							++FaceCounter;
						}
					}

					int SubmeshIndex = NumAddedMeshes;

					Ex.set_mesh_submeshes_count(SubmeshIndex, NumFaces * 3);
					Ex.set_submesh_triangles(SubmeshIndex, 0, Indices, NumFaces * 3);

					int FrontVol = fragment.front_volume;
					int BackVol = fragment.back_volume;

					Ex.set_mesh_name(SubmeshIndex,
						std::string("mesh_Brush") + std::to_string(BrushIdx) +"_Face" + std::to_string(FaceIdx) + "_Frag" + std::to_string(FragIdx)
						+ "_VFront" + std::to_string(FrontVol) + "_VBack" + std::to_string(BackVol)
					);

					++FragIdx;
					++NumAddedMeshes;
				}

				++FaceIdx;
			}

			++BrushIdx;

			brush = world->next(brush);
		}


		Ex.export_obj("out", "debug");
	}

	static InteropMesh* __Cpp_GetMesh_AllocsMesh()
	{
		if (world == nullptr)
		{
			return nullptr;
		}

		auto rebuilt = world->rebuild();

		//ExportBrushesToObj();

		std::vector<Face3D> Faces;
		Faces.reserve(1024); // todo

		brush_t* brush = world->first();
		while (brush)
		{
			auto faces = brush->get_faces();
			for (face_t& face : faces)
			{
				for (fragment_t& fragment : face.fragments)
				{

					if (bIntersectionMode)
					{
						bool bCond1 =
							(
								(fragment.back_volume == (AIR) && (fragment.front_volume != AIR)) &&
								((fragment.back_volume != fragment.front_volume))
								);

						bool bCond2 = (fragment.back_volume == AIR || fragment.front_volume == AIR) && (fragment.back_volume > SOLID || fragment.front_volume > SOLID);

						bool bCond3 = (fragment.back_volume == GlobalMeshCounter) && (fragment.front_volume != GlobalMeshCounter);

						if (bCond3)
						{
							const vector<vertex_t>& vertices = fragment.vertices;

							vec3 normal = face.plane->normal;
							//bool flip_face = fragment.back_volume == GlobalMeshCounter;
							bool flip_face = false;

							if (bCond1)
							{
								flip_face = fragment.back_volume != GlobalMeshCounter;
							}
							else
							{
								flip_face = fragment.back_volume == GlobalMeshCounter;
							}

							if (flip_face)
							{
								normal = -normal;
							}

							//glNormal3fv(value_ptr(normal));

							vector<triangle_t> tris = triangulate(fragment);
							for (triangle_t& tri : tris)
							{
								if (flip_face)
								{
									Faces.push_back({ GlmToVec(vertices[tri.i].position),GlmToVec(vertices[tri.k].position) ,GlmToVec(vertices[tri.j].position) });
								}
								else
								{
									Faces.push_back({ GlmToVec(vertices[tri.i].position),GlmToVec(vertices[tri.j].position) ,GlmToVec(vertices[tri.k].position) });
								}
							}
						}
					}
					else
					{
						if (fragment.back_volume == SOLID && fragment.front_volume == AIR ||
							fragment.back_volume == AIR && fragment.front_volume == SOLID)
						{
							const vector<vertex_t>& vertices = fragment.vertices;

							vec3 normal = face.plane->normal;
							bool flip_face = fragment.back_volume == SOLID;
							if (flip_face)
							{
								normal = -normal;
							}

							//glNormal3fv(value_ptr(normal));

							vector<triangle_t> tris = triangulate(fragment);
							for (triangle_t& tri : tris)
							{
								if (flip_face)
								{
									Faces.push_back({ GlmToVec(vertices[tri.i].position),GlmToVec(vertices[tri.k].position) ,GlmToVec(vertices[tri.j].position) });
								}
								else
								{
									Faces.push_back({ GlmToVec(vertices[tri.i].position),GlmToVec(vertices[tri.j].position) ,GlmToVec(vertices[tri.k].position) });
								}
							}
						}
					}



				}


			}

			brush = world->next(brush);
		}

		// TODO must be freed!
		//InteropMesh* IM = MeshConverter::GeoMeshToInteropMesh_AllocsMesh(AccumulatedMesh, CurrentCfg.bMeshValidation);


		// todo cleanup!
		InteropMesh* M = new InteropMesh();

		auto NumVerts = Faces.size() * 3;
		auto NumTris = Faces.size();

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
				const auto& IterFace = Faces[i];

				M->TriangleIndicesInt32[i * 3 + 0] = i * 3 + 0;
				M->TriangleIndicesInt32[i * 3 + 1] = i * 3 + 1;
				M->TriangleIndicesInt32[i * 3 + 2] = i * 3 + 2;

				Verts[i * 3 + 0] = IterFace.A;
				Verts[i * 3 + 1] = IterFace.B;
				Verts[i * 3 + 2] = IterFace.C;
			}
		}

		return M;

	}

	static void __Cpp_Cleanup()
	{
#if ROPO_LOG
		if (bInit)
		{
			GEO::Logger::out("") << "__Cpp_Cleanup 1" << std::endl;
		}
#endif

		if (world != nullptr)
		{
			delete world;
		}
		world = nullptr;

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


		brush_t* brush = world->add();
		brush->set_time(GlobalMeshCounter);

		if (BrushMode == BrushType::Positive)
		{
			brush->set_volume_operation(make_fill_operation(SOLID));
		}
		else if (BrushMode == BrushType::Subtractive)
		{
			brush->set_volume_operation(make_fill_operation(AIR));
		}
		else if (BrushMode == BrushType::Intersective)
		{
			bIntersectionMode = true;
			Helpers::ShowErrorMsg("Interesction is not supported in tinycsg, could probably be implemented but is cancer, can't get it to work");
			//brush->set_volume_operation(make_intersection_operation());
		}
		else
		{
			Helpers::ShowErrorMsg("unknown csg mode");
		}

		// get/set the bounding planes
		vector<plane_t> planes = InteropMeshToPlanes(Mesh);

		vector<plane_t> UniquePlanes;
		UniquePlanes.reserve(planes.size());

		for (int i = 0; i < planes.size(); ++i)
		{
			const auto& A = planes[i];
			bool bUnique = true;

			// Loop against already added ones, not original array!
			for (int j = 0; j < UniquePlanes.size(); ++j)
			{
				const auto& B = UniquePlanes[j];

				if (ArePlanesEqual(A, B, CurrentCfg.PlaneComp_DotEps, CurrentCfg.PlaneComp_PosEps))
				{
					bUnique = false;
				}
			}

			if (bUnique)
			{
				UniquePlanes.push_back(A);
			}
		}

		brush->set_planes(UniquePlanes);


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

world_t* CSGGeogramBoolean::world = nullptr;
int CSGGeogramBoolean::GlobalMeshCounter = 0;
bool CSGGeogramBoolean::bIntersectionMode = false;

ConfigBoolean CSGGeogramBoolean::CurrentCfg = ConfigBoolean();


static void __Cpp_Bool_ResetAndStart(
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


InteropMesh* __Cpp_Bool_GetMesh_AllocsMesh()
{
	return CSGGeogramBoolean::__Cpp_GetMesh_AllocsMesh();
}

void __Cpp_Bool_Cleanup()
{
	CSGGeogramBoolean::__Cpp_Cleanup();
}


void __Cpp_Bool_AddBrush(InteropMesh* Mesh, BrushType BrushMode)
{
	return CSGGeogramBoolean::__Cpp_AddBrush(Mesh, BrushMode);

}


void __Cpp_Bool_FreeMesh(InteropMesh* Mesh)
{
	return CSGGeogramBoolean::__Cpp_FreeMesh(Mesh);
}