#pragma once



#include <cstdlib>
#include <stdint.h>



#include <Windows.h>
#include <string>
#include <filesystem>
#include <iostream>


class Helpers
{
public:

	static void ShowErrorMsg(std::string Msg)
	{
		int msgboxID = MessageBoxA(
			NULL,
			(LPCSTR)Msg.c_str(),
			(LPCSTR)L"CSG Geogram Error",
			MB_OK
		);
	}


};

struct InteropVector3Float
{
	float x;
	float y;
	float z;
};

struct Face3D
{
	InteropVector3Float A;
	InteropVector3Float B;
	InteropVector3Float C;

};

struct InteropMesh
{
public:

	float_t* Float3DVertices;
	int32_t* TriangleIndicesInt32;
	int NumVertices;
	int NumIndices;

public:

	void Cleanup()
	{
		delete[] Float3DVertices;
		delete[] TriangleIndicesInt32;

		Float3DVertices = nullptr;
		TriangleIndicesInt32 = nullptr;

		NumVertices = 0;
		NumIndices = 0;
	}
};

enum class BrushType
{
	// DO NOT CHANGE ORDER! MAPS TO C# <-> C++!
	Positive = 0,
	Subtractive,
	Intersective,
	Unknown = 10000

};
