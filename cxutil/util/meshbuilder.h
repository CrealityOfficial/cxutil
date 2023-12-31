#ifndef CX_MESHBUILDER_1600070448786_H
#define CX_MESHBUILDER_1600070448786_H
#include "cxutil/math/Point3.h"
#include "cxutil/math/floatpoint.h"
#include <vector>

namespace cxutil
{
	class MeshObject;
	class Mesh;
	MeshObject* meshFromTriangleSoup(const std::vector<Point3>& soup);
	MeshObject* meshFromTriangleSoup(float* vertex, int vertexNum);
	MeshObject* meshFromTrimesh(float* vertex, int faceNum, int* faceIndex);
	Mesh* meshFromTrimesh2(int vertexNum, float* vertex, int faceNum, int* faceIndex);

	MeshObject* loadSTLBinaryMesh(const char* fileName);

	void loadSTLBinarySoup(const char* fileName, std::vector<Point3>& soups);
	void loadSTLBinarySoup(const char* fileName, std::vector<FPoint3>& soups);

	void convertf2i(const std::vector<FPoint3>& fpoints, std::vector<Point3>& points);
}

#endif // CX_MESHBUILDER_1600070448786_H