#ifndef SLICE_SLICEHELPER_1598712992630_H
#define SLICE_SLICEHELPER_1598712992630_H
#include "cxutil/input/meshobject.h"
#include "cxutil/math/point2.h"
#include <unordered_map>

namespace cxutil
{
	class SliceHelper
	{
	public:
		SliceHelper();
		~SliceHelper();

		void prepare(MeshObject* mesh);
		void sliceOneLayer(int z,
			std::vector<SlicerSegment>& segments, std::unordered_map<int, int>& face_idx_to_segment_idx);
		
		static SlicerSegment project2D(const Point3& p0, const Point3& p1, const Point3& p2, const coord_t z);
		static void buildMeshFaceHeightsRange(const MeshObject* mesh, std::vector<Point2>& heightRanges);
	protected:
		MeshObject* mesh;
		std::vector<Point2> faceRanges;
	};
}

#endif // SLICE_SLICEHELPER_1598712992630_H