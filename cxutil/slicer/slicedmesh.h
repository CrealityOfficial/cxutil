#ifndef CX_SLICEDMESH_1599726180169_H
#define CX_SLICEDMESH_1599726180169_H
#include <vector>
#include "cxutil/math/polygon.h"

namespace cxutil
{
	class SlicedMeshLayer
	{
	public:
		SlicedMeshLayer();
		~SlicedMeshLayer();

		Polygons polygons;
		Polygons openPolylines;
		ClipperLib::cInt z;
	};

	class SlicedMesh
	{
	public:
		SlicedMesh();
		~SlicedMesh();

		void save(int index, const std::string& prefix);

		std::vector<SlicedMeshLayer> m_layers;
	};

	class SlicedResult
	{
	public:
		SlicedResult();
		~SlicedResult();

		void save(const std::string& prefix);
		void connect();
		void simplify(coord_t resolution, coord_t deviation,
			float xy_offset, bool enable_xy_offset);
		int meshCount();
		int layerCount();

		std::vector<SlicedMesh> slicedMeshes;
	};
}

#endif // CX_SLICEDMESH_1599726180169_H