#include "cxutil/slicer/meshslice.h"
#include "cxutil/slicer/slicepolygonbuilder.h"

#include "cxutil/slicer/slicehelper.h"
namespace cxutil
{
	void sliceMeshes(std::vector<MeshObject*>& meshes, std::vector<SlicedMesh>& slicedMeshes, std::vector<int>& z)
	{
		size_t meshCount = meshes.size();
		if (meshCount > 0)
		{
			if (slicedMeshes.size() != meshCount)
				slicedMeshes.resize(meshCount);

			for (size_t i = 0; i < meshCount; ++i)
				sliceMesh(meshes.at(i), slicedMeshes.at(i), z);
		}
	}

	void sliceMesh(MeshObject* mesh, SlicedMesh& slicedMesh, std::vector<int>& z)
	{
		size_t size = z.size();
		if (size > 0)
		{
			slicedMesh.m_layers.resize(size);
			
			SliceHelper helper;
			helper.prepare(mesh);

#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				sliceMesh(mesh, slicedMesh.m_layers.at(i), z.at(i), &helper);
		}
	}

	void sliceMesh(MeshObject* mesh, SlicedMeshLayer& slicedMeshLayer, int z, SliceHelper* helper)
	{
		SlicePolygonBuilder builder;
		helper->sliceOneLayer(z, builder.segments, builder.face_idx_to_segment_idx);

		builder.makePolygon(&slicedMeshLayer.polygons, &slicedMeshLayer.openPolylines);
	}
}