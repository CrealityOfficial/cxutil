#include "cxutil/slicer/meshslice.h"
#include "cxutil/slicer/slicepolygonbuilder.h"
#include "trimesh2/TriMesh.h"

#include "cxutil/slicer/slicedmesh.h"
#include "cxutil/slicer/meshslice.h"
#include "cxutil/slicer/preslice.h"
#include "cxutil/slicer/slicehelper.h"

#include "ccglobal/tracer.h"

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
		slicedMeshLayer.z = z;
	}

	void sliceMeshes_src(std::vector<trimesh::TriMesh*>& meshes, std::vector<SlicedMesh>& slicedMeshes, std::vector<int>& z)
	{
		size_t meshCount = meshes.size();
		if (meshCount > 0)
		{
			if (slicedMeshes.size() != meshCount)
				slicedMeshes.resize(meshCount);

			for (size_t i = 0; i < meshCount; ++i)
				sliceMesh_src(meshes.at(i), slicedMeshes.at(i), z);
		}
	}


	void sliceMesh_src(trimesh::TriMesh* mesh, SlicedMesh& slicedMesh, std::vector<int>& z)
	{
		size_t size = z.size();
		if (size > 0)
		{
			slicedMesh.m_layers.resize(size);

			SliceHelper helper;
			helper.prepare(mesh);

#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				sliceMesh_src(slicedMesh.m_layers.at(i), z.at(i), &helper);
		}
	}


	void sliceMesh_src(SlicedMeshLayer& slicedMeshLayer, int z, SliceHelper* helper)
	{
		SlicePolygonBuilder builder;
		builder.sliceOneLayer_dst(helper, z, &slicedMeshLayer.polygons, &slicedMeshLayer.openPolylines);
	}

	bool sliceInput(const DLPInput& input, SlicedResult& result, ccglobal::Tracer* tracer)
	{
		int meshCount = input.meshCount();
		if (meshCount == 0)
			return false;

		std::vector<int> z;
		buildSliceInfos(input, z);

		size_t layerCount = z.size();
		if (layerCount == 0)
			return false;

		if (tracer)
		{
			tracer->progress(0.3f);
			if (tracer->interrupt())
			{
				tracer->progress(1.0f);
				return false;
			}
		}

		const std::vector<MeshObjectPtr>& meshptres = input.Meshes;
		std::vector<MeshObject*> meshes;
		for (size_t i = 0; i < meshCount; ++i)
		{
			meshes.push_back(&*meshptres.at(i));
		}

		sliceMeshes(meshes, result.slicedMeshes, z);
		result.save("sliced");

		if (tracer)
		{
			tracer->progress(0.6f);
			if (tracer->interrupt())
			{
				tracer->progress(1.0f);
				return false;
			}
		}

		const DLPParam& param = input.Param;
		//////////
		result.connect();
		result.save("connected");

		if (tracer)
		{
			tracer->progress(0.7f);
			if (tracer->interrupt())
			{
				tracer->progress(1.0f);
				return false;
			}
		}

		const coord_t line_segment_resolution = param.line_segment_resolution;
		const coord_t line_segment_deviation = param.line_segment_deviation;
		float xy_offset = param.xy_offset;
		bool enable_xy_offset = param.enable_xy_offset;
		result.simplify(line_segment_resolution, line_segment_deviation, xy_offset, enable_xy_offset);

		if (tracer)
		{
			tracer->progress(0.9f);
			if (tracer->interrupt())
			{
				tracer->progress(1.0f);
				return false;
			}
		}

		result.save("result");
		return true;
	}
}

