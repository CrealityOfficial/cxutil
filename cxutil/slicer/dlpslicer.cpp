#include "cxutil/slicer/dlpslicer.h"
#include "cxutil/input/dlpinput.h"
#include "cxutil/slicer/slicedmesh.h"

#include "cxutil/slicer/meshslice.h"
#include "cxutil/slicer/slicepolygonbuilder.h"

#include <assert.h>

#ifdef _OPENMP
#include "omp.h"
#endif

namespace cxutil
{
	DLPSlicer::DLPSlicer()
	{

	}

	DLPSlicer::~DLPSlicer()
	{

	}

	bool DLPSlicer::compute(const DLPInput& input, DLPData& data, ccglobal::Tracer* tracer)
	{
#ifdef _OPENMP
		omp_set_num_threads(omp_get_num_procs());
#endif

		SlicedResult result;
		if (!sliceInput(input, result, tracer))
			return false;

		int layerCount = result.layerCount();
		if (layerCount <= 0)
			return false;

		data.layersData.resize(layerCount);
		for (SlicedMesh& slicedMesh : result.slicedMeshes)
		{
#pragma omp parallel for
			for (int layer_nr = 0; layer_nr < static_cast<int>(layerCount); layer_nr++)
			{
				SlicedMeshLayer& layer = slicedMesh.m_layers.at(layer_nr);
				DLPLayer& dlpplayer = data.layersData.at(layer_nr);
				dlpplayer.printZ = layer.z;
				dlpplayer.polygons = dlpplayer.polygons.unionPolygons(layer.polygons);
			}
		}

		return true;
	}

    bool DLPSlicer::compute(const DLPInput& input, float z, DLPDebugger* debugger)
    {
        int meshCount = input.meshCount();
        if (meshCount == 0)
            return false;

        SlicedMesh slicedMesh;
        std::vector<int> zs;
        zs.push_back(z*1000);
        sliceMesh(input.Meshes[0].get(), slicedMesh, zs);

        if (debugger)
        {
            if (slicedMesh.m_layers.size()>0)
            {
                debugger->onConnected(slicedMesh.m_layers[0].polygons, slicedMesh.m_layers[0].openPolylines);
            }
        }

        return  true;
    }

    OneLayerSlicer::OneLayerSlicer(MeshObjectPtr mesh)
        :m_mesh(mesh)
    {
        m_helper.prepare(mesh.get());
    }

    OneLayerSlicer::~OneLayerSlicer()
    {

    }

    bool OneLayerSlicer::compute(float z, DLPDebugger* debugger)
    {
        cxutil::coord_t iz = MM2INT(z);
        SlicePolygonBuilder builder;
        m_helper.sliceOneLayer(iz, builder.segments, builder.face_idx_to_segment_idx);

        Polygons polygons;
        Polygons openPolygons;
        builder.makePolygon(&polygons, &openPolygons);

        if (debugger)
        {
            debugger->onConnected(polygons, openPolygons);
        }

        return true;
    }
}