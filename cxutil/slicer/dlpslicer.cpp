#include "cxutil/slicer/dlpslicer.h"
#include "cxutil/input/dlpinput.h"
#include "cxutil/slicer/slicedmesh.h"

#include "cxutil/slicer/meshslice.h"

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
}