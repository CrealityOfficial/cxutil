#include "cxutil/input/dlpdata.h"

namespace cxutil
{
	DLPData::DLPData()
	{
	}

	DLPData::~DLPData()
	{
	}

	int DLPData::layers() const
	{
		return (int)layersData.size();
	}

	ClipperLib::PolyTree* DLPData::layerData(int layer) const
	{
		ClipperLib::PolyTree* tree = nullptr;
		if (layer >= 0 && layer < layers())
			tree = layersData.at(layer).polygons.splitIntoPolyTree(true);

		return tree;
	}
}
