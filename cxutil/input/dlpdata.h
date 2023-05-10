#ifndef _CXSW_DLPDATA_1593762618888_H
#define _CXSW_DLPDATA_1593762618888_H
#include "clipper/clipper.hpp"
#include "cxutil/math/polygon.h"

namespace cxutil
{
	struct DLPLayer
	{
		ClipperLib::cInt printZ;     //!< The height at which this layer needs to be printed. Can differ from sliceZ due to the raft.
		cxutil::Polygons polygons;
	};

	class DLPData
	{
	public:
		DLPData();
		virtual ~DLPData();

		int layers() const;
		ClipperLib::PolyTree* layerData(int layer) const;
		bool isValid() const;

		std::vector<DLPLayer> layersData;
	};
}


#endif // _CXSW_DLPDATA_1593762618888_H
