#ifndef CX_DLPSLICER_1602646711722_H
#define CX_DLPSLICER_1602646711722_H
#include "cxutil/input/dlpdata.h"
#include "cxutil/input/dlpinput.h"
#include "cxutil/slicer/slicehelper.h"
#include "cxutil/math/polygon.h"

namespace ccglobal
{
	class Tracer;
}

namespace cxutil
{
    class DLPDebugger
    {
    public:
        virtual ~DLPDebugger() {}

        virtual void onSegments() = 0;
        virtual void onConnected(const Polygons& polygons, const Polygons& openPolygons, const ClipperLib::Path& intersectionPoints) = 0;
    };

	class DLPSlicer
	{
	public:
		DLPSlicer();
		~DLPSlicer();

		bool compute(const DLPInput& input, DLPData& data, ccglobal::Tracer* tracer);

	protected:
	};

    class OneLayerSlicer
    {
    public:
        OneLayerSlicer(MeshObjectPtr mesh);
        ~OneLayerSlicer();

        bool compute(float z, DLPDebugger* debugger = nullptr);
    protected:
        MeshObjectPtr m_mesh;
        SliceHelper m_helper;
    };
}

#endif // CX_DLPSLICER_1602646711722_H