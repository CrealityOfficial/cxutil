#ifndef CX_DLPSLICER_1602646711722_H
#define CX_DLPSLICER_1602646711722_H
#include "cxutil/input/dlpdata.h"
#include "cxutil/input/dlpinput.h"

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
        virtual void onConnected(const Polygons& polygons, const Polygons& openPolygons) = 0;
    };

	class DLPSlicer
	{
	public:
		DLPSlicer();
		~DLPSlicer();

		bool compute(const DLPInput& input, DLPData& data, ccglobal::Tracer* tracer);

        //debug
        bool compute(const DLPInput& input, float z, DLPDebugger* debugger = nullptr);
	protected:
	};
}

#endif // CX_DLPSLICER_1602646711722_H