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
	class DLPSlicer
	{
	public:
		DLPSlicer();
		~DLPSlicer();

		DLPData* compute(DLPInput* input, ccglobal::Tracer* tracer);
	protected:
	};
}

#endif // CX_DLPSLICER_1602646711722_H