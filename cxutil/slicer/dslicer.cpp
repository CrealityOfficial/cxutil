#include "dslicer.h"
#include "dlpslicer.h"

#include "cxutil/input/dlpinput.h"
#include "cxutil/util/meshbuilder.h"
#include "ccglobal/log.h"

namespace cxutil
{
	bool dlp_slice(const std::vector<DLPSliceObject>& objects, const DLPSliceParam& inParam,
		/*out*/DLPData& data, ccglobal::Tracer* tracer)
	{
		DLPInput input;

		cxutil::DLPParam& param = input.Param;

		param.initial_layer_thickness = DLP_MM2_S(inParam.initial_layer_thickness);
		param.layer_thickness = DLP_MM2_S(inParam.layer_thickness);

		param.xy_offset = inParam.xy_offset;
		//float hole_xy_offset = setting->floatParam("hole_xy_offset");
		param.enable_xy_offset = abs(inParam.xy_offset) > 0.000001;

		param.enable_b_offset = inParam.enable_b_offset;
		param.sThickness = inParam.sThickness * 1000.0f;
		param.sh = inParam.sh * 1000.0f;
		param.sLength = inParam.sLength * 1000.0f;
		param.sWidth = inParam.sWidth * 1000.0f;
		param.pixX = inParam.pixX;
		param.pixY = inParam.pixY;

		if (param.layer_thickness == 0)
		{
			LOGW("produceSliceInput error.  layer_thickness = [%d]", (int)param.layer_thickness);
			return false;
		}

		int count = (int)objects.size();
		if (count <= 0)
			return false;

		for (const DLPSliceObject& object : objects)
		{
			cxutil::MeshObjectPtr meshPtr(cxutil::meshFromTrimesh(object.mesh.get()));
			input.Meshes.push_back(meshPtr);
		}

		cxutil::DLPSlicer slicer;

		return slicer.compute(input, data, tracer);
	}
}