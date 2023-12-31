#include "openpolygonprocessor.h"
#include "cxutil/slicer/slicer.h"

namespace cxutil
{
	void connectOpenPolygons(Polygons& openPolygons, Polygons& closedPolygons)
	{
		SlicerLayer objSlicer;
		objSlicer.polygons = closedPolygons;
		objSlicer.connectOpenPolylines(openPolygons);

		closedPolygons = objSlicer.polygons;
	}

	void stitch(Polygons& openPolygons, Polygons& closedPolygons)
	{
		SlicerLayer objSlicer;
		objSlicer.polygons = closedPolygons;
		objSlicer.stitch(openPolygons);

		closedPolygons = objSlicer.polygons;
	}

	void stitchExtensive(Polygons& openPolygons, Polygons& closedPolygons)
	{
		SlicerLayer objSlicer;
		objSlicer.polygons = closedPolygons;
		objSlicer.stitch_extensive(openPolygons);

		closedPolygons = objSlicer.polygons;
	}
}