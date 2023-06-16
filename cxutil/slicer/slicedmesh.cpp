#include "cxutil/slicer/slicedmesh.h"
#include "cxutil/processor/openpolygonprocessor.h"
#include "cxutil/slicer/slicepolygonbuilder.h"

namespace cxutil
{
	SlicedMeshLayer::SlicedMeshLayer()
		:z(0)
	{

	}

	SlicedMeshLayer::~SlicedMeshLayer()
	{

	}

	SlicedMesh::SlicedMesh()
	{

	}

	SlicedMesh::~SlicedMesh()
	{

	}

	void SlicedMesh::save(int index, const std::string& prefix)
	{
#if _DEBUG
		int layerCount = (int)m_layers.size();
		for (int layer_nr = 0; layer_nr < layerCount; layer_nr++)
		{
			SlicedMeshLayer& layer = m_layers.at(layer_nr);
			char name[256];
			sprintf(name, "%s_%d_%d", prefix.c_str(), (int)index, (int)layer_nr);
			layer.polygons.save(name);
		}
#endif
	}

	SlicedResult::SlicedResult()
	{

	}

	SlicedResult::~SlicedResult()
	{

	}

	void SlicedResult::save(const std::string& prefix)
	{
#if _DEBUG
		int count = meshCount();
		for (int i = 0; i < count; ++i)
			slicedMeshes.at(i).save(i, prefix);
#endif
	}

	void SlicedResult::connect()
	{
		for (SlicedMesh& slicedMesh : slicedMeshes)
		{
			int layerCount = (int)slicedMesh.m_layers.size();
#pragma omp parallel for
			for (int j = 0; j < (int)layerCount; ++j)
			{
				SlicedMeshLayer& layer = slicedMesh.m_layers.at(j);

				Polygons closedPolygons;
				if(layer.openPolylines.size() > 0)
					connectOpenPolygons(layer.openPolylines, closedPolygons);

				for (size_t k = 0; k < closedPolygons.size(); ++k)
				{
					layer.polygons.add(closedPolygons[k]);
				}

				{
					Polygons stitchClosedPolygons;
					stitch(layer.openPolylines, stitchClosedPolygons);

					for (size_t k = 0; k < stitchClosedPolygons.size(); ++k)
					{
						layer.polygons.add(stitchClosedPolygons[k]);
					}
				}

				if (false)
				{
					stitchExtensive(layer.openPolylines, layer.polygons);
				}

				if (false)
				{
					for (PolygonRef polyline : layer.openPolylines)
					{
						if (polyline.size() > 0)
							layer.polygons.add(polyline);
					}
				}

				Polygons resultPolygons;
				for (PolygonRef polyline : layer.openPolylines)
				{
					if (polyline.size() > 0)
					{
						resultPolygons.add(polyline);
					}
				}
				layer.openPolylines = resultPolygons;
			}
		}
	}

	void SlicedResult::simplify(coord_t line_segment_resolution, coord_t line_segment_deviation,
		float xy_offset, bool enable_xy_offset)
	{
		for (SlicedMesh& slicedMesh : slicedMeshes)
		{
			SlicePolygonBuilder builder;
			int layerCount = (int)slicedMesh.m_layers.size();
#pragma omp parallel for
			for (int j = 0; j < (int)layerCount; ++j)
			{
				SlicedMeshLayer& layer = slicedMesh.m_layers.at(j);
				//const coord_t snap_distance = std::max(param.minimum_polygon_circumference, static_cast<coord_t>(1));
				//auto it = std::remove_if(layer.polygons.begin(), layer.polygons.end(), [snap_distance](PolygonRef poly) {
				//	return poly.shorterThan(snap_distance);
				//	});
				//layer.polygons.erase(it, layer.polygons.end());

				//Finally optimize all the polygons. Every point removed saves time in the long run.
				layer.polygons.simplify(line_segment_resolution, line_segment_deviation);
				layer.polygons.removeDegenerateVerts(); // remove verts connected to overlapping line segments

				if (enable_xy_offset && abs(xy_offset) > 0.000001)
				{
					layer.polygons = layer.polygons.offset(DLP_MM2_S(xy_offset));
				}

                //ºÏ²¢±ÕºÏÂÖÀª
                layer.polygons = layer.polygons.unionPolygons();

				ClipperLib::Path intersectPoints;
				if (layer.openPolylines.paths.size())
				{
					builder.connectOpenPolylines(layer.polygons, layer.openPolylines, intersectPoints);
				}
			}
		}
	}

	int SlicedResult::meshCount()
	{
		return (int)slicedMeshes.size();
	}

	int SlicedResult::layerCount()
	{
		if (meshCount() > 0)
			return (int)slicedMeshes.at(0).m_layers.size();
		return 0;
	}
}