#ifndef CX_DLPINPUT_1602651441705_H
#define CX_DLPINPUT_1602651441705_H
#include "cxutil/input/meshobject.h"

namespace cxutil
{
	struct DLPParam
	{
		coord_t initial_layer_thickness;
		coord_t layer_thickness;

		coord_t minimum_polygon_circumference;
		coord_t line_segment_resolution;
		coord_t line_segment_deviation;

		DLPParam()
		{
			initial_layer_thickness = 300;
			layer_thickness = 100;

			minimum_polygon_circumference = 1000;
			line_segment_resolution = 50;
			line_segment_deviation = 50;
		}
	};

	class DLPInput
	{
	public:
		DLPInput();
		~DLPInput();

		void addMeshObject(MeshObjectPtr object);
		const std::vector< MeshObjectPtr>& meshes() const;
		std::vector<MeshObjectPtr>& meshes();

		DLPParam& param();
		AABB3D box();
	protected:
		std::vector<MeshObjectPtr> m_meshes;
		DLPParam m_param;
	};
}

#endif // CX_DLPINPUT_1602651441705_H