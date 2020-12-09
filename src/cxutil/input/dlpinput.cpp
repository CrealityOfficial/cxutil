#include "cxutil/input/dlpinput.h"

namespace cxutil
{
	DLPInput::DLPInput()
	{

	}

	DLPInput::~DLPInput()
	{
	}

	void DLPInput::addMeshObject(MeshObjectPtr object)
	{
		m_meshes.push_back(object);
	}

	const std::vector< MeshObjectPtr>& DLPInput::meshes() const
	{
		return m_meshes;
	}

	std::vector<MeshObjectPtr>& DLPInput::meshes()
	{
		return m_meshes;
	}

	AABB3D DLPInput::box()
	{
		AABB3D box;
		for (MeshObjectPtr& ptr : m_meshes)
			box.include(ptr->box());
		return box;
	}

	DLPParam& DLPInput::param()
	{
		return m_param;
	}
}