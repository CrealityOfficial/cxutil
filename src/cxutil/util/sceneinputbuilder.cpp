#include "cxutil/util/sceneinputbuilder.h"
#include "cxutil/input/sceneinput.h"

#include "cxutil/util/jsonsettingsloader.h"
#include "cxutil/util/meshbuilder.h"

namespace cxutil
{
	SceneInputBuilder::SceneInputBuilder()
		:m_input(nullptr)
		, m_currentGroup(nullptr)
	{
		m_input = new SceneInput();
		m_input->addExtruderInput(ExtruderInputPtr(new ExtruderInput()));
		m_input->addExtruderInput(ExtruderInputPtr(new ExtruderInput()));
	}

	SceneInputBuilder::~SceneInputBuilder()
	{

	}

	SceneInput* SceneInputBuilder::build()
	{
		return m_input;
	}

	void SceneInputBuilder::loadGlobalSettings(const std::string& settingFile)
	{
		loadJsonSetting(settingFile.c_str(), m_input->settings());
	}

	void SceneInputBuilder::loadMesh(const std::string& fileName)
	{
		MeshObjectPtr mesh(loadSTLBinaryMesh(fileName.c_str()));
		if (mesh)
		{
			MeshInputPtr meshInput(new MeshInput());
			meshInput->setMeshObject(mesh);

			if (!m_currentGroup) addGroup("");
			m_currentGroup->addMeshInput(meshInput);
		}
	}

	void SceneInputBuilder::addGroup(const std::string&)
	{
		GroupInputPtr groupInput(new GroupInput());
		m_input->addGroupInput(groupInput);

		m_currentGroup = &*groupInput;
	}
}