#ifndef CX_SCENEINPUTBUILDER_1600068255227_H
#define CX_SCENEINPUTBUILDER_1600068255227_H
#include <string>

namespace cxutil
{
	class SceneInput;
	class GroupInput;
	class SceneInputBuilder
	{
	public:
		SceneInputBuilder();
		~SceneInputBuilder();

		SceneInput* build();

		void loadGlobalSettings(const std::string& settingFile);
		void loadMesh(const std::string& fileName);
		void addGroup(const std::string& name);
	protected:
		SceneInput* m_input;
		GroupInput* m_currentGroup;
	};
}

#endif // CX_SCENEINPUTBUILDER_1600068255227_H