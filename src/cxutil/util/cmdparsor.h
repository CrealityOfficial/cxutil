#ifndef CX_CMDPARSOR_1600066272729_H
#define CX_CMDPARSOR_1600066272729_H
#include <functional>
#include <unordered_map>

namespace cxutil
{
	typedef std::function<void(const std::string& param)> processFunc;
	class CmdParsor
	{
	public:
		CmdParsor(int argc, char* argv[]);
		~CmdParsor();

		void registerFunc(const char* cmd, processFunc func);
		void process();
	protected:
		int m_argc;
		char** m_argv;
		std::unordered_map<std::string, processFunc> m_funcs;
	};
}

#endif // CX_CMDPARSOR_1600066272729_H