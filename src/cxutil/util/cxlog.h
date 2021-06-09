#ifndef _CXLOG_H
#define _CXLOG_H
#include<string>
#include<memory>

namespace spdlog
{
	class logger;
}

namespace cxutil
{

#define MAX_LOG_LEN 128
#define CXLogInfo(...)  cxutil::CXLog::Instance().info(__VA_ARGS__)
#define CXLogDebug(...)  cxutil::CXLog::Instance().debug(__VA_ARGS__)
#define CXLogWarn(...) cxutil::CXLog::Instance().warn(__VA_ARGS__)
#define CXLogError(...)  cxutil::CXLog::Instance().error(__VA_ARGS__)
//#define CXLogCritical(...)  cxutil::CXLog::Instance().GetLogger()->critical(__VA_ARGS__)
//#define CXLogIns

	class CXLog
	{
	public:
		static CXLog& Instance();

		void InitCXLog(std::string file_name="cxslice.log", int maxSize=1024*10, int maxFiles=4, int log_level=0);

		void info(const char* fmt,...);
		void error(const char* fmt,...);
		void debug(const char* fmt, ...);
		void warn(const char* fmt, ...);

		void EndLog();

		void SetLevel(int level = 0);

		auto GetLogger()
		{
			return mp_logger_;
		}

	protected:
		CXLog() {};
		~CXLog() {};
		CXLog(const CXLog& other) = delete;
		CXLog& operator=(const CXLog& other) = delete;
	private:
		std::shared_ptr<spdlog::logger> mp_logger_;
		bool hasInitLog = false;
	};

}
#endif