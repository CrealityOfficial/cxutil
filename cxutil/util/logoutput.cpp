#include "logoutput.h"

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include"ccglobal/log.h"
//#include "spdlog/cxlog_macro.h"

#define MAX_LOG_LEN 256
namespace cxutil
{
    static int verbose_level;
    static bool progressLogging;

    void increaseVerboseLevel()
    {
        verbose_level++;
    }

    void enableProgressLogging()
    {
        progressLogging = true;
    }

//   void logError(const char* fmt, ...)
//   {
//       char buf[MAX_LOG_LEN] = { 0 };
//       va_list args;
//       va_start(args, fmt);
//       vsprintf(buf, fmt, args);
//       CXLogError(buf);
///#pragma omp critical
///        {
///            fprintf(stderr, "[ERROR] ");
///            vfprintf(stderr, fmt, args);
///            fflush(stderr);
///        }
//       va_end(args);
//   }
//
//   void logError(const long long logSortId,const char* fmt, ...)
//   {
//       char buf[MAX_LOG_LEN] = { 0 };
//       va_list args;
//       va_start(args, fmt);
//       vsprintf(buf, fmt, args);
//       if (logSortId > 0)
//       {
//           std::stringstream ss;
//           ss << "[" << logSortId << "] ";
//           ss << buf;
//           strncpy(buf, ss.str().c_str(), ss.str().length());
//       }
//       CXLogError(buf);
//       va_end(args);
//   }

//    void logWarning(const char* fmt, ...)
//    {
//        char buf[MAX_LOG_LEN] = { 0 };
//        va_list args;
//        va_start(args, fmt);
//        vsprintf(buf, fmt, args);
//        CXLogWarn(buf);
////#pragma omp critical
////        {
////            fprintf(stderr, "[WARNING] ");
////            vfprintf(stderr, fmt, args);
////            fflush(stderr);
////        }
//        va_end(args);
//    }
//
//    void logWarning(const long long logSortId, const char* fmt, ...)
//    {
//        char buf[MAX_LOG_LEN] = { 0 };
//        va_list args;
//        va_start(args, fmt);
//        vsprintf(buf, fmt, args);
//        if (logSortId > 0)
//        {
//            std::stringstream ss;
//            ss << "[" << logSortId << "] ";
//            ss << buf;
//            strncpy(buf, ss.str().c_str(), ss.str().length());
//        }
//        CXLogWarn(buf);
//        va_end(args);
//    }

    void logAlways(const char* fmt, ...)
    {
        va_list args;
        va_start(args, fmt);
//#pragma omp critical
//        {
//            vfprintf(stderr, fmt, args);
//            fflush(stderr);
//        }
        va_end(args);
    }

    void log(const char* fmt, ...)
    {
        va_list args;
        if (verbose_level < 1)
            return;

        va_start(args, fmt);
//#pragma omp critical
//        {
//            vfprintf(stderr, fmt, args);
//            fflush(stderr);
//        }
        va_end(args);
    }

    //void logDebug(const char* fmt, ...)
    //{
    //    va_list args;
    //    if (verbose_level < 2)
    //    {
    //        return;
    //    }
    //    va_start(args, fmt);
    //    char buf[MAX_LOG_LEN] = { 0 };
    //    vsprintf(buf, fmt, args);
    //    CXLogDebug(buf);
//#pragma omp critical
//        {
//            fprintf(stderr, "[DEBUG] ");
//            vfprintf(stderr, fmt, args);
//            fflush(stderr);
//        }
    //    va_end(args);
    //}

    void logProgress(const char* type, int value, int maxValue, float percent)
    {
        if (!progressLogging)
            return;

//#pragma omp critical
//        {
//            fprintf(stderr, "Progress:%s:%i:%i \t%f%%\n", type, value, maxValue, percent);
//            fflush(stderr);
//        }
    }

}