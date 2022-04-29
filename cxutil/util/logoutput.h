#ifndef CXUTIL_LOGOUTPUT_1606213492619_H
#define CXUTIL_LOGOUTPUT_1606213492619_H
#include "ccglobal/log.h"

#define logError(...) LOGE(__VA_ARGS__)
#define logError(logsortId,...) LOGE(logsortId,##__VA_ARGS__)
#define logWarning(...) LOGW(__VA_ARGS__)
#define logWarning(logsortId,...) LOGW(logsortId,##__VA_ARGS__)
#define logDebug(...) LOGD(__VA_ARGS__)
#define logDebug(logsortId,...) LOGD(logsortId,##__VA_ARGS__)


namespace cxutil
{
    /*
     * \brief Increase verbosity level by 1.
     */
    void increaseVerboseLevel();

    /*
     * \brief Enable logging the current slicing progress to the log.
     */
    void enableProgressLogging();

    ///*
    // * \brief Report an error message.
    // *
    // * This is always reported, regardless of verbosity level.
    // */
    //void logError(const char* fmt, ...);
    //void logError(const long long logSortId, const char* fmt, ...);
    /*
     * \brief Report a warning message.
     *
     * Always reported, regardless of verbosity level.
     */
    //void logWarning(const char* fmt, ...);
    //void logWarning(const long long logSortId, const char* fmt, ...);

    /*
     * \brief Report a message if the verbosity level is 1 or higher.
     */
    void log(const char* fmt, ...);

    /*
     * \brief Log a message, regardless of verbosity level.
     */
    void logAlways(const char* fmt, ...);

    /*
     * \brief Log a debugging message.
     *
     * The message is only logged if the verbosity level is 2 or higher.
     */
    //void logDebug(const char* fmt, ...);

    /*
     * \brief Report the progress in the log.
     *
     * Only works if ``enableProgressLogging()`` has been called.
     */
    void logProgress(const char* type, int value, int maxValue, float percent);
}

#endif // CXUTIL_LOGOUTPUT_1606213492619_H