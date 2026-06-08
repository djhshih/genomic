#ifndef cna_logging_h
#define cna_logging_h

#include <cstdarg>

enum LogLevel {
	LOG_ERROR_LEVEL,
	LOG_WARN_LEVEL,
	LOG_INFO_LEVEL,
	LOG_DEBUG_LEVEL,
	LOG_TRACE_LEVEL
};

bool log_enabled(LogLevel level);
void log_set_level(LogLevel level);
LogLevel log_level();
void log_write(LogLevel level, const char* file, int line, const char* function, const char* fmt, ...);
void log_trace(const char* file, int line, const char* function, const char* fmt, ...);
void log_debug(const char* file, int line, const char* function, const char* fmt, ...);
void log_info(const char* file, int line, const char* function, const char* fmt, ...);
void log_warn(const char* file, int line, const char* function, const char* fmt, ...);
void log_error(const char* file, int line, const char* function, const char* fmt, ...);

#endif
