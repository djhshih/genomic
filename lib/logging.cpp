#include "logging.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
	LogLevel current_level =
#if genomic_DEBUG == 1
		LOG_DEBUG_LEVEL;
#else
		LOG_WARN_LEVEL;
#endif

	LogLevel parse_level(const char* value) {
		if (value == NULL) return current_level;
		if (std::strcmp(value, "error") == 0) return LOG_ERROR_LEVEL;
		if (std::strcmp(value, "warn") == 0) return LOG_WARN_LEVEL;
		if (std::strcmp(value, "info") == 0) return LOG_INFO_LEVEL;
		if (std::strcmp(value, "debug") == 0) return LOG_DEBUG_LEVEL;
		if (std::strcmp(value, "trace") == 0) return LOG_TRACE_LEVEL;
		return current_level;
	}

	const char* level_name(LogLevel level) {
		switch (level) {
			case LOG_ERROR_LEVEL: return "ERROR";
			case LOG_WARN_LEVEL: return "WARN";
			case LOG_INFO_LEVEL: return "INFO";
			case LOG_DEBUG_LEVEL: return "DEBUG";
			case LOG_TRACE_LEVEL: return "TRACE";
			default: return "UNKNOWN";
		}
	}

	struct EnvConfig {
		EnvConfig() {
			current_level = parse_level(std::getenv("GENOMIC_LOG_LEVEL"));
		}
	};

	EnvConfig env_config;

	void log_vwrite(LogLevel message_level, const char* file, int line, const char* function, const char* fmt, va_list args) {
		if (!log_enabled(message_level)) return;
		std::fprintf(stderr, "[%s] %s:%d %s: ", level_name(message_level), file, line, function);
		std::vfprintf(stderr, fmt, args);
		std::fprintf(stderr, "\n");
	}
}

bool log_enabled(LogLevel level) {
	return level <= current_level;
}

void log_set_level(LogLevel level) {
	current_level = level;
}

LogLevel log_level() {
	return current_level;
}

void log_write(LogLevel message_level, const char* file, int line, const char* function, const char* fmt, ...) {
	if (!log_enabled(message_level)) return;

	std::fprintf(stderr, "[%s] %s:%d %s: ", level_name(message_level), file, line, function);
	va_list args;
	va_start(args, fmt);
	std::vfprintf(stderr, fmt, args);
	va_end(args);
	std::fprintf(stderr, "\n");
}

void log_trace(const char* file, int line, const char* function, const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_vwrite(LOG_TRACE_LEVEL, file, line, function, fmt, args);
	va_end(args);
}

void log_debug(const char* file, int line, const char* function, const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_vwrite(LOG_DEBUG_LEVEL, file, line, function, fmt, args);
	va_end(args);
}

void log_info(const char* file, int line, const char* function, const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_vwrite(LOG_INFO_LEVEL, file, line, function, fmt, args);
	va_end(args);
}

void log_warn(const char* file, int line, const char* function, const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_vwrite(LOG_WARN_LEVEL, file, line, function, fmt, args);
	va_end(args);
}

void log_error(const char* file, int line, const char* function, const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_vwrite(LOG_ERROR_LEVEL, file, line, function, fmt, args);
	va_end(args);
}
