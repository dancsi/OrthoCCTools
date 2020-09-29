#pragma once

#include <utility>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

namespace logging {
	using console_logger_t = std::shared_ptr<spdlog::logger>;
	extern console_logger_t console;

	void setup(console_logger_t& console) {
		namespace spd = spdlog;
		auto&& real_console = spd::stdout_logger_st("console");
		console.swap(real_console);
		console->set_pattern("[%T] [%l] %v");
		console->flush_on(spd::level::trace);
		console->set_level(spd::level::trace);
		console->info("Logging set up");
	}

	void set_process_name(std::string name) {
		std::string fmt = "[%T] [" + name + "] [%l] %v";
		console->set_pattern(fmt);
	}

	template<typename... Args>
	void log_and_throw(const char* fmt, const Args&... args) {
		auto&& msg = fmt::format(fmt, args...);
		console->error(msg);
		throw std::runtime_error(msg);
	}
}