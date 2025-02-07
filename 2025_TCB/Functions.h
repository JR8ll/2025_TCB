#pragma once

#include "Log.h"

namespace TCB {
	extern Logger logger;
}

class ExcSched {
private:
	std::string message;
public:
	ExcSched(const std::string& s) : message(s) {}
	const std::string& getMessage() const { return message; }
};