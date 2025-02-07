#pragma once

#include<fstream>
#include<string>

enum LogPriority { Info, Warning, Error };

class Logger {
private:
	static LogPriority verbosity;
	static const char* filepath;

public:
	static void SetVerbosity(LogPriority newPriority);
	static void Log(LogPriority prio, std::string message);
};

