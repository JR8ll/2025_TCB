#include "Log.h"

void Logger::SetVerbosity(LogPriority newPriority) {
	verbosity = newPriority;
}

void Logger::Log(LogPriority prio, std::string message) {
	if (prio >= verbosity) {
		std::ofstream FILE(filepath, std::ios_base::app);
		switch (prio) {
		case Info: FILE << "Info:\t"; break;
		case Warning: FILE << "Warning:\t"; break;
		case Error: FILE << "ERROR:\t"; break;
		}
		FILE << message << "\n";
		FILE.close();
	}
}