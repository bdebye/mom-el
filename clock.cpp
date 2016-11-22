
#include "clock.h"
#include <ctime>

Clock::Clock() {
	this->start = clock();
}

void Clock::reset() {
	this->start = clock();
}

int Clock::second() {
	return (clock() - start) / CLOCKS_PER_SEC;
}

int Clock::minute() {
	return second() / 60;
}

int Clock::hour() {
	return minute() / 60;
}

int Clock::millisecond() {
	return (clock() - start);
}

time_ Clock::time_elapse() {
	time_ elp;
	elp.second = this->second() % 60;
	elp.minute = this->minute() % 60;
	elp.hour = this->hour() % 60;
	elp.millisecond = this->millisecond() % 1000;
	return elp;
}

std::string Clock::time_elapse_format() {
	time_ elp = this->time_elapse();
	char format[100];
#ifdef _WIN32
	sprintf_s(format, "%d hour %d min %d second %d", elp.hour, elp.minute, elp.second, elp.millisecond);
#else
	sprintf(format, "%d hour %d min %d second %d", elp.hour, elp.minute, elp.second, elp.millisecond);
#endif
	return std::string(format);
}
