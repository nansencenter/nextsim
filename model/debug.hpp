/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   debug.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Wed Apr 27 10:53:14 2016
 */

#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_ 1

#include <iostream>
#include <sstream>


namespace Nextsim
{
	enum LogLevel
	{
		ERROR   = 0,
		WARNING = 1,
		INFO    = 2,
		DEBUG   = 3
	};


static const char* STRENUM[] = { "INFO", "WARNING", "DEBUG", "ERROR" };

class Log
{
public:
	Log(LogLevel loglevel = INFO)
	{
		_buffer << "["
		        << STRENUM[loglevel]
		        << "] :"
		        << std::string(1, ' ');
	}

	template <typename T>
	Log & operator<<(T const & value)
	{
		_buffer << value;
		return *this;
	}

	~Log()
	{
		//_buffer << std::endl;
		std::cerr << _buffer.str();
	}

private:
	std::ostringstream _buffer;
};

extern LogLevel M_log_level;

#define LOG(level) \
	if (level > M_log_level) ; \
	else Log(level)
#endif

} // Nextsim
