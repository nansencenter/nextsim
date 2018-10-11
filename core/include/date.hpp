/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   date.hpp
 * @author Abdoulaye Samake <abdama@beijing.wifi.ad.nersc.no>
 * @date   Wed Oct 14 16:25:10 2015
 */

#ifndef __Date_HPP
#define __Date_HPP 1

#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace gregorian = boost::gregorian;
namespace date_time = boost::date_time;
namespace posix_time = boost::posix_time;

namespace Nextsim
{
inline double dateStr2Num(std::string const& datestr)
{
	gregorian::date epoch = date_time::parse_date<gregorian::date>( "1900-01-01", date_time::ymd_order_iso);
	gregorian::date date = date_time::parse_date<gregorian::date>( datestr, date_time::ymd_order_iso);

	if (date.year() < 1900 )
	{
		std::cout << "bad year: year should be greater or equal 1900"<<"\n";
		throw std::logic_error("bad year: year should be greater or equal 1900");
	}

	gregorian::date_duration diff = date - epoch;
	return diff.days();
}

inline boost::gregorian::date parse_date( double date_time )
{
    boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1900-01-01", boost::date_time::ymd_order_iso );
    dt += boost::gregorian::date_duration( static_cast<long>( floor(date_time) ) );
    return dt;
}

inline boost::posix_time::time_duration parse_time( double date_time )
{
    double fractionalDay = date_time - floor(date_time);
    long milliseconds = static_cast<long>( floor( fractionalDay * 24.0 * 60.0 * 60.0 * 1000.0 + 0.5) );
    return boost::posix_time::milliseconds( milliseconds );
}

inline double from_date_string( const std::string& datestr )
{
    boost::gregorian::date epoch = boost::date_time::parse_date<boost::gregorian::date>( "1900-01-01", boost::date_time::ymd_order_iso);
    boost::gregorian::date date = boost::date_time::parse_date<boost::gregorian::date>( datestr, boost::date_time::ymd_order_iso);

    if (date.year() < 1900 )
	{
		std::cout << "bad year: year should be greater or equal 1900"<<"\n";
		throw std::logic_error("bad year: year should be greater or equal 1900");
	}

    boost::gregorian::date_duration diff = date - epoch;
    return diff.days();
}

inline std::string to_date_string( double date_time )
{
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    return (boost::format( "%1%-%2%-%3%" )
            % dt.year()
            % dt.month().as_number()
            % dt.day().as_number()
            ).str();
}

inline std::string to_date_string_yd( double date_time )
{
    // yyyymmdd
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    return (boost::format( "%1%%2%%3%" )
            % dt.year()
            % boost::io::group(std::setw(2), std::setfill('0'), dt.month().as_number())
            % boost::io::group(std::setw(2), std::setfill('0'), dt.day().as_number())
            ).str();
}

inline std::string to_date_string_ym( double date_time )
{
    // yyyymm
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    return (boost::format( "%1%%2%" )
            % dt.year()
            % boost::io::group(std::setw(2), std::setfill('0'), dt.month().as_number())
            ).str();
}

inline std::string to_date_string_y( double date_time )
{
    // yyyy
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    return (boost::format( "%1%" ) % dt.year()).str();
}

inline std::string to_date_string_md( double date_time )
{
    // mmdd
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    return (boost::format( "%1%%2%" )
            % boost::io::group(std::setw(2), std::setfill('0'), dt.month().as_number())
            % boost::io::group(std::setw(2), std::setfill('0'), dt.day().as_number())
            ).str();
}

inline double from_date_time_string( const std::string& datestr )
{
    double date = from_date_string( datestr );

    if (datestr.find(" ") != std::string::npos)
    {
        //std::string date_time_str = datestr + " 00:00:00";
        boost::posix_time::ptime t = boost::posix_time::time_from_string( datestr );
        double seconds = static_cast<double>(t.time_of_day().total_seconds());
        date += (seconds / 24.0 / 60.0 / 60.0);
    }

    // if milliseconds are needed
    // double milliseconds = static_cast<double>(t.time_of_day().total_milliseconds());
    // date += (milliseconds / 24.0 / 60.0 / 60.0 / 1000.0);

    return date;
}

inline std::string to_date_time_string( double date_time )
{
    // yyyy-mm-dd hh:mm:ss.fss
    boost::gregorian::date date_part = Nextsim::parse_date( date_time );
    boost::posix_time::time_duration time_part = Nextsim::parse_time( date_time );

    long long fractional_seconds = time_part.fractional_seconds();
    boost::date_time::time_resolutions resolution = time_part.resolution();
    if ( resolution == boost::date_time::micro )
    {
        fractional_seconds /= 1000;
    }
    else
    {
        if (resolution != boost::date_time::milli)
            throw std::logic_error( "Unexpected time resolution" );
    }

    return (boost::format( "%d-%02d-%02d %02d:%02d:%02d.%03d" )
            % date_part.year()
            % date_part.month().as_number()
            % date_part.day().as_number()
            % time_part.hours()
            % time_part.minutes()
            % time_part.seconds()
            % fractional_seconds
            ).str();
}

inline std::string to_date_time_string_for_filename( double date_time )
{
    // yyyymmddThhmmssZ
    boost::gregorian::date date_part = Nextsim::parse_date( date_time );
    boost::posix_time::time_duration time_part = Nextsim::parse_time( date_time );

    return (boost::format( "%d%02d%02dT%02d%02d%02dZ" )
            % date_part.year()
            % date_part.month().as_number()
            % date_part.day().as_number()
            % time_part.hours()
            % time_part.minutes()
            % time_part.seconds()
            ).str();
}

inline std::string current_time_local()
{
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    return posix_time::to_simple_string(today_local);
}

inline std::string current_time_UTC()
{
    posix_time::ptime today_utc(gregorian::day_clock::universal_day(), posix_time::second_clock::universal_time().time_of_day());
    return posix_time::to_simple_string(today_utc);
}

inline std::string time_spent( const std::string& value )
{
    posix_time::ptime epoch = posix_time::time_from_string( value );
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    posix_time::time_duration diff = today_local - epoch;
    return posix_time::to_simple_string(diff);
}

inline std::string model_time_str( const std::string& datestr, double time_in_seconds = 0)
{
    std::string date_time_str = datestr;
    if (date_time_str.find(" ") == std::string::npos)
    {
        date_time_str += " 00:00:00";
    }

    boost::posix_time::ptime posixtime = boost::posix_time::time_from_string( date_time_str );
    posixtime += boost::posix_time::time_duration(0,0,time_in_seconds);
    return to_simple_string(posixtime);
}


} // Nextsim
#endif
