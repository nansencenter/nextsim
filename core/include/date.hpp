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

//! get the nextsim reference time as a boost::gregorian::date object
inline boost::gregorian::date getEpoch()
{
    return boost::gregorian::date(1900, 1, 1);
}


//! replace a substring in a string with another substring
inline std::string replaceSubstring(std::string input, std::string const& str_from,  std::string const& str_to)
{
    if(str_from.empty())
        return input;
    size_t start_pos = 0;
    while((start_pos = input.find(str_from, start_pos)) != std::string::npos)
    {
        input.replace(start_pos, str_from.length(), str_to);
        start_pos += str_to.length(); // In case 'str_to' contains 'str_from', like replacing 'x' with 'yx'
    }
    return input;
}

 
//! convert boost::posix_time::ptime object to a date-time string according to a given input
//! format
//! * for available formats see:
//!   https://www.boost.org/doc/libs/1_68_0/doc/html/date_time/date_time_io.html
//! * eg 2005-12-06 would come from format="%Y-%m-%d" (this is the default)
inline std::string posixTimeToString( boost::posix_time::ptime const& p_time,
        std::string const& format="%Y-%m-%d %H:%M:%S" )
{
    std::stringstream ss;
    std::cout<<format<<"\n";
    std::locale loc(ss.getloc(),
            new boost::posix_time::time_facet(format.c_str()));
        //NB pointer is deleted when loc is deleted (at the end of the function)
    ss.imbue(loc);
    ss << p_time;
    return ss.str();
}


//! convert date string to boost::posix_time::ptime object
//! * accepts %Y-%m-%d %H:%M:%S format or %Y-%m-%d
inline boost::posix_time::ptime stringToPosixTime( std::string datestr)
{
    if (datestr.find(" ") != std::string::npos)
        return boost::posix_time::time_from_string(datestr);
    else
        return boost::posix_time::ptime(
                boost::date_time::parse_date<boost::gregorian::date>(
                    datestr, boost::date_time::ymd_order_iso));
}


//! convert boost::posix_time::ptime to nextsim time (decimal days since 1900-1-1 0:00)
inline double posixTimeToDatenum( boost::posix_time::ptime const& p_time)
{
    if(p_time.date().year()<1900)
    {
        std::string msg = "bad year (" + std::to_string(p_time.date().year())
            + "): year should be >= 1900";
		throw std::logic_error(msg);
    }

    auto diff = p_time.date() - Nextsim::getEpoch();
    double datenum = diff.days();//get number of days
    datenum += p_time.time_of_day().total_seconds()/(24.*3600.);//add fractional day
    return datenum;
}


//! convert string to nextsim time (decimal days since 1900-1-1 0:00)
inline double stringToDatenum( std::string const& datestr )
{
    return posixTimeToDatenum(
            Nextsim::stringToPosixTime(datestr));
}


//! convert nextsim time (decimal days since 1900-1-1 0:00) to boost::gregorian::date
//! * \note this just gets the date and ignores the time
inline boost::gregorian::date parse_date( double date_time )
{
    boost::gregorian::date dt = Nextsim::getEpoch();
    dt += boost::gregorian::date_duration( static_cast<long>(date_time) );
    return dt;
}


//! convert nextsim time (decimal days since 1900-1-1 0:00) to boost::posix_time::time_duration
//! * \note this just gets the time and ignores the date
inline boost::posix_time::time_duration parse_time( double date_time )
{
    double fractionalDay = date_time - floor(date_time);
    long milliseconds = static_cast<long>( floor( fractionalDay * 24.0 * 60.0 * 60.0 * 1000.0 + 0.5) );
    return boost::posix_time::milliseconds( milliseconds );
}


//! convert nextsim time (decimal days since 1900-1-1 0:00) to boost::posix_time::ptime
inline boost::posix_time::ptime datenumToPosixTime( double const& datenum)
{
    return boost::posix_time::ptime( Nextsim::parse_date(datenum),
            Nextsim::parse_time(datenum) );
}


//! convert nextsim time (decimal days since 1900-1-1 0:00) to string
inline std::string datenumToString( double const& datenum, std::string const& format="%Y-%m-%d %H:%M:%S")
{
    auto p_time = datenumToPosixTime(datenum);
    return posixTimeToString(p_time, format);
}


//! convert year, month, day, hour, minute, second info directly to boost::posix_time::ptime object
inline boost::posix_time::ptime getPosixTime(int const& year, int const& month, int const& day,
        int const& hour=0, int const& minute=0, int const& seconds=0)
{
    boost::gregorian::date pdate(year, month, day);
    boost::posix_time::time_duration dtime = boost::posix_time::hours(hour)
        + boost::posix_time::minutes(minute) + boost::posix_time::seconds(seconds);
    return boost::posix_time::ptime(pdate, dtime);
}


//! convert year, month, day, hour, minute, second info directly to nextsim time
//! (decimal days since 1900-1-1 0:00)
inline double getDatenum(int const& year, int const& month, int const& day,
        int const& hour=0, int const& minute=0, int const& seconds=0)
{
    auto p_time = Nextsim::getPosixTime(year, month, day, hour, minute, seconds);
    return Nextsim::posixTimeToDatenum(p_time);
}


//! convert nextsim time (double) to a standate date-time string for use
//! in output filenames
//! * eg 20180101T123143Z
inline std::string datenumToFilenameString( double date_time )
{
    return Nextsim::datenumToString(date_time, "%Y%m%dT%H%M%SZ");
}


//! get a string representation of the current local time
inline std::string current_time_local()
{
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    return posix_time::to_simple_string(today_local);
}


//! get a string representation of the current UTC time
inline std::string current_time_UTC()
{
    posix_time::ptime today_utc(gregorian::day_clock::universal_day(), posix_time::second_clock::universal_time().time_of_day());
    return posix_time::to_simple_string(today_utc);
}


//! get a string representation of the time difference between the current time and the reference time
//! given as input (as a string)
inline std::string time_spent( const std::string& value )
{
    posix_time::ptime epoch = posix_time::time_from_string( value );
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    posix_time::time_duration diff = today_local - epoch;
    return posix_time::to_simple_string(diff);
}


} // Nextsim
#endif
