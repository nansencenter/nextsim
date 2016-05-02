/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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

double dateStr2Num(std::string const& datestr)
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

boost::gregorian::date parse_date( double date_time )
{
    //boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1899-12-30", boost::date_time::ymd_order_iso );
    boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( "1900-01-01", boost::date_time::ymd_order_iso );
    dt += boost::gregorian::date_duration( static_cast<long>( floor(date_time) ) );
    return dt;
}

boost::posix_time::time_duration parse_time( double date_time )
{
    double fractionalDay = date_time - floor(date_time);
    long milliseconds = static_cast<long>( floor( fractionalDay * 24.0 * 60.0 * 60.0 * 1000.0 + 0.5) );
    return boost::posix_time::milliseconds( milliseconds );
}

std::string to_date_string( double date_time )
{
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    //return (boost::format( "%4-%02d-%02d" ) % dt.year() % dt.month().as_number() % dt.day().as_number()).str();
    //return (boost::format( "%1%-%2%-%3%" ) % dt.year() % dt.month().as_number() % dt.day().as_number()).str();
    return (boost::format( "%1%-%2%-%3%" ) % dt.year() % dt.month().as_number() % dt.day().as_number()).str();
}

std::string to_date_string_ym( double date_time )
{
    boost::gregorian::date dt = Nextsim::parse_date( date_time );
    //return (boost::format( "%4-%02d-%02d" ) % dt.year() % dt.month().as_number() % dt.day().as_number()).str();
    //return (boost::format( "%1%-%2%-%3%" ) % dt.year() % dt.month().as_number() % dt.day().as_number()).str();
    //return (boost::format( "%1%%2%" ) % dt.year() % dt.month().as_number()).str();
    return (boost::format( "%1%%2%" ) % dt.year() % boost::io::group(std::setw(2), std::setfill('0'), dt.month().as_number())).str();
}

double from_date_string( const std::string& value )
{
    //boost::gregorian::date epoch = boost::date_time::parse_date<boost::gregorian::date>( "1899-12-30", boost::date_time::ymd_order_iso);
    boost::gregorian::date epoch = boost::date_time::parse_date<boost::gregorian::date>( "1900-01-01", boost::date_time::ymd_order_iso);
    boost::gregorian::date dt = boost::date_time::parse_date<boost::gregorian::date>( value, boost::date_time::ymd_order_iso);

    boost::gregorian::date_duration diff = dt - epoch;
    return diff.days();
}

std::string to_date_time_string( double date_time )
{
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
            % date_part.year() % date_part.month().as_number() % date_part.day().as_number()
            % time_part.hours() % time_part.minutes() % time_part.seconds() % fractional_seconds ).str();
}

double from_date_time_string( const std::string& value )
{
    double date = from_date_string( value );

    boost::posix_time::ptime t = boost::posix_time::time_from_string( value );
    double milliseconds = static_cast<double>(t.time_of_day().total_milliseconds());

    return date + (milliseconds / 24.0 / 60.0 / 60.0 / 1000.0);
}

std::string current_time_local()
{
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    return posix_time::to_simple_string(today_local);
}

std::string current_time_UTC()
{
    posix_time::ptime today_utc(gregorian::day_clock::universal_day(), posix_time::second_clock::universal_time().time_of_day());
    return posix_time::to_simple_string(today_utc);
}

std::string time_spent( const std::string& value )
{
    posix_time::ptime epoch = posix_time::time_from_string( value );
    posix_time::ptime today_local(gregorian::day_clock::local_day(), posix_time::second_clock::local_time().time_of_day());
    posix_time::time_duration diff = today_local - epoch;
    return posix_time::to_simple_string(diff);
}

} // Nextsim
#endif
