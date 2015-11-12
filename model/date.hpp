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

namespace Nextsim
{
double
dateStr2Num(std::string const& datestr)
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
} // Nextsim
#endif
