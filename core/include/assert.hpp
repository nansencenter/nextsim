/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

/**
 * @file   assert.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Fri Jul  3 16:28:01 2015
 */

#ifndef NDEBUG
#define ASSERT(condition, message) \
	do { \
		if ( ! (condition)) { \
			std::cerr << std::endl \
			          << "ERROR in " \
			          <<  __FILE__ \
			          << " line " << __LINE__ << ": " << message << std::endl \
			          << std::endl; \
			std::exit(EXIT_FAILURE); \
		} \
	} while (false)
#else
#define ASSERT(condition, message) do { } while (false)
#endif
