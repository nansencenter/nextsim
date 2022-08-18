/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   optionhandler.hpp
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __OptionHandler_HPP
#define __OptionHandler_HPP 1

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>

namespace Nextsim
{

    class OptionHandler
    {
        public:

            //------------------------------------------------------------------------------------------------------
            //! given a map eg [("ec2", setup::AtmosphereType::EC2), ...]
            //! and an option name opt_name eg "setup.atmosphere-type" with
            //! vm[opt_name].as<std::string>() = "ec2", opt_val is set to setup::AtmosphereType::EC2
            //! Called by initOptAndParam()
            template<typename option_type>
            static option_type
            getOptionFromMap(boost::program_options::variables_map const &vm,
                    std::string const &opt_name, boost::unordered_map<const std::string, option_type> map)
            {

                if(vm.count(opt_name)==0)
                    throw std::runtime_error(
                            "OptionHandler::getOptionFromMap: Unknown option name: "
                            + opt_name+"\n");

                std::string const option_str = vm[opt_name].as<std::string>();
                if ( map.count(option_str) == 0 )
                {
                    // TODO: Get LOG() working instead of std::cerr
                    std::cerr << "OptionHandler::getOptionFromMap: Unknown option for "
                            << opt_name << ": " << option_str<< "\n";
                    std::cerr <<"Valid options are:\n";
                    for (auto ptr=map.begin(); ptr!=map.end(); ptr++)
                        std::cerr <<"  "<< ptr->first <<"\n";
                    throw std::runtime_error("Invalid option for "
                            +opt_name + ": " + option_str+"\n");
                }
                return map[option_str];
            };

            //------------------------------------------------------------------------------------------------------
            //! return an option, checking if it is allowed
            //! Called by initOptAndParam()
            template<typename option_type>
            static option_type
            getAllowedOption(boost::program_options::variables_map const &vm,
                    std::string const &opt_name, const std::vector<option_type> &options)
            {

                if(vm.count(opt_name)==0)
                    throw std::runtime_error(
                            "OptionHandler::getOptionFromMap: Unknown option name: "
                            + opt_name+"\n");

                option_type const option = vm[opt_name].as<option_type>();
                if(std::count(options.begin(), options.end(), option) == 0)
                {
                    // TODO: Get LOG() working instead of std::cerr
                    std::cerr << "OptionHandler::getAllowedOption: Unknown option for "
                            << opt_name << ": " << option << "\n";
                    std::cerr <<"Valid options are:\n";
                    for (auto opt : options)
                        std::cerr <<"  "<< opt <<"\n";
                    throw std::runtime_error("Invalid option for "
                            +opt_name + ": " + option+"\n");
                }
                return option;
            }; //getAllowedOption
    };
}
#endif
