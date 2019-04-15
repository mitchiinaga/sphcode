#pragma once

#include <exception>
#include <string>
#include <sstream>

#include "logger.hpp"

#define THROW_ERROR(...)\
    do {\
        std::ostringstream os;\
        sph::make_error_message(os, __VA_ARGS__);\
        throw sph::SPHException(os.str(), __FILE__, __func__, __LINE__);\
    } while(0)

namespace sph
{

inline void make_error_message(std::ostringstream & os) {}

template <typename Head, typename... Tail>
inline void make_error_message(std::ostringstream & os, Head && head, Tail &&... tail)
{
    os << head;
    make_error_message(os, std::move(tail)...);
}

class SPHException : public std::exception {
    std::string m_message;
    const char * m_file;
    const char * m_func;
    int m_line;
public:
    SPHException(const std::string & message) : m_message(message) {};
    SPHException(const std::string & message, const char * file, const char * func, const int line) :
        m_message(message), m_file(file), m_func(func), m_line(line) {};

    const char * get_file_name() const
    {
        return m_file;
    }

    const char * get_func_name() const
    {
        return m_func;
    }

    int get_line_number() const
    {
        return m_line;
    }
    const char * what() const throw()
    {
        return m_message.c_str();
    }
};

template <typename F>
inline void exception_handler(F && func)
{
    try {
        func();
    }
    catch(sph::SPHException & e) {
        if(sph::Logger::is_open()) {
            WRITE_LOG << "catch exception.";
            WRITE_LOG << e.what();
            WRITE_LOG << "file: " << e.get_file_name();
            WRITE_LOG << "func: " << e.get_func_name();
            WRITE_LOG << "line: " << e.get_line_number();
        } else {
            std::cerr << "catch exception." << std::endl;
            std::cerr << e.what() << std::endl;
            std::cerr << "file: " << e.get_file_name() << std::endl;
            std::cerr << "func: " << e.get_func_name() << std::endl;
            std::cerr << "line: " << e.get_line_number() << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
    catch(std::exception & e) {
        if(sph::Logger::is_open()) {
            WRITE_LOG << e.what();
        } else {
            std::cerr << e.what() << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
}

}