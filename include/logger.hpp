#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace sph
{

class Logger {
    static std::string dir_name;
    static std::ofstream log_io;
    static bool open_flag;
    std::ostringstream output_msg;
    bool output_log_file_only;
public:
    Logger(bool log_only = false) : output_log_file_only(log_only) {}
    ~Logger() {
        log_io << output_msg.str() << std::endl;
        if(!output_log_file_only) {
            std::cout << output_msg.str() << std::endl;
        }
    }
    static void open(const std::string & output_dir);
    static void open(const char * output_dir);
    static std::string get_dir_name() { return dir_name; }
    static bool is_open() { return open_flag; }

    template<typename T>
    Logger & operator<<(const T & msg) {
        output_msg << msg;
        return *this;
    }
};

#define WRITE_LOG      sph::Logger()
#define WRITE_LOG_ONLY sph::Logger(true)

}