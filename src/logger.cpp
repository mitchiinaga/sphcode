#include <ctime>

#include <boost/format.hpp>

#include "logger.hpp"
#include "defines.hpp"
#include "exception.hpp"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace sph
{
std::string Logger::dir_name;
std::ofstream Logger::log_io;
bool Logger::open_flag = false;

void Logger::open(const std::string & output_dir)
{
    open(output_dir.c_str());
}

void Logger::open(const char * output_dir)
{
    struct stat st;
    bool is_mkdir = false;
    if(stat(output_dir, &st)) {
#ifdef _WIN32
        if(_mkdir(output_dir) == -1) {
#else
        if(mkdir(output_dir, 0775) == -1) {
#endif
            THROW_ERROR("cannot open directory");
        }
        is_mkdir = true;
    }

    std::time_t now = std::time(nullptr);
    std::tm * pnow = std::localtime(&now);
    std::string logfile = output_dir + (boost::format("/%04d%02d%02d%02d%02d%02d.log") % (pnow->tm_year + 1900) % (pnow->tm_mon + 1) % pnow->tm_mday % pnow->tm_hour % pnow->tm_min % pnow->tm_sec).str();
    log_io.open(logfile);
    if(is_mkdir) {
        log_io << "mkdir " << output_dir << std::endl;
    }
    dir_name = output_dir;
    open_flag = true;
}

}