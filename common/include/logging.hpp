/* -*-c++-*-
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <cstdio>
#include <filesystem>
#include <stdio.h>
#include <unistd.h>

#define LOG_TRACE(...) \
    LOG_CALL(LOG_LEVEL_TRACE, LOG_NAME_TRACE, LOG_COLOR_TRACE, __VA_ARGS__)

#define LOG_DEBUG(...) \
    LOG_CALL(LOG_LEVEL_DEBUG, LOG_NAME_DEBUG, LOG_COLOR_DEBUG, __VA_ARGS__)

#define LOG_INFO(...) \
    LOG_CALL(LOG_LEVEL_INFO, LOG_NAME_INFO, LOG_COLOR_INFO, __VA_ARGS__)

#define LOG_WARN(...) \
    LOG_CALL(LOG_LEVEL_WARN, LOG_NAME_WARN, LOG_COLOR_WARN, __VA_ARGS__)

#define LOG_ERROR(...) \
    LOG_CALL(LOG_LEVEL_ERROR, LOG_NAME_ERROR, LOG_COLOR_ERROR, __VA_ARGS__)

#define LOG_CRITICAL(...) \
    LOG_CALL(LOG_LEVEL_CRITICAL, LOG_NAME_CRITICAL, LOG_COLOR_CRITICAL, __VA_ARGS__)

#define LOG_LEVEL_TRACE 0
#define LOG_LEVEL_DEBUG 1
#define LOG_LEVEL_INFO 2
#define LOG_LEVEL_WARN 3
#define LOG_LEVEL_ERROR 4
#define LOG_LEVEL_CRITICAL 5
#define LOG_LEVEL_OFF 6

#ifndef LOG_ACTIVE_LEVEL
#define LOG_ACTIVE_LEVEL LOG_LEVEL_INFO
#endif

#define LOG_NAME_TRACE "trace"
#define LOG_NAME_DEBUG "debug"
#define LOG_NAME_INFO "info"
#define LOG_NAME_WARN "warning"
#define LOG_NAME_ERROR "error"
#define LOG_NAME_CRITICAL "critical"

#define LOG_COLOR_TRACE "\033[37m"
#define LOG_COLOR_DEBUG "\033[36m"
#define LOG_COLOR_INFO "\033[32m"
#define LOG_COLOR_WARN "\033[33m\033[1m"
#define LOG_COLOR_ERROR "\033[31m\033[1m"
#define LOG_COLOR_CRITICAL "\033[1m\033[41m"
#define LOG_COLOR_OFF "\033[m"

#define LOG_CALL(level, name, color, format, ...)                         \
    do                                                                    \
    {                                                                     \
        if constexpr (level >= LOG_ACTIVE_LEVEL)                          \
        {                                                                 \
            if (isatty(fileno(stdout)))                                   \
            {                                                             \
                LOG_PRINT(color name LOG_COLOR_OFF, format, __VA_ARGS__); \
            }                                                             \
            else                                                          \
            {                                                             \
                LOG_PRINT(name, format, __VA_ARGS__);                     \
            }                                                             \
        }                                                                 \
    } while (0)

#define LOG_PRINT(name, format, ...)                                \
    std::printf("[" name "] [%s:%d] " format "\n",                  \
                std::filesystem::path(__FILE__).filename().c_str(), \
                __LINE__ __VA_OPT__(, ) __VA_ARGS__)

#endif