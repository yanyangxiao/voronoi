
/************************************************************************/
/* author: Yanyang Xiao                                                 */
/* email : yanyangxiaoxyy@gmail.com                                     */
/************************************************************************/

#ifndef XLOG_H
#define XLOG_H

#include <cstdio>
#include <fstream>
#include <time.h>

#if defined(linux) || defined(__LYNX)
#define LOG_FONT_COLOR_DEFAULT "\033[0m"
#define LOG_FONT_COLOR_RED     "\033[31m"
#define LOG_FONT_COLOR_YELLOW  "\033[33m"
#define LOG_FONT_COLOR_BLUE    "\033[34m"
#else
#define LOG_FONT_COLOR_DEFAULT
#define LOG_FONT_COLOR_RED
#define LOG_FONT_COLOR_YELLOW
#define LOG_FONT_COLOR_BLUE
#endif

static char* get_current_time()
{
	static char s[64];
	time_t t;
	struct tm* ltime;

	time(&t);
	ltime = localtime(&t);
	strftime(s, 64, "%Y-%m-%d %H:%M:%S", ltime);

	return s;
}

#define xlog(format, ...)         printf("[%s] %s: " format "\n", get_current_time(), __FUNCTION__, ##__VA_ARGS__)
#define xlog_warning(format, ...) xlog(LOG_FONT_COLOR_YELLOW format LOG_FONT_COLOR_DEFAULT, ##__VA_ARGS__)
#define xlog_error(format, ...)   xlog(LOG_FONT_COLOR_RED format LOG_FONT_COLOR_DEFAULT, ##__VA_ARGS__)

#ifdef _DEBUG
#define xlog_debug(format, ...)   xlog(LOG_FONT_COLOR_BLUE format LOG_FONT_COLOR_DEFAULT, ##__VA_ARGS__)
#else
#define xlog_debug(format, ...)
#endif

static void xlog_file(const char *func, const char *what)
{
	std::ofstream logfile = std::ofstream("log.txt", std::ios::app);

	if (!logfile)
		return;

	logfile << "[" << get_current_time() << "]\t" << func << "\t" << what << std::endl;
	logfile.close();
}

#endif
