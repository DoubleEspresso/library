#pragma once

#include <iostream>
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#include<time.h>
//std::string dir_char = "\\";
#ifndef UNICODE  
typedef std::string String;
#else
typedef std::wstring String;
#endif
#else
typedef std::string String;
std::string dir_char = "/";
#endif

namespace StringUtils
{
	bool only_ws(std::wstring test);
	bool only_ws(std::string test);
	bool null_or_ws(std::wstring test);
	bool null_or_ws(std::string test);
	std::string yyyyMMdd_HHmmSS(std::string sep = "");
	std::string HHmmSS(std::string sep = "");
};