#include "../utils/string_utils.h"

bool StringUtils::only_ws(std::wstring test)
{
	for (int j = 0; j < test.length(); ++j) if (test[j] != L' ') return false;
	return true;
}

bool StringUtils::only_ws(std::string test)
{
	for (int j = 0; j < test.length(); ++j) if (test[j] != ' ') return false;
	return true;
}

bool StringUtils::null_or_ws(std::wstring test)
{
	return (test.length() <= 0 || only_ws(test));
}

bool StringUtils::null_or_ws(std::string test)
{
	return (test.length() <= 0 || only_ws(test));
}

std::string StringUtils::yyyyMMdd_HHmmSS(std::string sep)
{
#ifdef _WIN32
	SYSTEMTIME st;
	GetLocalTime(&st);
	char buf[256];
	sprintf(buf, "%d%s%d%s%d_%d%s%d%s%d",
		st.wYear, sep.c_str(),
		st.wMonth, sep.c_str(),
		st.wDay,
		st.wHour, sep.c_str(),
		st.wMinute, sep.c_str(),
		st.wSecond);
#endif
	return std::string(buf);
}

std::string StringUtils::HHmmSS(std::string sep)
{
#ifdef _WIN32
	SYSTEMTIME st;
	GetLocalTime(&st);
	char buf[256];
	sprintf(buf, "%d%s%d%s%d",
		st.wHour, sep.c_str(),
		st.wMinute, sep.c_str(),
		st.wSecond);
#endif
	return std::string(buf);
}