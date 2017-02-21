#include "../utils/log.h"

using namespace StringUtils;

namespace Log
{
	Mutex log_mutex;
	FILE * log_ptr = 0;
}

bool Log::write(const char * msg, ...)
{
	log_mutex.lock();
	if (!log_ptr) return false;

	va_list args;
	va_start(args, msg);
	char dt_msg[1024];

	vsnprintf(dt_msg, 1024, msg, args);
	fprintf(log_ptr, "[%s]: %s", StringUtils::HHmmSS(".").c_str(), dt_msg);
	fflush(log_ptr);
	va_end(args);

	log_mutex.unlock();
	return true;
}

void Log::close()
{
	log_mutex.lock();
	if (log_ptr) fclose(log_ptr);
	log_mutex.unlock();
}

//#ifdef _WIN32
void Log::init(TCHAR * app_folder)
{
	TCHAR local_app_path[MAX_PATH];

	if (SUCCEEDED(SHGetFolderPath(NULL,
		CSIDL_LOCAL_APPDATA | CSIDL_FLAG_CREATE,
		NULL,
		0,
		local_app_path)))
	{
		if (null_or_ws(local_app_path))
		{
			printf("..invalid path passed for logfile!\n");
			return;
		}

		// append the requested app_folder (application defined) to the local-app path
		wcscat_s(local_app_path, MAX_PATH, app_folder);

		// check if we need to create the application defined folder (if so create it)
		int errorcode = GetFileAttributes(local_app_path);
		if (errorcode == INVALID_FILE_ATTRIBUTES)
		{
			printf("..attributes check for logging filepath failed (error-code %d)\n", errorcode);
			printf("..attempt to create logging filepath from scratch\n");
			errorcode = SHCreateDirectoryEx(NULL, local_app_path, NULL);
			if (errorcode != ERROR_SUCCESS)
			{
				printf("..attempt to create log folder failed (error-code %d)\n", errorcode);
				return;
			}
			else printf("..created log directory path successfully\n");
		}
		// convert wchar array "local_app_path" to standard char* for fopen
		char app_path[256];
		wcstombs(app_path, local_app_path, wcslen(local_app_path) + 1);

		char fname[256];
		sprintf(fname, "%s\\%s-log.txt", app_path, yyyyMMdd_HHmmSS().c_str());

		log_ptr = fopen(fname, "a");
		if (log_ptr) printf("..created log file @ (%s)\n", fname);
		else
		{
			printf("..failed to open log file!\n");
			return;
		}
	}
}
//#endif