#pragma once
#ifndef UTILS_FILEIO_H
#define UTILS_FILEIO_H

#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#include "Shlobj.h"
#include<time.h>
#ifndef UNICODE  
typedef std::string String;
#else
typedef std::wstring String;
#endif
#else
typedef std::string String;
std::string dir_char = "/";
#endif

namespace FileIO
{
	template<typename T> bool read_csv(std::string fname, int ncols, std::string format_str, std::vector<T>& data)
	{
		T element;
		FILE * file_ptr = fopen(fname.c_str(), "r");
		char sep_c;
		
		std::string fmt = format_str;
		if (format_str.find("%c") == std::string::npos) fmt += "%c";

		if (!file_ptr)
		{
			printf("..error: failed to open file @(%s)\n", fname.c_str());
			return false;
		}

		do
		{
			for (int j = 0; j < ncols - 1; ++j)
			{
				int scan_res = fscanf_s(file_ptr, fmt.c_str(), &element, &sep_c);

				if (scan_res == 2)
				{
					data.push_back(element);
				}
				else
				{
					fclose(file_ptr);
					if (file_ptr) { delete file_ptr; file_ptr = 0; }
					return false;
				}

			}

			fscanf_s(file_ptr, format_str.c_str(), &element); 
			data.push_back(element);

		} while (file_ptr && !feof(file_ptr));

		fclose(file_ptr);
		if (file_ptr) { delete file_ptr; file_ptr = 0; }

		return true;
	}

};

#endif