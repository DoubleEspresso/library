#pragma once
#ifndef UTILS_DEFINITIONS_H
#define UTILS_DEFINITIONS_H

// dll export definitions for win32
#ifdef _WIN32
#define DllExport __declspec(dllexport)
#define DECL __cdecl
#else
#define DllExport
#define LIB
#define DECL
#endif

#endif