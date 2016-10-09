#pragma once

#ifndef SYSTEM_TYPES_H
#define SYSTEM_TYPES_H

#ifdef __APPLE__
#include <stdint.h>
#elif __linux
#include <stdint.h>
#elif __unix 
#include <stdint.h>
#elif __posix
#include <stdint.h>
#elif _WIN32
#include <cstdint>

typedef unsigned int uint;

typedef __int8 int8t; // int8_t collides with signed char
typedef unsigned __int8 uint8_t;
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

typedef int8t int8;
typedef uint8_t uint8;
typedef const uint8_t cuint8;
typedef int16_t int16;
typedef const uint16_t cuint16;
typedef uint16_t uint16;
typedef int32_t int32;
typedef const uint32_t cuint32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef const uint64_t cuint64;
typedef uint64_t uint64;
#endif
#endif