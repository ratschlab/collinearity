//
// Created by Sayan Goswami on 27.11.2024.
//

#ifndef COLLINEARITY_PRELUDE_H
#define COLLINEARITY_PRELUDE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <stdexcept>
#include <list>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <cassert>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define RESET "\x1B[0m"

#define HIGH BLU
#define MED  MAG
#define LOW  CYN
#define LOGLVL LOW

static char *time_str(){
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    static char buf[9];
    strftime (buf, 9,"%T",timeinfo);
    return buf;
}

#define error(fmt, ...) do { \
    fprintf(stderr, "[%s]" RED "ERROR: " fmt RESET " at %s:%i\n",       \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
    assert(false);                                                                      \
} while(0)

#define info(fmt, ...) do { \
    fprintf(stderr, "[%s]" GRN "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
} while(0)

#define debug(lvl, fmt, ...) do { if (lvl[3] <= LOGLVL[3]) { \
    fprintf(stderr, "[%s]" lvl "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
}} while(0)

#define sitrep(fmt, ...) do { \
    fprintf(stderr, "\r" MAG "STATUS: " fmt RESET "", ##__VA_ARGS__); \
    fflush(stderr); \
} while(0)

#define warn(fmt, ...) do { \
    fprintf(stderr, "[%s]" YEL "ERROR: " fmt RESET " at %s:%i\n",       \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
} while(0)

#define alignup(n, a) (((n) + (a)-1) & ~((a)-1))

#define expect(expression) if (!(expression)) error("Expected " #expression "")

#ifdef SANITY_CHECK
#define verify(expression) if (!(expression)) error("Expected " #expression "")
#else
#define verify(expression)
#endif

/// aliases and typedefs
typedef uint8_t u1;
typedef uint32_t u4;
typedef uint64_t u8;

typedef u4 KeyT;
typedef u8 ValT;

/** memory allocation/deallcoation utils */

#define KiB <<10u
#define MiB <<20u
#define GiB <<30u

#endif //COLLINEARITY_PRELUDE_H
