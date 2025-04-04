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
#include <string>
#include <cstring>

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

#define log_error(fmt, ...) do { \
    fprintf(stderr, "[%s]" RED "ERROR: " fmt RESET " at %s:%i\n",       \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
    exit(1);                                                                      \
} while(0)

#define log_info(fmt, ...) do { \
    fprintf(stderr, "[%s]" GRN "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
} while(0)

#define log_debug(lvl, fmt, ...) do { if (lvl[3] <= LOGLVL[3]) { \
    fprintf(stderr, "[%s]" lvl "INFO: " fmt RESET "\n", time_str(), ##__VA_ARGS__); \
}} while(0)

#define sitrep(fmt, ...) do { \
    fprintf(stderr, "\r" MAG "STATUS: " fmt RESET "", ##__VA_ARGS__); \
    fflush(stderr); \
} while(0)

#define stderrflush fprintf(stderr, "\n")

#define log_warn(fmt, ...) do { \
    fprintf(stderr, "[%s]" YEL "WARN: " fmt RESET " at %s:%i\n",       \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);     \
} while(0)

#define prompt(fmt, ...) ({ \
    fprintf(stderr, "" GRN "" fmt RESET ". Hit enter to continue..", ##__VA_ARGS__); \
    while( getchar() != '\n' ); \
})

#define alignup(n, a) (((n) + (a)-1) & ~((a)-1))

#ifdef NDEBUG
#define verify(expression)
#else
#define verify(expression) if (SANITY_CHECKS && !(expression)) log_error("Expected " #expression "")
#endif

#define expect(expression) if (!(expression)) log_error("Expected " #expression "")

/// aliases and typedefs
typedef uint8_t u1;
typedef uint16_t u2;
typedef uint32_t u4;
typedef uint64_t u8;

/** memory allocation/deallcoation utils */
#define KiB <<10u
#define MiB <<20u
#define GiB <<30u

/** generic macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define alignup(n, a) (((n) + (a)-1) & ~((a)-1))

/** string macros */
#define streq(s1, s2) (!std::strcmp((s1), (s2)))
#define str_endswith(str, suffix) (streq((suffix), (str) + (strlen(str) - strlen(suffix))))
#define str_startswith(str, prefix)

#define BATCH_SZ 4096

#endif //COLLINEARITY_PRELUDE_H
