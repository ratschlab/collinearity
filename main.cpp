#include "collinearity.h"

#include "streamvbyte.h"

#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wdeclaration-after-statement"
#pragma clang diagnostic ignored "-Wunused-variable"
#endif

parlay::sequence<long> prime_sieve(long n) {
    if (n < 2) return parlay::sequence<long>();
    else {
        long sqrt = std::sqrt(n);
        auto primes_sqrt = prime_sieve(sqrt);
        parlay::sequence<bool> flags(n+1, true);  // flags to mark the primes
        flags[0] = flags[1] = false;              // 0 and 1 are not prime
        parlay::parallel_for(0, primes_sqrt.size(), [&] (size_t i) {
            long prime = primes_sqrt[i];
            parlay::parallel_for(2, n/prime + 1, [&] (size_t j) {
                flags[prime * j] = false;
            });
        }, 1);
        return parlay::filter(parlay::iota<long>(n+1),
                              [&](size_t i) { return flags[i]; });
    }
}

int main() {
    process_fasta("/scratch/Zymo/reads-tiny.fasta", 10, 4);
    return 0;
}
