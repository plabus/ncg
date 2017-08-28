#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include "Random.h"


/* PCG RANDOM NUMBER GENERTOR FUNCTIONS: */

inline void pcg_setseq_64_step_r(struct pcg32_random_t* rng)
{
  rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_64 + rng->inc;
}

inline uint32_t pcg_output_xsh_rr_64_32(uint64_t state)
{
  return pcg_rotr_32(((state >> 18u) ^ state) >> 27u, state >> 59u);
}

inline uint32_t pcg_rotr_32(uint32_t value, unsigned int rot)
{
  return (value >> rot) | (value << ((- rot) & 31));
}

inline void pcg32_srandom_r(struct pcg32_random_t* rng,
    uint64_t initstate, uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg_setseq_64_step_r(rng);
  rng->state += initstate;
  pcg_setseq_64_step_r(rng);
}

inline uint32_t pcg32_random_r(struct pcg32_random_t* rng)
{
  uint64_t oldstate = rng->state;
  pcg_setseq_64_step_r(rng);
  return pcg_output_xsh_rr_64_32(oldstate);
}

inline uint32_t	pcg32_boundedrand_r(struct pcg32_random_t* rng, uint32_t bound)
{
  uint32_t threshold = -bound % bound;
  for (;;) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold)
      return r % bound;
  }
}

// Return random unsigned double in range [0, 1) uniformly distributed
inline double uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    )
{
  return ldexp( pcg32_random_r(rng), -32 );
}

// Return random uniformly distributed int in range [0, N)
inline int uniform_int(
    struct pcg32_random_t* rng, // pointer to random number generator
    const int upper
    )
{
  return pcg32_boundedrand_r( rng, upper );
}

// Return random signed double in range (-1, 1) uniformly distributed
inline double signed_uniform(
    struct pcg32_random_t* rng // pointer to random number generator
    )
{
  return ( uniform_int(rng, 2) ? -1 : 1 ) * uniform( rng );
}
