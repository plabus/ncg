#ifndef _RANDOM_H_
#define _RANDOM_H_

#define PCG_DEFAULT_MULTIPLIER_64  6364136223846793005ULL
struct pcg32_random_t {
    uint64_t state;
    uint64_t inc;
};
 
void pcg_setseq_64_step_r(struct pcg32_random_t* rng);
uint32_t pcg_output_xsh_rr_64_32(uint64_t state);
uint32_t pcg_rotr_32(uint32_t value, unsigned int rot);

void pcg32_srandom_r(struct pcg32_random_t* rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_random_r(struct pcg32_random_t* rng);
uint32_t pcg32_boundedrand_r(struct pcg32_random_t* rng, uint32_t bound);

#endif
