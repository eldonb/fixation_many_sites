/* 2N diploid individuals in the population at the start; so 4N allele copies */
/* c++-10.2 -Wall -Wextra -Wsign-compare -march=native -m64 -O3 diploid_excursions_random_bottlenecks.cpp -lm -lgsl -lgslcblas */
/* https://www.nature.com/articles/s41467-021-21379-x */
/* https://academic.oup.com/mbe/article/17/1/32/975527 */
/* https://bioinformatics.cvr.ac.uk/calculating-dnds-for-ngs-datasets/ */
/* https://academic.oup.com/mbe/article/37/8/2450/5804989 */
/* https://github.com/danny-wilson/genomegaMap/releases */
/* https://academic.oup.com/bib/article-abstract/22/5/bbaa431/6105943 */
const double GLOBAL_CONST_N = 5000. ;
/* maximum number of alleles at one chromosome */
const double GLOBAL_CONST_A = 4. * GLOBAL_CONST_N ;
/* starting number  of diploid individuals */
const double GLOBAL_CONST_I = 2. * GLOBAL_CONST_N ;
const unsigned GLOBAL_CONST_II = 2 * static_cast<unsigned>( GLOBAL_CONST_N);
const int GLOBAL_NUMBER_SITES = 2 ;
const double GLOBAL_NUMBER_SITES_d = static_cast<double>( GLOBAL_NUMBER_SITES) ;
const int GLOBAL_NUMBER_TYPES = static_cast< int>( pow(3, GLOBAL_NUMBER_SITES_d) );
const unsigned int GLOBAL_CONST_CUTOFF_ONE = static_cast<unsigned int>( 2*GLOBAL_CONST_N) ;
const unsigned int GLOBAL_CONST_CUTOFF_TWO = static_cast<unsigned int>( 2*GLOBAL_CONST_N) ; 
const double GLOBAL_CONST_PSI_ONE = static_cast<double>(GLOBAL_CONST_CUTOFF_ONE) ;
const double GLOBAL_CONST_PSI_TWO = static_cast<double>(GLOBAL_CONST_CUTOFF_TWO) ;
const double GLOBAL_CONST_ALPHA_ONE = 0.75 ;
const double GLOBAL_CONST_ALPHA_TWO = 3.0 ;
/* epsilon is probability of alpha one */
const double GLOBAL_CONST_EPSILON = 0.1 ;
const double GLOBAL_CONST_SELECTION  = 1. ;
const unsigned int GLOBAL_CONST_BOTTLENECK = 100 ;
const double GLOBAL_CONST_PROBABILITY_BOTTLENECK = 0.1 ;
const int GLOBAL_CONST_NUMBER_EXPERIMENTS = 25000 ;
const char GLOBAL_A = 'G' ;
