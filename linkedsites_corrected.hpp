const unsigned long GLOBAL_CONST_N = 5000 ;
const unsigned long GLOBAL_CONST_CARRYING_CAPACITY = 2*GLOBAL_CONST_N ;
const unsigned long GLOBAL_CONST_CUTOFF = 2*GLOBAL_CONST_N;
const double GLOBAL_CONST_CUTOFFd = static_cast<double>(GLOBAL_CONST_CUTOFF) ;
const double GLOBAL_CONST_ALPHA_ONE = 0.75 ;
const double GLOBAL_CONST_ALPHA_TWO = 3. ;
const double GLOBAL_CONST_RECOMBINATION = 0.75 ; 
const int GLOBAL_CONST_NUMBER_SITES = 2 ;
//const std::vector<double> GLOBAL_CONST_RECOMBINATION = {0.75, 1}; 
const double GLOBAL_CONST_SELECTION = 10.0 ;
/* ********* */
const double GLOBAL_CONST_EPSILON = -0.1 ;
const double GLOBAL_CONST_PROBABILITY_BOTTLENECK = -0.1 ;
/* ********* */
//const unsigned long GLOBAL_CONST_BOTTLENECK = 1000 ;
const double GLOBAL_CONST_NUMBER_SITESd = static_cast<double>(GLOBAL_CONST_NUMBER_SITES) ;
/* make sure the length of  string "11...1"  and number of sites match */
const std::string GLOBAL_CONST_ONES (GLOBAL_CONST_NUMBER_SITES, '1') ;
const unsigned long GLOBAL_CONST_MAX_INDEX =  std::bitset<GLOBAL_CONST_NUMBER_SITES>(GLOBAL_CONST_ONES).to_ulong() ;
/* total number of possible phased L-site types */
const unsigned long GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES = (GLOBAL_CONST_MAX_INDEX + 1)*(GLOBAL_CONST_MAX_INDEX + 2) / 2 ; 
const int GLOBAL_CONST_NUMBER_EXPERIMENTS = 1 ;
