#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <bitset>
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include "Lsitesalllinked.hpp"


gsl_rng * rngtype ; 
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


/* generate a lookup table where entry with index i contains the corresponding  haplotype indexes */
static void generate_lookup_table( std::vector< std::pair<unsigned long, unsigned long> >& table)
{
  table.clear();
  table.shrink_to_fit();
  assert( table.size() < 1 );
  
    for( unsigned long i = 0 ; i <= GLOBAL_CONST_MAX_INDEX ; ++i){
    for( unsigned long j = i ; j <= GLOBAL_CONST_MAX_INDEX ; ++j){
      table.push_back( std::make_pair(i,j) );
    }}
}

/* check if parent haplotypes recombine */
static int recombination(  )
{
  /* if recombination return a discrete uniform x  between  1 and  L - 1 */
  /* recombination then occurs between x -1 and x */
  return ( gsl_rng_uniform(rngtype) <=  (GLOBAL_CONST_RECOMBINATION * (GLOBAL_CONST_NUMBER_SITESd - 1)) ? gsl_rng_uniform_int(rngtype, GLOBAL_CONST_NUMBER_SITES -1) + 1 :  GLOBAL_CONST_NUMBER_SITES + 2 ) ; 

  /* gsl_rng_uniform_int( rngtype, n) returns a random int between 0 and n-1  */
  /* shifting by one returns an integer between 1 and L - 1 */
}

/* generate recombinant haplotypes given recombination occurs  */
/* return one of the two recombined haplotypes picked with equal probability */
static unsigned long recombine( const int lrec, const unsigned long indexhapone, const unsigned long indexhaptwo )
{
  /* indexhapone and indexhaptwo are the indexes for the two parent haplotypes  */
  /* recombination happens between site lrec - 1 and lrec */
  /* for lrec between 1 and L - 1 ; L is number of sites */

  std::string hapone =  std::bitset<GLOBAL_CONST_NUMBER_SITES>(indexhapone).to_string() ;
  std::string haptwo =  std::bitset<GLOBAL_CONST_NUMBER_SITES>(indexhaptwo).to_string() ;
  std::string tmp = hapone ;
  
  /* generate recombinant hapone [0, ..., l-1] - [l, ... , L-1] hap two  */

  hapone.replace( lrec , GLOBAL_CONST_NUMBER_SITES - lrec , haptwo.substr( lrec,  GLOBAL_CONST_NUMBER_SITES - lrec ) ) ;
  haptwo.replace( lrec , GLOBAL_CONST_NUMBER_SITES - lrec , tmp.substr( lrec, GLOBAL_CONST_NUMBER_SITES - lrec )) ;  
  
  /* return the index of the recombinant haplotype picked for the juvenile */
  return(  gsl_rng_uniform(rngtype) < 0.5 ? std::bitset<GLOBAL_CONST_NUMBER_SITES>(hapone).to_ulong() : std::bitset<GLOBAL_CONST_NUMBER_SITES>(haptwo).to_ulong() ) ;  
}

/* sample haplotype index from a parent given haplotype indexes of the parent */
static unsigned long sample_haplotype_index( const unsigned long hone, const unsigned long htwo )
{
  /* hone and htwo are the haplotype indexes of the two parents */
  const int x = recombination() ;
  /* if x < number of sites then  recombination occurs  */
  /* otherwise  return the index of  one or the other of the parent haplotypes  */
  return( x < GLOBAL_CONST_NUMBER_SITES ?  recombine(x, hone, htwo) : (gsl_rng_uniform(rngtype) < 0.5 ? hone : htwo) ) ;
}


/* going along rows in the array of number of individuals of each */
/* diploid type; */
static unsigned long lookup( const unsigned long i, const unsigned long j)
{
  /* i and j are the two haplotype indexes */
  /*  i <= j */
  /* max index is the index of "1...1" */
  return ( j  +  (i*GLOBAL_CONST_MAX_INDEX) -  (i*(i-1)/2) ) ;
}


/* initializing the array of number of diploid individuals */
/* of each phased  L-site  type */
static void initializearray( std::vector<unsigned long>& p)
{
  /* p is the population; each entry is the number of individuals of the respective phased L-site type */
  std::string s (GLOBAL_CONST_NUMBER_SITES, '0') ;

  std::fill( std::begin(p), std::end(p), 0);
  assert( std::accumulate(std::begin(p), std::end(p) , 0) < 1 );

  
  for( int i = 0 ; i < GLOBAL_CONST_NUMBER_SITES ; ++i){
    s[i] = '1' ;
    // std::cout << s << '\n';
    p[ lookup(0, std::bitset<GLOBAL_CONST_NUMBER_SITES>(s).to_ulong() )] =  1 ;
    s[i] = '0' ; }

  /* Nsites number of diploid individuals carry a fit type at one site, i.e. are heterozygous at one site only */
  /* and at all other sites homozygous for the wild type;  */
  /* all other diploid individuals are homozygous for the wild type at all sites */
  p[0] = GLOBAL_CONST_CARRYING_CAPACITY - GLOBAL_CONST_NUMBER_SITES ;
}

static unsigned long current_number_individuals( const std::vector<unsigned long>& p)
{
  return  std::accumulate( std::begin(p), std::end(p), 0 ) ;
}

/* compute weight of juvenile given haplotype indexes */
static double weight( const std::vector< unsigned long>& h)
{
  /* h contains haplotype indexes  of juvenile */

  const std::string hapone = std::bitset<GLOBAL_CONST_NUMBER_SITES>( h[0]).to_string() ;
  const std::string haptwo = std::bitset<GLOBAL_CONST_NUMBER_SITES>( h[1]).to_string() ;
  double w = 0 ;
  /* need to be homozygous for fit type to increase weight by one */
  for( int i = 0 ; i < GLOBAL_CONST_NUMBER_SITES ; ++i){
    w += ( hapone[i] == '1' ? (haptwo[i] == '1' ? 1. : 0) : 0.) ; }

  /* adding a small deviation not necessary :  +  gsl_ran_gaussian_ziggurat( rngtype, 0.1) */
  return ( gsl_ran_exponential( rngtype, 1./( 1. +  (GLOBAL_CONST_SELECTION * w)  ) ) ) ;
  /* random uniform only temporary */
  // return gsl_rng_uniform( rngtype);
}

/* add juvenile to vector of juveniles */
static void add_juvenile( const std::vector< unsigned long>& pone, const std::vector<unsigned long>& ptwo, std::vector< std::pair< std::vector<unsigned long>, double>>& vj, std::vector<double>& v_juvenileweights )
{
  /* pone and ptwo are the pairs of haplotype indexes for the two parents */

  std::vector< unsigned long> h {} ;
  h.push_back( sample_haplotype_index( pone[0], pone[1] ) ) ;
  h.push_back( sample_haplotype_index( ptwo[0], ptwo[1] ) ) ;

  const double w = weight(h);
  /* record the weight of the juvenile for computing the nth smallest */
  v_juvenileweights.push_back( w) ;
  /* h now contains the indexes of the two haplotypes */
  /* compute the weight of the juvenile given h */
  vj.push_back( std::make_pair( h, w ));
}

static double masskernel( const double k, const double alpha)
{
  return ( pow( 1./k, alpha) - pow( 1./(k + 1.), alpha) ) ;
}

/* generate the CDF for sampling random number of juveniles */
static void cdf_number_juveniles( std::vector<double>& x , const double a )
{
  // const double cconst = 1./( pow(0.5, GLOBAL_CONST_ALPHA) -  pow( 1./( GLOBAL_CONST_CUTOFFd + 1.), GLOBAL_CONST_ALPHA) ) ;
  for( unsigned long k = 2; k <=  GLOBAL_CONST_CARRYING_CAPACITY ; ++k){
    x[k] =   ( masskernel( static_cast<double>(k), a) ) ; }

  const double cconst  =  std::accumulate( std::begin( x), std::end(x), 0.);

  for( unsigned long k = 2; k <=  GLOBAL_CONST_CARRYING_CAPACITY ; ++k){
    x[k] = x[k-1] + ( masskernel( static_cast<double>(k), a ) / cconst ) ; }

  /* the cdf must have last element one to guarantee a value within limits is sampled */
  x.back() = 1. ;

  // for (const auto &p : x){ std::cout << p << '\n';}
    
}

/* sample a random number of juveniles for one family */
static int sample_litter_size( const std::vector<double>& x)
{
  int j = 2 ;
  const double u = gsl_rng_uniform(rngtype);
  while( u > x[j]){ ++j ;}
  
  assert( j > 1);
  assert( j <= static_cast< int>( GLOBAL_CONST_CUTOFF) ) ;

  return j ;
}

/* add litter to pool of juveniles */
static void add_litter( const std::vector<double>& v_cdf, const std::vector<unsigned long>&  po, const std::vector<unsigned long>& pt, std::vector<  std::pair< std::vector<unsigned long>, double>>& pool, std::vector< double>& vw  )
{
  /* sample litter size */
  const int littersize = sample_litter_size(v_cdf) ;
  for( int j = 0 ; j < littersize; ++j){
    add_juvenile( po, pt, pool, vw);}
}

/* sample one  index of phased L-site type of diploid parent */ 
static unsigned long sample_genotype_parent( std::vector<unsigned long>& p )
{
  /* p is population */ 
  unsigned long i  = 0 ;
  unsigned int nothers = static_cast< unsigned int>( current_number_individuals(p) - p[i]) ;
  unsigned int x =  gsl_ran_hypergeometric( rngtype, p[0], nothers, 1);

  while( (x < 1) && (i <  GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES ) ){
    ++i ;
    nothers -= p[ i];
    x =  gsl_ran_hypergeometric(rngtype, p[i], nothers, 1);
  }
  /* check if an individual has been sampled */
  i += (x < 1 ? 1 : 0); 
  /* adjust the number of remaining parents */
  /* an individual of type with index i sampled, so subtract one from the number of remaining individuals with same type */
  --p[i] ;
  /* return the index  of the diploid  phased L-site type of the parent */
  /* index is between 0 and GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES */
  return i ;
}

static void clearpooljuveniles( std::vector< std::pair< std::vector< unsigned long>, double>>& x  )
{
  for( auto &y : x){
    std::get<0>(y).clear() ;
    std::get<0>(y).shrink_to_fit() ;
    std::vector< unsigned long>().swap( std::get<0>(y) ) ; }

  x.clear() ;
  x.shrink_to_fit() ;
  std::vector< std::pair< std::vector< unsigned long>, double>>().swap(x) ;
  assert( x.size() < 1 ) ; 
}

/* generate a new pool of juveniles */
static void new_pool_juveniles( const unsigned long Nt,  std::vector<  std::pair< std::vector<unsigned long>, double>>& pool, std::vector< double>& vw, std::vector<unsigned long>& population,  const std::vector<double>& vcdf, const std::vector< std::pair< unsigned long, unsigned long>>& table )
{
    /* Nt is current number of individuals */
  /* pool is the vector for the new pool of juveniles */
  clearpooljuveniles( pool ) ;

  vw.clear();
  vw.shrink_to_fit() ;
  assert( vw.size() < 1);

  unsigned long indexpone {} ;
  unsigned long indexptwo {} ;

  std::vector<unsigned long> pone (2,0);
  std::vector<unsigned long> ptwo (2,0);
  
  /* Nt is current number of individuals */
  /* each time sample two parents */
  for( unsigned long i = 0 ; i < Nt/2 ; ++i){
    /* sample index of a phased Lsite-type of parent one */
    /* table[indexpone] contains the corresponding haplotype indexes */
    indexpone = sample_genotype_parent( population );
    pone[0] = std::get<0>(table[ indexpone]) ;
    pone[1] = std::get<1>(table[ indexpone]) ;
    /* table[indexptwo] contains the corresponding haplotype indexes */
    indexptwo = sample_genotype_parent( population ); 
    ptwo[0] = std::get<0>(table[ indexptwo]) ;
    ptwo[1] = std::get<1>(table[ indexptwo]) ;
    add_litter( vcdf, pone, ptwo, pool, vw ) ;
  }
}

/* given a site for which to check, count how many diploid individuals are homozygous for the wild type at the site */
static unsigned long checknullsite( const unsigned long nsite,  const std::vector< unsigned long>& population, const std::vector< std::pair<unsigned long, unsigned long>>& tafla )
{
  /* nsite records the number of individuals homozygous for the wild type for each site */
  unsigned long sumof = 0 ;
  for( unsigned long i = 0 ; i < GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES ; ++i){
    /* if the type is  homozygous for the wild type at site nsite */
    /* add the  number of individuals with the phased L-site type */
    sumof +=  (std::bitset<GLOBAL_CONST_NUMBER_SITES>( std::get<0>(tafla[i]) ).to_string()[nsite] == '0' ? (std::bitset<GLOBAL_CONST_NUMBER_SITES>( std::get<1>(tafla[i]) ).to_string()[nsite] == '0' ?  population[i] : 0) : 0 ) ; 
  }
  /* return the number of individuals homozygous for the wild type at site nsite */
  return sumof ; 
}

static bool comp( const double a, const double b )
{
  return a < b ;
}

static double nthelm(  std::vector< double>& weights )
{
  std::nth_element( weights.begin(), weights.begin() + (GLOBAL_CONST_CARRYING_CAPACITY - 1), weights.end(), comp ) ;
  return weights[ GLOBAL_CONST_CARRYING_CAPACITY - 1] ;
}

/* sort juveniles and update population  */
static void  select_juveniles_according_to_weight( std::vector<double>& v_weights,  std::vector< unsigned long>& v_population, const std::vector< std::pair< std::vector< unsigned long>, double>>& v_pool )
{
  /* compute the nth smallest viability  weight */
  /* only sample juveniles according to weight if the total number of juveniles exceeds the carrying capacity */
  /* otherwise  all juveniles survive */
  const double wnth =  ( static_cast<unsigned long>( v_pool.size()) > GLOBAL_CONST_CARRYING_CAPACITY ? nthelm( v_weights) :  std::max_element( std::begin(v_weights), std::end(v_weights) )[0] + 1. )    ;

  assert( wnth > 0.) ;
  /* set number of all phased L-site types to zero */
  std::fill( std::begin( v_population), std::end(v_population), 0); 
  assert( std::accumulate( std::begin(v_population), std::end(v_population), 0.) < 1);
  /* add to count of phased L-site type of juvenile if juvenile survives */
  for( const auto &j : v_pool){
    v_population[ lookup( std::get<0>(j)[0], std::get<0>(j)[1]) ] += (  std::get<1>(j) <= wnth ? 1 : 0) ; }

  /* without bottlenekcs */
  // assert( current_number_individuals( v_population) == GLOBAL_CONST_CARRYING_CAPACITY) ;
}

/* count the number of individuals homozygous for the wild type at each site  (ns) */
/* if at any site s we have  ns = Nt where Nt is the current number of individuals then  lost type at the site */
static bool not_lost_type( std::vector<unsigned long>& nhomozygous ,  const std::vector< unsigned long>& population, const std::vector< std::pair<unsigned long, unsigned long>>& mtafla )
{

  std::fill( std::begin(nhomozygous), std::end(nhomozygous), 0);
  const unsigned long Nt =  current_number_individuals( population) ; 
  for( unsigned long nullsite = 0 ; nullsite < GLOBAL_CONST_NUMBER_SITES ; ++nullsite){
    nhomozygous[nullsite] =  checknullsite( nullsite, population, mtafla) ; }

  return std::all_of( std::begin(nhomozygous), std::end(nhomozygous), [Nt](unsigned long n){return n < Nt;} );
}

/* remove individuals not surviving a bottleneck */
static unsigned long remove_not_surviving_bottleneck(std::vector< unsigned long>& v_p)
{
  // std::vector<unsigned long> surviving ( GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES, 0);

  unsigned long x {} ;
  const unsigned long currentNt =   current_number_individuals( v_p ); 
  for ( unsigned long i = 0 ; i < currentNt - GLOBAL_CONST_BOTTLENECK ; ++i){
    x = sample_genotype_parent( v_p) ; }

  assert( current_number_individuals( v_p) == GLOBAL_CONST_BOTTLENECK );

  return current_number_individuals( v_p)  ;
  // v_p = surviving ;
}



static void bbb()
{
  std::vector< std::pair< std::vector<unsigned long>, double>> l {} ;
  std::vector< unsigned long> v (2,0); 
  l.push_back( std::make_pair(v,0) ) ;
  std::cout << std::get<0>( l[0])[0] << '\n' ;
}

static void prufa()
{
  std::vector< unsigned long> pop ( GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES, 0);

  initializearray( pop );

  std::cout << pop[ lookup( GLOBAL_CONST_MAX_INDEX, GLOBAL_CONST_MAX_INDEX) ] << '\n' ;

  const unsigned long i = 45 ;

  std::cout << i / 2 << '\n' ;

  for( unsigned long i = 0 ; i <= GLOBAL_CONST_MAX_INDEX ; ++i){
    for( unsigned long j = i ; j <= GLOBAL_CONST_MAX_INDEX ; ++j){
      std::cout <<   std::bitset<GLOBAL_CONST_NUMBER_SITES>(i).to_string() << ' ' << std::bitset<GLOBAL_CONST_NUMBER_SITES>(j).to_string() <<  ' ' << i << ' ' << j  << ' '  << lookup(i,j) << ' '   << pop[lookup(i,j)] <<    '\n' ; 
    }}

  std::vector< unsigned long> nhomozygous (GLOBAL_CONST_NUMBER_SITES, 0); 
    std::vector< std::pair< unsigned long, unsigned long>> m {} ;
  m.clear() ;
  generate_lookup_table(m);
  std::cout << not_lost_type(nhomozygous,  pop, m ) << '\n'  ;

  /* *************
  std::cout << sample_genotype_parent( pop  ) << '\n';
  for( unsigned long i = 0 ; i <= GLOBAL_CONST_MAX_INDEX ; ++i){
    for( unsigned long j = i ; j <= GLOBAL_CONST_MAX_INDEX ; ++j){
      std::cout << std::bitset<GLOBAL_CONST_NUMBER_SITES>(i).to_string() << ' ' << std::bitset<GLOBAL_CONST_NUMBER_SITES>(j).to_string() <<  ' ' << pop[lookup(i,j)] <<    '\n' ; 
    }}
    ******************* */
}

static void run( )
{

  /* define the population */
   std::vector< unsigned long> pop ( GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES, 0);

   std::vector< unsigned long> nhomozygouswt (GLOBAL_CONST_NUMBER_SITES, 0); 
   std::vector< std::pair< unsigned long, unsigned long>> m {} ;
   m.clear() ;
   generate_lookup_table(m);


   std::vector< std::pair< std::vector< unsigned long>, double>> pooljuveniles {} ;
   pooljuveniles.clear() ;
   std::vector< double> viabilityweights {} ;
   viabilityweights.clear() ;

   std::vector< double> v_cdf_number_juveniles_one (GLOBAL_CONST_CUTOFF + 1, 0) ;
   std::vector< double> v_cdf_number_juveniles_two (GLOBAL_CONST_CUTOFF + 1, 0) ;
   cdf_number_juveniles( v_cdf_number_juveniles_one, GLOBAL_CONST_ALPHA_ONE );
   cdf_number_juveniles( v_cdf_number_juveniles_two, GLOBAL_CONST_ALPHA_TWO );
   std::vector<double> excursion {} ;
   unsigned long currentNt  {} ;
   int trials = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1 ;
   int timi {} ;
   while( --trials > 0){ 
     initializearray( pop );
     timi = 0 ;
     excursion.clear() ;
     assert( excursion.size() < 1) ;
     pooljuveniles.clear() ;
     pooljuveniles.shrink_to_fit() ;
     assert( pooljuveniles.size() < 1) ;
     while( (pop.back() < current_number_individuals( pop)) & not_lost_type(nhomozygouswt, pop, m) ){
       ++ timi ;
       currentNt =  current_number_individuals( pop) ;
       excursion.push_back( static_cast<double>( pop.back() ) / static_cast<double>(currentNt) ) ;

       if( gsl_rng_uniform(rngtype) < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
	 currentNt = remove_not_surviving_bottleneck( pop);}
       
       new_pool_juveniles( currentNt,  pooljuveniles,  viabilityweights, pop, (gsl_rng_uniform(rngtype) < GLOBAL_CONST_EPSILON ? v_cdf_number_juveniles_one : v_cdf_number_juveniles_two), m);
       select_juveniles_according_to_weight( viabilityweights, pop, pooljuveniles ) ; }
     std::cout << (pop.back() < current_number_individuals(pop) ? 0 : 1) << ' ' << timi << '\n' ;
     if ( pop.back() == current_number_individuals(pop) ){
       /* fixation occurs so record the excursion */
       std::ofstream outfile ("tmpexcurs", std::ios_base::app) ;
       for( const auto &x : excursion){
	 outfile << x << ' ' ; }
       outfile << '\n' ;
       outfile.close() ;
     }
   }
}

int main(int argc, char * argv[])
{
  
  setup_rng( static_cast<unsigned long>(atoi(argv[argc - 1])) );
  run () ;
  gsl_rng_free( rngtype) ;
  return 0 ;
}
