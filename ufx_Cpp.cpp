/*        ***
simulate L unlinked sites in a diploid population evolving according to a model of random sweepstakes,
randomly occurring bottlenecks, and viability selection
Copyright (C) 2022  Bjarki Eldon
the R code for producing the table of L-site types
**********************
#!/usr/bin/R
m <- cbind( 0:26,   rep( 0:2, rep(9,3)), rep( 0:2, rep(3,3)), rep(0:2, 9)  )
print(m)
fall <- function(i,j,k){ (i*(3*3))  + 3*j + k}
rm(m)
m <- cbind(  rep( 0:2, rep( 27, 3)), rep( 0:2, rep(9, 3)), rep(  rep( 0:2, rep(3,3) ), 3), rep( 0:2, 27) )
g <- function(l, LL)
{
    x <- numeric()
    while( length(x) < (3**LL) ){
        x <- c(x,  rep(0,3**(LL-l)),  rep(1, 3**(LL-l)), rep(2, 3**(LL-l))) }
    return (x)
}
w <- 5
rm(m)
m <- matrix( 0, 3**w, w)
for(l in 1:w){
    m[, l]  <- t(g(l, w))}
print(m)
write(t(m), "types.txt", ncolumns =  w)
rm(m)
****************** 
w above is number of sites, produces all possible w-site types ordered in the
dictionary order
   */
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include "Lsites.hpp"


gsl_rng * rngtype ; 
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}

/* look up entry in population array given genotypes at each site */
static size_t lookup(const std::vector<short>& types )
  {
    size_t n = 0 ;
    for( int ell = 1 ; ell <= GLOBAL_NUMBER_SITES ; ++ell){
      n +=  types[ell - 1]*static_cast<short>( gsl_sf_pow_int(3, GLOBAL_NUMBER_SITES - ell) ) ; 
    }
    return (n) ;
  }



static void crep( const double LL,  const double l, const short n, std::vector<short>& y )
{
  y.clear();
  y.resize( pow(3, LL-l)) ;
  y.assign( pow(3, LL - l), n) ;

}


static void onesite( const double LL, const double l, std::vector< short>& x )
{
  x.clear();
  std::vector<short> tmp {};
  const size_t m =   static_cast<size_t>( pow( 3., LL)) ;
  while( x.size() < m){
    for( short i = 0 ; i < 3 ; ++i){
      crep( LL, l, i, tmp);
      x.insert( x.end(),  tmp.begin(), tmp.end() ) ;}
  }
}

static void writearray( const double LLw,  gsl_matrix_short * M)
{
  std::vector< short> d {} ;
  for( size_t k = 1; k <= static_cast<size_t>(LLw) ; ++k){
    onesite( LLw, k, d);
    /* copy into array */
    for( size_t u = 0 ; u < d.size() ; ++u){
      gsl_matrix_short_set(M, u, k-1, d[u]);}
  }
}



static void init_containers( std::vector< unsigned>& population, std::vector<double>& cdf_one, std::vector<double>& cdf_two)
  {

    population.clear() ;
    population.assign( GLOBAL_NUMBER_TYPES, 0);
 
    population[0] = GLOBAL_CONST_II - GLOBAL_NUMBER_SITES ;
   

    std::vector< short> types (GLOBAL_NUMBER_SITES, 0 );
    for( unsigned i = 0 ; i < GLOBAL_NUMBER_SITES ; ++i){
      std::fill( types.begin(), types.end(), 0);
      types[i] = 1 ;
      
      population[ lookup( types) ] = 1 ;}
    
    cdf_one.clear() ;
    cdf_two.clear() ;
    cdf_one.reserve( GLOBAL_CONST_CUTOFF_ONE + 2) ;
    cdf_two.reserve( GLOBAL_CONST_CUTOFF_TWO + 2) ;
    cdf_one.push_back( 0.);
    cdf_one.push_back( 0.);
    cdf_two.push_back(0.);
    cdf_two.push_back(0.);
    assert( cdf_one.size() == 2 );
    assert( cdf_two.size() == 2 );
  }


static void init_for_trajectory( std::vector< unsigned>& population)
  {
    std::fill( population.begin(), population.end(), 0);
    population[0] = GLOBAL_CONST_II - GLOBAL_NUMBER_SITES ;
    
     std::vector< short> types (GLOBAL_NUMBER_SITES, 0 );
    for( unsigned i = 0 ; i < GLOBAL_NUMBER_SITES ; ++i){
      std::fill( types.begin(), types.end(), 0);
      /* set one individual as heterozygous at site i and homozygous for wild type at all other sites */
      types[i] = 1 ;
      population[ lookup( types) ] = 1 ;}
  }


/* return the current total number of individuals in the population */
static  unsigned int current_number_individuals(const std::vector<unsigned>& population)
  {
    return std::accumulate( std::begin(population), std::end(population), 0); 
  }


/* sample index of L-site type of parent */ 
static  int sample_genotype_parent( std::vector<unsigned>& p,   gsl_rng *r )
{
  /* p is population */ 
  int i  = 0 ;
  unsigned int nothers = current_number_individuals(p) - p[i] ;
  unsigned int x =  gsl_ran_hypergeometric( r, p[0], nothers, 1);

  while( (x < 1) && (i < GLOBAL_NUMBER_TYPES ) ){
    ++i ;
    nothers -= p[ i];
    x =  gsl_ran_hypergeometric(r, p[i], nothers, 1);
  }
  /* check if an individual has been sampled */
  i += (x < 1 ? 1 : 0); 
  /* adjust the number of remaining parents */
  /* an individual of type with index i sampled, so subtract one from the number of remaining individuals with same type */
  --p[i] ;
  /* return the index  of the genotype of the parent */
  /* index is between 0 and GLOBAL_NUMBER_TYPES */
  return i ;
}


/* all juveniles survive */
static void update_population_all_juveniles(const std::vector< std::pair< size_t, double>>& juveniles, std::vector<unsigned>& population)
  {
    std::fill(population.begin(), population.end(), 0 ) ;
    assert( current_number_individuals(population) < 1 ); 
    for( const auto &j : juveniles){
      /* j[0] is type index of juvenile j */
      population[ std::get<0>( j ) ] += 1;}
  }



/* assign single site  genotype to juvenile given  genotypes in parents */
static short assign_type_one_site( const short gone, const short gtwo, gsl_rng *r)
{

  short int g {} ;
  const double u = gsl_rng_uniform(r) ;
  switch(gone){
  case 0 : {
    g = (gtwo < 1 ? 0 : (gtwo < 2 ? (u < 0.5 ? 0 : 1) : 1) );
    break ;}
  case 1 : {
    g = (gtwo < 1 ? (u < .5 ? 0 : 1) : (gtwo < 2 ? (u < 0.25 ? 0 : ( u < 0.75 ? 1 : 2)) : (u < 0.5 ? 1 : 2) ) ) ;
    break ; }
  case 2 : {
    g = (gtwo < 1 ? 1 : (gtwo < 2 ? (u < .5 ? 1 : 2) : 2) ) ;
    break ; }
  default : break ; } 

  assert( g == 0 || g == 1 || g == 2 ) ;
  
  return g ;
}


/* read in file with all types for given number of sites */
static void read_types( const std::vector<unsigned>& p ,  gsl_matrix_short * M  )
{
  std::ifstream f("types_five_sites.txt");

  short x {} ;
 
  for( int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( int j = 0 ; j < GLOBAL_NUMBER_SITES ; ++j){
      f >> x ; 
      gsl_matrix_short_set(M, i, j, x) ;}}
  f.close() ;

  /* print matrix for check */
  int z = 0 ;
   for( int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( int j = 0 ; j < GLOBAL_NUMBER_SITES ; ++j){
      std::cout << gsl_matrix_short_get(M, i,j) << ' ' ;}
    std::cout <<  p[z] <<   '\n' ;
    ++z ; }
}



/* sample a random number  of juveniles */
static size_t sample_random_number_juveniles( const size_t c_twoone, const std::vector<double>& cdfone, const std::vector<double>& cdftwo,   gsl_rng *r)
{
  const double u = gsl_rng_uniform(r);
  size_t j = 2 ;
  if( c_twoone < 2 ){
    while( u > cdfone[j] ){ ++j ;}}
  else{
    while( u > cdftwo[j] ){ ++j ;}}

  assert( j > 1) ;
  return j ;
}

/* new type confers advantage */
static double underdominant( const short z)
{
   /* z  is type at given site */
  return( (z < 1 ? 0.5 : (z < 2 ? .0 : 1.) ) /GLOBAL_NUMBER_SITES_d ) ;
}


static double recessive( const short z)
{
  /* z is type at given site */
  return( (z < 2 ? 0. : 1.) /GLOBAL_NUMBER_SITES_d ) ;
}


/* new type confers advantage */
static double semidominant( const short z)
{
   /* z  is type at given site */
  return( (z < 1 ? 0. : (z < 2 ? .5 : 1.) ) /GLOBAL_NUMBER_SITES_d ) ;
}


/* mutation confers advantage */
static double dominant( const short z )
{
  /* z  is type at given site */
  return( (z < 1 ? 0. : 1. )/GLOBAL_NUMBER_SITES_d ) ;
}


/* mutation confers advantage */
static double overdominant( const short z)
{
    /* z  is type at given site */ 
  return( (z < 1 ? 0. : ( z < 2 ? 1. : .5 ) ) / GLOBAL_NUMBER_SITES_d ) ;
}


/* mutation confers disadvantage */
static double negativedominant( const short z)
{
    /* z  is type at given site */
  return(  (z < 1 ? 1. : 0.)  / GLOBAL_NUMBER_SITES_d ) ;
}


/* mutation confers disadvantage */
static double negativeoverdominant( const short z)
{
    /* z  is type at given site */
  return ( (z < 1 ? .5 : ( z < 2 ? 1. : 0.))  / GLOBAL_NUMBER_SITES_d ) ;
}


static double negativerecessive(const short z)
{
  /* */
  return ( (z < 2 ? 0.5 : 0.) / GLOBAL_NUMBER_SITES_d );
}

/* trying out epistasis for two sites */
static double epistasis( const std::vector<short>& types, gsl_rng *r)
{
  /* types is the vector of types at each site of juvenile */
  double g = 0;
  assert( (types[0] == 0) || ((types[0] == 1) || (types[0] == 2)) ) ;
  
  switch( types[0] ){
  case 0 : {
    g = ( types[1] < 1 ? 1. : (types[1] < 2 ? 1.5 : 2.0) ) ;
      break ;}
  case 1 : {
    g =  (types[1] < 1 ? 0. : (types[1] < 2 ? .5 : 1.) );
      break ; }
  case 2 : {
    g =   (types[1] < 1 ? 0. : 1.) ;
      break ;}
  default : break ;}

  return ( gsl_ran_exponential( r,  1./( 1. + (GLOBAL_CONST_SELECTION*g) ) ) ) ;
}

static double half( const double x,  gsl_rng * rr)
{ 
  return ( gsl_rng_uniform(rr) < x ? GLOBAL_CONST_SELECTION : 0.) ;
}

static double computeweight(const std::vector<short>& types,  gsl_rng *r)
{
  double g = 1 ;
  /* p = 0.25  initially */
  const double p = 0.01 ; 

  for( int s = 0 ; s < GLOBAL_NUMBER_SITES ; ++s){
    switch(s){
    case 0 : {
      // g +=  (types[1] < 1 ? recessive( types[s]) : dominant(types[s])) +   (types[2] < -1 ? recessive( types[s]) : dominant(types[s])) +  (types[3] < -1 ? recessive( types[s]) : dominant(types[s])) ;
      // g+=  (types[1] < 1 ? overdominant(types[s]) : underdominant( types[s])) ;
      // g+=  (types[1] < 1 ? ( gsl_rng_uniform(r) < p ? underdominant(types[s]) : overdominant( types[s])) : (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s]) ) ) ;
      // g += dominant( types[s]) ;
      // g += types[1] < 2 ?  (gsl_rng_uniform(r) < 1. - p ? recessive( types[s]) : dominant(types[s])) :  (gsl_rng_uniform(r) < p ? recessive( types[s]) : dominant(types[s])) ;
      // g *= 1. +  (types[1] < 2 ?  (types[s] < 2 ? 0 :  GLOBAL_CONST_SELECTION) : (types[s] < 1 ? 0 : GLOBAL_CONST_SELECTION)) ;
      switch( types[1]){
      case 0 : {
	// g *=  1.  +  (types[1] < 2 ? (gsl_rng_uniform(r) < p ? GLOBAL_CONST_SELECTION : 0.) : (gsl_rng_uniform(r) < 1. - p ? GLOBAL_CONST_SELECTION : 0.)) ;
	//g *=  1.  +  (types[s] < 2 ? half(p,r) : half(1-p,r)) ;
	g *=  1. + (  GLOBAL_CONST_SELECTION *  (gsl_rng_uniform(r) < p ? underdominant( types[s]) : overdominant( types[s])) ) ;
	break ; }
      case 1 : {
	// g *=  1.  +   (types[s] < 2 ? half(p,r) : half(1-p,r)) ;
	g *=  1. +  (GLOBAL_CONST_SELECTION *  (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s])) ) ;
	break ; }
	case 2 : {
	  // g *=  1.  +  (types[s] < 1 ? half(p,r) : half(1.-p,r)) ;
	  g *= 1. +  (GLOBAL_CONST_SELECTION *   (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s])) ) ;
	  break ; }
	  default : break ;} 
      break ;}
    case 1 : {
      //  g +=  (types[0] < 1 ? recessive( types[s]) : dominant(types[s])) +   (types[2] < -1 ? recessive( types[s]) : dominant(types[s])) +  (types[3] < -1 ? recessive( types[s]) : dominant(types[s])) ;
      // g+=  (types[0] < 1 ? ( gsl_rng_uniform(r) < p ? underdominant(types[s]) : overdominant( types[s])) : (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s]) ) ) ;
      // g += recessive( types[s]) ;
      //  g += types[0] < 2 ?  (gsl_rng_uniform(r) < 1.-p ? recessive( types[s]) : dominant(types[s])) :  (gsl_rng_uniform(r) < p ? recessive( types[s]) : dominant(types[s])) ;
      // g *= 1. +  (types[0] < 2 ?  ( 0 ?  GLOBAL_CONST_SELECTION : (types[s] < 1 ? 0 : GLOBAL_CONST_SELECTION)) ;
      /* overdominance */
      switch( types[0]){
      case 0 : {
	// g *=  1.  +  (types[s] < 2 ? half(p,r) : half(1-p,r)) ;
	g *=  1. +   ( GLOBAL_CONST_SELECTION *  (gsl_rng_uniform(r) < p ? underdominant( types[s]) : overdominant( types[s])) ) ;
	break ; }
      case 1 : {
	// g *=  1.  +  (types[s] < 2 ? half(p,r) : half(1-p,r))  ;
        g *=  1.  +  ( GLOBAL_CONST_SELECTION *  (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s])) ) ;
	break ; }
	case 2 : {
	  // g *=  1.  +  (types[s] < 1 ? half(p,r) : half(1.-p,r))  ;
	  g *=  1.  +  ( GLOBAL_CONST_SELECTION *  (gsl_rng_uniform(r) < 1. - p ? underdominant(types[s]) : overdominant(types[s])) ) ;
	  break ; }
	  default : break ;}
      break ;}
    case 2 : {
      g +=  (types[1] < 1 ? recessive( types[s]) : dominant(types[s])) +   (types[0] < -1 ? recessive( types[s]) : dominant(types[s])) +  (types[3] < -1 ? recessive( types[s]) : dominant(types[s])) ;
      break ;}
       case 3 : {
	 g +=  (types[1] < 1 ? recessive( types[s]) : dominant(types[s])) +   (types[2] < -1 ? recessive( types[s]) : dominant(types[s])) +  (types[0] < -1 ? recessive( types[s]) : dominant(types[s])) ;
      break ;}
    default : break;}
  }
  //  return ( gsl_ran_exponential( r,  1./( 1. + (GLOBAL_CONST_SELECTION*g) )  ) ) ;
    return ( gsl_ran_exponential( r,  1./( g )  ) ) ;
}


/* compute weight two sites with epistasis */
static double computeweighttwosites( const std::vector<short>& types,  gsl_rng *r)
{
  return( dominant(types[0]) +  (types[0] > 0 ? recessive(types[1]) : negativerecessive(types[1])) ) ;
}

/* g is average of two diploid genotypes 0,1,2 */
static double computeweight_threesites(  const std::vector<short>& types,   gsl_rng *r)
{
  /* types is the vector of type at each site of juvenile  */

  double g = 0 ;
  
  /* compute the average contribution  over the sites */
  // GLOBAL_NUMBER_SITES; ++s){
  for (int s = 0; s < GLOBAL_NUMBER_SITES; ++s){
    /* there are L sites so if homoz of good type at all sites then returns one */
    /* if unfit type at all sites then returns  2L/L = 2 so weight is Exp( exp(-4s) ) */ 
    // g +=  (types[s] < 1 ? 0.5 :  (types[s] < 2 ? 0.0 : 1.) )/GLOBAL_NUMBER_SITES_d;
    // g +=  (types[s] < 2 ? 0. : 1 )/GLOBAL_NUMBER_SITES_d;
    // g += ( s < 1 ? recessive( types[s]) : ( types[s-1] < 2 ? ( gsl_rng_uniform(r) < .5 ? recessive(types[s]) : negativedominant(types[s]))  : dominant( types[s]) ) ) ;
    switch(s){
    case 0 : {
      g += (types[1] < 2 ? (types[2] < 2 ? (gsl_rng_uniform(r) < .5 ? recessive(types[0]) : dominant(types[0])) : dominant(types[0])) : (types[2] < 2 ? recessive(types[0]) : dominant(types[0])));
      break ;}
    case 1 : {
      g += (types[0] < 2 ? (types[2] < 2 ? (gsl_rng_uniform(r) < .5 ? recessive(types[s]) : dominant(types[s])) : dominant(types[s])) : (types[2] < 2 ? recessive(types[s]) : dominant(types[s])));
      break ;}
    case 2 : {
      g += (types[0] < 2 ? (types[1] < 2 ? (gsl_rng_uniform(r) < .5 ? recessive(types[s]) : dominant(types[s])) : dominant(types[s])) : (types[1] < 2 ? recessive(types[s]) : dominant(types[s])));
      break ; }
    default : break ;}
    // g += negativedominant( types[s] );
  }
  //  g = epistasis( types, r); 
  /* computing rate as Exp( exp(-sf(g)) ) 
   return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( g , 2.) ) ) ) ;
   **** 
   ***  compute rate as  1 + sf(g) if mutant confers advantage */
   return ( gsl_ran_exponential( r,  1./( 1. + (GLOBAL_CONST_SELECTION*g) )  )   ) ;
   // return( g );
}




static void add_juvenile(const int type_index_one, const int type_index_two, gsl_matrix_short * Mtypes, std::vector< std::pair<size_t, double>>& jvs,   gsl_rng * r)
{
  std::vector<short> types_juvenile {} ;
  types_juvenile.clear() ;

  /* get the types of the juvenile */
  for( int s = 0; s < GLOBAL_NUMBER_SITES ; ++s){
    /* sample and record the type at site s */
    types_juvenile.push_back( assign_type_one_site( gsl_matrix_short_get( Mtypes, type_index_one, s), gsl_matrix_short_get( Mtypes, type_index_two, s), r) ) ;}

  /* compute the weight and add the juvenile to the vector of juves */
  jvs.push_back( std::make_pair( lookup( types_juvenile), computeweight( types_juvenile, r) ) );
}


/* gone and gtwo are the type indexes for the two parents */
static void add_juveniles_for_given_parent_pair( const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,   const int gone, const int gtwo, const size_t conetwo,  gsl_matrix_short * Mtypes,   gsl_rng *r)
{
  /* first sample a random number of juveniles */
  const size_t numberj = sample_random_number_juveniles(  conetwo, cdfone, cdftwo, r);
  assert( numberj > 1 ) ;

  /* add the sampled number of juveniles to the pool */
  for( size_t j = 0; j < numberj ; ++j){
    add_juvenile( gone, gtwo,  Mtypes, jvs, r); }
}



static void generate_pool_juveniles( std::vector< std::pair<size_t, double>>& jvs, std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo,  gsl_matrix_short * Mtypes,  gsl_rng *r)
{
  jvs.clear() ;
  jvs.shrink_to_fit();
  assert( jvs.size() < 1 );
  int gone {} ;
  int gtwo {} ;
  /* sample distribution of number of juveniles */
  const size_t conetwo = (gsl_rng_uniform(r) < GLOBAL_CONST_EPSILON ? 1 : 2) ; 
  /* i runs over number of pairs that can be formed from the current number of individuals */
  const double currenti = current_number_individuals(p) ;  
  assert( currenti  < GLOBAL_CONST_I + 1) ;

  for ( double i = 0 ; i < floor( currenti / 2. ) ; ++i){
    gone = sample_genotype_parent(p, r);
    gtwo = sample_genotype_parent(p, r);
    /* gone and gtwo are the type indexes of the two parents */
    assert(gone > -1);
    assert( gtwo > -1);
    add_juveniles_for_given_parent_pair(cdfone, cdftwo, jvs, gone, gtwo, conetwo, Mtypes, r) ;}

  assert( jvs.size()  >= static_cast<size_t>( currenti ) ) ;
}


/* sample diploid individuals surviving a bottleneck */
static void sample_surviving_bottleneck( std::vector<unsigned>& p,   gsl_rng * r)
{
  int i = 0 ;

  unsigned int nothers = current_number_individuals(p) - p[i] ;
  unsigned newn = gsl_ran_hypergeometric( r,  p[i],  nothers,  GLOBAL_CONST_BOTTLENECK);
  unsigned int remaining = GLOBAL_CONST_BOTTLENECK - newn ;
  
  /* update count of individuals of type index i surviving bottleneck */
  p[i] = newn ;
  while( i < GLOBAL_NUMBER_TYPES - 2 ){
    ++i ;
    nothers -= p[i] ;
    newn = (remaining > 0 ? gsl_ran_hypergeometric( r, p[i], nothers, remaining) : 0) ;
    p[i] = newn ;
    remaining -= newn ;
  }
  /* update for index 8 with the remaining to sample */
  assert( GLOBAL_NUMBER_TYPES - 1 < p.size() ) ;
  p[ GLOBAL_NUMBER_TYPES - 1 ] = (remaining < GLOBAL_CONST_BOTTLENECK ? remaining : GLOBAL_CONST_II );

  assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
}

static bool not_lost_type( const std::vector<unsigned>& p, gsl_matrix_short * M )
{
  std::vector<unsigned> x (GLOBAL_NUMBER_SITES, 0 ) ;
  for (int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( size_t s = 0; s < GLOBAL_NUMBER_SITES ; ++s){
      x[s] += gsl_matrix_short_get(M, i, s) < 1 ? p[i] : 0 ; }}

  /* GLOBAL_CONST_II is 2N the maximum number of diploid individuals */
 return std::all_of( x.begin(), x.end(), []( unsigned n ){ return n < GLOBAL_CONST_II ; } );
}

  
static  void sample_juveniles_according_to_weight(std::vector<unsigned>& population, const std::vector< std::pair< size_t, double>>& juveniles,   const double c_nth)
  {
    assert( c_nth > 0. );
    std::fill( population.begin(), population.end(), 0) ;
    assert( current_number_individuals(population) < 1);

    /* check number of juveniles and  nth element */
    assert( juveniles.size() >= GLOBAL_CONST_II ) ;
    // countjuvenilesurvivingselection( juveniles, c_nth) ;
    
    size_t j = 0 ;
    while( j < GLOBAL_CONST_II ){
      assert( j < GLOBAL_CONST_II) ;
      population[ std::get<0>(juveniles[j]) ] += std::get<1>(juveniles[j]) <= c_nth ? 1 : 0 ;
      ++j ;
    }
    assert( current_number_individuals(population) == GLOBAL_CONST_II );
  }


static bool comp( std::pair<size_t, double> a, std::pair<size_t, double> b) 
  { 
    return ( std::get<1>(a) < std::get<1>(b) ); 
  }



static  double nthelm( std::vector< std::pair< size_t, double>>& juveniles )
  {
    std::nth_element( juveniles.begin(), juveniles.begin() + (GLOBAL_CONST_II - 1), juveniles.end(), comp);
    
    return( std::get<1>(juveniles[ GLOBAL_CONST_II - 1]) );
  }



static bool onestep( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs, gsl_matrix_short * M, gsl_rng *r )
{
  double nth {} ;
  double u = gsl_rng_uniform(r) ; 
  /* check if bottleneck */
  if( u < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
    /* bottleneck occurs ; sample surviving types and update population */
   sample_surviving_bottleneck(p, r) ; }
  /* first check if lost type at either loci */
  if ( not_lost_type(p, M) ){
    /* not lost type; check if fixed at all sites */
    if( p.back() < GLOBAL_CONST_II ){
        /* not all  individuals of type 2 at all sites, so sample juveniles */
      generate_pool_juveniles(jvs, p, cdfone, cdftwo, M, r);
      if( jvs.size() <= GLOBAL_CONST_II )
          {
            /* total number of juveniles not over capacity so all survive */
            update_population_all_juveniles(jvs, p) ;
            assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
          }
        else{
          /* need to sort juveniles and sample according to weight */
          nth = nthelm(jvs) ;
          sample_juveniles_according_to_weight(p, jvs, nth) ;
          assert( current_number_individuals(p) >= GLOBAL_CONST_II ); 
        }
      }
      /* mutation has fixed at both loci */
    }
    /* mutation has been lost */


  return (u < GLOBAL_CONST_PROBABILITY_BOTTLENECK) ;
}



static void trajectory( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,  const int numer,  gsl_matrix_short * M,  gsl_rng *r)
{
  init_for_trajectory(p) ;
  
  const std::string skra = "twosites_uoadd_enullpnulls1Ce4_" + std::to_string(numer) ;

  std::vector< double > excursion_to_fixation {} ;
  std::vector< int > time_bottleneck {} ;
  time_bottleneck.clear() ;
  bool b { } ;
  int timi = 0 ;
  while( ( not_lost_type(p, M) ) && ( p.back() < GLOBAL_CONST_II ) ){
    /* record the number of diploid individuals  homozygous  1/1 at all sites over current number of diploid individuals  */
    assert( GLOBAL_CONST_PROBABILITY_BOTTLENECK > 0. ? (current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK) : (1 > 0) ) ;
    excursion_to_fixation.push_back( static_cast< double>( p.back() ) / static_cast<double>( current_number_individuals(p) ) ) ;
    ++ timi ;
    // std::cout << excursion_to_fixation.back() << ' ' ;
    b = onestep(p, cdfone, cdftwo, jvs, M, r) ;
    if( b ){
      time_bottleneck.push_back( timi ); }
  }
  //  std::cout << '\n' ;

  assert( p.back() == p[GLOBAL_NUMBER_TYPES - 1] ) ;
  if( p.back() == GLOBAL_CONST_II ){

    /* fixation occurs so  print excursion to file */
    // FILE *fptr= fopen("dfdfd"  std::to_string(33), "a");
    std::ofstream outfile (skra, std::ios_base::app) ;
      for( const auto& y: excursion_to_fixation){
        // fprintf(fptr, "%g ", y) ;
        outfile << y << ' ' ;}
      outfile << '\n' ;
      // fprintf(fptr, "\n");
      // fclose( fptr) ;
      outfile.close() ;

      /* print time of bottleneck to file */
      std::ofstream bfile("timesbottle_.txt", std::ios_base::app ) ;
      for( const auto &y : time_bottleneck){
	bfile << y << ' ' ;}
      bfile << '\n';
      bfile.close() ;
  }
  
  std::cout << ( p[GLOBAL_NUMBER_TYPES - 1] < GLOBAL_CONST_II ? 0 : 1) << ' ' << timi << '\n' ;
}


/* the mass function for number of juveniles */
static double px(const double k, const double calpha, const double ccutoff)
{
  return ( (pow( 1./k, calpha) - pow( 1./(k + 1.), calpha) )/( pow( .5, calpha) - pow( 1./(ccutoff + 1.), calpha) ) ) ;
}


static void initialise_cdf( std::vector<double>& cdfo, std::vector<double>& cdft )
{
  
  for( double i = 2; i <=  GLOBAL_CONST_PSI_ONE ; ++i){
    cdfo.push_back( cdfo.back() +   px( i, GLOBAL_CONST_ALPHA_ONE, GLOBAL_CONST_PSI_ONE) ) ;}

  for( double j = 2; j <= GLOBAL_CONST_PSI_TWO; ++j){
    cdft.push_back( cdft.back() +  px( j, GLOBAL_CONST_ALPHA_TWO, GLOBAL_CONST_PSI_TWO) ) ; }

}


static void runsims( const int x,   gsl_rng *r)
{

  std::vector< unsigned int> population (GLOBAL_NUMBER_TYPES, 0) ;
  std::vector< std::pair< size_t, double>> juveniles {} ;

  std::vector<double> cdf_one {} ;
  std::vector<double> cdf_two {} ;

  init_containers(population, cdf_one, cdf_two) ;
  initialise_cdf(cdf_one, cdf_two) ;

  gsl_matrix_short * M = gsl_matrix_short_calloc( GLOBAL_NUMBER_TYPES, GLOBAL_NUMBER_SITES) ;
  writearray( static_cast<double>( GLOBAL_NUMBER_SITES),  M);
  /* *** 
     std::ifstream f("types.txt");
  short y {} ;
  
  for( int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( int j = 0 ; j < GLOBAL_NUMBER_SITES ; ++j){
      f >> y ; 
      gsl_matrix_short_set(M, i, j, y) ;}}
  f.close() ;
  ****** */
    
  int z = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1;
  while( --z > 0){
    
    trajectory(population, cdf_one, cdf_two, juveniles, x, M,  r) ;
  } 
  
  //freememory(population, juveniles, cdf_one, cdf_two) ;
  gsl_matrix_short_free(M);
}


/* a few basic  checks */
static void wprofa(  )
{

  std::vector<unsigned> p {} ;
  std::vector<double> cdfo {};
  std::vector<double> cdft {};
  gsl_matrix_short * M = gsl_matrix_short_calloc( GLOBAL_NUMBER_TYPES, GLOBAL_NUMBER_SITES);
  
  init_containers(p, cdfo, cdft);
  read_types( p , M ) ;

  std::cout << (not_lost_type(p, M) ? 1 : 0) << '\n' ;

  gsl_matrix_short_free(M);
}


static void  profa()
{
  double g = 1 ;
  for (double i = 1 ; i < 5 ; ++i){
    g *= i ;}

  std::cout << 'g' << g << '\n' ;
}



int main(int argc, char * argv[] )
{

  setup_rng( static_cast<unsigned long>( atoi(argv[1]) ) );
  
   runsims(  atoi(argv[1]), rngtype ) ;
  //profa() ;
   gsl_rng_free( rngtype) ;
  
   return GSL_SUCCESS;
}
