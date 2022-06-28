\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\IN}{\ensuremath{\mathds{N}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
}
\makeatother
\title{Fixation at many sites}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany} \footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
 through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by the Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge,   Wolfgang Stephan, and BE. BE also acknowledges Start-up module grants through SPP 1819  with Jere Koskela and Maite Wilke-Berenguer, and  with Iulia Dahmer. \\ \today} \orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
 This code generates excursions of the evolution of a diploid
      population partitioned into an arbitrary number of unlinked
      sites and two genetic types at each site, with viability weight
      determined by $W = e^{-sf(g)}$, where $g = (g_1, \ldots, g_L)$
      are the genotypes of a given individual at $L$ sites, $f$ is a
      function determining how the genotypes affect the trait value,
      and $s > 0$ is the strength of selection on the trait.  The
      population evolves according to a model of random sweepstakes
      and randomly occurring bottlenecks and viability selection.  The code can be used to 
      estimate the probability of fixation of the ($L$-site) type
      conferring advantage, and the expected time to fixation
      conditional on fixation of the type conferring advantage. 
\end{abstract}


\tableofcontents



@* {\bf Copyright}. 

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline


This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf introduction}. 
\label{SEC:intro}

We are interested in how random sweepstakes, randomly occurring
bottlenecks, and dominance mechanisms affect the probability of
fixation, and the time to fixation conditional on fixation of the type
conferring selective advantage.  We consider a diploid population
partitioned into an arbitrary number of unlinked sites and two types
at each site. The population is evolving according to random
sweepstakes and randomly occurring bottlenecks. In between bottlenecks
the population experiences viability selection provided the number of
juveniles produced each time exceeds the given carrying capacity.



Consider a diploid population starting with $2N$ diploid individuals. Let
$X, X_1^N, \ldots, X_N^N$ be i.i.d.\ discrete random variables
taking values in $\{2, \ldots, u(N)\}$; the $X_1^N, \ldots, X_M^N$
denote the random number of juveniles independently produced in a
given generation by $M$ parent pairs  according to
\begin{equation}
\label{eq:1}
   \prb{X = k} = \frac{ 1 }{  2^{-\alpha} -  ( u(N) + 1)^{-\alpha}  } \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^{\alpha}} \right),
\quad 2 \leq k \leq u(N).
\end{equation}
The mass in Eq~\eqref{eq:1} is normalised so that
$\prb{ 1 \leq X \leq u(N)} =1 $, and $\prb{X = k} \ge \prb{X =
k+1}$. Given a pool of at least $2M$ juveniles, we sample $2M$
juveniles for the next generation.  Leaving out an atom at zero and
one gives $X_1^N + \cdots + X_N^N \ge 2M$ almost surely, guaranteeing
that we always have at least $2M$ juveniles to choose from in each
generation.  We randomize on $\alpha$ so that with probability
$1 - \varepsilon_N$ all parent pairs produce juveniles according to
Eq~\eqref{eq:1} with $\alpha = \alpha_2 \ge 2$, and with $\alpha = \alpha_1  \in (1,2)$
with probability   $\varepsilon_N$. 


We admit randomly occurring bottlenecks; in the beginning of every
generation a bottleneck occurs with a fixed probability, and when it
happens we sample uniformly at random and without replacement a fixed
number of diploid individuals to survive the bottleneck. The surviving
individuals then continue to produce juveniles.  If the total number
of juveniles does not exceed the carrying capacity $(C)$ all the
juveniles survive, otherwise we sample $C$ of them without replacement
and with weights. Let $A$ denote the event a bottleneck occurs, the
population size at time $N_{t + 1}$ is then given by
$N_{t + 1} = \one{S_{\lfloor M/2 \rfloor} \le C}S_{\lfloor M/2
\rfloor} + \one{S_{\lfloor M/2 \rfloor} > C}C$ where
$S_{\lfloor M/2 \rfloor}$ is the total number of juveniles produced by
$\lfloor M/2 \rfloor$ parent pairs, and
$M = \one{A}B + \one{A^\complement}N_t$ where $B$ is the bottleneck
`size', i.e.\ the number of individuals surviving a
bottleneck. Viability selection therefore does not kick in unless the
number of juveniles at any time exceeds the carrying capacity; given
recurrent bottlenecks the population may be evolving neutrally for several generations.


The code can be used to estimate the probability of fixation at all
sites for the type conferring advantage, and the time to fixation
conditional on fixation.  In this context `fixation' means fixation at
all sites of the type conferring advantage; if the type is lost at any
one site fixation cannot occur.

At the time of writing  weight of a juvenile required for viability selection is taken as a random exponential with rate
$\exp\left(-sf(g) \right)$ where $s > 0$ is the strength of selection
and $g = (g_1, \ldots, g_L)$ are the genotypes at the $L$ sites (\S~\ref{sec:weight}, \S~\ref{sec:compute}).  One
can assume various dominance mechanisms, and different dominance
mechanisms for different sites. 

@* {\bf code}. 
\label{sec:code}




The R code \S~\ref{sec:rcodetypes} produces all possible $L$-site
types,  for a given $L$ there's $3^L$ of them. For two sites $(L = 2)$ with nine possible types the matrix is
\begin{center}
\begin{tabular}[center]{lll}
\hline
type index Eq~\eqref{eq:2}&  two-site type \\ 
\hline
0 & (0,0) \\
1 & (0,1) \\
2 & (0,2) \\
3 & (1,0) \\
4 & (1,1) \\
5 & (1,2) \\
6 & (2,0) \\
7 & (2,1) \\
8 & (2,2) \\
\hline
\end{tabular}
\end{center}

Sections \S~\ref{sec:includes}--\S~\ref{sec:runsims} describe the individual modules. 
The algorithm is formulated in the following pseudocode.
Let $Y(t) := \left( Y_{1}(t),
\ldots, Y_{L}(t) \right) \in [0,1]^{L}$ be the L-site type  frequency
process, where $Y_{\ell}(t)$ is the number of copies of the type
conferring selective   advantage at site $\ell$  relative to the
population size at time $t$.   Write $[n] := \{1, 2, \ldots, n\}$ for
$n \in \N$, $N_{t}$ for the population size at time $t$, $B$ for the
population size right after a bottleneck, $p$ for the probability of a
bottleneck in a given generation,  $g = (g_{1}, \ldots, g_{L})$ an 
L-site type, and $f(g)$ a function for how $g$ determines the
trait value
\begin{enumerate}
\item $\left( Y_{1}(0), \ldots, Y_{L}(0) \right)\leftarrow \left(
\tfrac {1}{2N}, \ldots, \tfrac {1}{2N}  \right) $ with population set
in the starting configuration as described above
\item while( $\left\{ \prod_{\ell \in [L]} Y_{\ell}(t) > 0 \right\}
\bigcap  \left\{ \sum_{\ell \in [L]}Y_{\ell}(t) < L   \right\}  $  ) \\
then for generation $t + 1$ started with $N_{t}$ individuals:
\begin{enumerate}
\item $U \leftarrow$ random uniform on $[0,1]$
\item $M \leftarrow \one{U \le p}B + \one{U > p}N_{t} $ where $p$ is the probability of a bottleneck
\item  the $M$ individuals in the case of a bottleneck are
sampled uniformly at random and without replacement from the $N_t$ individuals
\item sample $X_{1}, \ldots, X_{\lfloor M/2 \rfloor}$ random number of
juveniles produced independently  by $\lfloor M/2 \rfloor$ parent
pairs and assign types according to Mendel's laws
\item $S_{ \lfloor M/2 \rfloor } \leftarrow X_{1} + \cdots + X_{\lfloor M/2 \rfloor}$ the
total number of juveniles 
\item   sample a random exponential $E  \leftarrow {\rm Exp}\left(  \exp\left( -sf(g) \right) \right)$ for each
juvenile with rate $   \exp\left( -sf(g) \right)  $
\item $\mathfrak{E} \leftarrow \one{S_{\lfloor M/2 \rfloor} > C }E_{(C)}  +
\one{ S_{ \lfloor M/2 \rfloor } \le C  }E_{ \left( S_{ \lfloor M/2
\rfloor } \right) }$ where $E_{(j)}$ is the $j$th smallest exponential
out of the $S_{\lfloor M/2 \rfloor }$ exponentials from (d)
\item  a juvenile with sampled exponential $E$ survives if  $ E \le \mathfrak{E}  $
\end{enumerate}
\item record the result of the trajectory according to  $\left\{
\prod_{\ell \in [L]} Y_{\ell}(t) = 0 \right\}$ or   $ \left\{ \sum_{\ell \in [L]}Y_{\ell}(t) = L   \right\}  $
\end{enumerate}








@*1 {\bf R code for types}. 
\label{sec:rcodetypes}

\begin{verbatim}
# given L sites there are 3**L possible L-site types
g <- function(l, LL)
{
    x <- numeric()
    while( length(x) < (3**LL) ){
        x <- c(x,  rep(0,3**(LL-l)),  rep(1, 3**(LL-l)), rep(2, 3**(LL-l))) }
    return (x)
}
# example:  set the number of sites w to five
w <- 5
# define the matrix for the types
m <- matrix( 0, 3**w, w)
# compute the  matrix of types at w  sites
for(l in 1:w){
    m[, l]  <- t(g(l, w))}
# write the matrix to a file
write(t(m), "types.txt", ncolumns =  w)
\end{verbatim}





@*1 {\bf the included libraries}. 
\label{sec:includes}


@<Includes@>=@#
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
#include "cweblsites.hpp"




@*1 {\bf GSL random number generator}. 
\label{sec:gslrng}

initialise a gsl random number generator; called in \S~\ref{sec:main}


@<gsl random number generator@>=@#
gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf look up an index}.
\label{sec:lookup}

Compute an index given an $L$-site type $(i_1, \ldots, i_L) \in \{0,1,2\}^L$ as
\begin{equation}
\label{eq:2}
f( (i_1, \ldots, i_L)) =  \sum_{\ell = 1} ^L  i_\ell 3^{L - \ell}.
\end{equation}
This enables us to look up the number of individuals in the population with the given $L$-site type


@<lookup@>=@#
static size_t lookup(const std::vector<short>& types )
  {
    size_t n = 0 ;
    for( int ell = 1 ; ell <= GLOBAL_NUMBER_SITES ; ++ell){
      n +=  types[ell - 1]*static_cast<short>( gsl_sf_pow_int(3, GLOBAL_NUMBER_SITES - ell) ) ; 
    }
    return (n) ;
  }



@*1 {\bf replicate a value}. 
\label{sec:crep}

replicate  a given value a given number of times as in R function {\tt rep}; required for enumerating all the possible L-site types

@<as R rep@>=@#
static void crep( const double LL,  const double l, const short n, std::vector<short>& y )
{
  y.clear();
  y.resize( pow(3, LL-l)) ;
  y.assign( pow(3, LL - l), n) ;

}

@*1 {\bf one column of type array}.
\label{sec:onecolumn}


write out one column of the type array, i.e.\ the expanded types at one site; required for producing all he possible L-site types


@< L-site type array one column@>=@#
static void g( const double LL, const double l, std::vector< short>& x )
{
  x.clear();
  std::vector<short> tmp {};
  const size_t m =   static_cast<size_t>( pow( 3., LL)) ;
  while( x.size() < m){
    for( short i = 0 ; i < 3 ; ++i){
    /* \newline \S~\ref{sec:crep} */
      crep( LL, l, i, tmp);
      x.insert( x.end(),  tmp.begin(), tmp.end() ) ;}
  }
}


@*1 {\bf write array of L-site types}.
\label{sec:allLsitetypes}

write the array of all L-site types


@<all L-site types@>=@#
static void writearray( const double LLw,  gsl_matrix_short * M)
{
  std::vector< short> d {} ;
  for( size_t k = 1; k <= static_cast<size_t>(LLw) ; ++k){
  /* \newline \S~\ref{sec:onecolumn} */
    g( LLw, k, d);
    /* \newline copy into array */
    for( size_t u = 0 ; u < d.size() ; ++u){
      gsl_matrix_short_set(M, u, k-1, d[u]);}
  }
}






@*1 {\bf initialize containers}.
\label{sec:initconts}

initialize the main containers; we start with $2N - L$ individuals
homozygous for the wild type at all sites, and $L$ distinct
individuals each heterozygous at a distinct site.

@<container initialization@>=@#
static void init_containers( std::vector< unsigned>& population, std::vector<double>& cdf_one, std::vector<double>& cdf_two)
  {

    population.clear() ;
    population.assign( GLOBAL_NUMBER_TYPES, 0);

/* \newline set $2N - L$ individuals homozygous for the wild type at all sites */
    population[0] = GLOBAL_CONST_II - GLOBAL_NUMBER_SITES ;
   
/* \newline  set an $L$ site type  heterozygous at one site and homozygous for wild type at all other sites */
    std::vector< short> types (GLOBAL_NUMBER_SITES, 0 );
    for( unsigned i = 0 ; i < GLOBAL_NUMBER_SITES ; ++i){
      std::fill( types.begin(), types.end(), 0);
      types[i] = 1 ;

/* \newline set the number of individuals heterozygous at site $i$ to one \S~\ref{sec:lookup} */
      population[ lookup( types) ] = 1 ;}


/* \newline  initialize the containers for  the CDFs for the number of juveniles */
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



@*1 {\bf  initialize for an excursion}. 
\label{sec:inittraject}

initialize the population array for a new trajectory

@<newtrajectory@>=@#
static void init_for_trajectory( std::vector< unsigned>& population)
  {
    std::fill( population.begin(), population.end(), 0);
    /* \newline set $2N - L$ individuals as homozygous for wild type at all sites */
    population[0] = GLOBAL_CONST_II - GLOBAL_NUMBER_SITES ;
    
     std::vector< short> types (GLOBAL_NUMBER_SITES, 0 );
    for( unsigned i = 0 ; i < GLOBAL_NUMBER_SITES ; ++i){
      std::fill( types.begin(), types.end(), 0);
      /* \newline  set one individual as heterozygous at site |i| and homozygous for wild type at all other sites */
      types[i] = 1 ;
      population[ lookup( types) ] = 1 ;}
  }


@*1 {\bf total number of individuals in population}. 
\label{sec:Nt}

 return the current total number of individuals in the population

@<Nt@>=@#
static  unsigned int current_number_individuals(const std::vector<unsigned>& population)
  {
    return std::accumulate( std::begin(population), std::end(population), 0); 
  }



@*1 {\bf sample index of parent }. 
\label{sec:sampleparentindex}


 sample index of $L$-site type of parent without replacement  and update the number of remaining individuals accordingly

@<getparentindex@>=@#
static  int sample_genotype_parent( std::vector<unsigned>& p,   gsl_rng *r )
{
  /* \newline  |p| is population */ 
  int i  = 0 ;
  unsigned int nothers = current_number_individuals(p) - p[i] ;
  unsigned int x =  gsl_ran_hypergeometric( r, p[0], nothers, 1);

  while( (x < 1) && (i < GLOBAL_NUMBER_TYPES ) ){
    ++i ;
    nothers -= p[ i];
    x =  gsl_ran_hypergeometric(r, p[i], nothers, 1);
  }
  /* \newline  check if an individual has been sampled */
  i += (x < 1 ? 1 : 0); 
  /* \newline adjust the number of remaining parents */
  /* \newline  an individual of type with index |i| sampled, so subtract one from the number of remaining individuals with same type */
  --p[i] ;
  /* return the index  of the genotype of the parent */
  /* \newline  index is between |0| and |GLOBAL_NUMBER_TYPES| */
  return i ;
}


@*1 {\bf all juveniles survive }. 
\label{sec:alljuvssurvive}


 all juveniles survive

@<all juveniles@>=@#
static void update_population_all_juveniles(const std::vector< std::pair< size_t, double>>& juveniles, std::vector<unsigned>& population)
  {
  /* \newline set number of all types to zero */
    std::fill(population.begin(), population.end(), 0 ) ;
    /* \newline  \S~\ref{sec:Nt} */
    assert( current_number_individuals(population) < 1 ); 
    for( const auto &j : juveniles){
      /* \newline  |j[0]| is type index of juvenile |j| */
      population[ std::get<0>( j ) ] += 1;}
  }


@*1  {\bf sample a juvenile type}.
\label{sec:sampletype}

assign single site genotype to juvenile given genotypes in parents
following Mendel's laws; in our coding $0$ corresponds to homozygous
$0/0$ for the wild type, $1$ corresponds to the heterozygote,  and $2$ to the homozygote  for the mutation

@<sample single site type@>=@#
static short assign_type_one_site( const short gone, const short gtwo, gsl_rng *r)
{

/* \newline |gone| and |gtwo| are the two parent types */
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


@*1 {\bf read in file with types}.
\label{sec:readtypesfile}

 read in file with all types for given number of sites produced by the R code in \S~\ref{sec:rcodetypes}

@<file of types@>=@#
static void read_types( const std::vector<unsigned>& p ,  gsl_matrix_short * M  )
{
  std::ifstream f("types_five_sites.txt");

  short x {} ;
 
  for( int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( int j = 0 ; j < GLOBAL_NUMBER_SITES ; ++j){
      f >> x ; 
      gsl_matrix_short_set(M, i, j, x) ;}}
  f.close() ;

  /* \newline  print matrix for check */
  int z = 0 ;
   for( int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( int j = 0 ; j < GLOBAL_NUMBER_SITES ; ++j){
      std::cout << gsl_matrix_short_get(M, i,j) << ' ' ;}
    std::cout <<  p[z] <<   '\n' ;
    ++z ; }
}


@*1 {\bf a random number of juveniles}.
\label{sec:randomnojuvs}

 sample a random number  of juveniles with a distribution based on Eq~\eqref{eq:1}, \newline  returning $\min\{j \in \IN :  F(j) \ge u\}$ where $u$ is a given
 random uniform and $F$ the CDF

@<randomnumberjuvs@>=@#
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


@*1 {\bf weight of a site}.
\label{sec:weight}


compute the contribution of a site to the weight

@<weight@>=@#
static double weight( const short x )
{
 return ( (x < 1 ? 2. : 0.) / GLOBAL_NUMBER_SITES_d ) ;
}

@*1 {\bf compute weight }.
\label{sec:compute}

compute weight of a juvenile given the type at all sites

@<compute the weight@>=@#
static double computeweight(  const std::vector<short>& types,   gsl_rng *r)
{
  /* \newline  |types| is the vector of type at each site of juvenile  */

  double g = 0 ;
  
  /* \newline  compute the average contribution  over the sites  */
  for (const auto &s : types){
  /* compute the contribution of site |s| to the weight \S~\ref{sec:weight} */ 
    g += weight(s) ; }


/* \newline return the weight as a random exponential with rate $\exp\left( -sg^2 \right)$ where $g$ is from \S~\ref{sec:weight} */  
  return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( g , 2.) ) ) ) ;
}


@*1 {\bf add one juvenile}.
\label{sec:addjuv}

add one juvenile to  the pool of juveniles 

@<add one juvenile@>=@#
static void add_juvenile(const int type_index_one, const int type_index_two, gsl_matrix_short * Mtypes, std::vector< std::pair<size_t, double>>& jvs,   gsl_rng * r)
{
  std::vector<short> types_juvenile {} ;
  types_juvenile.clear() ;

  /* \newline  get the types of the juvenile */
  for( int s = 0; s < GLOBAL_NUMBER_SITES ; ++s){
    /* \newline  sample and record the type at site |s| \S~\ref{sec:sampletype} */
    types_juvenile.push_back( assign_type_one_site( gsl_matrix_short_get( Mtypes, type_index_one, s), gsl_matrix_short_get( Mtypes, type_index_two, s), r) ) ;}

  /* \newline  compute the weight \S~\ref{sec:compute} and type index \S~\ref{sec:lookup}  and add the juvenile to the vector of juves */
  jvs.push_back( std::make_pair( lookup( types_juvenile), computeweight( types_juvenile, r) ) );
}


@*1 {\bf juveniles for a parent pair}. 
\label{sec:juvsforparentpair}

produce a set of juveniles for one parent pair with given type indexes

@<produce juvs for parent pair@>=@#
static void add_juveniles_for_given_parent_pair( const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,   const int gone, 
const int gtwo, const size_t conetwo,  gsl_matrix_short * Mtypes,   gsl_rng *r)
{
        /* \newline  |gone| and |gtwo| are the type indexes for the two parents; |Mtypes| the matrix of types  */
  /* \newline  first sample a random number of juveniles \S~\ref{sec:randomnojuvs} */
  const size_t numberj = sample_random_number_juveniles(  conetwo, cdfone, cdftwo, r);
  assert( numberj > 1 ) ;

  /* \newline  add the sampled number of juveniles to the pool \S~\ref{sec:addjuv} */
  for( size_t j = 0; j < numberj ; ++j){
    add_juvenile( gone, gtwo,  Mtypes, jvs, r); }
}


@*1 {\bf a new pool of juveniles}.
\label{sec:newpooljuvs}

generate a new pool of juveniles 

@<new pool juvs@>=@#
static void generate_pool_juveniles( std::vector< std::pair<size_t, double>>& jvs, std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo,  gsl_matrix_short * Mtypes,  gsl_rng *r)
{
  jvs.clear() ;
  jvs.shrink_to_fit();
  assert( jvs.size() < 1 );
  int gone {} ;
  int gtwo {} ;
  /* \newline sample distribution of number of juveniles */
  const size_t conetwo = (gsl_rng_uniform(r) < GLOBAL_CONST_EPSILON ? 1 : 2) ; 
  /* \newline  |i| runs over number of pairs that can be formed from the current number of individuals \S~\ref{sec:Nt} */
  const double currenti = current_number_individuals(p) ;  
  assert( currenti  < GLOBAL_CONST_I + 1) ;

  for ( double i = 0 ; i < floor( currenti / 2. ) ; ++i){
  /* \newline \S~\ref{sec:sampleparentindex} */
    gone = sample_genotype_parent(p, r);
    gtwo = sample_genotype_parent(p, r);
    /* \newline  |gone| and |gtwo| are the type indexes of the two parents */
    assert(gone > -1);
    assert( gtwo > -1);
    /* \newline \S~\ref{sec:juvsforparentpair} */
    add_juveniles_for_given_parent_pair(cdfone, cdftwo, jvs, gone, gtwo, conetwo, Mtypes, r) ;}

  assert( jvs.size()  >= static_cast<size_t>( currenti ) ) ;
}


@*1 {\bf surviving a bottleneck}. 
\label{sec:survivebottle}

sample diploid individuals surviving a bottleneck; we sample uniformly at random without replacement

@<survive bottleneck@>=@#
static void sample_surviving_bottleneck( std::vector<unsigned>& p,   gsl_rng * r)
{
  int i = 0 ;

  unsigned int nothers = current_number_individuals(p) - p[i] ;
  unsigned newn = gsl_ran_hypergeometric( r,  p[i],  nothers,  GLOBAL_CONST_BOTTLENECK);
  unsigned int remaining = GLOBAL_CONST_BOTTLENECK - newn ;
  
  /* \newline  update count of individuals of type index |i| surviving bottleneck */
  p[i] = newn ;
  while( i < GLOBAL_NUMBER_TYPES - 2 ){
    ++i ;
    nothers -= p[i] ;
    newn = (remaining > 0 ? gsl_ran_hypergeometric( r, p[i], nothers, remaining) : 0) ;
    p[i] = newn ;
    remaining -= newn ;
  }
  
  assert( GLOBAL_NUMBER_TYPES - 1 < p.size() ) ;
  p[ GLOBAL_NUMBER_TYPES - 1 ] = (remaining < GLOBAL_CONST_BOTTLENECK ? remaining : GLOBAL_CONST_II );

  assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
}


@*1 {\bf check if lost a type}. 
\label{sec:checklosttype}


return True if not lost the type conferring advantage at any site, otherwise False


@<mutation is still around@>=@#
static bool not_lost_type( const std::vector<unsigned>& p, gsl_matrix_short * M )
{
  std::vector<unsigned> x (GLOBAL_NUMBER_SITES, 0 ) ;
  for (int i = 0 ; i < GLOBAL_NUMBER_TYPES ; ++i){
    for( size_t s = 0; s < GLOBAL_NUMBER_SITES ; ++s){
      x[s] += gsl_matrix_short_get(M, i, s) < 1 ? p[i] : 0 ; }}

  /* \newline  |GLOBAL_CONST_II| is $2N$  the maximum number of diploid individuals */
 return std::all_of( x.begin(), x.end(), []( unsigned n ){ return n < GLOBAL_CONST_II ; } );
}


@*1 {\bf sample juveniles by weight }.
\label{sec:samplebyweight}

sample juveniles by weight; the surviving juveniles in number equal the carrying capacity

@<pick by weight@>=@#
static  void sample_juveniles_according_to_weight(std::vector<unsigned>& population, const std::vector< std::pair< size_t, double>>& juveniles,   const double c_nth)
  {
    assert( c_nth > 0. );
    std::fill( population.begin(), population.end(), 0) ;
    /* \S~\ref{sec:Nt} */
    assert( current_number_individuals(population) < 1);

    /* \newline check number of juveniles and  |nth| element */
    assert( juveniles.size() >= GLOBAL_CONST_II ) ;
    
    size_t j = 0 ;
    while( j < GLOBAL_CONST_II ){
      assert( j < GLOBAL_CONST_II) ;
      population[ std::get<0>(juveniles[j]) ] += std::get<1>(juveniles[j]) <= c_nth ? 1 : 0 ;
      ++j ;
    }
     /* \S~\ref{sec:Nt} */
    assert( current_number_individuals(population) == GLOBAL_CONST_II );
  }


@*1 {\bf compare two juveniles }. 
\label{sec:compare}

compare the weight of two juveniles, needed for computing the $2N$th smallest weight among the weight of juveniles

@<compare@>=@#
static bool comp( std::pair<size_t, double> a, std::pair<size_t, double> b) 
  { 
    return ( std::get<1>(a) < std::get<1>(b) ); 
  }


@*1 {\bf  computing the $2N$th smallest weight }. 
\label{sec:nth}

compute the $2N$th smallest weight among the weight of juveniles,
needed for sampling the juveniles according to weight when the total number of juveniles exceeds the carrying capacity

@<nth@>=@#
static  double nthelm( std::vector< std::pair< size_t, double>>& juveniles )
  {
  /* \S~\ref{sec:compare} */
    std::nth_element( juveniles.begin(), juveniles.begin() + (GLOBAL_CONST_II - 1), juveniles.end(), comp);
    
    return( std::get<1>(juveniles[ GLOBAL_CONST_II - 1]) );
  }


@*1 {\bf take one step }.
\label{sec:onestep}

step through one generation


@<one step@>=@#
static void onestep( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs, gsl_matrix_short * M, gsl_rng *r )
{
  double nth {} ;
  /* \newline  check if bottleneck */
  if( gsl_rng_uniform( r) < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
    /* \newline  bottleneck occurs ; sample surviving types and update population \S~\ref{sec:survivebottle} */
   sample_surviving_bottleneck(p, r) ; }
  /* \newline  first check if lost type at any site \S~\ref{sec:checklosttype} */
  if ( not_lost_type(p, M) ){
    /* \newline not lost type; check if fixed at all sites */
    if( p.back() < GLOBAL_CONST_II ){
        /* \newline  not all  individuals of type 2 at all sites, so sample juveniles \S~\ref{sec:newpooljuvs} */
      generate_pool_juveniles(jvs, p, cdfone, cdftwo, M, r);
      if( jvs.size() <= GLOBAL_CONST_II )
          {
            /* \newline  total number of juveniles not over capacity so all survive \S~\ref{sec:alljuvssurvive} */
            update_population_all_juveniles(jvs, p) ;
            /* \newline \S~\ref{sec:Nt} */
            assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
          }
        else{
          /* \newline  need to sort juveniles and sample according to weight \S~\ref{sec:nth}  */
          nth = nthelm(jvs) ;
          /* \newline \S~\ref{sec:samplebyweight} */
          sample_juveniles_according_to_weight(p, jvs, nth) ;
          assert( current_number_individuals(p) >= GLOBAL_CONST_II ); 
        }
      }
      /* \newline mutation has fixed at both loci */
    }
    /* \newline  mutation has been lost */ 
}


@*1 {\bf generate one excursion}.
\label{sec:traj}

generate one excursion and record the result

@<one excursion@>=@#
static void trajectory( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,  const int numer,  gsl_matrix_short * M,  gsl_rng *r)
{
  init_for_trajectory(p) ;
  
  const std::string skra = "excursions_sites_" + std::to_string(numer) ;

  std::vector< double > excursion_to_fixation {} ;
  int timi = 0 ;
  while( ( not_lost_type(p, M) ) && ( p.back() < GLOBAL_CONST_II ) ){
    /* \newline record the number of diploid individuals  homozygous  1/1 at all sites over current number of diploid individuals  */
    assert( GLOBAL_CONST_PROBABILITY_BOTTLENECK > 0. ? (current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK) : (1 > 0) ) ;
    excursion_to_fixation.push_back( static_cast< double>( p.back() ) / static_cast<double>( current_number_individuals(p) ) ) ;
    ++ timi ;
 
    onestep(p, cdfone, cdftwo, jvs, M, r) ; }
 


  if( p.back() == GLOBAL_CONST_II ){

    /* \newline  fixation occurs so  print excursion to file */
    
    std::ofstream outfile (skra, std::ios_base::app) ;
      for( const auto& y: excursion_to_fixation){
        outfile << y << ' ' ;}
      outfile << '\n' ;

      outfile.close() ;
  }
  
  std::cout << ( p[GLOBAL_NUMBER_TYPES - 1] < GLOBAL_CONST_II ? 0 : 1) << ' ' << timi << '\n' ;
}



@*1 {\bf mass function for number of juveniles}. 
\label{sec:mass}


 the mass function for number of juveniles as in Eq~\eqref{eq:1}

@<mass@>=@#
static double px(const double k, const double calpha, const double ccutoff)
{
  return ( (pow( 1./k, calpha) - pow( 1./(k + 1.), calpha) )/( pow( .5, calpha) - pow( 1./(ccutoff + 1.), calpha) ) ) ;
}



@*1 {\bf initialise the CDF }.
\label{sec:initcdf}

initialise the CDF  corresponding to Eq~\eqref{eq:1} for the distribution of the  number of juveniles

@<init cdf@>=@#
static void initialise_cdf( std::vector<double>& cdfo, std::vector<double>& cdft )
{
  
  for( double i = 2; i <=  GLOBAL_CONST_PSI_ONE ; ++i){
    cdfo.push_back( cdfo.back() +   px( i, GLOBAL_CONST_ALPHA_ONE, GLOBAL_CONST_PSI_ONE) ) ;}

  for( double j = 2; j <= GLOBAL_CONST_PSI_TWO; ++j){
    cdft.push_back( cdft.back() +  px( j, GLOBAL_CONST_ALPHA_TWO, GLOBAL_CONST_PSI_TWO) ) ; }
}


@*1 {\bf run a given number of experiments}.
\label{sec:runsims}

run a given number of experiments and record the results

@<generate many excursions@>=@#
static void runsims( const int x,   gsl_rng *r)
{

  std::vector< unsigned int> population (GLOBAL_NUMBER_TYPES, 0) ;
  std::vector< std::pair< size_t, double>> juveniles {} ;

  std::vector<double> cdf_one {} ;
  std::vector<double> cdf_two {} ;

/* \S~\ref{sec:initconts} */
  init_containers(population, cdf_one, cdf_two) ;
  /* \S~\ref{sec:initcdf} */
  initialise_cdf(cdf_one, cdf_two) ;

  /* \newline write the matrix of all L-site types \S~\ref{sec:allLsitetypes} */
   gsl_matrix_short * M = gsl_matrix_short_calloc( GLOBAL_NUMBER_TYPES, GLOBAL_NUMBER_SITES ) ; 
  writearray( static_cast<double>(GLOBAL_NUMBER_SITES) ,  M) ;

  int z = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1;
  while( --z > 0){
    /* \newline  \S~\ref{sec:traj} */
    trajectory(population, cdf_one, cdf_two, juveniles, x, M,  r) ;
  } 
  
  gsl_matrix_short_free(M);
}



@*1 {\bf main module}. 
\label{sec:main}

the {\tt main} function calling 


@C


/* \newline \S \ref{sec:includes} */
@<Includes@>@#
/* \newline \S~\ref{sec:gslrng} */
@<gsl random number generator@>@#
/* \newline \S~\ref{sec:lookup} */
@<lookup@>@#
/* \newline \S~\ref{sec:initconts} */
@<container initialization@>@#
/* \newline \S~\ref{sec:inittraject} */
@<newtrajectory@>@#
/* \newline \S~\ref{sec:Nt} */
@<Nt@>@#
/* \newline \S~\ref{sec:sampleparentindex} */
@<getparentindex@>@#
/* \newline \S~\ref{sec:alljuvssurvive} */
@<all juveniles@>@#
/* \newline \S~\ref{sec:sampletype} */
@<sample single site type@>@#
/* \newline \S~\ref{sec:readtypesfile} */
@<file of types@>@#
/* \newline \S~\ref{sec:randomnojuvs} */
@<randomnumberjuvs@>@#
/* \newline \S~\ref{sec:weight} */
@<weight@>@#
/* \newline \S~\ref{sec:compute} */
@<compute the weight@>@#
/* \newline \S~\ref{sec:addjuv} */
@<add one juvenile@>@#
/* \newline \S~\ref{sec:juvsforparentpair} */
@<produce juvs for parent pair@>@#
/* \newline \S~\ref{sec:newpooljuvs} */
@<new pool juvs@>@#
/* \newline \S~\ref{sec:survivebottle} */
@<survive bottleneck@>@#
/* \newline \S~\ref{sec:checklosttype} */
@<mutation is still around@>@#
/* \newline \S~\ref{sec:samplebyweight} */
@<pick by weight@>@#
/* \newline \S~\ref{sec:compare} */
@<compare@>@#
/* \newline \S~\ref{sec:nth} */
@<nth@>@#
/* \newline \S~\ref{sec:onestep} */
@<one step@>@#
/* \newline \S~\ref{sec:traj} */
@<one excursion@>@#
/* \newline \S~\ref{sec:mass} */
@<mass@>@#
/* \newline \S~\ref{sec:initcdf} */
@<init cdf@>@#
@<as R rep@>@#
/* \newline \S~\ref{sec:onecolumn} */
@< L-site type array one column@>@#
/* \newline \S~\ref{sec:allLsitetypes} */
@<all L-site types@>@#
/* \newline \S~\ref{sec:runsims} */
@<generate many excursions@>@#
/* \newline \S~\ref{sec:crep} */



int main(int argc, char * argv[] ){

/* \newline \S~\ref{sec:gslrng} */
setup_rng( static_cast<unsigned long>( atoi(argv[1]) ) );

/* \newline \S~\ref{sec:runsims} */
runsims( atoi(argv[1]), rngtype);

gsl_rng_free( rngtype);

return GSL_SUCCESS ;

}



@
\end{document}