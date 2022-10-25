\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{graphicx}
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
\usepackage{refcheck}
\usepackage{environ}
\usepackage{orcidlink}
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}
\newcommand{\im}{\ensuremath{\imath} }
\newcommand{\jm}{\ensuremath{\jmath} }
\newcommand{\be}{\begin{equation}}
\newcommand{\ee}{\end{equation}}
\newcommand{\bd}{\begin{displaymath}}
\newcommand{\ed}{\end{displaymath}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\C}{\ensuremath{\mathfrak{C}} }
\newcommand{\F}{\ensuremath{\mathbb{F}}}
\newcommand{\G}{\ensuremath{\mathds{G}}}
\newcommand{\NN}{\ensuremath{\mathds{N}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\ao}{\ensuremath{{\alpha_1}}}
\newcommand{\at}{\ensuremath{{\alpha_2}}}
\newcommand{\OO}[1]{\ensuremath{\mathcal{O}\left( #1\right)}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
} \makeatother \title{Fixation at linked sites}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany\\
Supported by Deutsche Forschungsgemeinschaft (DFG) - Projektnummer
273887127 through DFG SPP 1819: Rapid Evolutionary Adaptation grant
STE 325/17 to Wolfgang Stephan; BE acknowledges funding by Icelandic
Centre of Research through an Icelandic Research Fund Grant of
Excellence no.\ 185151-051 to Einar \'Arnason, Katr\'in
Halld\'orsd\'ottir, Wolfgang Stephan, Alison M.\ Etheridge, and BE;
and SPP 1819 Start-up module grants with Jere Koskela and Maite Wilke
Berenguer, and with Iulia Dahmer \\ \today}\orcidlink{0000-0001-9354-2391}}

\begin{document}
\maketitle


\rule{\textwidth}{.8pt}


\begin{abstract}
A simulator for the evolution of a diploid population evolving
according to random sweepstakes, recurrent bottlenecks, and  viability
selection acting on linked sites. At each site there are two types,
the wild and the fit type,  with the latter conferring selective
advantage.    We record excursions to fixation of the fit type jointly
at all the sites under  all kinds of  dominance and epistatic
mechanisms. 
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



@* {\bf intro}. 
\label{SEC:intro}


Consider a diploid population evolving in discrete (non-overlapping)
generations according to random sweepstakes and recurrent bottlenecks.  In each
generation the current individuals randomly form pairs, and each pair
independently produces a random number of potential offspring (juveniles) according to a given law. If the
total number of juveniles so produced exceeds a fixed carrying
capacity \C we sample \C of them without replacement based on their
viability weight, otherwise all the juveniles survive.


Let $X, X_{1},\ldots, X_{M}$ for $M \in \NN := \{1,2, \ldots\}$  denote iid copies of positive random
variables with
\begin{equation}
\label{eq:1}
\prb{X=k} = \one{2 \le k \le u(N)} C\left( \frac{1}{k^{\alpha}}  -  \frac{1}{(k+1)^{\alpha}}  \right)
\end{equation}
where $\alpha, C > 0$ are constants so that
\begin{displaymath}
\begin{split}
& \prb{X = k}   \ge \prb{X = k+1}, \quad   k \ge 0; \\
& \prb{2 \le X \le u(N) }  = 1
\end{split}
\end{displaymath}
We choose the law so that  $\prb{X \ge 2M} = 1$.


We randomize on $\alpha$ in Eq~\eqref{eq:1}.  Fix $0 < \alpha_{1} < 2$
and $\alpha_{2} \ge 2$.  With probability
$\varepsilon$ all current  parent pairs independently  produce
juveniles with $\alpha = \alpha_{1}$, and with probability $1 -
\varepsilon$  we have $\alpha = \alpha_{2}$.   Such a mixture can be
shown to have better properties than fixing $\alpha$.  


Recurrent bottlenecks are modelled by tossing a coin at the start of a
given generation; if a bottleneck occurs we sample a fixed number $B$
of diploid individuals to survive a bottleneck and produce juveniles
who all survive if the total number of juveniles does not exceed the
carrying capacity {\C}.  Let $N_{t}$ denote the population size, the
number of diploid individuals at time $t$.  Then
\begin{equation}
\label{eq:2}
N_{t+1} =  \C \one{S_{ \lfloor M/2 \rfloor } > \C} +  S_{M} \one{ S_{ \lfloor M/2 \rfloor} \le \C  }
\end{equation}
where  $S_{N}$ denotes the  total number of juveniles produced by $N$
parent pairs and 
\begin{equation}
\label{eq:3}
M =  B\one{\text{bottleneck occurs}} +  N_{t}\one{\text{bottleneck does not occur}}
\end{equation}
When the total number of juveniles exceeds \C we sample for each
juvenile a random exponential with rate the viability weight of the
juvenile. The juveniles with the \C smallest exponentials survive.
This is a way of applying viability selection to the evolution of the
population.  The viability weight of each juvenile is determined by
the genotypes at all the sites. 



@*1 {\bf pseudocode}. 
\label{sec:pseudocode}



@*  {\bf code}. 


@*1 {\bf the included libraries}. 
\label{sec:includes}

the included libraries and the header file

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
#include <bitset>
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include "Lsitesalllinked.hpp"



@*1 {\bf GSL random number generator}. 
\label{sec:gslrng}

initialise a GSL   random number generator; see  \S~\ref{sec:main}


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


@*1 {\bf lookup table}. 
\label{sec:lookup}

 generate a lookup table where entry with index |i|  contains the
 corresponding  haplotype indexes; see \S~\ref{sec:run}

@<generate lookup table@>=@#
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


@*1 {\bf the lookup function}. 
\label{sec:lookupfunction}

going along rows in the array of number of individuals of each diploid
phased $L$-site  type; see \S~\ref{sec:initarray} 

@<lookup function@>=@#
static unsigned long lookup( const unsigned long i, const unsigned long j)
{
  /* |i| and |j| with $i \le j$ are the two haplotype indexes \newline  */
  /* |GLOBAL_CONST_MAX_INDEX | is the index of haplotype   $1\ldots 1$ \newline  */
  return ( j  +  (i*GLOBAL_CONST_MAX_INDEX) -  (i*(i-1)/2) ) ;
}



@*1 {\bf check if recombination occurs}.
\label{sec:checkforrecombination}


check if parent haplotypes recombine and produce recombinant
haplotypes; we check this regardless of the configuration of the
haplotypes, i.e.\ if the recombinants would be different from the
originals or not. If a random uniform does not exceed  $r(L-1)$ where
$r$ is the recombination probability between a pair of sites
recombination occurs 

@<check for recombination@>=@#
static int recombination(  )
{
  /* if recombination return a discrete uniform |x| between  1 and  |L
  - 1| \newline */
  /* recombination then occurs between |x -1| and |x| \newline  */
  return ( gsl_rng_uniform(rngtype) <=  (GLOBAL_CONST_RECOMBINATION * (GLOBAL_CONST_NUMBER_SITESd - 1)) ? gsl_rng_uniform_int(rngtype, GLOBAL_CONST_NUMBER_SITES -1) + 1 :  GLOBAL_CONST_NUMBER_SITES + 2 ) ; 

  /* |gsl_rng_uniform_int( rngtype, n)| returns a random int between zero and |n-1| \newline */
  /* shifting by one returns an integer between one  and |L - 1|
  \newline   */
}



@*1 {\bf generate recombinant haplotypes}. 
\label{sec:recombinanthaplotypes}

generate recombinant haplotypes given recombination occurs; return one
of the two recombined haplotypes picked with equal probability; called
in \S~\ref{sec:samplehaplotypeindex};

if the parent haplotypes are $(a_{1}, a_{2}, \ldots, a_{L})$ and
$(b_{1}, \ldots, b_{L})$ and recombination occurs between sites
$\ell-1$ and $\ell$   then the two recombinant haplotypes are
$(a_{1}, \ldots, a_{\ell-1}, b_{\ell}, b_{\ell + 1}, \ldots, b_{L})$
and  $(b_{1}, \ldots, b_{\ell-1}, a_{\ell}, a_{\ell + 1}, \ldots, a_{L})$

@<recombine haplotypes@>=@#
static unsigned long recombine( const int lrec, const unsigned long indexhapone, const unsigned long indexhaptwo )
{
  /* |indexhapone| and |indexhaptwo| are the indexes for the two parent haplotypes \newline */
  /* recombination happens between site |lrec - 1| and |lrec| \newline
  */
  /* for |lrec| between one and |L - 1|; |L| is number of sites \newline */

  std::string hapone =  std::bitset<GLOBAL_CONST_NUMBER_SITES>(indexhapone).to_string() ;
  std::string haptwo =  std::bitset<GLOBAL_CONST_NUMBER_SITES>(indexhaptwo).to_string() ;
  std::string tmp = hapone ;

  /* generate the recombinat haplotypes by replacing the respective
  ends \newline */
  hapone.replace( lrec , GLOBAL_CONST_NUMBER_SITES - lrec , haptwo.substr( lrec,  GLOBAL_CONST_NUMBER_SITES - lrec ) ) ;
  haptwo.replace( lrec , GLOBAL_CONST_NUMBER_SITES - lrec , tmp.substr( lrec, GLOBAL_CONST_NUMBER_SITES - lrec )) ;  
  
  /* return the index of the recombinant haplotype picked for the
  juvenile \newline */
  return(  gsl_rng_uniform(rngtype) < 0.5 ? std::bitset<GLOBAL_CONST_NUMBER_SITES>(hapone).to_ulong() : std::bitset<GLOBAL_CONST_NUMBER_SITES>(haptwo).to_ulong() ) ;  
}




@*1 {\bf sample a haplotype index}. 
\label{sec:samplehaplotypeindex}


sample haplotype index from a parent given haplotype indexes of the
parent; the index may represent a recombinant 

@<sample haplotype index@>=@#
static unsigned long sample_haplotype_index( const unsigned long hone, const unsigned long htwo )
{
  /* |hone| and |htwo| are the haplotype indexes of the two parents; 
     \newline  */
     /* check if recombination occurs; see 
     \S~\ref{sec:checkforrecombination} for |recombination| \newline */
  const int x = recombination() ;
  /* if |x < | number of sites then  recombination occurs  \newline  */
  /* otherwise  return the index of  one or the other of the parent
  haplotypes;  see \S~\ref{sec:recombinanthaplotypes} for |recombine| \newline   */
  return( x < GLOBAL_CONST_NUMBER_SITES ?  recombine(x, hone, htwo) : (gsl_rng_uniform(rngtype) < 0.5 ? hone : htwo) ) ;
}



@*1 {\bf initializing the population array}. 
\label{sec:initarray}

 initializing the array of number of diploid individuals 
 of each phased  |L|-site  type

@<initialize population array@>=@#
static void initializearray( std::vector<unsigned long>& p)
{
  /* |p| is the population; each entry is the number of individuals of
  the respective phased |L|-site type \newline */
  std::string s (GLOBAL_CONST_NUMBER_SITES, '0') ;

  std::fill( std::begin(p), std::end(p), 0);
  assert( std::accumulate(std::begin(p), std::end(p) , 0) < 1 );

  
  for( int i = 0 ; i < GLOBAL_CONST_NUMBER_SITES ; ++i){
    s[i] = '1' ;
    /* \newline \S~\ref{sec:lookupfunction} for  |lookup| */
    p[ lookup(0, std::bitset<GLOBAL_CONST_NUMBER_SITES>(s).to_ulong() )] =  1 ;
    s[i] = '0' ; }

  /* |L| number of diploid individuals carry a fit type at one site,
  i.e. are heterozygous at one site only \newline */
  /* and at all other sites homozygous for the wild type; \newline  */
  /* all other diploid individuals are homozygous for the wild type at
  all sites \newline */
  p[0] = GLOBAL_CONST_CARRYING_CAPACITY - GLOBAL_CONST_NUMBER_SITES ;
}



@*1 {\bf current number of diploid individuals}. 
\label{sec:Nt}

compute the current number of diploid individuals

@<current number |Nt|@>=@#
static unsigned long current_number_individuals( const std::vector<unsigned long>& p)
{
  return  std::accumulate( std::begin(p), std::end(p), 0 ) ;
}




@*1 {\bf viability weight}. 
\label{sec:viabilityweight}

 compute  viability   weight of a  juvenile given haplotype indexes
 @<viability weight@>=@#
static double weight( const std::vector< unsigned long>& h)
{
  /* |h| contains haplotype indexes  of juvenile \newline  */

 /* convert the indexes to a string of zeros and ones; zero for wild
 type and one for fit type \newline */
  const std::string hapone = std::bitset<GLOBAL_CONST_NUMBER_SITES>( h[0]).to_string() ;
  const std::string haptwo = std::bitset<GLOBAL_CONST_NUMBER_SITES>( h[1]).to_string() ;
  double w = 0 ;
  /* need to be homozygous for fit type to increase weight by one
  \newline */
 /* $(g_{1}, g_{2}) \mapsto \one{g_{1} = 2} +  \one{g_{2} = 2}$ and no
 epistasis \newline  */
  for( int i = 0 ; i < GLOBAL_CONST_NUMBER_SITES ; ++i){
    w += ( hapone[i] == '1' ? (haptwo[i] == '1' ? 1. : 0) : 0.) ; }

  /* adding a small deviation not necessary :
  |gsl_ran_gaussian_ziggurat( rngtype, 0.1)| \newline  */
  return ( gsl_ran_exponential( rngtype, 1./( 1. +  (GLOBAL_CONST_SELECTION * w)  ) ) ) ;
}



@*1 {\bf probability kernel for number of juveniles}. 
\label{sec:kernel}

compute a kernel according to Eq~\eqref{eq:1} for sampling a random
number of juveniles


@<kernel@>=@#
static double masskernel( const double k, const double alpha)
{
  return ( pow( 1./k, alpha) - pow( 1./(k + 1.), alpha) ) ;
}




@*1 {\bf the CDF for distribution of number of juveniles}. 
\label{sec:cdf}

 generate the CDF for sampling random number of juveniles

@<cdf@>=@#
static void cdf_number_juveniles( std::vector<double>& x , const double a )
{
/* |a| is $\alpha$ in Eq~\eqref{eq:1} \newline  */  
  for( unsigned long k = 2; k <=  GLOBAL_CONST_CARRYING_CAPACITY ;
  ++k){
  /* \newline \S~\ref{sec:kernel} */
    x[k] =   ( masskernel( static_cast<double>(k), a) ) ; }

  const double cconst  =  std::accumulate( std::begin( x), std::end(x), 0.);

  for( unsigned long k = 2; k <=  GLOBAL_CONST_CARRYING_CAPACITY ; ++k){
    x[k] = x[k-1] + ( masskernel( static_cast<double>(k), a ) / cconst ) ; }

  /* the cdf must have last element one to guarantee a value within
  limits is sampled \newline  */
  x.back() = 1. ;
}


@*1 {\bf litter size}. 
\label{sec:littersize}

 sample a random number of juveniles for one family 

@<sample litter size@>=@#
static int sample_litter_size( const std::vector<double>& x)
{
  int j = 2 ;
  const double u = gsl_rng_uniform(rngtype);
  /* with $F$ denoting the CDF  return $\min\{j \in \{2,3, \ldots, u(N) \} :  u
  \le  F(j) \}$ \newline */
  while( u > x[j]){ ++j ;}
  
  assert( j > 1);
  assert( j <= static_cast< int>( GLOBAL_CONST_CUTOFF) ) ;

  return j ;
}


@*1 {\bf adding juvenile}. 
\label{sec:addingjuv}

 add a juvenile to the  litter 

@<addjuv@>=@#
static void add_juvenile( const std::vector< unsigned long>& pone, const std::vector<unsigned long>& ptwo, std::vector< std::pair< std::vector<unsigned long>, double>>& vj, std::vector<double>& v_juvenileweights )
{
  /* |pone| and |ptwo| are the pairs of haplotype indexes for the two
  parents \newline  */

  std::vector< unsigned long> h {} ;
  /* sample haplotype from parent one (arbitrarily enumerated)  with haplotype indexes in |pone|
  \newline */
  h.push_back( sample_haplotype_index( pone[0], pone[1] ) ) ;
    /* sample haplotype from parent two (aribtrarily enumerated)  with haplotype indexes in |ptwo|
  \newline */
  h.push_back( sample_haplotype_index( ptwo[0], ptwo[1] ) ) ;
/* |h| now contains the indexes of the two haplotypes carried by the
juvenile; the genome of the juvenile   \newline  */
/* compute the weight of the juvenile given the indexes \S~\ref{sec:viabilityweight}   \newline  */
  const double w = weight(h);
  /* record the weight of the juvenile for computing the |n|th
  smallest \newline */
  v_juvenileweights.push_back( w) ;
  
  vj.push_back( std::make_pair( h, w) );
}


@*1 {\bf add sibship to pool of juveniles}. 
\label{sec:sibship}

 add litter or sibship   to pool of juveniles

@<add sibship@>=@#
static void add_litter( const std::vector<double>& v_cdf, const std::vector<unsigned long>&  po, const std::vector<unsigned long>& pt, std::vector<  std::pair< std::vector<unsigned long>, double>>& pool, std::vector< double>& vw  )
{
  /* sample litter size \S~\ref{sec:littersize} \newline */
  const int littersize = sample_litter_size(v_cdf) ;
  for( int j = 0 ; j < littersize; ++j){
  /* given sibship size add juveniles one by one \S~\ref{sec:addingjuv} \newline */
    add_juvenile( po, pt, pool, vw);}
}



@*1 {\bf sample one diploid individual}. 
\label{sec:sampleoneindividual}

 sample one  index of phased |L|-site type of diploid parent and
 update the population; used in \S~\ref{sec:bottleneck} for removing
 individuals not surviving a bottleneck

@<onehypergeometric@>=@#
static unsigned long sample_genotype_parent( std::vector<unsigned long>& p )
{
  /* |p| is population \newline */ 
  unsigned long i  = 0 ;
  unsigned int nothers = static_cast< unsigned int>( current_number_individuals(p) - p[i]) ;
  unsigned int x =  gsl_ran_hypergeometric( rngtype, p[0], nothers, 1);

  while( (x < 1) && (i <  GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES ) ){
    ++i ;
    nothers -= p[ i];
    x =  gsl_ran_hypergeometric(rngtype, p[i], nothers, 1);
  }
  /* check if an individual has been sampled \newline  */
  i += (x < 1 ? 1 : 0); 
  /* adjust the number of remaining parents \newline  */
  /* an individual of type with index |i| sampled, so subtract one
  from the number of remaining individuals with same type \newline */
  --p[i] ;
  /* return the sampled  index  of the diploid  phased |L|-site type
  of the parent \newline */
  /* index is between zero and
  |GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES| \newline */
  return i ;
}



@*1 {\bf clear container for juveniles}. 
\label{sec:clearpool}

clear the container for the pool of juveniles; may not need this but
is here nevertheless

@<clear the pool@>=@#
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



@*1 {\bf a new pool of juveniles}. 
\label{sec:newpool}

 generate a new pool of juveniles 

@<a new pool@>=@#
static void new_pool_juveniles( const unsigned long Nt,  std::vector<  std::pair< std::vector<unsigned long>, double>>& pool, std::vector< double>& vw, std::vector<unsigned long>& population,  const std::vector<double>& vcdf, const std::vector< std::pair< unsigned long, unsigned long>>& table )
{
    /* |Nt| is current number of individuals \newline  */
  /* |pool| is the vector for the new pool of juveniles \S~\ref{sec:clearpool}  \newline  */
  clearpooljuveniles( pool ) ;

  vw.clear();
  vw.shrink_to_fit() ;
  assert( vw.size() < 1);

  unsigned long indexpone {} ;
  unsigned long indexptwo {} ;

  std::vector<unsigned long> pone (2,0);
  std::vector<unsigned long> ptwo (2,0);
  
  /* |Nt| is current number of individuals \newline */
  /* each time sample two parents \newline  */
  for( unsigned long i = 0 ; i < Nt/2 ; ++i){
    /* sample index of a phased |L|-site type of parent one \newline  */
    /* |table[indexpone]| contains the corresponding haplotype indexes \S~\ref{sec:sampleoneindividual}
    \newline  */
    indexpone = sample_genotype_parent( population );
    pone[0] = std::get<0>(table[ indexpone]) ;
    pone[1] = std::get<1>(table[ indexpone]) ;
    /* |table[indexptwo]|  contains the corresponding haplotype
    indexes \newline  */
    indexptwo = sample_genotype_parent( population ); 
    ptwo[0] = std::get<0>(table[ indexptwo]) ;
    ptwo[1] = std::get<1>(table[ indexptwo]) ;
    /* given parent genotypes sample a sibship  \newline \S~\ref{sec:sibship} */
    add_litter( vcdf, pone, ptwo, pool, vw ) ;
  }
}


@*1 {\bf count homozygous}. 
\label{sec:counthomoz}

 given a site for which to check, count how many diploid individuals
 are homozygous for the wild type at the site

@<count homozygous at site@>=@#
static unsigned long checknullsite( const unsigned long nsite,  const std::vector< unsigned long>& population, const std::vector< std::pair<unsigned long, unsigned long>>& tafla )
{
  /* |nsite| is the site of interest; we record how many individuals
  are  homozygous for the wild
  type at  site |nsite|  \newline */
  unsigned long sumof = 0 ;
  for( unsigned long i = 0 ; i < GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES ; ++i){
    /* if the type is  homozygous for the wild type at site |nsite|
    add the  number of individuals with the phased |L|-site type \newline */
    sumof +=  (std::bitset<GLOBAL_CONST_NUMBER_SITES>( std::get<0>(tafla[i]) ).to_string()[nsite] == '0' ? (std::bitset<GLOBAL_CONST_NUMBER_SITES>( std::get<1>(tafla[i]) ).to_string()[nsite] == '0' ?  population[i] : 0) : 0 ) ; 
  }
  /* return the number of individuals homozygous for the wild type at
  site |nsite| \newline */
  return sumof ;
}


@*1 {\bf compare two values}. 
\label{sec:comp}

compare two values for computing the $n$th smallest element \S~\ref{sec:nth}

@<compare@>=@#
static bool comp( const double a, const double b )
{
  return a < b ;
}

@*1 {\bf $n$th smallest element}. 
\label{sec:nth}

compute the $n$th element for  sampling juveniles when the total
number of juveniles exceeds the carrying capacity 

@<$n$th element@>=@#
static double nthelm(  std::vector< double>& weights )
{
/* \S~\ref{sec:comp} \newline */
  std::nth_element( weights.begin(), weights.begin() + (GLOBAL_CONST_CARRYING_CAPACITY - 1), weights.end(), comp ) ;
  return weights[ GLOBAL_CONST_CARRYING_CAPACITY - 1] ;
}



@*1 {\bf sample juveniles}. 
\label{sec:sortjuveniles}

 sort juveniles and update population 

@<sorting@>=@#
static void  select_juveniles_according_to_weight( std::vector<double>& v_weights,  std::vector< unsigned long>& v_population, const std::vector< std::pair< std::vector< unsigned long>, double>>& v_pool )
{
  /* compute the |n|th smallest viability  weight \newline */
  /* only sample juveniles according to weight if the total number of
  juveniles exceeds the carrying capacity \newline */
  /* otherwise  all juveniles survive \newline */
  const double wnth =  ( static_cast<unsigned long>( v_pool.size()) > GLOBAL_CONST_CARRYING_CAPACITY ? nthelm( v_weights) :  std::max_element( std::begin(v_weights), std::end(v_weights) )[0] + 1. )    ;

  assert( wnth > 0.) ;
  /* set number of all phased |L|-site types to zero \newline */
  std::fill( std::begin( v_population), std::end(v_population), 0); 
  assert( std::accumulate( std::begin(v_population), std::end(v_population), 0.) < 1);
  /* add to count of phased |L|-site type of juvenile if juvenile
  survives \newline  */
  for( const auto &j : v_pool){
  /* \newline  see \S~\ref{sec:lookupfunction} for |lookup| */
    v_population[ lookup( std::get<0>(j)[0], std::get<0>(j)[1]) ] += (  std::get<1>(j) <= wnth ? 1 : 0) ; }
}


@*1 {\bf check if lost a fit type}. 
\label{sec:checklosttype}

count the number of individuals homozygous for the wild type at each
site (|n|); if at any site s we have |n = Nt| where |Nt| is the
current number of individuals then lost type at the site

@<not lost a type@>=@#
static bool not_lost_type( std::vector<unsigned long>& nhomozygous ,  const std::vector< unsigned long>& population, const std::vector< std::pair<unsigned long, unsigned long>>& mtafla )
{

  std::fill( std::begin(nhomozygous), std::end(nhomozygous), 0);
  /* \newline  see \S~\ref{sec:Nt} for |current_number_individuals|  */
  const unsigned long Nt =  current_number_individuals( population) ; 
  for( unsigned long nullsite = 0 ; nullsite <
  GLOBAL_CONST_NUMBER_SITES ; ++nullsite){
  /* \newline see \S~\ref{sec:counthomoz} for |checknullsite| */
    nhomozygous[nullsite] =  checknullsite( nullsite, population, mtafla) ; }

  return std::all_of( std::begin(nhomozygous), std::end(nhomozygous), [Nt](unsigned long n){return n < Nt;} );
}



@*1 {\bf bottleneck}. 
\label{sec:bottleneck}

 remove individuals not surviving a bottleneck using
 |sample_genotype_parent| in \S~\ref{sec:sampleoneindividual}

@<generate a bottleneck@>=@#
static unsigned long remove_not_surviving_bottleneck(std::vector< unsigned long>& v_p)
{
  
  unsigned long x {} ;
  /* \newline see \S~\ref{sec:Nt} for |current_number_individuals| */
  const unsigned long currentNt =   current_number_individuals( v_p ); 
  for ( unsigned long i = 0 ; i < currentNt - GLOBAL_CONST_BOTTLENECK ; ++i){
       /* see \S~\ref{sec:sampleoneindividual} for
       |sample_genotype_parent|  */
        x = sample_genotype_parent( v_p) ; }

  assert( current_number_individuals( v_p) == GLOBAL_CONST_BOTTLENECK );

  return current_number_individuals( v_p)  ;
}


@*1 {\bf run experiments}.
\label{sec:run}

run experiments and record excursions to fixation at all sites if any occur


@<run experiments and record result@>=@#
static void run( )
{

  /* define the population \newline */
   std::vector< unsigned long> pop ( GLOBAL_CONST_TOTAL_NUMBER_PHASED_TYPES, 0);

   std::vector< unsigned long> nhomozygouswt (GLOBAL_CONST_NUMBER_SITES, 0); 
   std::vector< std::pair< unsigned long, unsigned long>> m {} ;
   m.clear() ;
   /* \newline see \S~\ref{sec:lookup} for  |generate_lookup_table| */
   generate_lookup_table(m);


   std::vector< std::pair< std::vector< unsigned long>, double>> pooljuveniles {} ;
   pooljuveniles.clear() ;
   std::vector< double> viabilityweights {} ;
   viabilityweights.clear() ;

   std::vector< double> v_cdf_number_juveniles_one (GLOBAL_CONST_CUTOFF + 1, 0) ;
   std::vector< double> v_cdf_number_juveniles_two
   (GLOBAL_CONST_CUTOFF + 1, 0) ;
   /* \newline compute the CDF for sampling random number of
   juveniles; see \S~\ref{sec:cdf} for  |cdf_number_juveniles| */
   cdf_number_juveniles( v_cdf_number_juveniles_one, GLOBAL_CONST_ALPHA_ONE );
   cdf_number_juveniles( v_cdf_number_juveniles_two, GLOBAL_CONST_ALPHA_TWO );
   std::vector<double> excursion {} ;
   unsigned long currentNt  {} ;
   int trials = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1 ;
   int timi {} ;
   while( --trials > 0){
   /* \newline see \S~\ref{sec:initarray} for  |initializearray| */
     initializearray( pop );
     timi = 0 ;
     excursion.clear() ;
     assert( excursion.size() < 1) ;

     /* \newline see \S~\ref{sec:clearpool} for |clearpooljuveniles| */
     clearpooljuveniles( pooljuveniles) ;

/* \newline the population evolves until fixed for the fit type at all
sites or lost a fit type at a site; see \S~\ref{sec:checklosttype} for  |not_lost_type|
*/
     while( (pop.back() < current_number_individuals( pop)) & not_lost_type(nhomozygouswt, pop, m) ){
       ++ timi ;
       /* \newline see \S~\ref{sec:Nt} for
       |current_number_individuals| */
       currentNt =  current_number_individuals( pop) ;
       excursion.push_back( static_cast<double>( pop.back() ) / static_cast<double>(currentNt) ) ;

/* \newline  check if a  bottleneck occurs and if so  remove individuals not
surviving a bottleneck; see \S~\ref{sec:bottleneck} for
|remove_not_surviving_bottleneck|  */
       if( gsl_rng_uniform(rngtype) <
       GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
       /* record the number of individuals surviving a bottleneck \newline */
	 currentNt = remove_not_surviving_bottleneck( pop);}
       /* generate a new pool of juveniles; see \S~\ref{sec:newpool}
       for   |new_pool_juveniles| */
       new_pool_juveniles( currentNt,  pooljuveniles,
       viabilityweights, pop, (gsl_rng_uniform(rngtype) <
       GLOBAL_CONST_EPSILON ? v_cdf_number_juveniles_one :
       v_cdf_number_juveniles_two), m);
       /* \newline sort juveniles; see \S~\ref{sec:sortjuveniles} for
       |select_juveniles_according_to_weight| */
       select_juveniles_according_to_weight( viabilityweights, pop,
       pooljuveniles ) ; }
       /* \newline  generated one experiment; record the result;  see \S~\ref{sec:Nt} for
       |current_number_individuals| */
     std::cout << (pop.back() < current_number_individuals(pop) ? 0 : 1) << ' ' << timi << '\n' ;
     if ( pop.back() == current_number_individuals(pop) ){
       /* fixation occurs so record the excursion \newline */
       std::ofstream outfile ("tmpexcurs", std::ios_base::app) ;
       for( const auto &x : excursion){
	 outfile << x << ' ' ; }
       outfile << '\n' ;
       outfile.close() ;
     }
   }
}


@*1 {\bf the main module}. 
\label{sec:main}


the {\tt main} function requires a random seed; so either give a
specific number or   call with
\begin{center}
{\tt ./outfile  \$(shuf -i <range> -n1)}
\end{center}
where {\tt range} is a range of numbers, e.g. {\tt 4343-232383}. 

@C

/* \newline \S~\ref{sec:includes} */
@<Includes@>@#
/* \newline \S~\ref{sec:gslrng} */
@<gsl random number generator@>@#
/* \newline \S~\ref{sec:lookup} */
@<generate lookup table@>@#
/* \newline \S~\ref{sec:lookupfunction} */
@<lookup function@>@#
/* \newline \S~\ref{sec:checkforrecombination} */
@<check for recombination@>@#
/* \newline \S~\ref{sec:recombinanthaplotypes} */
@<recombine haplotypes@>@#
/* \newline \S~\ref{sec:samplehaplotypeindex} */
@<sample haplotype index@>@#
/* \newline \S~\ref{sec:initarray} */
@<initialize population array@>@#
/* \newline \S~\ref{sec:Nt} */
@<current number |Nt|@>@#
/* \newline \S~\ref{sec:viabilityweight} */
@<viability weight@>@#
/* \newline \S~\ref{sec:addingjuv} */
@<addjuv@>@#
/* \newline \S~\ref{sec:kernel} */
@<kernel@>@#
/* \newline \S~\ref{sec:cdf} */
@<cdf@>@#
/* \newline \S~\ref{sec:littersize} */
@<sample litter size@>@#
/* \newline \S~\ref{sec:sibship} */
@<add sibship@>@#
/* \newline \S~\ref{sec:sampleoneindividual} */
@<onehypergeometric@>@#
/* \newline \S~\ref{sec:clearpool} */
@<clear the pool@>@#
/* \newline \S~\ref{sec:newpool} */
@<a new pool@>@#
/* \newline \S~\ref{sec:counthomoz} */
@<count homozygous at site@>@#
/* \newline \S~\ref{sec:comp} */
@<compare@>@#
/* \newline \S~\ref{sec:nth} */
@<$n$th element@>@#
/* \newline \S~\ref{sec:sortjuveniles} */
@<sorting@>@#
/* \newline \S~\ref{sec:checklosttype} */
@<not lost a type@>@#
/* \newline \S~\ref{sec:bottleneck} */
@<generate a bottleneck@>@#
/* \newline \S~\ref{sec:run} */
@<run experiments and record result@>@#



int main(int argc, char * argv[] ){

/* initialise the GSL random number generator \S~\ref{sec:gslrng}
\newline */
setup_rng( static_cast<unsigned long>( atoi(argv[1]) ) );
 

/* free the GSL random number generator \newline */
gsl_rng_free( rngtype);

return GSL_SUCCESS ;

}



@
\end{document}