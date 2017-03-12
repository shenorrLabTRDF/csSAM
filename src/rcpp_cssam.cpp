// csSAM is an R package that performs partial gene expression deconvolution
// and estimate cell-specific gene differential expression.
// Copyright (C) 2011-2012  Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
// Copyright (C) 2012 Renaud Gaujoux (implementation in C++ of some of the functions)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "rcpp_cssam.h"
#include <R.h>
#include <Rdefines.h>
#include <iostream>

/**
 * Compute significant genes from csSAM FDR results
 */
SEXP findSigGenes(SEXP rhat, SEXP cutp, SEXP fdr){
    using namespace Rcpp ;

    // rhat is ncell x ngenes
    NumericMatrix m_rhat(rhat);
    int numgene = m_rhat.ncol();
    int numcell = m_rhat.nrow();

    NumericMatrix m_cutp(cutp);
    NumericMatrix m_fdr(fdr);
    int thresholdLen = m_fdr.ncol();

    // initialise result matrix full of 1s
    NumericMatrix sigGene(numcell, numgene);
    std::fill(sigGene.begin(), sigGene.end(), 1);

    // find genes
    for(int curThresh=0; curThresh<thresholdLen; ++curThresh){
    		for (int curcell=0; curcell<numcell; ++curcell) {
    			for (int curgene=0; curgene<numgene; ++curgene) {
    				if ( m_rhat(curcell, curgene) >= m_cutp(curcell,curThresh) ) {
    					sigGene(curcell, curgene) = m_fdr(curcell,curThresh);
    				}
    			}
    		}
    	}

    return(wrap(sigGene));
}

/*
 csSAM_fdr <- function(rhat, rhatperm, alternative, numcell = nrow(rhat), nperms = dim(rhatperm)[1L]){

    vmessage("Alternative '", alternative,"' ... ", appendLF=FALSE)
    cutp.g=matrix(NA,nrow=numcell,ncol=100)
    numcut = ncol(cutp.g)

    fdr.g=ncall.g=nperm.g<-array(dim = c(numcell, numcut))

 for(j in 1:numcell)
	cutp.g[j,]=seq(0,max(abs(rhat[j,])),length=100)
	for (i in 1:numcut) {
		for (curcell in 1:numcell) {
			if(alternative == 'two.sided') {
				fdr.g[curcell,i]=sum(abs(rhatperm[,curcell,])>cutp.g[curcell,i])/nperms /sum(abs(rhat[curcell,])>cutp.g[curcell,i])
                				ncall.g[curcell,i]=sum(abs(rhat[curcell,])>cutp.g[curcell,i])
			}
			if(alternative == 'greater') {
				fdr.g[curcell,i]=sum(rhatperm[,curcell,]>cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
                				ncall.g[curcell,i]=sum(rhat[curcell,]>cutp.g[curcell,i])
			}
			if(alternative == 'less') {
# [RG] BUG FIX: should be < - cutp.g[curcell,i]
#			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
				fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,] < -cutp.g[curcell,i])
                				ncall.g[curcell,i]=sum(rhat[curcell,]< -cutp.g[curcell,i])
			}
		}
	}
*/

#define PERMUTATION_COUNTS(comparison) \
{ \
	if( n_obs > 0 ){\
		int n_in_perms = 0; \
		const double* perm_val = p_rhatperm + curcell * Nperms; \
		for(int l=0; l<Ngene; ++l){ \
			for(int j=0; j<_nperms; ++j){ \
				if( comparison ){\
					++n_in_perms; \
				} \
			} \
			perm_val += perm_step; \
		}\
		fdr_g(curcell, i) = ((double) n_in_perms) / _nperms / n_obs;\
	} else fdr_g(curcell, i) = INFINITY; \
}

/**
 * Computes FDR from csSAM permutations
 */
SEXP csSAM_fdr(SEXP rhat, SEXP rhatperm, SEXP alternative
				, SEXP numcell, SEXP nperms, SEXP nperms_total
				, SEXP cutp_g){
    using namespace Rcpp ;
    BEGIN_RCPP

    List ret;

	// input
    NumericMatrix _rhat(rhat);
    int Ngene = _rhat.ncol();
    int Ncell = _rhat.nrow();
    int Nperms = NumericVector(nperms_total)[0];
    const int perm_step = Nperms * Ncell;
    // direct pointer to 3D-array that contains the permutation results
    const double* p_rhatperm = NUMERIC_POINTER(rhatperm);
    const std::string _alternative = (const char*) CharacterVector(alternative)[0];
    int _numcell = NumericVector(numcell)[0];
    int _nperms = NumericVector(nperms)[0];

    //    	Rprintf("numcell = %i/%i, numcut= %i, Ngene = %i, nperm=%i/%i\n"
    //    			, _numcell, Ncell, numcut, Ngene, _nperms, Nperms);

    // init data
    const NumericMatrix _cutp_g(cutp_g);
    int numcut = _cutp_g.ncol();
    NumericMatrix fdr_g(_numcell, numcut);
    std::fill(fdr_g.begin(), fdr_g.end(), NA_REAL);
    IntegerMatrix ncall_g(_numcell, numcut);
    std::fill(ncall_g.begin(), ncall_g.end(), NA_INTEGER);

    // main loops
    if( _alternative == "two.sided" ){

    	for ( int curcell = 0; curcell < _numcell; ++curcell) {
    		const NumericVector v_obs = abs( _rhat(curcell, _ ) );
    		for (int i=0; i<numcut; ++i) {
    			int n_obs = sum( v_obs > _cutp_g(curcell, i) );
    			// compute hits from permutations
    			PERMUTATION_COUNTS( fabs(perm_val[j]) > _cutp_g(curcell, i) )
    			ncall_g(curcell, i) = n_obs;
    		}
    	}
    }else if( _alternative == "greater" ){
    	for ( int curcell = 0; curcell < _numcell; ++curcell) {
    		const NumericVector v_obs = _rhat(curcell, _ );
    		for (int i=0; i<numcut; ++i) {
    			int n_obs = sum( v_obs > _cutp_g(curcell, i) );
    			// compute hits from permutations
    			PERMUTATION_COUNTS( perm_val[j] > _cutp_g(curcell, i) )
    			ncall_g(curcell, i) = n_obs;
    		}
    	}
    }else if( _alternative == "less" ){
    	for ( int curcell = 0; curcell < _numcell; ++curcell) {
    		const NumericVector v_obs = _rhat(curcell, _ );
    		for (int i=0; i<numcut; ++i) {
    			int n_obs = sum( v_obs < - _cutp_g(curcell, i) );
    			// compute hits from permutations
    			PERMUTATION_COUNTS( perm_val[j] < - _cutp_g(curcell, i) )
    			ncall_g(curcell, i) = n_obs;
    		}
    	}
    }else{
    	std::stringstream s;
    	s << "Invalid alternative '" << _alternative << "'";
    	stop(s.str());
    }

    ret["fdr"] = fdr_g; ret["ncall"] = ncall_g;
    return( ret );

    END_RCPP
}

