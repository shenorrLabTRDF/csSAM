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

#ifndef _CellMix_RCPP_H
#define _CellMix_RCPP_H

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP findSigGenes(SEXP rhat, SEXP cutp, SEXP fdr) ;

RcppExport SEXP csSAM_fdr(SEXP rhat, SEXP rhatperm, SEXP alternative
							, SEXP numcell, SEXP nperms, SEXP nperms_total
							, SEXP cutp_g);

#endif
