/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-5   The R Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#ifndef R_CTEST_H
#define R_CTEST_H

#include <R.h>

void leader_cluster ( double * delta, double * points, double * weights,
                      int * cluster_id, int * nrow, int * ncol, int * type,
                      double * p );


#endif
