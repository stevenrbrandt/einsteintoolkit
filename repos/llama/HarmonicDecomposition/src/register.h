
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef _HDECOMP_
#define _HDECOMP_

#include <vector>

namespace HDecomp {

using namespace std;

// the slice ids for each slice and variable
extern vector<vector<int> > sid;

// the slice ids for each slice and precalculated sYlm
//extern vector<vector<int> > sid_sYlm_re;
//extern vector<vector<int> > sid_sYlm_im;



}


#endif
