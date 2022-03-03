
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

// This file has been generated automatically.  Do not change it manually.
J[0][0]=0.;
J[0][1]=-(zp/(yp*yp+zp*zp));
J[0][2]=yp/(yp*yp+zp*zp);
J[1][0]=yp/(xp*xp+yp*yp);
J[1][1]=-(xp/(xp*xp+yp*yp));
J[1][2]=0.;
J[2][0]=(-1.*(ri-rm)*rm*xp*yp*(ri+yp)*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-\
2.))/sqrt(xp*xp+yp*yp+zp*zp);
J[2][1]=(-1.*(ri-rm)*(rm*pow(yp,3.)+ri*(-(rm*(xp*xp+zp*zp))+pow(xp*xp+yp*yp+\
zp*zp,1.5)))*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-2.))/sqrt(xp*xp+yp*yp+zp\
*zp);
J[2][2]=(-1.*(ri-rm)*rm*yp*(ri+yp)*zp*pow(rm*yp+ri*sqrt(xp*xp+yp*yp+zp*zp),-\
2.))/sqrt(xp*xp+yp*yp+zp*zp);
