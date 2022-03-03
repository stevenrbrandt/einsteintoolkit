        o1 = 5.0000000000000d-1*eta(i,j,k)
        o2 = exp(o1)
        o3 = psi2d(i,j,k)
        o4 = 1/o3
        o5 = cos(phi(i,j,k))
        o6 = o5**2
        o7 = cos(q(i,j,k))
        o8 = o7**2
        o9 = detapsi2d(i,j,k)
        o10 = -2.50000000000000d0*eta(i,j,k)
        o11 = exp(o10)
        o12 = dqqpsi2d(i,j,k)
        o13 = sin(phi(i,j,k))
        o14 = o13**2
        o15 = detaqpsi2d(i,j,k)
        o16 = sin(q(i,j,k))
        o17 = dqpsi2d(i,j,k)
        o18 = detaetapsi2d(i,j,k)
        o19 = o16**2
        o20 = tan(q(i,j,k))
        o21 = 1/o20
        psixx(i,j,k) = o2*o4*(1.00000000000000d0*o11*o14*o17*o21 - 5.000
     &  0000000000d-1*o11*o14*o3 + o11*o18*o19*o6 + 7.5000000000000d-1*o
     &  11*o19*o3*o6 + 2.00000000000000d0*o11*o15*o16*o6*o7 - 3.00000000
     &  000000d0*o11*o16*o17*o6*o7 + o11*o12*o6*o8 - 5.0000000000000d-1*
     &  o11*o3*o6*o8 + o11*o14*o9 - 2.00000000000000d0*o11*o19*o6*o9 + o
     &  11*o6*o8*o9)
        psixy(i,j,k) = o2*o4*(o11*o13*o18*o19*o5 - o11*o13*o17*o21*o5 + 
     &  5.0000000000000d-1*o11*o13*o3*o5 + 7.5000000000000d-1*o11*o13*o1
     &  9*o3*o5 + 2.00000000000000d0*o11*o13*o15*o16*o5*o7 - 3.000000000
     &  00000d0*o11*o13*o16*o17*o5*o7 + o11*o12*o13*o5*o8 - 5.0000000000
     &  000d-1*o11*o13*o3*o5*o8 - o11*o13*o5*o9 - 2.00000000000000d0*o11
     &  *o13*o19*o5*o9 + o11*o13*o5*o8*o9)
        psixz(i,j,k) = o2*o4*(-(o11*o15*o19*o5) + 1.50000000000000d0*o11
     &  *o17*o19*o5 - o11*o12*o16*o5*o7 + o11*o16*o18*o5*o7 + 1.25000000
     &  000000d0*o11*o16*o3*o5*o7 + o11*o15*o5*o8 - 1.50000000000000d0*o
     &  11*o17*o5*o8 - 3.00000000000000d0*o11*o16*o5*o7*o9)
        psiyy(i,j,k) = o2*o4*(o11*o14*o18*o19 + 7.5000000000000d-1*o11*o
     &  14*o19*o3 + 1.00000000000000d0*o11*o17*o21*o6 - 5.0000000000000d
     &  -1*o11*o3*o6 + 2.00000000000000d0*o11*o14*o15*o16*o7 - 3.0000000
     &  0000000d0*o11*o14*o16*o17*o7 + o11*o12*o14*o8 - 5.0000000000000d
     &  -1*o11*o14*o3*o8 - 2.00000000000000d0*o11*o14*o19*o9 + o11*o6*o9
     &   + o11*o14*o8*o9)
        psiyz(i,j,k) = o2*o4*(-(o11*o13*o15*o19) + 1.50000000000000d0*o1
     &  1*o13*o17*o19 - o11*o12*o13*o16*o7 + o11*o13*o16*o18*o7 + 1.2500
     &  0000000000d0*o11*o13*o16*o3*o7 + o11*o13*o15*o8 - 1.500000000000
     &  00d0*o11*o13*o17*o8 - 3.00000000000000d0*o11*o13*o16*o7*o9)
        psizz(i,j,k) = o2*o4*(o11*o12*o19 - 5.0000000000000d-1*o11*o19*o
     &  3 - 2.00000000000000d0*o11*o15*o16*o7 + 3.00000000000000d0*o11*o
     &  16*o17*o7 + o11*o18*o8 + 7.5000000000000d-1*o11*o3*o8 + o11*o19*
     &  o9 - 2.00000000000000d0*o11*o8*o9)