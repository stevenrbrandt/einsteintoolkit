        o1 = cos(qgrd(j))
        o2 = 1.0
        o3 = o2**2
        o4 = 1/o3
        o5 = sin(qgrd(j))
        o6 = 2.0 + n
        o7 = o5**2
        o8 = -(o6*o7)
        o9 = 1.0 + n + o8
        o10 = o3**2
        o11 = o10**2
        o12 = dngtil**2
        o13 = o1**2
        o14 = gtil**2
        o15 = o9**2
c        exc33(i,j) = 0
c        exc32(i,j) = -(dngtil*o1*o4*o5**n)
c        exc31(i,j) = gtil*o4*o5**(-1.0 + n)*o9
c        exc22(i,j) = 0
c        exc21(i,j) = 0
c        exc11(i,j) = 0
        ksq(i,j) = (2.0*o5**(2*(
     $	-2 + n))*(o14*o15 + o12*o13*o7))/(o10*o11)