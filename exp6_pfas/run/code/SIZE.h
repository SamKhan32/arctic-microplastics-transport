#include "CPP_OPTIONS.h"
CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
CEOP
        INTEGER sNx
        INTEGER sNy
        INTEGER OLx
        INTEGER OLy
        INTEGER nSx
        INTEGER nSy
        INTEGER nPx
        INTEGER nPy
        INTEGER Nx
        INTEGER Ny
        INTEGER Nr
        INTEGER MAX_OLX
        INTEGER MAX_OLY
        PARAMETER (
     &            sNx =   27,
     &            sNy =   30,
     &            OLx =    2,
     &            OLy =    2,
     &            nSx =    1,
     &            nSy =    1,
     &            nPx =   10,
     &            nPy =    9,
     &            Nr  =    1 )
        PARAMETER (
     &            Nx  = sNx*nSx*nPx,
     &            Ny  = sNy*nSy*nPy )
        PARAMETER ( MAX_OLX = OLx )
        PARAMETER ( MAX_OLY = OLy )
