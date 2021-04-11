c comdiag.h last update 2/12/01 C. Molling
c
c This holds the parameters and common blocks for diagnostics
c
c
      integer nvars, nfiles
      parameter (nvars=119, nfiles=10)   ! nvars = # diagnostic variables
c
      integer ldiag(nvars,nfiles)      ! chosen diagnostic variables
c
      common /diag1/ ldiag
c
c
      integer diagstart(nfiles),      ! year diagnostic output begins
     >        diagend(nfiles),        ! year diagnostic output stops
     >        ndiagpt(nfiles),        ! point in an npoi array
     >        nfreq(nfiles)           ! frequency of diagnostic output
c
      common /diag2/diagstart,diagend,ndiagpt,nfreq
c
