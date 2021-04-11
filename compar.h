c ------
c compar
c ------
c
c ----------------------------
c main parameters for TRIPLEX
c ----------------------------
c
      integer nlon,    ! longitude dimension of domain
     >        nlat,    ! latitude dimension of domain
     >        npoi,    ! total number of land points
     >        nband,   ! number of solar radiation wavebands
     >        nsoilay, ! number of soil layers
     >        nsnolay, ! number of snow layers
     >        npft,    ! number of plant functional types (PFT)
     >        npftu,    ! number of upper canopy pfts
     >        ncrop    ! number of crops types: only C3 and C4 two types
c
      real xres,       ! longitude resolution (in degrees)
     >     yres,       ! latitude resolution (in degrees)
     >     pi          ! you know, that constant thingy
c
c ------------------------------
c typical TRIPLEX configurations
c ------------------------------
c
c global 0.5 x 0.5 deg Geographic grid: 720 by 360 array, with 58920 land points:
c
      parameter (nlon = 720,
     >           nlat = 360,
     >           npoi = 58920,
     >           xres = 0.50,
     >           yres = 0.50)
c
c global 1.0 x 1.0 deg Geographic grid: 360 by 180 array, with 14545 land points:
c
c      parameter (nlon = 360,
c     >           nlat = 180,
c     >           npoi = 14545,  
c     >           xres = 1.00,
c     >           yres = 1.00)
c
c site-based simualtion Geographic grid: 1 by 1 array, with 1 land point:
c
c      parameter (nlon = 1,
c     >           nlat = 1,
c     >           npoi = 1,
c     >           xres = 0.50,
c     >           yres = 0.50)
c
c -------------------------------
c state description configuration
c -------------------------------
c
      parameter (nband   = 2,
     >           nsoilay = 6,    ! 8
     >           nsnolay = 3,
     >           npft    = 13,           !!13
     >           npftu   = 8,
     >           ncrop   = 2)
c
c
c --------------
c some constants
c --------------
c
      parameter (pi = 3.1415927)
c
c ----------------------
c required common blocks
c ----------------------
c
      real
     >  epsilon,       ! small quantity to avoid zero-divides and other
     >                 ! truncation or machine-limit troubles with small
     >                 ! values. should be slightly greater than o(1)
     >                 ! machine precision
     >  dtime,         ! model timestep (seconds)
     >  zweight,
     >  stef,          ! stefan-boltzmann constant (W m-2 K-4)
     >  vonk,          ! von karman constant (dimensionless)
     >  grav,          ! gravitational acceleration (m s-2)
     >  tmelt,         ! freezing point of water (K)
     >  hfus,          ! latent heat of fusion of water (J kg-1)
     >  hvap,          ! latent heat of vaporization of water (J kg-1)
     >  hsub,          ! latent heat of sublimation of ice (J kg-1)
     >  ch2o,          ! specific heat of liquid water (J deg-1 kg-1)
     >  cice,          ! specific heat of ice (J deg-1 kg-1)
     >  csno,          ! specific heat of snow (J deg-1 kg-1)
     >  cair,          ! specific heat of dry air at constant pressure (J deg-1 kg-1)
     >  cvap,          ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
     >  rair,          ! gas constant for dry air (J deg-1 kg-1)
     >  rvap,          ! gas constant for water vapor (J deg-1 kg-1)
     >  cappa,         ! rair/cair
     >  rhow           ! density of liquid water (all types) (kg m-3)
c
      common /compar1/ epsilon, dtime, zweight, stef, vonk, grav, tmelt, hfus,
     >   hvap, hsub, ch2o,  cice, csno, cair, cvap, rair, rvap, cappa, rhow
c
      real
     >  garea(npoi),   ! area of each gridcell (m**2)
     >  vzero(npoi),   ! a real array of zeros, of length npoi
     >  landcl(npoi)   ! a real array to hold landclass info.c

      common /compar2/ garea, vzero, landcl
c
      integer
     >  ndaypy,        ! number of days per year
     >  nlonsub,       ! number of longitude points for subsetting
     >  nlatsub        ! number of latitude points for subsetting
c
      common /compar3/ ndaypy, nlonsub, nlatsub
c
      integer
     >  ndaypm(12)     ! number of days per month
c
      common /compar4/ ndaypm
c
