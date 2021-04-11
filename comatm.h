c
c ------
c comatm
c ------
c      
      real 
     >  coszen(npoi),      ! cosine of solar zenith angle
     >  fira(npoi)         ! incoming ir flux (W m-2)
c
      common /comatm1/ coszen, fira
c
      real 
     >  solad(npoi,nband), ! direct downward solar flux (W m-2)
     >  solai(npoi,nband), ! diffuse downward solar flux (W m-2)
     >  asurd(npoi,nband), ! direct albedo of surface system
     >  asuri(npoi,nband)  ! diffuse albedo of surface system 
c
      common /comatm2/ solad, solai, asurd, asuri
c
      real 
     >  ua(npoi),          ! wind speed (m s-1)
     >  ta(npoi),          ! air temperature (K)
     >  qa(npoi),          ! specific humidity (kg_h2o/kg_air)
     >  raina(npoi),       ! rainfall rate (mm/s or kg m-2 s-1)
     >  rh(npoi),          ! relative humidity(%)
     >  snowa(npoi)        ! snowfall rate (mm/s or kg m-2 s-1 of water)
c
      common /comatm3/ ua, ta, qa, raina, rh, snowa
c
      real 
     >  psurf(npoi),       ! surface pressure (Pa)
     >  cloud(npoi),       ! cloud fraction
     >  td(npoi),          ! daily average temperature (K)
     >  tmax(npoi),        ! maximum daily temperature (K)
     >  tmin(npoi),        ! maximum daily temperature (K)
     >  qd(npoi),          ! daily average specific humidity (kg_h2o/kg_air)
     >  ud(npoi),          ! daily average wind speed (m/sec)
     >  precip(npoi),      ! daily precitation (mm/day)
     >  precipday(npoi,31),       
     >  precipdaysum(npoi)       
c
      common /comatm4/ psurf, cloud, td, tmax, tmin, qd, ud, precip,
     >                   precipday, precipdaysum
c
      real 
     >  xstore(npoi,3)     ! weather generator 'memory' matrix
c
      common /comatm5/ xstore
c
      integer 
     >  iwet(npoi),        ! wet day / dry day flag
     >  iwetday(npoi,31), 
     >  iwetdaysum(npoi) 
c
      common /comatm6/ iwet, iwetday, iwetdaysum
c
      real 
     >  co2conc,           ! co2 concentration (mol/mol)
     >  o2conc             ! o2 concentration (mol/mol)
c
      common /comatm7/ co2conc, o2conc
c
