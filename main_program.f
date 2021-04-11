c #    #    ##       #    #    #
c ##  ##   #  #      #    ##   #
c # ## #  #    #     #    # #  #
c #    #  ######     #    #  # #
c #    #  #    #     #    #   ##
c #    #  #    #     #    #    #
c
c ---------------------------------------------------------------
      program main
c ---------------------------------------------------------------
c
c common blocks
c 
      include 'implicit.h'
      include 'compar.h'
      include 'comatm.h'
      include 'comdiag.h'
      include 'comagrp.h'
c
c
c local variables
c
      integer icount,        ! number of times ran2 is called (for restart)
     >        icountdum      !
c
      integer istep,         ! timestep counter (per day)
     >        iday,          ! daily loop counter
     >        imonth,        ! monthly loop counter
     >        istyrd,        ! first yr daily means exist (from precip file)
     >        istyrm,        ! first yr monthly anoms exist (from file)
     >        iwest,         ! 1st lon index for subset
     >        iy1,           ! first year for year loop
     >        iy2,           ! last year for year loop
     >        iyear,         ! yearly loop counter
     >        iyear0,        ! first year of simulation
     >        iyranom,       ! year to start reading monthly anomalies
     >        iyrdaily,      ! year to start reading daily means
     >        iyrlast,       ! last year of previous run (for restart)
     >        idiag,         ! number of diagnostic files requested
     >        idailyout,     ! 0: no daily output 1: daily output
     >        imonthout,     ! 0: no monthly output 1: monthly output
     >        iyearout,      ! 0: no yearly output 1: yearly output
     >        irestart,      ! 0: normal mode 1: restart mode
     >        isimveg,       ! 0: static veg 1: dynam veg initialized w/ fixed
     >                       ! 2: dynam veg initialized w/ cold start
     >        isimfire,      ! 0: fixed fire  1: dynamic fire
     >        isimco2,       ! 0: fixed co2   1: changing co2
c     
     >        isim_ac,       ! 0: single age-class   1: multiple age-classes
     >        jday,          ! julian day of the simulation
     >        jnorth,        ! 1st lat index for subset
     >        nanom,         ! # of years in the anomaly files
     >        niter,         ! total number of time iterations per day
     >        nday,          ! number of days since the start of the simulation
     >        nrun,          ! # of years in this run
     >        nspinsoil,     ! year of simulation that soil c reaches equilibrium 
     >        plen,          ! length of precipitation event in timesteps (see plens)
     >        plenmax,       ! upper bound on plen
     >        plenmin,       ! lower bound on plen
     >        seed,          ! value used to initialize the random # generator
     >        spin,          ! counter for iterations used to spin up soilbgc
     >        spinmax,       ! maximum iterations used to spin up soilbgc
     >        eqyears,       !
     >        soilcspin,     ! 0: no spinup procedure for soil c  1: acceleration procedure used
     >        lun, ifile,    ! file indices
     >        istat,         ! NetCDF Error status code
     >        i, n, ivar     ! loop indices
c
      real    co2init,       ! initial co2 concentration in mol/mol
     >        o2init,        ! initial o2 concentration in mol/mol
     >        dran,          ! random number from generator
     >        plens,         ! length of precipitation event in seconds
     >        startp,        ! time to start precipitation event (seconds since midnight)
     >        endp,          ! time to end precipitation event (seconds since midnight)
     >        spincons,      ! # times soil gbc called each day during spinup
     >        spinfrac,      ! fraction of nspinsoil time with max iteration spinmax
     >        slope,         ! rate of decrease of number of iteration for soil C spinup 
     >        snorth,        ! north latitude for subsetting std grid
     >        ssouth,        ! south latitude for subsetting std grid
     >        swest,         ! west longitude for subsetting std grid
     >        seast,         ! east longitude for subsetting std grid
     >        test,          ! test on timestep 
     >        time           ! time of day since midnight (in seconds)
c
c External
c
      real ran2              ! Function : random number generator
c
c number of days per month
c
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data snorth, ssouth, swest, seast / 90., -90., -180., 180. /

c ---------------------------------------------------------------
c       read in parameters from external parameter files
c ---------------------------------------------------------------
c
      call rd_param    
c
c ---------------------------------------------------------------
c                       j o b   c o n t r o l 
c ---------------------------------------------------------------
c
c open local input file called 'TRIPLEX.infile'
c and read in pertinent information
c
      lun = 12
      open (lun, status='old', file='TRIPLEX.infile')
c
      read (lun,*) irestart
      read (lun,*) iyear0
      read (lun,*) nrun
      read (lun,*) iyranom
      read (lun,*) nanom
      read (lun,*) iyrdaily
      read (lun,*) soilcspin
      read (lun,*) iyearout
      read (lun,*) imonthout
      read (lun,*) idailyout
      read (lun,*) isimveg
      read (lun,*) isimfire
c
      read (lun,*) isim_ac
      read (lun,*) isimco2
      read (lun,*) co2init
      read (lun,*) o2init
      read (lun,*) dtime
      read (lun,*) idiag
      read (lun,*,end=99) snorth
      read (lun,*,end=99) ssouth
      read (lun,*,end=99) swest
      read (lun,*,end=99) seast
c
 99   close (lun)
c
c tell user about the simulation
c
      write (*,*) ' '
      write (*,*) '**************************************************'
      write (*,*) '* TRIPLEX model for green-house gas emission GHG *'
c      write (*,*) '* Version 1.0                                    *'
      write (*,*) '* Version 2.0                                    *'
      write (*,*) '* with Trace Green House Gas emission: CH4,N2O,NO*'
      write (*,*) '* Last modified March 2020                       *'
      write (*,*) '**************************************************'
      write (*,*) ' '
c
      if (irestart .eq. 1) then
        write (*,*) 'running in restart mode'
      end if
c
      write (*,*) ' '
      write (*,9000) nrun
      write (*,9005) iyranom
      write (*,*) ' '
      write (*,9010) xres, yres
      write (*,9020) nlon, nlat
      write (*,*) ' '
c
 9000 format (1x,'length of this simulation (years)   : ',i8)
 9005 format (1x,'year to begin using anomalies       : ',i8)
 9006 format (1x,'number of iterations per day        : ',i8)
 9010 format (1x,'model lon, lat resolution (degrees) : ',2f8.2)
 9020 format (1x,'model domain (nlon x nlat)          : ',i5,
     >        ' by ',i5)
 9030 format (1x,'last year run in this sequence      : ',i8)
c
c ---------------------------------------------------------------
c                     t i m e   c o n t r o l 
c ---------------------------------------------------------------
c
c determine the number of timesteps per day
c
      niter = int (86400.0 / dtime)
c
      zweight = exp(-1. / (10.0 * niter))
c
c test on length of the time step and total number of iteration per day
c
      test = amod (86400.0, dtime)
c
      if (test .gt. 1.e-20) then
        write (*,*) 'dtime ', dtime, ' should be divisible into 86400'  
        stop
      else
        write (*,9006) niter
      end if
c
c check for restart conditions and set "iyrlast" accordingly
c
      if (irestart .eq. 1) then
        open (13, status='old', file='yearsrun.dat', err=9050)
        read (13,*) iyrlast
        close (13)
        goto 9059
c
        write (*,*) 'beginning from ', iyear0
        iyrlast = iyear0 - 1
      else
        iyrlast = iyear0 - 1
      end if
c
      write (*,9030) nrun + iyrlast
      write (*,*) ' '
c
c ---------------------------------------------------------------
c                   i n i t i a l i z a t i o n
c ---------------------------------------------------------------
c
c initial concentration of co2, o2
c
      co2conc = co2init
      o2conc  = o2init
      if (isimco2.eq.1) call co2 (co2init, co2conc, iyear0)
c
c read global boundary condition datasets
c
      call readit (isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
c
c check if diagnostic output is requested, if so read info from 'diag.infile'
c
      if (idiag .ne. 0) call inidiag(idiag)
c
c preliminary analysis of climate data
c
      if (irestart .eq. 0) then
        call climanl
      end if
c
c initialize monthly anomalies, daily means if needed

      istyrm = 1961 
      istyrd = 1901
c
c initialize the model
c
      call initial (isimveg, irestart, iyrlast)
c
c initialize random number generator seed
c
      seed = -1
      dran = ran2 (seed)
c
c ---------------------------------------------------------------
c              s t a r t   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
      write (*,*) ' '
      write (*,*) '***************************************'
      write (*,*) '* start of the TRIPLEX-GHG simulation *'
      write (*,*) '***************************************'
      write (*,*) ' '
c
c reset elapsed number of days, accumulating past days if restart mode
c
c start of yearly loop
c
      iy1 = iyrlast + 1
      iy2 = iyrlast + nrun
c
      do 200 iyear = iy1, iy2

      write (*,9101) iyear, co2conc * 1e6
 9101 format (1x, 'Year: ', I4, ' atmospheric CO2 concentration: ',
     >                   F8.3, ' ppm')

c reset julian date
c
        jday = 0
c
c determine the calendar date for this year
c
        ndaypm(2) = 28
        ndaypy = 365
c
c read in manure application and vegetation fraction information
          if (iyear .ge. 1860) then
           call rdmanureN(iyear,iyear0,iwest,jnorth)    !!read Nmanure(1860-2014)
           call rdyear3(iyear,iyear0,iwest,jnorth)      !!read F_crop,F_forest,F_nonveg,F_grass,yndep(1860-2015)
          else
c           call rdmanureN(1860,iyear0,iwest,jnorth)
           call rdyear3(1860,iyear0,iwest,jnorth)
          end if
c
c start of monthly loop
c 
         do 210 imonth = 1, 12
c          call rdanom (imonth, iyear, iyranom, nanom, iy2, istyrm, iwest, jnorth)
c
c read in chemical N fertilizer application         
           if (iyear .ge. 1960) then
               call rdfertilizer (imonth,iyear,iyranom,nanom,iy2,istyrm,iwest,jnorth)
           endif
c
c start of daily loop
c
         do 220 iday = 1, ndaypm(imonth)
c
c update elapsed number of days and julian date
c
            nday = nday + 1
            jday = jday + 1
c
c get daily means (if requested) and determine today's climatic conditions
c
            if (iyear.ge.iyrdaily) then
              call rdday(jday,imonth,iyear,istyrd,iwest,jnorth)
              call daily (iyear, imonth, iday, seed, 1, iyranom, iyrlast, nrun)
            else
              call rdday(jday,imonth,1901,istyrd,iwest,jnorth)
              call daily (1901, imonth, iday, seed, 1, iyranom, iyrlast, nrun)
            end if
c
c determine the daily vegetation cover characteristics
c
            call pheno
c
c call soil biogeochemistry model
c
c soil spin up procedure calls soil bgc model spincons times
c for each time step to spin up the process (accelerate C & N pools)
c
            spinfrac  = 0.75          ! fraction of nspinsoil time used to
c                                     ! spin up soil at maximum spin up rate
c
            spincons  = 40.0          ! number of times soilbgc subroutine is
c                                     ! called for each day in simulation during
c                                     ! acceleration procedure
c
            eqyears   = 50            ! years used past spin up to slowly bring
c                                     ! entire soil to a steady equilibrium
c
            nspinsoil = iyear0 + 150  ! year of simulation that soil c reaches equilibrium
c
            if (soilcspin .eq. 1) then

              if ((iyear - iyear0) .le.
     >            (spinfrac * (nspinsoil - iyear0 - eqyears))) then
                 spinmax = int(spincons)
c
              else if ((iyear - iyear0) .lt.
     >               (nspinsoil - iyear0 -  eqyears)) then
c
                slope   = spincons / ((nspinsoil - iyear0 - eqyears) -
     >                    (spinfrac * (nspinsoil - iyear0 - eqyears)))
c
                spinmax = int (spincons - (slope * ((iyear - iyear0) -
     >            (spinfrac * (nspinsoil - iyear0 - eqyears)))))
c
                spinmax = max(spinmax,1)
c
              else
c
                spinmax = 1
c
              endif
c
            else 
c
                spinmax = 1
c
            endif
c
            do 230 spin = 1, spinmax
              call soilbgc (iyear, iyear0, imonth, iday, nspinsoil,
     >                      spin, spinmax)
 230        continue
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
            plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
            plenmax = max (int (24.0 * 3600. / dtime), plenmin)
c
            plen    = min (plenmax, int (
     >                     plenmin + ran2(seed) * (plenmax-plenmin+1) ))
c
            plens   = dtime * plen
            startp  = dtime * min (niter-plen, 
     >                             int(ran2(seed)*(niter-plen+1)))
            endp    = startp + plens
c
c start of hourly loop
c
            do 240 istep = 1, niter
c
c calculate the time of day since midnight (in seconds)
c
              time = (istep - 1) * dtime
c
c determine climatic conditions for the given timestep
c
              call diurnal (time, jday, plens, startp, endp, seed)
c
c call the land surface model
c
              call lsxmain
              
c call new irrigation submodule (if necessary) to simulate GHG mostly
c attention! this module is in hourly loop and we simply assume that irrigation happens at 12:00 only.
              if (istep .eq. 12) then
                 call irrigation (iyear, imonth, iday, nday, iyear0)
              end if
c           
c accumulate some variables every timestep
c
              call sumnow
c
              call sumday   (istep,plens)
c call new tillage practice subsubmodule (if necessary) to simulate GHG mostly
              if (istep .eq. niter) then
                 call tillage (iyear,imonth,iday,nday,iyear0)     
              end if
c
c call green house gas emission
              if (istep .eq. niter) then
                 call GreenHouseGassEmission (iyear,imonth,iday,nday,iyear0)       
              end if
c calculate green house gas emission and begin summary for different timesteps
c          
              call summonth (istep, iday, imonth)
              call sumyear  (istep, iday, imonth)
c
c write out diagnostics
c
              if (idiag .ne. 0) then
c
                do 250 i = 1, idiag
                  if (iyear.ge.diagstart(i) .and.
     >                iyear.le.diagend(i)   .and. 
     >                ndiagpt(i).ge.1)  then  
                    if (mod(istep,nfreq(i)).eq.0) then
                      call wdiag (i, iyear, imonth, iday, istep)
                    end if
                  end if
 250            continue
c
              end if
c
c ---------------------------------------------------------------
c               e n d   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
c end of the hourly loop
c
 240        continue
c
c write out daily output
            if (idailyout.eq.1) then
              write (*,*) ' '
              write (*,*) 'writing daily output'
              write (*,*) ' '
              call wdaily (nday, iyear, iyear0) 
            endif
c
c end of the daily loop
c

c        call vdiag (iyear, iyear0,imonth,iday)

 220      continue
c
c write out monthly output (use 1st day of this month because of GrADS bug)
c
          if (imonthout.eq.1) then
            write (*,*) ' '
            write (*,*) 'writing monthly output for ', iyear, '/', imonth
            write (*,*) ' '
            call wmonthly (nday-ndaypm(imonth)+1, imonth, iyear, 
     >                     iyear0)
          endif
c
c end of the monthly loop
c
 210    continue
c
c recalculate bioclimatic parameters if using anomalies
c
        if (iyear .ge. iyranom-1 .or. iyear.ge.iyrdaily) call climanl2
c
c get new o2 and co2 concentrations for this year
c
        if (isimco2.eq.1) call co2 (co2init, co2conc, iyear)
        if (isimco2.eq.2) call co2_constafter2000 (co2init, co2conc, iyear)
c
c perform vegetation dynamics
c 
**** DTP 2002/05/12: Get disturbance data from external data file
c        if ((isim_ac .EQ. 1) .AND. ( mod(iyear, ageint) .EQ. 0)) then
c           call getdist (iyear, iwest, jnorth)
c        end if
        if (isimveg.ne.0) call dynaveg (isimfire, isim_ac, iyear)
c
c calculate simple annual diagnostic variables for the globe
c
c        call gdiag (iyear, iyear0)
c
c calculate simple annual diagnostic variables for each vegetation type
c
        call vdiag (iyear, iyear0)
c
c write out annual output
c
        if (iyearout.eq.1) then
          write (*,*) ' '
          write (*,*) 'writing annual output',iyear
          write (*,*) ' '
          call wyearly (nday, iyear, iyear0)
        endif
c
c update iyrlast value and file yearsrun.dat
c
        iyrlast = iyrlast + 1
        open(12,file='yearsrun.dat',status='unknown')
        write(12,*) iyrlast
        close(12)
c
c end of the yearly loop
c
 200  continue
c
c end of the simulation
c
      write (*,*) ' '
      write (*,*) '*** end of run ***'
      write (*,*) ' '
c
      stop
      end
c