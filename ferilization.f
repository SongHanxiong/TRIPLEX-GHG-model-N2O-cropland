c  -------------------------------------------------------------
      subroutine fertilization (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to simulate the process associated with chemical
c fertilizer application to soil. This subroutine is called in the GHG submodule (GHGasE.f).
c The major framworks are adapted from DNDC model.

      include 'implicit.h'
      include 'compar.h'
      include 'comatm.h'
      include 'comdiag.h'
      include 'comagrp.h'

c
      integer i,        ! gridpoint for which
     >        iyear,    ! current calendar year
     >        iyear0,   ! first year
     >        imonth,   ! current month
     >        iday,     ! current date in the month
     >        k,
c
     >        nday    !nday for debug
     
c###songhx### 2019/4/29 add these lines because the fertilizer dataset starts from 1961.
       if (iyear. ge. 1961) then
          experimentF(i) = 1.0
       else
          experimentF(i) = 0.0
       end if
c###songh###
c###songh### 2019/5/1 decide which day the fertlizer is applied 
c            (from published paper information or fertlization dataset).
c
        if (ferdate(i,imonth) .ge. 1 .and. ferdate(i,imonth) .le. 31) then
           ferday(i) = ferdate(i,imonth)
        else if (ferdate(i,imonth) .ge. 32 .and. ferdate(i,imonth) .le. 59) then
           ferday(i) = ferdate(i,imonth) - 31
        else if (ferdate(i,imonth) .ge. 60 .and. ferdate(i,imonth) .le. 90) then
           ferday(i) = ferdate(i,imonth) - 59
        else if (ferdate(i,imonth) .ge. 91 .and. ferdate(i,imonth) .le. 120) then
           ferday(i) = ferdate(i,imonth) - 90
        else if (ferdate(i,imonth) .ge. 121 .and. ferdate(i,imonth) .le. 151) then
           ferday(i) = ferdate(i,imonth) - 120
        else if (ferdate(i,imonth) .ge. 152 .and. ferdate(i,imonth) .le. 181) then
           ferday(i) = ferdate(i,imonth) - 151
        else if (ferdate(i,imonth) .ge. 182 .and. ferdate(i,imonth) .le. 212) then
           ferday(i) = ferdate(i,imonth) - 181
        else if (ferdate(i,imonth) .ge. 213 .and. ferdate(i,imonth) .le. 243) then
           ferday(i) = ferdate(i,imonth) - 212
        else if (ferdate(i,imonth) .ge. 244 .and. ferdate(i,imonth) .le. 273) then
           ferday(i) = ferdate(i,imonth) - 243
        else if (ferdate(i,imonth) .ge. 274 .and. ferdate(i,imonth) .le. 304) then
           ferday(i) = ferdate(i,imonth) - 273
        else if (ferdate(i,imonth) .ge. 305 .and. ferdate(i,imonth) .le. 334) then
           ferday(i) = ferdate(i,imonth) - 304
        else if (ferdate(i,imonth) .ge. 335 .and. ferdate(i,imonth) .le. 365) then
           ferday(i) = ferdate(i,imonth) - 334
        end if
c set the application indicator
c
        if (iday .eq. ferday(i)) then
           applyfer(i) = 1.0
        else 
           applyfer(i) = 0.0
        end if
c check the amount of N fertilizer to prevent unrealistic value
c also, values could be set according to published papers. 
         NH4fer(i,imonth) = max(NH4fer(i,imonth),0.0)
         NO3fer(i,imonth) = max(NO3fer(i,imonth),0.0)
c
c add the applied fertilizer into the N pools of the top soil layer (1st layer only)
c
         NH4_add(i) = NH4fer(i,imonth)*applyfer(i)*experimentF(i)
         NO3_add(i) = NO3fer(i,imonth)*applyfer(i)*experimentF(i)
c
         NH4(i,1) = NH4(i,1) + NH4_add(i)
         NO3(i,1) = NO3(i,1) + NO3_add(i)
c
     end
c
c 
c read in chemical fertilizer application data, this data is in monthly time step including applied amount of 
c NH4 and NO3 (kgN ha-1), application date and NH4/NO3 ratio. This dataset is provided by Nishina et al. 2017.
c ---------------------------------------------------------------------
      subroutine rdfertilizer(imonth,iyear,iyranom,nanom,iy2,istyr,iwest,jnorth)
c ---------------------------------------------------------------------
c###songhx### 2019/3/13 today is end of an area!!! Prof.Zhu helped me out with this subroutine and realized that
c             it's because of miss use 'nanom' in ibisinfile which should = 50, the number of years in anom file
c             should be more careful with these.
      include 'implicit.h'
c
      include 'compar.h'
      include 'combcs.h'
      include 'comveg.h'
      include 'comwork.h'
c
c Arguments
c
      integer imonth, ! month
     >        iyear,  ! year
     >        iyranom,! year to start reading anoms
     >        iy2,    ! final year of the run
     >        istyr,  ! 1st year in data files
     >        iwest,  ! 1st lon index for subset
     >        jnorth, ! 1st lat index for subset
     >        nanom   ! # of years in the anomaly files
c
c Local variables
c
      real anom(npoi)
      equivalence(anom(1),cdummy(1))
c
      integer imon,   ! month + 1
     >        iyr,    ! number of years after begin of data file(istyr)
     >        istat,  ! error flag for netcdf
     >        i,      ! loop indice on land points
     >        jyear   ! iyr divided by nanom (for looping through anomaly files)
c
      integer istart(3), icount(3) ! forreading rdanom vars

      character*80 filen
c
c ---------------------------------------------------------------------
c
      data istart / 1,1,1 /, icount / nlon,nlat,1 /
      istart(1) = iwest
      istart(2) = jnorth
      icount(1) = nlonsub
      icount(2) = nlatsub
c
c determine which (if any) month to read
c If prior to November of the year before an anomaly year, then return.
c Else read anomalies for the following month.
c
      if (iyear .lt. iyranom-1) then
         return
      else if (iyear .eq. iyranom-1 .and. imonth .lt. 11) then
         return
c
c If timestep equals December of final year in anomaly file and also equals
c the final year of the run, then set January anomalies to zero and return.
c If not the final year of the run, then loop back to start of anomaly file
c (see jyear).
c
      else if (iyear .eq. istyr+nanom-1 .and. imonth .eq. 12) then  
        if (iyear .eq. iy2) then
c
          print *, 'WARNING: last month of run; no anomalies for January'
          print *, 'Using fertilization for month year ='
          print *, imonth+1,iyear+1
          do 4 i = 1, npoi
            NH4fer(i,1) = 0.0
            NO3fer(i,1) = 0.0
 4        continue
          return
        end if
      end if
c
      iyr = iyear-istyr
      imon = imonth + 1
      if (imon .eq. 13) then
         imon = 1
         iyr = iyr + 1
      end if
      jyear = iyr/nanom
      istart(3) = (iyr - nanom*jyear)*12 + imon

      if (iyr.gt.0 .and. (iyr - nanom*jyear).eq.0) then
        print *, 'WARNING: Attempted to read past last month in anomaly file'
        print *, 'Looping back to the beginning of the file'
      end if

      if (istart(3) .gt. 0) then
        print *, 'rdanom reading month year ',imon,iyr+istyr-nanom*jyear,
     >           '. Step = ', istart(3)
      else
        print *, 'WARNING, anomalies begin in year ',istyr
        print *, 'Not reading in anomalies for month year ='
        print *, imon,iyr+istyr
        return
      end if
      
c###songhx### read in the applied monthly amount of NH4-N

      aname = 'NH4fer'
      filen = '/type in local file path/NH4_fertilization2015.nc'
      call readvar(filen,aname,'',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdfertilizer, NH4_fertilization_try.nc'
         stop 1
      end if
      call arr2vec( work, anom)
c      NH4fer(i,imonth) = max (anom(i), 0.0)
      do 10 i = 1, npoi
         NH4fer(i,imon) = max (anom(i), 0.0)
10    continue
c      write(*,*) 'NH4fer =', NH4fer(i,imonth)
 
c###songhx### read in the applied monthly amount of NO3-N
c      
      aname = 'NO3fer'
      filen = '/type in local file path/NO3_fertilization2015.nc'
      call readvar(filen,aname,'',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdfertilizer, NO3_fertilization_try.nc'
         stop 1
      end if
      call arr2vec( work, anom)
c      NO3fer(i) = max (NO3fer(i), 0.0)
      do 11 i = 1, npoi
         NO3fer(i,imon) = max (anom(i), 0.0)
11      continue
c      write(*,*) 'NO3fer =', NO3fer(i,imonth)
      
c###songhx### read in the fertilizer application date in monthly run
      aname = 'ferdate'
      filen = '/type in local file path/FerApplyDate_month.nc'
      call readvar(filen,aname,'',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > cdummy(3*nlonsub+1),istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in FerApplyDate_month_new.nc'
         stop 1
      end if
      call arr2vec( work, anom)
c
      do 12 i = 1, npoi
         ferdate(i,imon) = max (anom(i), 0.0)
12      continue
c      write(*,*) 'ferdate =', ferdate(i,imonth)

      return
      end
c
c return to main program
c
      return
      end
c             