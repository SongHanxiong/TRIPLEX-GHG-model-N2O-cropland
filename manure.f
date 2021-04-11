c  -------------------------------------------------------------
      subroutine manureapp (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to simulate the process associated with manure
c  application to soil. This subroutine is called in the GHG submodule (GHGasE.f).
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
     
c###songhx### 2019/4/29 add these lines because the manure dataset starts from 1861.
       if (iyear. ge. 1861) then
             experimentM(i) = 1.0
       else
             experimentM(i) = 0.0
       end if
c###songh###
c calculate the manure chemical properties (i.e., NH4, NO3 and organic carbon and nitrogen).
c becasue the applied N amount is settle, we re-calculate the manure sourced C using the C:N ratios
c of different soil organic matter pools for further decomposition to release souable nurtient 
c within the biogeochmstry subroutine (biogeochem.f).
c
         manure_N(i)= max(0.0,manure_N(i))
c the amount of manure_N can also be set according to published information.
c###songh### the unit of input manure_N from dataset (Zhang et al. 2017) is kgN km-2 yr-1
c         
         NH4_manure(i)= max(0.0, COEMnh4 * manure_N(i))                   !manure sourced NH4 
         NO3_manure(i)= max(0.0, COEMno3 * manure_N(i))                  !manure sourced NO3
c         
c###songhx### notably the value of parameter COEMnh4 is read in from parameter file
c             meanwhile, it can be set according to published information when available.
c
        manureC2pas(i) = COEMpas * manure_N(i) * cnsoipas(i)/1000000.0  !manure_C to passive SOM 
        manureC2slop(i) = COEMslop * manure_N(i) * cnsoislop(i)/1000000.0  !manure_C to protected SOM
        manureC2slon(i) = COEMslon * manure_N(i) * cnsoislon(i)/1000000.0  !manure_C to non-protected SOM
        manure_C(i) = manureC2pas(i)+manureC2slop(i)+manureC2slon(i)  !total manure_C
        csoislon(i) = csoislon(i) + max(manureC2slon(i),0.0)
        csoipas(i) = csoipas(i) + max(manureC2pas(i),0.0)
        csoislop(i) = csoislop(i) + max(manureC2slop(i),0.0)
c        
c these results are passed to biogeochem subroutine
c
c set the application indicator
c
        if (imonth .eq. 'type in month' .and. iday .eq. 'type in day') then
           applyman(i) = 1.0
        else 
           applyman(i) = 0.0
        end if
c
c add the applied manure sourced NH4 and NO3 into the N pools of the top soil layer (1st layer only)
c
         NH4_add(i) = NH4_manure*applyman(i)*experimentM(i)
         NO3_add(i) = NO3_manure*applyman(i)*experimentM(i)
c
         NH4(i,1) = NH4(i,1) + NH4_add(i)
         NO3(i,1) = NO3(i,1) + NO3_add(i)
c
     end
c
c 
c read in manure application to cropland data, this data is in yearly time step.
c This dataset is provided by Zhang et al. 2018.
c ---------------------------------------------------------------------
      subroutine rdmanureN (iyear,iyear0,iwest,jnorth)
c ---------------------------------------------------------------------
c
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
     >        nanom,   ! # of years in the anomaly files
     >        iyear0
c
c Local variables
c
      real yyear(npoi)
      equivalence(yyear(1),cdummy(1))
c
      integer imon,   ! month + 1
     >        iyr,    ! number of years after begin of data file(istyr)
     >        istat,  ! error flag for netcdf
     >        i,      ! loop indice on land points
     >        jyear,   ! iyr divided by nanom (for looping through anomaly files)
     >        iiyear
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

      if(iyear.ge.2015) then
       istart(3)=iyear-1860
      else
       istart(3)=iyear-1860+1
      end if

      aname = 'Nmanure'
      filen = '/storage03/zhangkr/SongHX/input0.5/Nmanure.nc'
      call readvard(filen,aname,'',istart,icount,
     > work,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),
     > istat)
      if (istat .ne. 0) then
         write(*,*) 'ERROR in rdyear, nmanure.nc'
         stop 1
      end if
      call arr2vec( work, manure_N )
c      write(*,*) "manure_N=",manure_N

      return
      end
c
c