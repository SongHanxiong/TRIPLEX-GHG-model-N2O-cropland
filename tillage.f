c
c  -------------------------------------------------------------
      subroutine tillage (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to simulate the process of tillage process by 
c 1) decide if there is a tillage event by indicator which is provided by site-info
c or set timing and 2) confirm such information and 3) average the major soil variables of the 
c top 4 layers of soil profile and send those updated value in related submodule for calculation.
c NOTE: This subroutine is called in hourly loop of the main program.

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
     >        tillevent,       ! tillage event sender (we mark it here and the definiation is in comsoi.h) ###songhx###
     >        tillyr,
     >        tillmon,
     >        tillday,
c
     >        nday    !nday for debug
     
c     initial variables in initial.f at the begining of the program
c
c     initial some relative soil parameters
c      call initial_1stYear(i,k,iyear)    !put in the initial file

c
c     day_wfps: in DNDC daily soil moisture water filled prorosity (0~1)
c     not the same as the varaible "wsoi" in IBIS,
c     wsoi is relative to pore space not occupied by ice and water
c     thus must include the ice fraction in the calculation
c     wfps = (1.0 - wisoi(i,1)) * wsoi(i,1)
c     the unit of tsoi(i,k) in IBIS is K, to centigrade degree, minus 273.13

      do i= 1, npoi !  global iteration

      call Zero_day()

      do k = 1, nsoilay
          day_wfps(i,k) = (1.0 - wisoi(i,k)) * wsoi(i,k) ! 0~1
          SoilT(i,k) = tsoi(i,k) - 273.13
          Soil_Ice(i,k) = wisoi(i,k)
          SoilPoros(i,k) = poros(i,k)
          SoilDensi(i,k) = (rhosoi(i,k) * (1.0-poros(i,k)) * 0.001)
          SoilWater(i,k) = wsoi(i,k)
          SoilClay(i,k) = clay(i,k)/100.0   ! in DNDC, clay is fraction (0~1)
      end do
c
c initialize tillage event   
          tillevent(i) = 0.0          
c decide when the tillage happens
c the tillage dates are generated from published papers.
          tillyr = 'type in year here'
          tillmon = 'type in month here'
          tillday = 'type in day here'
             if (iyear .eq. tillyr .and. imonth .eq. tillmon .and. iday .eq. tillday) then
                tillevent(i) = 1.0
             else 
                tillevent(i) = 0.0
             end if
          end if
c average the values of major soil variables of the top 4 layers of soil profile
c NOTE: soil physical properties (i.e., sand, clay, slit fraction) will not change.
c first, initialize the cumulative values of soil variables
          cum_NO3(i) = 0.0
          cum_NO2(i) = 0.0
          cum_NH4(i) = 0.0
          cum_DOC(i) = 0.0
          cum_wsoi(i) = 0.0
          cum_wisoi(i) = 0.0

c average them when there is a tillevent         
          do k = 1, 4
             cum_NO3(i) = cum_NO3(i) + NO3(i,k)
             cum_NO2(i) = cum_NO2(i) + NO2(i,k)
             cum_NH4(i) = cum_NH4(i) + NH4(i,k)
             cum_DOC(i) = cum_DOC(i) + DOC(i,k)
             cum_wsoi(i) = cum_wsoi(i) + wsoi(i,k)
             cum_wisoi(i) = cum_wisoi(i) + wisoi(i,k)
             
          end do
c        
          if (tillevent(i) = 1.0) then
             do k = 1,4
                NO3(i,k) = cum_NO3(i) / 4.0
                NO2(i,k) = cum_NO2(i) / 4.0
                NH4(i,k) = cum_NH4(i) / 4.0
                DOC(i,k) = cum_DOC(i) / 4.0
                wsoi(i,k) = cum_wsoi(i) / 4.0
                wisoi(i,k) = cum_wisoi(i) / 4.0
                
             end do   
           else
             do k = 1,4
                NO3(i,k) = max(NO3(i,k), 0.0)
                NO2(i,k) = max(NO2(i,k), 0.0)
                NH4(i,k) = max(NH4(i,k), 0.0) 
                DOC(i,k) = max(DOC(i,k), 0.0)
                wsoi(i,k) = max(wsoi(i,k), 0.0)
                wisoi(i,k) = max(wisoi(i,k), 0.0)
             end do
           end if
c           
c more exposure to air with same amount of O2 concentration for top soil layer and atmos
           O2_0layer = 0.209
c           
           if (tillevent(i) = 1.0) then
             do k = 1,4    
                if (k .eq. 1) then
                  O2(i,k) = O2_0layer
                else
                  O2(i,k) = max(0.5*(O2(i,k)+O2_0layer), 0.0)
                end if
              end do
            end if       

c return to main program
c
      return
      end
c             