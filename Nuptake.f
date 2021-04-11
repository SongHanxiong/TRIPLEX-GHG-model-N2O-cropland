c
c  -------------------------------------------------------------
      subroutine Nplantuptake (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to improve the simulation of the process of plant
c uptake soil inorganic N.  The major framworks are adapted from original TRIPLEX-GHG and DNDC model.

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

      real NPlantUpLayer_NO3(nsoilay),
     >     NPlantUpLayer_NH4(nsoilay),
     >     NPlantUpLayer(nsoilay),
     >     NPlantUpin2layer(npoi)
c
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

c      SoilPoros(i,1) = 0.9
c
c     update after plant uptake
     
      do k = 1, nsoilay 
         TotalAN(i,k) = NH4(i,k)+NO3(i,k)
      end do
      do k = 1,nsoilay           
          NPlantUpLayer_NO3(k) = NO3pref*NPlantUp(i)* NO3(i,k) /((nsoilay)*TotalAN(i,k))* 10000.0 ! change the unit from kg N /m2 to kg N /ha
          NPlantUpLayer_NH4(k) = NPlantUp(i) - NPlantUpLayer_NO3(k)
      end do
      do k = 1, nsoilay
          if (k .LT. nsoilay) then
            NO3(i,k) = NO3(i,k) - NPlantUpLayer_NO3(k)! * (NH4(i,k)/(NH4(i,k) + NO3(i,k)))
            if (NO3(i,k) .LT. 0.0) then
                NH4(i,k) = NH4(i,k) - abs(NO3(i,k))
                NO3(i,k)=0.0
                if (NH4(i,k) .LT. 0.0) then
                  NPlantUpLayer_NO3(k+1)= NPlantUpLayer_NO3(k+1) + abs(NH4(i,k))
                  NH4(i,k) = 0.0
                end if
            endif
          else
            NO3(i,k) = NO3(i,k) - NPlantUpLayer_NO3(k)! * (NH4(i,k)/(NH4(i,k) + NO3(i,k)))
            if (NO3(i,k) .LT. 0.0) then
               NH4(i,k) = NH4(i,k) - abs(NO3(i,k))
               NO3(i,k)=0.0
               if (NH4(i,k) .LT. 0.0) then
c                 NPlantUpLayer_NO3(k+1)= NPlantUpLayer_NO3(k+1) + abs(NH4(i,k))
                 NH4(i,k)=0.0
               endif
            endif
          endif
          
          if (k .LT. nsoilay) then
            NH4(i,k) = NH4(i,k) - NPlantUpLayer_NH4(k)! * (NH4(i,k)/(NH4(i,k) + NO3(i,k)))
            if (NH4(i,k) .LT. 0.0) then
               NO3(i,k) = NO3(i,k) - abs(NH4(i,k))
               NH4(i,k)=0.0
               if(NO3(i,k) .LT. 0.0) then
                  NPlantUpLayer_NH4(k+1)= NPlantUpLayer_NH4(k+1) + abs(NH4(i,k))
                  NO3(i,k) = 0.0
               endif
            endif
          else
            NH4(i,k) = NH4(i,k) - NPlantUpLayer_NH4(k)! * (NH4(i,k)/(NH4(i,k) + NO3(i,k)))
            if (NH4(i,k) .LT. 0.0) then
               NO3(i,k) = NO3(i,k) - abs(NH4(i,k))
               NH4(i,k)=0.0
               if (NO3(i,k) .LT. 0.0) then
c                 NPlantUpLayer_NH4(k+1)= NPlantUpLayer_NH4(k+1) + abs(NO3(i,k))
                 NO3(i,k)=0.0
               endif
            endif
          endif
c        if ((NH4(i,k)+ NO3(i,k)) .GT. 0.0) then
cc          write (*,*) NPlantUp, "NH4", NH4(i,k)/10000 ,"NO3", NO3(i,k)/10000
c      	  NH4(i,k) = NH4(i,k) - NPlantUp/(nsoilay +  OrgLayers)* 10000.0 *(NH4(i,k)/(NH4(i,k) + NO3(i,k))) ! change the unit from kg N /m2 to kg N /ha
c      	  NO3(i,k) = NO3(i,k) - NPlantUp/(nsoilay +  OrgLayers)* 10000.0 *(NO3(i,k)/(NH4(i,k) + NO3(i,k)))
c      	end if
      end do
c
      end 
c      