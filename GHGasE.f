c
c  -------------------------------------------------------------
      subroutine GreenHouseGassEmission (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to simulate the process of nitification
c and denitrification. The gases like N2O, NO, NH3 and N2 are generated
c with these two processes. The major framworks are adapted from DNDC model.

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
c###songhx###
      call Decomposition_up(i,iday,imonth)

      call fertilization (iyear,imonth,iday,nday,iyear0)
      call manureapp (iyear,imonth,iday,nday,iyear0)
      
         NH4(i,1)=NH4(i,1)+NH4_add(i)
         NO3(i,1)=NO3(i,1)+NO3_add(i)
c               
      call WaterTable(i, iday, imonth)
c
      call Nplantuptake (iyear,imonth,iday,nday,iyear0)
      
      call Denitrification_up(i,iday)
      
c     update after leaching
c     In biogeochem.f : nleach(i)
      call DeLeaching(i)
c
c###songhx###
c
c     update nm(i) after Gas emission

      totnsoi(i) = totnsoi(i) + nm(i)/1000 ! change the unit from kg N /ha to kg N /m^2
      totneco(i) = totneco(i) + nm(i)/1000

c     updata MicCO2 emission
      tco2mic(i) = DailyMicCO2(i)/(86400. * 0.012)
c
c     update Total soil carbon
      totcsoi(i) = TotCSoil(i) + totcmic(i)      !!!!!!

      call DailyEmissionSum(i)

      end do !do i= 1, npoi !  global iteration
c
c
      end

c