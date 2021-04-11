c ------
c comagrp
c ------
c
      real

     >  xinrfCropFra(npoi),  ! C3 crops fraction for each pixel
     >  xinC4CropFra(npoi),  ! C4 crops fraction for each pixel
     >  parnno(npoi),              !nitrification rate coefficient
     >  CNratio(npoi),             !soil C/N
     >  F_crops(npoi),             !###songhx### tested upland cropland fraction
     >  F_rf_norice(npoi),         !rainfed cropland fraction
     >  F_ir_norice(npoi),         !irrigated cropland fractioon
     >  COEdNO3(npoi),              ! parameter for coefficient of NO3- consumption
     >  experimentF(npoi),         !for sites validation only according to papers and should be compared with Ito(2016)
                                   !index for if there is a experiment fertilization senario becasue some of them are extremely higher than Ito(2016)
     >  experimentM(npoi)
c
      common /comagrp1/ xinrfCropFra,xinC4CropFra, parnno,CNratio,F_crops,F_rf_norice,F_ir_norice, 
     > COEdNO3,experimentF, experimentM                         

c
      integer
     >  lmask(nlon,nlat)  ! landmask 0=water, 1=land
c
      common /comagrp2/ lmask
c

c
      real
     >  manure_N(npoi),    ! manure_N content including organic and inorganic !###songhx###
     >  manure_C(npoi),    ! manure_C content calculated by manure_N and related C/N ratios !###songhx###
     >  NH4_manure(npoi),  ! NH4 released directly from manure_N !###songhx###
     >  NO3_manure(npoi),  ! NO3 released directly from manure_N !###songhx### the portion are global mean for now.
     >  NH4fer(npoi,12),     ! fertilizer NH4-N kg-N/ha per month but with fertilizer date it will become daily! ###songhx###2019/2/27
     >  NO3fer(npoi,12),     ! fertilizer NO3-N kg-N/ha per mouth
     >  ferdate(npoi,12),       !###songhx### 2019/10/29 fertilizer application date
     >  ferday(npoi),           !###songhx### 2019/10/29 becasue the ferdate ranges from (1 to 365) so that we have to recalculate into 1-31.
     >  yndep(npoi),
     >  applyfer(npoi),     !###songhx###
     >  applyman(npoi),     !###songhx###
     >  irriamount(npoi),   !###songhx### how much to irrigate (kg H2O m-2)
     >  cum_NO3(i),         !###songhx### cumulative NO3 of the first 4 soil layers for tillage effect
     >  cum_NO2(i),
     >  cum_NH4(i),
     >  cum_DOC(i),
     >  cum_wsoi(i),
     >  cum_wisoi(i)
c
      common /comagrp3/  manure_N,NH4_manure,NH4fer,NO3fer,ferdate,ferday,yndep,applyfer,  !###songhx###
     >                   irriamount, cum_NO3, cum_NO2, cum_NH4, cum_DOC, cum_wsoi, cum_wisoi
c
      real

     >  NH4_add(npoi),         !all NH4 input from fertilizer, manure(directly), deposition  !###songhx###
     >  NO3_add(npoi)
c

* DTP 2002/04/17: Declare scaling factors for compressed anomaly data
*                 Used in rdanom(), readit() and init_anoms()

      common /comagrp4/ NH4_add, NO3_add
      
      real

     >day_wfps(npoi,nsoilay),
c
     >SoilDensi(npoi,nsoilay),
     >SoilT(npoi,nsoilay),
     >Soil_Ice(npoi,nsoilay),
     >SoilDepth(nsoilay),
     >SoilPoros(npoi,nsoilay),
     >SoilWater(npoi,nsoilay),
     >SoilClay(npoi,nsoilay),
c
     > NH4(npoi,nsoilay),
     > NH3(npoi,nsoilay),
     > NO3(npoi,nsoilay),
     > an_NO3(npoi,nsoilay),
     > N2O(npoi,nsoilay),
     > an_N2O(npoi,nsoilay),
     > NO(npoi,nsoilay),
     > an_NO(npoi,nsoilay),
     > NO2(npoi,nsoilay),
     > N2(npoi,nsoilay),
     > DOC(npoi,nsoilay),
     > an_DOC(npoi,nsoilay),
     > O2(npoi,nsoilay),
     > SOC(npoi,nsoilay)


c
      common /comargp5/ day_wfps, SoilDensi, SoilT, Soil_Ice, SoilDepth,
     > SoilPoros, SoilWater, SoilClay, NH4, NH3, NO3, an_NO3, N2O,
     > an_N2O, NO, an_NO, NO2, N2, DOC, an_DOC, O2, SOC
c
      real
     >  NPlantUp(npoi),         !bio_1117
     >  N_Imm(npoi),            !bio_1025
     >  GHG_fleach(npoi),        !bio_1902
     >  GHG_runoff(npoi),         !#songhx 2019/4/16
     >  TotalAN(npoi,nsoilay)     !#songhx 2019/4/18
c     >  WT_Factors        !physiology_317   !effects of Water table on photosynthesis
      common /comagrp6/ NPlantUp,N_Imm,GHG_fleach !,WT_Factors
      
c
      real 
     >  wpud(npoi),              ! liquid content of puddles per soil area (kg m-2)
     >  wipud(npoi),             ! ice content of puddles per soil area (kg m-2)
     >  manureC2pas(npoi),       !manure_C to passive SOM
     >  manureC2slop(npoi),      !manure_C to protected SOM
     >  manureC2slon(npoi),       !manure_C to non-protected SOM
     >  csoipas(npoi),      ! carbon in soil - passive humus               (kg_C m-2)
     >  csoislo(npoi),      ! carbon in soil - slow humus                  (kg_C m-2)
     >  csoislon(npoi),     ! carbon in soil - slow nonprotected humus     (kg_C m-2)
     >  csoislop(npoi)     ! carbon in soil - slow protected humus        (kg_C m-2)
c
      common /comagrp7/ wpud, wipud, manureC2pas, manureC2slop, manureC2slon, csoipas, csoislo,
     >                 csoislon, csoislop
c
      real
     >  hsoi(nsoilay+1)          ! soil layer thickness (m)
c
      common /comagrp8/ hsoi
c
      real 
     >  tsoi(npoi,nsoilay),      ! soil temperature for each layer (K)
     >  wsoi(npoi,nsoilay),      ! fraction of soil pore space containing liquid water
     >  wisoi(npoi,nsoilay),     ! fraction of soil pore space containing ice
     >  rhosoi(npoi,nsoilay),    ! soil density (without pores, not bulk) (kg m-3)
     >  poros(npoi,nsoilay),     ! porosity (mass of h2o per unit vol at sat / rhow)
     >  porosflo(npoi,nsoilay),  ! porosity after reduction by ice content
     >  sand(npoi,nsoilay),      ! percent sand of soil
     >  clay(npoi,nsoilay),      ! percent clay of soil
     >  irrievent(npoi)
c
      common /comagrp9/ tsoi, wsoi, wisoi, rhosoi, poros, porosflo, 
     >     sand, clay, irrievent
c
