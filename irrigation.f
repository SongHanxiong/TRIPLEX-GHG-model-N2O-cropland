c
c  -------------------------------------------------------------
      subroutine irrigation (iyear,imonth,iday,nday,iyear0)
c  -------------------------------------------------------------
c This routine is constructed to simulate the process of irrigation process by 
c 1) decide if there is an irrigation event by indicator which is provided by site-info
c or set timing and 2) confirm such information and 3) add extra water to the variable 
c 'wpud' which is going back to soilcli module in soil.f file (subroutine lsxmain).
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
     >        irrievent,       ! irrigation event sender (we mark it here and the definiation is in comsoi.h) ###songhx###
     >        irriyr,         
     >        irrimon,
     >        irriday,
c
     >        nday    !nday for debug
c initialize irrigation event   
          irrievent(i) = 0.0          
c decide when the irrigation happens
c the irrigation dates are generated from published papers.
          irriyr = 'type in year here'
          irrimon = 'type in month here'
          irriday = 'type in day here'
             if (iyear .eq. irriyr .and. imonth .eq. irrimon .and. iday .eq. irriday) then
                irrievent(i) = 1.0
             else 
                irrievent(i) = 0.0
             end if
          end if
c decide how much to irrigate (kg H2O m-2)
c initialize irrigation amount
          irriamount(i) = 0.0
c add more water into soil surface puddle liquid         
          if (irrievent(i) .eq. 1.0) then
             irriamount(i) = 'type in irrigation amount here'
             wpud(i) = wpud(i) + irriamount(i)
          else
             irriamount(i) = 0.0
             wpud(i) = wpud(i)
          end if
c return to main program
c
      return
      end
             