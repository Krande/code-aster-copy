! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine dtmallo(sd_dtm_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmallo : Memory allocation for a dynamic simulation on a projected basis.
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvr8.h"
#include "asterfort/mdallo.h"
#include "asterfort/mdlibe.h"
#include "asterfort/nlget.h"

!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
!
!   -0.2- Local variables
    integer(kind=8)           :: nbsauv, nbmode, iret, nbnli, nbvint, nbsteps
    integer(kind=8)           :: jordr, jdisc, jptem, jdepl
    integer(kind=8)           :: jvite, jacce, jvint
    integer(kind=8)           :: adapt, iarch_sd, iret1, iret2, nltreat
    real(kind=8)      :: dt, dtmin, dtmax, deltadt, epsi, taille, taille2
    character(len=8)  :: sd_dtm, nomres, basemo, riggen, masgen
    character(len=8)  :: amogen, sd_nl
    character(len=16) :: schema
!
    character(len=8), pointer :: inticho(:) => null()
    character(len=8), pointer :: noeucho(:) => null()
    character(len=8), pointer :: fonred(:) => null()
    character(len=8), pointer :: fonrev(:) => null()
    character(len=8), target  :: blanc(1)
    integer(kind=8), pointer :: vindx(:) => null()
    integer(kind=8), pointer :: allocs(:) => null()
    integer(kind=8) :: ordr
    real(kind=8) :: disc, ptem
    real(kind=8), pointer :: v_depl(:) => null()
    real(kind=8), pointer :: v_vite(:) => null()
    real(kind=8), pointer :: v_acce(:) => null()
    real(kind=8), pointer :: v_vint(:) => null()

!
!   0 - Initializations
    sd_dtm = sd_dtm_

    epsi = 100.d0*r8prem()
    blanc(1) = ' '
    inticho => blanc
    noeucho => blanc
    fonred => blanc
    fonrev => blanc

!
!   1 - Retrieval of the necessary information
    call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
    call dtmget(sd_dtm, _ARCH_NB, iscal=nbsauv)
    call dtmget(sd_dtm, _SCHEMA, kscal=schema)
    call dtmget(sd_dtm, _BASE_MOD, kscal=basemo)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    call dtmget(sd_dtm, _RIGI_MAT, kscal=riggen)
    call dtmget(sd_dtm, _MASS_MAT, kscal=masgen)
    call dtmget(sd_dtm, _NL_TREAT, iscal=nltreat)

    amogen = ' '
    call dtmget(sd_dtm, _AMOR_MAT, lonvec=iret)
    if (iret .gt. 0) call dtmget(sd_dtm, _AMOR_MAT, kscal=amogen)

    call dtmget(sd_dtm, _DT, rscal=dt)
!
    sd_nl = ' '
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)
    if (nbnli .gt. 0) call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)
!

    call dtmget(sd_dtm, _IND_ALOC, lonvec=iret1, vi=allocs)
    if (iret1 .ne. 0) then
        call dtmget(sd_dtm, _ADAPT, iscal=adapt)
        call dtmget(sd_dtm, _IARCH_SD, iscal=iarch_sd)
    else
        adapt = 0
        iarch_sd = 0

        if ((schema(1:5) .eq. 'RUNGE') .or. (schema(1:5) .eq. 'ADAPT') .or. &
            (schema(1:6) .eq. 'DEVOGE')) then
            call getvr8('INCREMENT', 'PAS     ', iocc=1, scal=dt)
            call getvr8('SCHEMA_TEMPS', 'PAS_MINI', iocc=1, scal=dtmin, nbret=iret1)
            call getvr8('SCHEMA_TEMPS', 'PAS_MAXI', iocc=1, scal=dtmax, nbret=iret2)
            if (iret1 .ne. 1) dtmin = 1.d-6*dt
            if (iret2 .ne. 1) dtmax = 1.d6*dt
            deltadt = abs(dtmax/dtmin-1)
            if (deltadt .gt. epsi) then
                adapt = 1
                iarch_sd = 1
            end if
        end if

        if (nltreat .eq. 1) iarch_sd = 1

        call dtmsav(sd_dtm, _ADAPT, 1, iscal=adapt)
        call dtmsav(sd_dtm, _IARCH_SD, 1, iscal=iarch_sd)
    end if

    nbvint = 0
    if (nbnli .gt. 0) then
        call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx)
        nbvint = vindx(nbnli+1)-1
    end if

!   2 - Calculate archiving memory usage, based on nbsauv, nbmode, nbnli
!       > real(kind=8) => 8 bytes (64 bits) per value
!       > integer      => 4 bytes (32 bits) per value
!   ---------------------------------------------------------------------------
!   + ORDR object (integer)
    taille = nbsauv*4.d0

!   + DISC/PTEM objects (real)
    taille = taille+2*(nbsauv*8.d0)

!   + DEPL/VITE/ACCE objects (real)
    taille = taille+3*(nbsauv*nbmode*8.d0)

    write (*, *) '--------------------------------------------------------'
    write (*, *) 'D Y N A    V I B R A    S D    S I Z E   I N F O'
    write (*, *) '--------------------------------------------------------'
    if (nbnli .gt. 0) then
        write (*, *) '> RESULTS (EXCLUDING NLS):', taille/(1024*1024), 'MB'
!       + .NL.VINT object (real)
        taille = taille+nbsauv*nbvint*8.d0

!       + .NL.VIND object (integer)
        taille2 = taille+(nbnli+1)*4.d0

!       + .NL.TYPE object (integer)
        taille2 = taille2+(nbnli)*4.d0

!       + .NL.INTI object (K24)
        taille2 = taille2+(5*nbnli)*24.d0
        write (*, *) '> NL INTERNAL VARS       :', &
            nbsauv*nbvint*8.d0/(1024*1024), 'MB'
        write (*, *) '> NL FIXED SIZE          :', &
            (taille2-taille)/(1024*1024), 'MB'
        write (*, *) '> NL TOTAL               :', &
            (nbsauv*nbvint*8.d0+taille2-taille)/(1024*1024), 'MB'
        write (*, *) '--------------------------------------------------------'
        write (*, *) '> TOTAL                  :', &
            taille2/(1024*1024), 'MB'

    else
        write (*, *) '> TOTAL                  :', taille/(1024*1024), 'MB'
    end if
    write (*, *) '> MEMORY SIZE PER STEP   :', taille/(1024*nbsauv), 'KB'

    ! implementation en fonction de nbsauv
    if (nbsauv .le. 2) then
        ! 1 bloc
        nbsteps = nbsauv
    else if (nbsauv .lt. 400) then
        ! 2 blocs si < 2**2*100
        nbsteps = (nbsauv-1)/2+1
    else
        ! sinon nbsteps ~= 100*n_blocs
        nbsteps = int(10*sqrt(real(nbsauv)))+1
    end if
    ! on limite à 10000 blocs maximum
    if (nbsauv .ge. nbsteps*10000) then
        nbsteps = nbsauv/10000
    end if
    ! au moins 2 increments par bloc
    if (nbsteps .lt. 2) then
        nbsteps = 2
    end if
    write (*, *) '> NB STEPS PER BLOC      :', nbsteps, 'STEPS'
    write (*, *) '> MEMORY SIZE PER BLOC   :', taille/(1024*1024)*nbsteps/nbsauv, 'MB'

    if (nbsauv .gt. nbsteps) then
        if (iarch_sd .eq. 0) then
            iarch_sd = 1
            call dtmsav(sd_dtm, _ARCH_NB, 1, iscal=nbsauv)
            call dtmsav(sd_dtm, _ADAPT, 1, iscal=adapt)
            call dtmsav(sd_dtm, _IARCH_SD, 1, iscal=iarch_sd)
        end if
        nbsauv = nbsteps
    end if

    if (iarch_sd .gt. 1) then
        ! copy last value from previous bloc ...
        AS_ALLOCATE(vr=v_depl, size=nbmode)
        AS_ALLOCATE(vr=v_vite, size=nbmode)
        AS_ALLOCATE(vr=v_acce, size=nbmode)
        if (nbvint .gt. 0) then
            AS_ALLOCATE(vr=v_vint, size=nbvint)
        end if
        ordr = zi(allocs(1)-1+nbsauv)
        disc = zr(allocs(2)-1+nbsauv)
        ptem = zr(allocs(3)-1+nbsauv)
        v_depl = zr(allocs(4)+nbmode*(nbsauv-1):allocs(4)-1+nbmode*nbsauv)
        v_vite = zr(allocs(5)+nbmode*(nbsauv-1):allocs(5)-1+nbmode*nbsauv)
        v_acce = zr(allocs(6)+nbmode*(nbsauv-1):allocs(6)-1+nbmode*nbsauv)
        if (nbvint .gt. 0) then
            v_vint = zr(allocs(7)+nbvint*(nbsauv-1):allocs(7)-1+nbvint*nbsauv)
        end if
        call dtmsav(sd_dtm, _ARCH_STO, 4, ivect=[1, 0, 0, 0])
        call mdlibe(nomres(1:8), nbnli, iarch_sd-1)
    end if
    call mdallo(nomres(1:8), 'TRAN', nbsauv, sauve='GLOB', method=schema, &
                base=basemo, nbmodes=nbmode, rigi=riggen, mass=masgen, amor=amogen, &
                dt=dt, nbnli=nbnli, checkarg=.false._1, &
                jordr=jordr, jdisc=jdisc, jptem=jptem, jdepl=jdepl, jvite=jvite, &
                jacce=jacce, jvint=jvint, sd_nl_=sd_nl, sd_index=iarch_sd)

    if (iarch_sd .gt. 1) then
        ! ... paste to next bloc first value
        zi(jordr) = ordr
        zr(jdisc) = disc
        zr(jptem) = ptem
        zr(jdepl:jdepl+nbmode) = v_depl
        zr(jvite:jvite+nbmode) = v_vite
        zr(jacce:jacce+nbmode) = v_acce
        if (nbvint .gt. 0) then
            zr(jvint:jvint+nbvint) = v_vint
        end if
        AS_DEALLOCATE(vr=v_depl)
        AS_DEALLOCATE(vr=v_vite)
        AS_DEALLOCATE(vr=v_acce)
        if (nbvint .gt. 0) then
            AS_DEALLOCATE(vr=v_vint)
        end if
    end if
    call dtmsav(sd_dtm, _IND_ALOC, 7, ivect=[jordr, jdisc, jptem, jdepl, jvite, jacce, jvint])

end subroutine
