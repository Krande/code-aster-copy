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

subroutine dtminfo(sd_dtm_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtminfo : Prints out general information about the transient calculation to
!           be undertaken, before any integration taking place, but after all
!           preliminary verifications are carried out.
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8depi.h"
#include "asterc/r8prem.h"
#include "asterfort/dtmget.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/nlget.h"
#include "asterfort/utmess.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
!
!   -0.2- Local variables
    integer(kind=8)               :: nbnli, nbchoc, nbants, nbdiscret, nbflam
    integer(kind=8)               :: nbrede, nbrevi, nbpas, nbpas_min, nbpas_max
    integer(kind=8)               :: substruct, fsicase, nbpheq, adapt, nltreat
    integer(kind=8)               :: nbmode, append, nbsav_forc, iparch
    integer(kind=8)               :: i, iret, iret1, iret2, nbfreq
    integer(kind=8)               :: nbrfis, nbpal, etinit
    real(kind=8)          :: dt, dtmin, dtmax, deuxpi, epsi
    real(kind=8)          :: tinit, tfin, cdivi, toler, fsivgap
    real(kind=8)          :: fmin, fmax, freq, nbpas_max_r
    character(len=8)      :: sd_dtm, basemo, modgen, fsibase
    character(len=8)      :: dep0, vit0, contres, matmas, matrig
    character(len=8)      :: matamo, sd_nl
    character(len=16)     :: numddl, schema, schtyp, nltreat_k
    character(len=18)     :: nomres
!
    real(kind=8), pointer :: amogen(:) => null()
    real(kind=8), pointer :: puls(:) => null()
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    deuxpi = r8depi()
    epsi = 1000.d0*r8prem()

    call utmess('I', 'DYNAMIQUE_55', nk=3, valk=['D Y N A _ V I B R A', &
                                                 '    TRANsitoire    ', &
                                                 '    GENEralisee    '])
!
!   1 - Calculation type : standard, substructuring, fluid-structure
    call dtmget(sd_dtm, _SUB_STRU, iscal=substruct)
    call dtmget(sd_dtm, _FSI_CASE, iscal=fsicase)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    if ((substruct+fsicase) .eq. 0) then
        call dtmget(sd_dtm, _BASE_MOD, kscal=basemo)
        call dtmget(sd_dtm, _NB_PHYEQ, iscal=nbpheq)
        call utmess('I', 'DYNAMIQUE_56', nk=1, valk=[basemo], ni=1, vali=[nbpheq])
    end if

    if (substruct .eq. 1) then
        call dtmget(sd_dtm, _BASE_MOD, kscal=modgen)
        call dtmget(sd_dtm, _NUM_DDL, kscal=numddl)
        call utmess('I', 'DYNAMIQUE_57', nk=2, valk=[modgen, numddl])
    end if

    if (fsicase .eq. 1) then
        call dtmget(sd_dtm, _FSI_BASE, kscal=fsibase)
        call dtmget(sd_dtm, _FSI_VGAP, rscal=fsivgap)
        call utmess('I', 'DYNAMIQUE_58', nk=1, valk=[fsibase], nr=1, valr=[fsivgap])
    end if

!   2 - Minimum and max frequencies
    call dtmget(sd_dtm, _OMEGA, vr=puls, lonvec=nbfreq)
    fmin = 0.d0
    fmax = 0.d0
    do i = 1, nbfreq
        freq = puls(i)/(deuxpi)
        if (freq .lt. fmin) fmin = freq
        if (freq .gt. fmax) fmax = freq
    end do
    call utmess('I', 'DYNAMIQUE_59', ni=1, vali=[nbmode], nr=2, valr=[fmin, fmax])

!   3 - Dynamic matrices
    call utmess('I', 'DYNAMIQUE_60')

    if (fsicase .eq. 1) then
        call utmess('I', 'DYNAMIQUE_65')
    else
        call dtmget(sd_dtm, _MASS_MAT, kscal=matmas)
        call dtmget(sd_dtm, _RIGI_MAT, kscal=matrig)
        call utmess('I', 'DYNAMIQUE_61', nk=2, valk=[matmas, matrig])

        call dtmget(sd_dtm, _AMOR_MAT, lonvec=iret)
        if (iret .gt. 0) then
            call dtmget(sd_dtm, _AMOR_MAT, kscal=matamo)
            call utmess('I', 'DYNAMIQUE_62', nk=1, valk=[matamo])
        else
            call dtmget(sd_dtm, _AMOR_DIA, vr=amogen)
            do i = 1, nbmode
                if (abs(amogen(i)) .gt. epsi) then
                    call utmess('I', 'DYNAMIQUE_63')
                    goto 10
                end if
            end do
            call utmess('I', 'DYNAMIQUE_64')
10          continue
        end if
    end if

!   4 - Integration scheme, type, timestep, principal numerical parameters
    call dtmget(sd_dtm, _SCHEMA, kscal=schema)
    call dtmget(sd_dtm, _ADAPT, iscal=adapt)
    call dtmget(sd_dtm, _DT, rscal=dt)
    call dtmget(sd_dtm, _INST_INI, rscal=tinit)
    call dtmget(sd_dtm, _INST_FIN, rscal=tfin)

    schtyp = 'explicite'
    if (schema(1:7) .eq. 'NEWMARK' .or. schema(1:6) .eq. 'TRBDF2') then
        schtyp = 'implicite'
    end if
    if (adapt .eq. 1) then
        call dtmget(sd_dtm, _DT_MIN, rscal=dtmin)
        call dtmget(sd_dtm, _DT_MAX, rscal=dtmax)
        call utmess('I', 'DYNAMIQUE_66', nk=1, valk=[schema], &
                    nr=3, valr=[dt, dtmin, dtmax])
        if (schema(1:5) .eq. 'ADAPT') then
            call getvr8('SCHEMA_TEMPS', 'COEF_DIVI_PAS', iocc=1, scal=cdivi)
            call utmess('I', 'DYNAMIQUE_68', nr=1, valr=[cdivi])
        else
            call getvr8('SCHEMA_TEMPS', 'TOLERANCE', iocc=1, scal=toler)
            call utmess('I', 'DYNAMIQUE_67', nr=1, valr=[toler])
        end if
        nbpas_min = nint((tfin-tinit)/dtmax)
        nbpas_max = 1000000000
        if (dtmin .gt. epsi) then
            nbpas_max_r = min(1.d0*nbpas_max, (tfin-tinit-epsi)/dtmin)
            nbpas_max = int(nbpas_max_r)+1
        end if

        call utmess('I', 'DYNAMIQUE_69', ni=2, vali=[nbpas_min, nbpas_max])
    else
        nbpas = nint((tfin-tinit)/dt)
        call utmess('I', 'DYNAMIQUE_70', nk=2, valk=[schema, schtyp], &
                    nr=1, valr=[dt], &
                    ni=1, vali=[nbpas])
    end if

!   5 - Localized non linearities
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)
    if (nbnli .gt. 0) then
        call utmess('I', 'DYNAMIQUE_71')
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)
        call nlget(sd_nl, _NB_CHOC, iscal=nbchoc)
        call nlget(sd_nl, _NB_REL_FX, iscal=nbrede)
        call nlget(sd_nl, _NB_REL_FV, iscal=nbrevi)

        if (nbchoc .gt. 0) then
            call dtmget(sd_dtm, _NL_TREAT, iscal=nltreat)
            nltreat_k = 'explicite'
            if (nltreat .eq. 1) then
                nltreat_k = 'implicite'
            end if
            call utmess('I', 'DYNAMIQUE_72', nk=1, valk=[nltreat_k], ni=1, vali=[nbchoc])
        end if

        call nlget(sd_nl, _NB_ANTSI, iscal=nbants)
        if (nbants .gt. 0) call utmess('I', 'DYNAMIQUE_73', ni=1, vali=[nbants])

        call nlget(sd_nl, _NB_FLAMB, iscal=nbflam)
        if (nbflam .gt. 0) call utmess('I', 'DYNAMIQUE_74', ni=1, vali=[nbflam])

        if (nbrede .gt. 0) call utmess('I', 'DYNAMIQUE_75', ni=1, vali=[nbrede])

        if (nbrevi .gt. 0) call utmess('I', 'DYNAMIQUE_76', ni=1, vali=[nbrevi])

        call nlget(sd_nl, _NB_DIS_VISC, iscal=nbdiscret)
        if (nbdiscret .gt. 0) then
            call utmess('I', 'DYNAMIQUE_77', ni=1, vali=[nbdiscret], nk=1, valk=['DIS_VISC'])
        end if

        call nlget(sd_nl, _NB_DIS_ECRO_TRAC, iscal=nbdiscret)
        if (nbdiscret .gt. 0) then
            call utmess('I', 'DYNAMIQUE_77', ni=1, vali=[nbdiscret], nk=1, valk=['DIS_ECRO_TRAC'])
        end if

        call nlget(sd_nl, _NB_R_FIS, iscal=nbrfis)
        if (nbrfis .gt. 0) call utmess('I', 'DYNAMIQUE_98', ni=1, vali=[nbrfis])

        call nlget(sd_nl, _NB_PALIE, iscal=nbpal)
        if (nbpal .gt. 0) call utmess('I', 'DYNAMIQUE_99', ni=1, vali=[nbpal])

    end if

!   6 - Initial state, continue case ?
    call dtmget(sd_dtm, _APPND_SD, iscal=append)
    if (append .ne. 0) then
        call dtmget(sd_dtm, _RESU_SD, kscal=nomres)
        nomres(9:18) = ' (reprise)'
        call utmess('I', 'DYNAMIQUE_78', nk=1, valk=[nomres])
    else
        call getfac('ETAT_INIT', etinit)
        if (etinit .eq. 0) then
            dep0 = 'nul'
            vit0 = 'nul'
            call utmess('I', 'DYNAMIQUE_79', nk=2, valk=[dep0, vit0])
            call utmess('I', 'DYNAMIQUE_84', nk=1, valk=['calculee'])
        else
            call getvid('ETAT_INIT', 'RESULTAT', iocc=1, scal=contres, nbret=iret)
            if (iret .ne. 0) then
                call utmess('I', 'DYNAMIQUE_78', nk=1, valk=[contres])
            else
                call getvid('ETAT_INIT', 'DEPL', iocc=1, scal=dep0, nbret=iret1)
                call getvid('ETAT_INIT', 'VITE', iocc=1, scal=vit0, nbret=iret2)
                if (iret1 .eq. 0) dep0 = 'nul'
                if (iret2 .eq. 0) vit0 = 'nul'
                call utmess('I', 'DYNAMIQUE_79', nk=2, valk=[dep0, vit0])
                call utmess('I', 'DYNAMIQUE_84', nk=1, valk=['calculee'])
            end if
        end if
    end if

!   7 - Simulation duration
    call utmess('I', 'DYNAMIQUE_80', nr=2, valr=[tinit, tfin])

!   8 - Archivage
    call dtmget(sd_dtm, _ARCH_PER, iscal=iparch)
    call dtmget(sd_dtm, _AR_LINST, lonvec=nbsav_forc)
    call utmess('I', 'DYNAMIQUE_81', ni=2, vali=[iparch, nbsav_forc-1])

end subroutine
