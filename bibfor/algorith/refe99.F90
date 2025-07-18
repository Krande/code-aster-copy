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

subroutine refe99(nomres)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/refdaj.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres
!
!     BUT:
!       RECUPERER LES NOMS UTILISATEUR DES CONCEPTS ASSOCIES AUX
!       MATRICES ASSEMBLEES CONSIDEREES - EFFECTUER QUELQUES CONTROLES
!       CREER LE .REFD
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NOMRES   : NOM DE LA SD_RESULTAT
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: i, ioc1, ioc3, ioc4, ioc5, ier, ibid, ibmo, imint, inmax
    integer(kind=8) :: ltmome, nbg, nbmome, ltnbmo, ltnbmax, nbli, nbmax, vali(2)
    integer(kind=8) :: nbtot, nbold(1), nbmod1, nbmod2, nbmout, nbmodo(1)
    integer(kind=8) :: nbmm, nbbm, nbmi, ioccmi, ioccbase
!
    real(kind=8) :: rbid
    complex(kind=8) :: cbid
!
    character(len=8) :: k8b, resul1, resul2, momeca, mostat
    character(len=19) :: numddl, numbis
    character(len=24) :: raid, mass, intf, intfb, amor, concep(3), valk(4), kbid
!
    aster_logical :: noseul
!
!-----------------------------------------------------------------------
!
    call jemarq()
    numddl = ' '
    raid = ' '
    mass = ' '
    amor = ' '
!
! --- DETERMINATION DU TYPE DE BASE
!
    call getfac('CLASSIQUE', ioc1)
    call getfac('RITZ', ioc3)
    call getfac('DIAG_MASS', ioc4)
    call getfac('ORTHO_BASE', ioc5)
!
! --- CAS CLASSIQUE
!
    if (ioc1 .gt. 0) then
        numbis = ' '
        call getvid('CLASSIQUE', 'INTERF_DYNA', iocc=1, scal=intf, nbret=ier)
        call dismoi('NOM_NUME_DDL', intf, 'INTERF_DYNA', repk=numddl)
        call getvid('CLASSIQUE', 'MODE_MECA', iocc=1, nbval=0, nbret=nbmome)
        nbmome = -nbmome
!
        call wkvect('&&REFE99.LIST.MODE_MECA', 'V V K8', nbmome, ltmome)
        call wkvect('&&REFE99.LIST.NBMOD', 'V V I', nbmome, ltnbmo)
        call wkvect('&&REFE99.LIST.NBMODMAX', 'V V I', nbmome, ltnbmax)
!
        call getvid('CLASSIQUE', 'MODE_MECA', iocc=1, nbval=nbmome, vect=zk8(ltmome), &
                    nbret=ibid)
!
        call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, nbval=0, nbret=nbli)
        nbli = -nbli
!       if nbli = 0: one will take all modes in each MODE_MECA, otherwise:
        if (nbli .ge. 1) then
            if (nbli .eq. 1) then
!               Apply the single NMAX_MODE criterion to all of the modal base
                nbmax = 0
                call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, scal=nbmax, nbret=ier)
                do i = 1, nbmome
                    zi(ltnbmax+i-1) = nbmax
                end do
            else if (nbli .eq. nbmome) then
!               Use the NMAX_MODE criteria, defined for each modal base separately
                call getvis('CLASSIQUE', 'NMAX_MODE', iocc=1, nbval=nbmome, vect=zi(ltnbmax), &
                            nbret=ier)
            else
!               Incoherence in the input data
                vali(1) = nbmome
                vali(2) = nbli
                call utmess('F+', 'DEFIBASEMODALE1_31', ni=2, vali=vali)
                call utmess('F', 'ALGORITH14_32')
            end if
        end if
!
        do i = 1, nbmome
            momeca = zk8(ltmome-1+i)
            call dismoi('REF_RIGI_PREM', momeca, 'RESU_DYNA', repk=raid, arret='C')
            call dismoi('REF_MASS_PREM', momeca, 'RESU_DYNA', repk=mass, arret='C')
            call dismoi('REF_AMOR_PREM', momeca, 'RESU_DYNA', repk=amor, arret='C')
            call dismoi('NUME_DDL', momeca, 'RESU_DYNA', repk=numbis)
!           Check for nume_ddl coherence of each mode_meca with that of interf_dyna
            if (numbis(1:14) .ne. numddl(1:14)) then
                valk(1) = momeca
                valk(2) = numbis(1:8)
                valk(3) = intf
                valk(4) = numddl(1:8)
                call utmess('F', 'ALGORITH14_24', nk=4, valk=valk)
            end if
!           Determine the real number of modes to recuperate from each base
            call rsorac(momeca, 'LONUTI', 0, rbid, kbid, &
                        cbid, rbid, 'ABSOLU', nbmodo, 1, &
                        ibid)
            if (nbli .ge. 1) then
                nbmout = zi(ltnbmax+i-1)
                if (nbmodo(1) .lt. nbmout) then
                    valk = momeca
                    vali(1) = nbmout
                    vali(2) = nbmodo(1)
                    call utmess('I', 'ALGORITH15_92', sk=valk(1), ni=2, vali=vali)
                else
                    nbmodo(1) = nbmout
                end if
            end if
            zi(ltnbmo+i-1) = nbmodo(1)
!           Add a reference in the dynamic result data structure
            concep(1) = raid
            concep(2) = mass
            concep(3) = amor
            call refdaj('F', nomres, nbmodo(1), numddl, 'DYNAMIQUE', &
                        concep, ier)
            concep(1) = intf
            call refdaj('F', nomres, 0, numddl, 'INTERF_DYNA', &
                        concep, ier)
        end do
!
    end if
!
! --- CAS RITZ
!
    if (ioc3 .gt. 0) then
!
        ASSERT(ioc3 .le. 2)
        nbg = 0
        ioccmi = 1
        ioccbase = 1
        do i = 1, ioc3
            call getvid('RITZ', 'MODE_MECA', iocc=i, nbval=0, nbret=nbmm)
            call getvid('RITZ', 'BASE_MODALE', iocc=i, nbval=0, nbret=nbbm)
            call getvid('RITZ', 'MODE_INTF', iocc=i, nbval=0, nbret=nbmi)
            nbg = nbg-nbmi
            if (nbmi .ne. 0) then
                ioccmi = i
            end if
            if ((nbmm .ne. 0) .or. (nbbm .ne. 0)) then
                ioccbase = i
            end if
        end do
!
        if (ioc3 .eq. 2) then
            if (nbg .ne. 1) then
                call utmess('F', 'DEFIBASEMODALE1_51')
            end if
        else
            if (nbmm .eq. 0) then
                call utmess('F', 'DEFIBASEMODALE1_1')
            end if
        end if
!
        noseul = .false.
        call getvid('RITZ', 'MODE_MECA', iocc=ioccbase, nbval=0, nbret=nbg)
        nbg = -nbg
        call getvid('RITZ', 'MODE_INTF', iocc=ioccmi, nbval=0, nbret=ier)
        if ((ier .gt. 0) .or. (nbg .gt. 1)) noseul = .true.
!
!       Reference numbering is required in case more than one modal base is given
        call getvid('    ', 'NUME_REF', iocc=1, scal=numddl, nbret=ier)
        if ((ier .eq. 0) .and. noseul) then
            call utmess('F', 'DEFIBASEMODALE1_9')
        end if
!
        intf = ' '
        call getvid('  ', 'INTERF_DYNA', iocc=1, nbval=0, nbret=ier)
        if (ier .lt. 0) then
            call getvid('  ', 'INTERF_DYNA', iocc=1, scal=intf, nbret=ier)
            call dismoi('NOM_NUME_DDL', intf, 'INTERF_DYNA', repk=numddl)
        end if
!
        call getvid('RITZ', 'BASE_MODALE', iocc=ioccbase, scal=resul1, nbret=ibmo)
        call getvid('RITZ', 'MODE_INTF', iocc=ioccmi, scal=resul2, nbret=imint)
!
!       BASE_MODALE kw treatment (with INTERF_DYNA and MODE_INTF on the 2nd occurence)
        if (ibmo .ne. 0) then
            call refdcp(resul1, nomres)
            call dismoi('NB_MODES_TOT', resul1, 'RESULTAT', repi=nbmod1)
        else
!           MODE_MECA kw treatment, similar to what is done for the "classique" case
            call getvid('RITZ', 'MODE_MECA', iocc=ioccbase, nbval=0, nbret=nbmome)
            nbmome = -nbmome
!
            call wkvect('&&REFE99.LIST.MODE_MECA', 'V V K8', nbmome, ltmome)
            call wkvect('&&REFE99.LIST.NBMOD', 'V V I', nbmome, ltnbmo)
            call wkvect('&&REFE99.LIST.NBMODMAX', 'V V I', nbmome, ltnbmax)
!
            call getvid('RITZ', 'MODE_MECA', iocc=ioccbase, nbval=nbmome, vect=zk8(ltmome), &
                        nbret=ibid)
!
            call getvis('RITZ', 'NMAX_MODE', iocc=ioccbase, nbval=0, nbret=nbli)
            nbli = -nbli
            if (nbli .eq. 0) then
!               Select all modes from each modal base
                do i = 1, nbmome
                    zi(ltnbmax+i-1) = 9999
                end do
            else if (nbli .eq. 1) then
!               Apply the single NMAX_MODE criterion to all of the modal base
                call getvis('RITZ', 'NMAX_MODE', iocc=ioccbase, scal=nbmax, nbret=ibid)
                do i = 1, nbmome
                    zi(ltnbmax+i-1) = nbmax
                end do
            else if (nbli .eq. nbmome) then
!               Use the NMAX_MODE criteria, defined for each modal base separately
                call getvis('RITZ', 'NMAX_MODE', iocc=ioccbase, nbval=nbmome, vect=zi(ltnbmax), &
                            nbret=ibid)
            else
!               Incoherence in the input data
                vali(1) = nbmome
                vali(2) = nbli
                call utmess('F', 'DEFIBASEMODALE1_31', ni=2, vali=vali)
            end if
!
            nbmod1 = 0
            do i = 1, nbmome
                momeca = zk8(ltmome-1+i)
                call dismoi('REF_RIGI_PREM', momeca, 'RESU_DYNA', repk=raid, arret='C')
                call dismoi('REF_MASS_PREM', momeca, 'RESU_DYNA', repk=mass, arret='C')
                call dismoi('REF_AMOR_PREM', momeca, 'RESU_DYNA', repk=amor, arret='C')
                call dismoi('NUME_DDL', momeca, 'RESU_DYNA', repk=numbis)
                if (numddl .eq. ' ') numddl = numbis
!               Determine the real number of modes to recuperate from each base
                call rsorac(momeca, 'LONUTI', 0, rbid, kbid, &
                            cbid, rbid, 'ABSOLU', nbmodo, 1, &
                            ibid)
                nbmout = zi(ltnbmax+i-1)
                if (nbmodo(1) .lt. nbmout) then
                    valk = momeca
                    vali(1) = nbmout
                    vali(2) = nbmodo(1)
                    call utmess('I', 'ALGORITH15_92', sk=valk(1), ni=2, vali=vali)
                else
                    nbmodo(1) = nbmout
                end if
                zi(ltnbmo+i-1) = nbmodo(1)
                nbmod1 = nbmod1+nbmodo(1)
!               Add a reference in the dynamic result data structure
                concep(1) = raid
                concep(2) = mass
                concep(3) = amor
                call refdaj('F', nomres, nbmodo(1), numddl, 'DYNAMIQUE', &
                            concep, ier)
            end do
        end if
        if (imint .gt. 0) then
!           Treating the MODE_INTF kw (2nd RITZ entry) for the static modes
!           Maximum number of static modes to extract : nbmod2
            call getvis('RITZ', 'NMAX_MODE', iocc=ioccmi, scal=nbmod2, nbret=inmax)
!           Number of modes that actually exist in the static base : nbold
            call rsorac(resul2, 'LONUTI', 0, rbid, k8b, &
                        cbid, rbid, 'ABSOLU', nbold, 1, &
                        ibid)
            if (inmax .eq. 0) then
                nbmod2 = nbold(1)
            else
                nbmod2 = min(nbmod2, nbold(1))
            end if
!
            concep(1) = resul2
            call refdaj('F', nomres, nbmod2, numddl, 'INTERF_STAT', &
                        concep, ier)
        else
            nbmod2 = 0
        end if
!
        nbtot = nbmod1+nbmod2
        if (nbtot .le. 0) then
            call utmess('F', 'DEFIBASEMODALE1_50')
        end if
!
        call dismoi('REF_INTD_DERN', nomres, 'RESU_DYNA', repk=intfb, arret='C', &
                    ier=ier)
        if ((intf .ne. ' ') .and. (intf .ne. intfb)) then
            concep(1) = intf
            call refdaj('F', nomres, 0, numddl, 'INTERF_DYNA', &
                        concep, ier)
        end if
    end if
!
! --- DIAGONALISATION DE LA MATRICE DE MASSE
!
    if (ioc4 .gt. 0) then
        intf = ' '
        call getvid('DIAG_MASS', 'MODE_MECA', iocc=1, scal=momeca, nbret=ibid)
!
        call dismoi('REF_RIGI_PREM', momeca, 'RESU_DYNA', repk=raid, arret='C')
        call dismoi('REF_MASS_PREM', momeca, 'RESU_DYNA', repk=mass, arret='C')
        call dismoi('REF_AMOR_PREM', momeca, 'RESU_DYNA', repk=amor, arret='C')
        call dismoi('NOM_NUME_DDL', mass, 'MATR_ASSE', repk=numddl)
        call dismoi('NB_MODES_TOT', momeca, 'RESULTAT', repi=nbmod1)
        concep(1) = raid
        concep(2) = mass
        concep(3) = amor
        call refdaj('F', nomres, nbmod1, numddl, 'DYNAMIQUE', &
                    concep, ier)
!
        call getvid('DIAG_MASS', 'MODE_STAT', iocc=1, scal=mostat, nbret=ibid)
        call dismoi('NB_MODES_TOT', mostat, 'RESULTAT', repi=nbmod2)
        concep(1) = mostat
!       Note that it is volontary to save the numbering associated to the dynamic modes
!       because later on we call copmod upon the static modes, modifying their nume_ddl
        call refdaj('F', nomres, nbmod2, numddl, 'INTERF_STAT', &
                    concep, ier)
    end if
!
! --- CAS ORTHO_BASE
!
    if (ioc5 .gt. 0) then
        call getvid('ORTHO_BASE', 'BASE', iocc=1, scal=resul1, nbret=ibid)
        call refdcp(resul1, nomres)
    end if
!
    call jedema()
!
end subroutine
