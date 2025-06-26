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
!
subroutine op0149()
    implicit none
!
!     OPERATEUR:  MODI_BASE_MODALE
!
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/modiba.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=8) :: nomres, basemo, modefl, typflu
    character(len=8) :: kbid
    character(len=16) :: typres, nomcmd
    character(len=19) :: basefl
    character(len=24) :: numo, vite, refefl, fsic, fsvi
    aster_logical :: newres, lnuor, lamor, lamoru, nocopl, numok
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iamo1, iamo2, iamor, ibid, idec, ifm
    integer(kind=8) :: ifsic, ifsvi, imasse, imin, inumo, inuo1, inuo2
    integer(kind=8) :: inuor, ireffl, itypfl, j, jpara, k
    integer(kind=8) :: na1, nbamo1, nbamun, nbmfl, nbmod2, nbmode, nbnuo1
    integer(kind=8) :: nbnuo2, nbnuor, nbvite, niv, numode, numvit, nuomin
!
    real(kind=8) :: amorun, rtamp
    integer(kind=8), pointer :: ordr(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
!-----0.VERIFICATIONS AVANT EXECUTION
!
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nbnuo1)
    nbnuo1 = abs(nbnuo1)
    if (nbnuo1 .ne. 0) then
        call getvr8(' ', 'AMOR_REDUIT', nbval=0, nbret=na1)
        nbamo1 = abs(na1)
        if (nbamo1 .ne. 0) then
            if (nbamo1 .ne. nbnuo1) then
                call utmess('F', 'ALGORITH9_57')
            end if
        end if
    end if
!
!
!     ---RECUPERATION DU NIVEAU D'IMPRESSION---
!
    call infmaj()
    call infniv(ifm, niv)
!
!-----1.RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getres(nomres, typres, nomcmd)
!
    newres = .true.
    call getvid(' ', 'BASE', scal=basemo, nbret=ibid)
    if (basemo .eq. nomres) newres = .false.
!
    call getvid(' ', 'BASE_ELAS_FLUI', scal=basefl, nbret=ibid)
    call getvis(' ', 'NUME_VITE_FLUI', scal=numvit, nbret=ibid)
!
    lnuor = .false.
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nbnuo1)
    nbnuo1 = abs(nbnuo1)
    if (nbnuo1 .ne. 0) then
        lnuor = .true.
        call wkvect('&&OP0149.TEMP.NUO1', 'V V I', nbnuo1, inuo1)
        call getvis(' ', 'NUME_ORDRE', nbval=nbnuo1, vect=zi(inuo1), nbret=ibid)
    end if
!
    lamor = .false.
    lamoru = .false.
    call getvr8(' ', 'AMOR_REDUIT', nbval=0, nbret=na1)
    nbamo1 = abs(na1)
    if (nbamo1 .ne. 0) then
        lamor = .true.
        call wkvect('&&OP0149.TEMP.AMO1', 'V V R', nbamo1, iamo1)
        if (na1 .ne. 0) then
            call getvr8(' ', 'AMOR_REDUIT', nbval=nbamo1, vect=zr(iamo1), nbret=ibid)
        end if
    else
        call getvr8(' ', 'AMOR_UNIF', nbval=0, nbret=nbamun)
        if (nbamun .ne. 0) then
            lamoru = .true.
            call getvr8(' ', 'AMOR_UNIF', scal=amorun, nbret=ibid)
        end if
    end if
!
!
!-----2.VERIFICATIONS A L'EXECUTION
!
!-----2.1.ERREUR FATALE SI LE CONCEPT MODE_MECA D'ENTREE N'EST PAS
!         CELUI AYANT SERVI AU CALCUL DE COUPLAGE FLUIDE-STRUCTURE
!
    refefl = basefl//'.REMF'
    call jeveuo(refefl, 'L', ireffl)
    modefl = zk8(ireffl+1)
    if (basemo .ne. modefl) then
        call utmess('F', 'ALGORITH9_58')
    end if
!
!-----2.2.ERREUR FATALE SI NUME_VITE_FLUI INVALIDE
!
    vite = basefl//'.VITE'
    call jelira(vite, 'LONUTI', nbvite)
    ASSERT(numvit .gt. 0 .and. numvit .le. nbvite)
!
!-----2.3.ERREUR FATALE SI TOUS LES MODES NON COUPLES SONT RETENUS
!         (MOT-CLE <NUME_ORDRE> NON UTILISE) ET NOMBRE D'ARGUMENTS
!         INVALIDE POUR LE MOT-CLE <AMOR_REDUIT>
!
    call jelira(basemo//'           .ORDR', 'LONUTI', nbmode)
    call jeveuo(basemo//'           .ORDR', 'L', vi=ordr)
!
!
!--------------------------------------------------------------------
    numo = basefl//'.NUMO'
    call jelira(numo, 'LONUTI', nbmfl)
    call jeveuo(numo, 'L', inumo)
!
    nbmod2 = nbmode-nbmfl
    if (.not. lnuor .and. lamor .and. nbamo1 .ne. nbmod2) then
        call utmess('F', 'ALGORITH9_60')
    end if
!
!
!-----3.CONSTITUTION DE LA LISTE DES NUMEROS D'ORDRE DES MODES RETENUS
!       POUR LA RECONSTRUCTION DE LA BASE MODALE
!       (MODES NON PERTURBES + MODES PRIS EN COMPTE POUR LE COUPLAGE)
!       LE CAS ECHEANT ON CREE UNE LISTE D'AMORTISSEMENTS REDUITS QUI
!       SERONT AFFECTES AUX MODES NON PERTURBES
!
!     NUMOI = BASEMO//'           .NUMO'
!     CALL JEVEUO(NUMOI,'L',INUMOI)
!
!-----3.1.SI ON CREE UN NOUVEAU CONCEPT DE TYPE MODE_MECA EN SORTIE
!
    if (newres) then
!
!-------3.1.1.SI DONNEE D'UNE LISTE DE NUMEROS D'ORDRE PAR <NUME_ORDRE>
!
        if (lnuor) then
!
            call wkvect('&&OP0149.TEMP.NUO2', 'V V I', nbnuo1, inuo2)
            nbnuo2 = 0
            if (lamor) call wkvect('&&OP0149.TEMP.AMO2', 'V V R', nbnuo1, iamo2)
!
!---------ON NE RETIENT QUE LES NUMEROS D'ORDRE QUI CORRESPONDENT
!         EFFECTIVEMENT A DES MODES NON COUPLES ET ON NOTE LE CAS
!         ECHEANT LES VALEURS D'AMORTISSEMENTS FOURNIES EN REGARD
!
            do i = 1, nbnuo1
                nocopl = .true.
                numok = .false.
                numode = zi(inuo1+i-1)
                do j = 1, nbmfl
                    if (zi(inumo+j-1) .eq. numode) then
                        nocopl = .false.
                        goto 12
                    end if
                end do
12              continue
                do k = 1, nbmode
                    call rsadpa(basemo, 'L', 1, 'NUME_MODE', ordr(k), &
                                0, sjv=jpara, styp=kbid)
                    if (zi(jpara) .eq. numode) then
                        numok = .true.
                        goto 14
                    end if
                end do
14              continue
                if (nocopl .and. numok) then
                    nbnuo2 = nbnuo2+1
                    zi(inuo2+nbnuo2-1) = numode
                    if (lamor) zr(iamo2+nbnuo2-1) = zr(iamo1+i-1)
                end if
            end do
!
!---------CONSTITUTION DES LISTES
!
            if (nbnuo2 .eq. 0) then
                call utmess('F', 'ALGORITH9_61')
            else
                nbnuor = nbnuo2+nbmfl
                call wkvect('&&OP0149.TEMP.NUOR', 'V V I', nbnuor, inuor)
                call wkvect('&&OP0149.TEMP.AMOR', 'V V I', nbnuor, iamor)
                do i = 1, nbnuo2
                    zi(inuor+i-1) = zi(inuo2+i-1)
                    if (lamor) then
                        zr(iamor+i-1) = zr(iamo2+i-1)
                    else if (lamoru) then
                        zr(iamor+i-1) = amorun
                    end if
                end do
                do i = nbnuo2+1, nbnuor
                    zi(inuor+i-1) = zi(inumo+i-nbnuo2-1)
                end do
                do i = 1, nbnuor-1
                    nuomin = zi(inuor+i-1)
                    imin = i
                    do j = i+1, nbnuor
                        if (zi(inuor+j-1) .lt. nuomin) then
                            nuomin = zi(inuor+j-1)
                            imin = j
                        end if
                    end do
                    zi(inuor+imin-1) = zi(inuor+i-1)
                    zi(inuor+i-1) = nuomin
                    if (lamor .or. lamoru) then
                        rtamp = zr(iamor+imin-1)
                        zr(iamor+imin-1) = zr(iamor+i-1)
                        zr(iamor+i-1) = rtamp
                    end if
                end do
            end if
!
!-------3.1.2.SINON
!
        else
!
!---------SI DONNEE D'AMORTISSEMENTS REDUITS, ON RETIENT TOUS LES MODES
!
            if (lamor .or. lamoru) then
                nbnuor = nbmode
                call wkvect('&&OP0149.TEMP.NUOR', 'V V I', nbnuor, inuor)
                call wkvect('&&OP0149.TEMP.AMOR', 'V V I', nbnuor, iamor)
                do i = 1, nbnuor
                    call rsadpa(basemo, 'L', 1, 'NUME_MODE', ordr(i), &
                                0, sjv=jpara, styp=kbid)
                    zi(inuor+i-1) = zi(jpara)
                end do
                idec = 0
                do i = 1, nbnuor
                    nocopl = .true.
                    numode = zi(inuor+i-1)
                    do j = 1, nbmfl
                        if (zi(inumo+j-1) .eq. numode) then
                            nocopl = .false.
                            goto 33
                        end if
                    end do
33                  continue
                    if (nocopl) then
                        if (lamor) then
                            idec = idec+1
                            zr(iamor+i-1) = zr(iamo1+idec-1)
                        else if (lamoru) then
                            zr(iamor+i-1) = amorun
                        end if
                    end if
                end do
!
!---------SINON, SEULS LES MODES COUPLES SONT RETENUS
!
            else
                nbnuor = nbmfl
                call wkvect('&&OP0149.TEMP.NUOR', 'V V I', nbnuor, inuor)
                call wkvect('&&OP0149.TEMP.AMOR', 'V V I', nbnuor, iamor)
                do i = 1, nbmfl
                    zi(inuor+i-1) = zi(inumo+i-1)
                end do
            end if
!
        end if
!
!-----3.2.SINON (ON MODIFIE LE CONCEPT D'ENTREE DE TYPE MODE_MECA)
!         => TOUS LES MODES SONT RETENUS
!
    else
!
        nbnuor = nbmode
        call wkvect('&&OP0149.TEMP.NUOR', 'V V I', nbnuor, inuor)
        call wkvect('&&OP0149.TEMP.AMOR', 'V V I', nbnuor, iamor)
        do i = 1, nbnuor
            call rsadpa(basemo, 'L', 1, 'NUME_MODE', ordr(i), &
                        0, sjv=jpara, styp=kbid)
            zi(inuor+i-1) = zi(jpara)
        end do
        if ((lnuor .and. lamor) .or. (lnuor .and. lamoru)) then
            do i = 1, nbnuo1
                nocopl = .true.
                numok = .false.
                numode = zi(inuo1+i-1)
                do j = 1, nbmfl
                    if (zi(inumo+j-1) .eq. numode) then
                        nocopl = .false.
                        goto 53
                    end if
                end do
53              continue
                do k = 1, nbmode
                    call rsadpa(basemo, 'L', 1, 'NUME_MODE', ordr(k), &
                                0, sjv=jpara, styp=kbid)
                    if (zi(jpara) .eq. numode) then
                        numok = .true.
                        goto 55
                    end if
                end do
55              continue
                if (nocopl .and. numok) then
                    if (lamor) zr(iamor+numode-1) = zr(iamo1+i-1)
                    if (lamoru) zr(iamor+numode-1) = amorun
                end if
            end do
        else if (lamor .or. lamoru) then
            idec = 0
            do i = 1, nbnuor
                nocopl = .true.
                numode = zi(inuor+i-1)
                do j = 1, nbmfl
                    if (zi(inumo+j-1) .eq. numode) then
                        nocopl = .false.
                        goto 58
                    end if
                end do
58              continue
                if (nocopl) then
                    if (lamor) then
                        idec = idec+1
                        zr(iamor+i-1) = zr(iamo1+idec-1)
                    else if (lamoru) then
                        zr(iamor+i-1) = amorun
                    end if
                end if
            end do
        end if
!
    end if
!
!
!-----4.RECUPERATION DU TYPE DE LA CONFIGURATION ETUDIEE
!
    typflu = zk8(ireffl)
    fsic = typflu//'           .FSIC'
    call jeveuo(fsic, 'L', ifsic)
    itypfl = zi(ifsic)
    imasse = -1
    if (itypfl .eq. 4) then
        fsvi = typflu//'           .FSVI'
        call jeveuo(fsvi, 'L', ifsvi)
        imasse = zi(ifsvi)
    end if
!
!
!-----5.RECONSTRUCTION OU MODIFICATION DE LA BASE MODALE EN FONCTION
!       DU TYPE DE LA CONFIGURATION ETUDIEE
!
    call modiba(nomres, basemo, basefl, numvit, newres, &
                itypfl, imasse, zi(inuor), nbnuor, zi(inumo), &
                nbmfl)
!
    call jedema()
end subroutine
