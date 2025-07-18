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
subroutine caldis(fremax, fremin, pas, frexci, nbptmd, &
                  nbmode, lismod, fremod, amomod, nindex, &
                  npdsc3, frefin)
    implicit none
#include "jeveux.h"
#include "asterfort/discrt.h"
#include "asterfort/disexc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ordr8.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbptmd, nbmode, nindex, npdsc3, lismod(*)
    real(kind=8) :: fremax, fremin, pas, fremod(*), amomod(*), frefin
    character(len=4) :: frexci
!
!  BUT: CALCUL DE LA DISCRETISATION FREQUENTIELLE POUR
!                 LE CALCUL DYNAMIQUE ALEATOIRE
!
! IN  : FREMAX : FREQ MAX DE LA DISCRETISATION
! IN  : FREMIN : FREQ MIN DE LA DISCRETISATION
! IN  : PAS    : PAS DE LA DISCRETISATION
! IN  : FREXCI : FREQUENCES EXCITATION: AVEC OU SANS
! IN  : NBPTMD : NOMBRE DE POINTS PAR PICS
! IN  : NBMODE : NOMBRE DE MODES DYNAMIQUES
! IN  : LISMOD : LISTE DES NUMEROS D'ORDRE DES MODES
! IN  : FREMOD : FREQUENCE DU MODE MECA
! IN  : AMOMOD : AMORTISSEMENTS MODAUX
! IN  : NINDEX : NOMBRE  D INDICES RECUPERES DANS L INTERPECTRE EXC
! OUT : NPDSC3 : NOMBRE DE VALEURS DE FREQUENCE
! OUT : FREFIN : DERNIERE FREQUENCE ( LA PLUS GRANDE)
!-----------------------------------------------------------------------
!
    real(kind=8) :: frema1, fredeb
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, iadii1, iadr1, iadsc1, iadsc2, iadsc3
    integer(kind=8) :: ibid1, icont1, ifreq1, igex1, ilfex, illex, ilong1
    integer(kind=8) :: imode, inajou, nbmax, npdsc0, npdsc2, mxval
    real(kind=8) :: amor, f0, f1, f2, freq
    real(kind=8) :: pasmin, r8ecar
!-----------------------------------------------------------------------
    call jemarq()
!
!---------CAS D UNE DISCRETISATION DEMANDEE PAR L UTILISATEUR
!
    call jeveuo('&&OP0131.LIADRFEX1', 'E', ilfex)
    call jeveuo('&&OP0131.LIADRLEX1', 'E', illex)
!
    if ((fremin .ne. -1.d0) .and. (fremax .ne. -1.d0)) then
        if (pas .eq. -1.d0) then
            pas = (fremax-fremin)/100.d0
        end if
        npdsc3 = int((fremax-fremin)/pas+1.d-6)+1
        call wkvect('&&OP0131.DISCR3', 'V V R8', npdsc3, iadsc3)
        do i1 = 1, npdsc3-1
            zr(iadsc3-1+i1) = fremin+(i1-1)*pas
        end do
        zr(iadsc3-1+npdsc3) = fremax
        fredeb = fremin
        frefin = fremax
        goto 332
    end if
!
!---------CALCUL DE LA DISCRETISATION MINI IMPOSEE PAR LES FREQUENCES
!         PROPRES
!
!---------ON ALLOUE UN TABLEAU DE TRAVAIL POUR STOCKER TOUTES LES
!         DISCRETISATIONS
!
    ilong1 = 0
    if (frexci .eq. 'AVEC') then
        do igex1 = 1, nindex*(nindex+1)/2
            ilong1 = ilong1+zi(illex)
        end do
    end if
    call wkvect('&&OP0131.DISCR1', 'V V R8', nbmode*nbptmd+ilong1+nbmode, iadsc1)
    frema1 = 0.d0
    do i1 = 1, nbmode
        imode = lismod(i1)
        f0 = fremod(imode)
        if (f0 .gt. frema1) frema1 = f0
        amor = amomod(i1)
        fredeb = 0.d0
        frefin = 2.d0*f0
        call discrt(f0, fredeb, frefin, nbptmd, amor, &
                    zr(iadsc1+(i1-1)*nbptmd))
    end do
    icont1 = 0
    mxval = nindex*(nindex+1)/2
    if (frexci .eq. 'AVEC') then
        do igex1 = 1, nindex*(nindex+1)/2
            iadr1 = zi(illex+mxval+1)
            do i2 = 1, zi(illex)
                icont1 = icont1+1
                zr(iadsc1+nbmode*nbptmd-1+icont1) = zr(iadr1-1+i2)
            end do
        end do
    end if
!
!
!----AJOUT DES FREQUENCES PROPRES DANS LE CAS D AMORTISSEMENT NON NUL
!
    do i2 = 1, nbmode
        imode = lismod(i2)
        amor = amomod(i2)
        freq = fremod(imode)
        if (amor .ne. 0.d0) then
            zr(iadsc1+nbmode*nbptmd-1+ilong1+i2) = freq
        end if
    end do
    npdsc0 = ilong1+nbmode
!
!---------TRI DES FREQUENCES OBTENUES-RANGEMENT DANS &&..DISCR2
!
    call wkvect('&&OP0131.DISCI1', 'V V I', nbmode*nbptmd+npdsc0, iadii1)
    call ordr8(zr(iadsc1), nbmode*nbptmd+npdsc0, zi(iadii1))
    call wkvect('&&OP0131.DISCR2', 'V V R8', nbmode*nbptmd+npdsc0+50, iadsc2)
    zr(iadsc2) = zr(iadsc1-1+zi(iadii1))
    i2 = 1
    do i1 = 2, nbmode*nbptmd+npdsc0
        if (zr(iadsc1-1+zi(iadii1+i1-1)) .gt. zr(iadsc1-1+zi(iadii1+i1-2))) then
            i2 = i2+1
            zr(iadsc2+i2-1) = zr(iadsc1-1+zi(iadii1+i1-1))
        end if
    end do
    npdsc2 = i2
!
!---------PRISE EN COMPTE DU PAS MINI
!
!  - ON TIENT COMPTE D'UN PAS MINI EXPRIME OU NON PAR L UTILISATEUR
!  --- IL EST PAR DEFAUT EGAL A (FREFIN -FREDEB)/100
!
    if ((fremin .eq. -1.d0) .and. (fremax .eq. -1.d0)) then
        fredeb = 0.d0
        frefin = 2.0d0*frema1
    else
        fredeb = fremin
        frefin = fremax
    end if
    if (pas .eq. -1.d0) then
        pasmin = (frefin-fredeb)/100.0d0
    else
        pasmin = pas
    end if
    nbmax = int((frefin-fredeb)/pasmin)
    call wkvect('&&OP0131.DISCR3', 'V V R8', npdsc2+nbmax, iadsc3)
    npdsc3 = npdsc2
    inajou = 0
    zr(iadsc3) = zr(iadsc2)
    do ifreq1 = 2, npdsc2
        f2 = zr(iadsc2-1+ifreq1)
        f1 = zr(iadsc2-1+ifreq1-1)
        r8ecar = f2-f1
        if ((r8ecar .gt. pasmin) .and. (f2 .le. frefin)) then
            inajou = int(r8ecar/pasmin)
            do ibid1 = 1, inajou
                zr(iadsc3-1+ifreq1-1+npdsc3-npdsc2+ibid1) = zr(iadsc2+ &
                                                               ifreq1-2)+ibid1*pasmin
            end do
        end if
        npdsc3 = npdsc3+inajou
        inajou = 0
        zr(iadsc3-1+ifreq1+npdsc3-npdsc2) = zr(iadsc2-1+ifreq1)
    end do
332 continue
!
!------CALCUL DES DSP EXCITS DANS LA DISCRETISATION REPONSE
!
    call disexc(nindex, ilfex, illex, npdsc3, iadsc3)
    call jedema()
end subroutine
