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
subroutine nmvcmx(mate, mailla, comref, comval)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmvcex.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=24) :: mate, comref
    character(len=19) :: comval
    character(len=8) :: mailla
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (UTILITAIRE)
!
! RECHERCHE DES MAXIMUM/ MIMINUM DES VARIABLES DE COMMANDES
!
! ----------------------------------------------------------------------
!
!
! IN  MATE   : CHAMP MATERIAU
! IN  COMREF : VARI_COM DE REFERENCE
! IN  COMVAL : VARI_COM
!
!
!
!
    integer(kind=8) :: nbcmp, nbcmp2
    character(len=8) :: valk(5)
    character(len=19) :: chsref, chscom
    character(len=24) :: vrcplu, vrcref
    integer(kind=8) :: jcesd, jcesl, nbma, nbpt, nbsp, icmp
    integer(kind=8) :: jcrsd, jcrsl, ima, ipt, isp, iad, iad2
    integer(kind=8) :: imamax, imamin, iref
    real(kind=8) :: valmin, valmax, valr(2)
    real(kind=8) :: valeur, valref
    character(len=8), pointer :: cvrcvarc(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: crsv(:) => null()
    character(len=8), pointer :: cvrcnom(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    chscom = '&&NMVCMX.COMVAL_SIM'
    chsref = '&&NMVCMX.COMREF_SIM'
!
! --- EXTRACTION DES VARIABLES DE COMMANDE
!
    call nmvcex('TOUT', comref, vrcref)
    call nmvcex('TOUT', comval, vrcplu)
!
! --- TRANSFO. EN CHAM_NO_S
!
    call celces(vrcplu, 'V', chscom)
    call celces(vrcref, 'V', chsref)
!
    call utmess('A+', 'MECANONLINE2_97')
!
!     CALCUL DU MIN / MAX
!
!     DESCRIPTEUR
    call jeveuo(chscom//'.CESD', 'L', jcesd)
    call jeveuo(chsref//'.CESD', 'L', jcrsd)
!     PRESENCE DES CMP (R)
    call jeveuo(chscom//'.CESL', 'L', jcesl)
    call jeveuo(chsref//'.CESL', 'L', jcrsl)
!     VALEUR DES CMP (R)
    call jeveuo(chscom//'.CESV', 'L', vr=cesv)
    call jeveuo(chsref//'.CESV', 'L', vr=crsv)
!
!     RECUPERATION DES NOMS DES VARC
    call jelira(mate(1:8)//'.CVRCNOM', 'LONMAX', ival=nbcmp2)
    call jeveuo(mate(1:8)//'.CVRCNOM', 'L', vk8=cvrcnom)
    call jeveuo(mate(1:8)//'.CVRCVARC', 'L', vk8=cvrcvarc)
!
    nbma = zi(jcesd-1+1)
!
    do icmp = 1, nbcmp2
        valmax = -r8maem()
        valmin = r8maem()
        imamin = 0
        imamax = 0
        iref = 0
        if (cvrcvarc(icmp) .eq. 'TEMP' .or. cvrcvarc(icmp) .eq. 'SECH') then
            iref = 1
        end if
!
        do ima = 1, nbma
            nbcmp = zi(jcesd-1+5+4*(ima-1)+3)
            if (nbcmp .eq. 0) goto 40
            call cesexi('C', jcrsd, jcrsl, ima, 1, &
                        1, icmp, iad2)
            if (iad2 .le. 0) goto 40
!
!
!           VALEURS DE REFERENCE
            if (iref .eq. 1) then
                call cesexi('C', jcrsd, jcrsl, ima, 1, &
                            1, icmp, iad2)
                valref = crsv(iad2)
            end if
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        valeur = cesv(iad)
                        if (isnan(valeur)) goto 20
!
                        if (iref .eq. 1) then
                            valeur = abs(valeur-valref)
                        end if
                        if (valeur .gt. valmax) then
                            imamax = ima
                            valmax = valeur
                        end if
                        if (valeur .lt. valmin) then
                            imamin = ima
                            valmin = valeur
                        end if
                    end if
20                  continue
                end do
            end do
40          continue
        end do
        if (imamax .gt. 0) then
            valk(2) = cvrcnom(icmp)
            valk(1) = cvrcvarc(icmp)
            valr(1) = valmax
            valr(2) = valmin
            valk(3) = int_to_char8(imamax)
            valk(4) = int_to_char8(imamin)
            if (iref .eq. 1) then
                valk(5) = valk(1)
                call utmess('A+', 'MECANONLINE2_95', nk=5, valk=valk, nr=2, &
                            valr=valr)
            else
                call utmess('A+', 'MECANONLINE2_94', nk=4, valk=valk, nr=2, &
                            valr=valr)
            end if
        end if
    end do
    call utmess('A', 'MECANONLINE2_93')
!
    call jedetr(chscom)
    call jedetr(chsref)
    call jedema()
end subroutine
