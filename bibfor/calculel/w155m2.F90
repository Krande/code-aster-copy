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

subroutine w155m2(chin, carele, ligrel, chextr, nomsym, &
                  nocmp, tymaxi)
! person_in_charge: jacques.pellet at edf.fr
! ======================================================================
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liglma.h"
#include "asterfort/utmess.h"
#include "asterfort/w155m3.h"
    character(len=8) :: carele, nocmp, tymaxi
    character(len=16) :: nomsym
    character(len=19) :: chin, chextr, ligrel
!
! ----------------------------------------------------------------------
! BUT : EXTRACTION DU CHAM_ELEM CORRESPONDANT AU MIN/MAX
!       SUR LES SOUS-POINTS
!
! IN/JXIN  CHIN  : CHAM_ELEM (PLUSIEURS SOUS-POINTS) DANS LEQUEL
!                  ON DOIT EXTRAIRE CHEXTR
! IN/JXIN  CARELE  : CARA_ELEM ASSOCIE A CHIN
! IN/JXIN  LIGREL  : LIGREL SUR LEQUEL CREER CHEXTR
! IN/JXVAR CHEXTR  : CHAM_ELEM (1 SEUL SOUS-POINT) A REMPLIR
! IN       NOCMP   : NOM DE LA COMPOSANTE DONT ON CHERCHE LE MIN/MAX
! IN       TYMAXI  : /'MAXI' /'MINI' /'MAXI_ABS', /'MINI_ABS'
    character(len=24) :: linuma, linute, valk(5)
    character(len=19) :: ces1, ces2, ces3, ces4
    character(len=16) :: option
    character(len=8) :: licmp(4), ma, nomgd, tsca, nompar
    character(len=3) :: exituy
    integer(kind=8) :: iret, nbma, nbmat, numa, kma
    integer(kind=8) :: nbpt, kpt, kcmp, nncp, jlite
    integer(kind=8) :: iad1, iad4, jlima, ncmp
    integer(kind=8) :: jce2l, jce2d, jce2v, jce3k, jce3l, jce3d
    integer(kind=8) :: jce4l, jce4d, nbspmx, nucmp
    integer(kind=8) :: isp, nbsp, ksp, nusec, nucou, nufib, posic, posis, kcmp2
    real(kind=8) :: val, vmima
    character(len=8), pointer :: ce3c(:) => null()
    character(len=8), pointer :: ce4c(:) => null()
    real(kind=8), pointer :: ce3v(:) => null()
    real(kind=8), pointer :: ce4v(:) => null()
!
! ----------------------------------------------------------------------
    call jemarq()
    call dismoi('NOM_MAILLA', chin, 'CHAM_ELEM', repk=ma)
    call dismoi('NOM_GD', chin, 'CHAM_ELEM', repk=nomgd)
    call dismoi('TYPE_SCA', chin, 'CHAM_ELEM', repk=tsca)
    call dismoi('MXNBSP', chin, 'CHAM_ELEM', repi=nbspmx)
    call dismoi('EXI_TUYAU', ligrel, 'LIGREL', repk=exituy)
    if (nbspmx .le. 1) then
        call utmess('F', 'CALCULEL2_15')
    end if
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbmat)
    ASSERT(tsca .eq. 'R')
    ASSERT(exituy .eq. 'OUI' .or. exituy .eq. 'NON')
!
!
!     1.  LISTE DES MAILLES A TRAITER :
!     ---------------------------------
    linuma = '&&W155M2.LIMA'
    linute = '&&W155M2.LITE'
    call liglma(ligrel, nbma, linuma, linute)
    ASSERT(nbma .gt. 0)
    call jeveuo(linuma, 'L', jlima)
    call jeveuo(linute, 'L', jlite)
!
!
!     2.  NOMBRE DE COUCHES, SECTEURS ET FIBRES  DES ELEMENTS :
!     -----------------------------------------------------------
    ces1 = '&&W155M2.CES1'
    ces2 = '&&W155M2.CES2'
    call celces(carele//'.CANBSP', 'V', ces1)
!
!     -- L'ORDRE DES CMPS EST IPORTANT (UTILISE DANS W155M3)
    licmp(1) = 'COQ_NCOU'
    licmp(2) = 'TUY_NCOU'
    licmp(3) = 'TUY_NSEC'
    licmp(4) = 'NBFIBR'
    call cesred(ces1, nbma, zi(jlima), 4, licmp, &
                'V', ces2)
    call detrsd('CHAM_ELEM_S', ces1)
    call jeveuo(ces2//'.CESD', 'L', jce2d)
    call jeveuo(ces2//'.CESV', 'L', jce2v)
    call jeveuo(ces2//'.CESL', 'L', jce2l)
!
!
!     3. CHIN -> CES3 :
!     ------------------
    ces3 = '&&W155M2.CES3'
    call celces(chin, 'V', ces3)
    call jeveuo(ces3//'.CESK', 'L', jce3k)
    call jeveuo(ces3//'.CESD', 'L', jce3d)
    call jeveuo(ces3//'.CESC', 'L', vk8=ce3c)
    call jeveuo(ces3//'.CESL', 'L', jce3l)
    call jeveuo(ces3//'.CESV', 'L', vr=ce3v)
    call jelira(ces3//'.CESC', 'LONMAX', ncmp)
!
!
!     4. REMPLISSAGE DU CHAM_ELEM_S RESULTAT CES4 :
!     ---------------------------------------------
    ces4 = '&&W155M2.CES4'
    call celces(chextr, 'V', ces4)
    call jeveuo(ces4//'.CESD', 'L', jce4d)
    call jeveuo(ces4//'.CESV', 'E', vr=ce4v)
    call jeveuo(ces4//'.CESL', 'L', jce4l)
    call jeveuo(ces4//'.CESC', 'L', vk8=ce4c)
    do kcmp = 1, ncmp
        if (ce3c(kcmp) .eq. nocmp) then
            nucmp = kcmp
            goto 20
!
        end if
    end do
    valk(1) = nocmp
    valk(2) = nomsym
    call utmess('F', 'CHAMPS_19', nk=2, valk=valk)
!
20  continue
    do kma = 1, nbma
        numa = zi(jlima-1+kma)
        ASSERT(numa .ge. 1 .and. numa .le. nbmat)
!       si l'option MINMAX_SP n'existe pas pour la maille on boucle
        nbpt = zi(jce4d-1+5+4*(numa-1)+1)
!       pour les champs elga on a pas de points de gauss dans ce cas
        if (nbpt .le. 0) goto 60
!       pour les champs elno, nbpt n'est pas à zéro quand l'option
!       n'existe pas sur la maille, on fait une vérification en plus
        call cesexi('C', jce4d, jce4l, numa, 1, &
                    1, 1, iad4)
        if (iad4 .eq. 0) goto 60

        nbpt = zi(jce3d-1+5+4*(numa-1)+1)
        nbsp = zi(jce3d-1+5+4*(numa-1)+2)
        if (nbsp .le. 0) goto 60
!
        do kpt = 1, nbpt
!         -- 4.1 CALCUL DE VMIMA ET ISP :
!            VMIMA : VALEUR MIN/MAX ATTEINTE SUR LES SOUS-POINTS
!            ISP   : NUMERO DU SOUS-POINT REALISANT LE MIN/MAX
            do ksp = 1, nbsp
                call cesexi('C', jce3d, jce3l, numa, kpt, &
                            ksp, nucmp, iad1)
                if (iad1 .gt. 0) then
                    val = ce3v(iad1)
                    if (tymaxi(5:8) .eq. '_ABS') val = abs(val)
                    if (ksp .eq. 1) then
                        vmima = val
                        isp = ksp
                    else
                        if (tymaxi(1:4) .eq. 'MAXI') then
                            if (val .gt. vmima) then
                                vmima = val
                                isp = ksp
                            end if
                        else if (tymaxi(1:4) .eq. 'MINI') then
                            if (val .lt. vmima) then
                                vmima = val
                                isp = ksp
                            end if
                        else
                            ASSERT(.false.)
                        end if
                    end if
                end if
            end do
!
!         -- 4.2  CALCUL DE NUCOU, NUSEC, ... A PARTIR DE ISP :
            call w155m3(numa, jce2d, jce2l, jce2v, isp, &
                        nucou, nusec, nufib, posic, posis)
!
!         -- 4.3 STOCKAGE DE VMIMA, NUCOU, NUSEC, ...
            do kcmp2 = 1, 6
                call cesexi('C', jce4d, jce4l, numa, kpt, &
                            1, kcmp2, iad4)
                ASSERT(iad4 .gt. 0)
                if (kcmp2 .eq. 1) then
                    ASSERT(ce4c(kcmp2) .eq. 'VAL')
                    ce4v(iad4) = vmima
                else if (kcmp2 .eq. 2) then
                    ASSERT(ce4c(kcmp2) .eq. 'NUCOU')
                    ce4v(iad4) = dble(nucou)
                else if (kcmp2 .eq. 3) then
                    ASSERT(ce4c(kcmp2) .eq. 'NUSECT')
                    ce4v(iad4) = dble(nusec)
                else if (kcmp2 .eq. 4) then
                    ASSERT(ce4c(kcmp2) .eq. 'NUFIBR')
                    ce4v(iad4) = dble(nufib)
                else if (kcmp2 .eq. 5) then
                    ASSERT(ce4c(kcmp2) .eq. 'POSIC')
                    ce4v(iad4) = dble(posic)
                else if (kcmp2 .eq. 6) then
                    ASSERT(ce4c(kcmp2) .eq. 'POSIS')
                    ce4v(iad4) = dble(posis)
                else
                    ASSERT(.false.)
                end if
            end do
        end do
60      continue
    end do
!
!
!     5 CES4 -> CHEXTR :
!     ------------------------------------
    call dismoi('NOM_OPTION', chextr, 'CHAM_ELEM', repk=option)
    call dismoi('NOM_PARAM', chextr, 'CHAM_ELEM', repk=nompar)
    call detrsd('CHAM_ELEM', chextr)
    call cescel(ces4, ligrel, option, nompar, 'OUI', &
                nncp, 'G', chextr, 'F', iret)
    ASSERT(nncp .eq. 0)
!
!
!     6. MENAGE :
!     ------------
    call detrsd('CHAM_ELEM_S', ces2)
    call detrsd('CHAM_ELEM_S', ces3)
    call detrsd('CHAM_ELEM_S', ces4)
    call jedetr(linuma)
    call jedetr(linute)
!
    call jedema()
end subroutine
