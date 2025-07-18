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

subroutine regres(nomres, mailsk, result, pfchn2)
    implicit none
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/rsexch.h"
#include "asterfort/profchno_crsd.h"
#include "asterfort/nueq_chck.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, mailsk, result
!  BUT : < PROJECTION D'UNE RESTITUTION SUR UN SQUELETTE ENRICHI >
!
! SUITE A UNE RESTITUTION GLOBALE GENERALISEE ON PROJETE LES CHAMPS SUR
! UN SQUELETTE DONT ON A FUSIONNE LES NOEUDS D'INTERFACE DYNAMIQUE
!
! NOMRES  /I/ : NOM K8 DU CONCEPT MODE_MECA RESULTAT
! MAILSK  /I/ : NOM K8 DU MAILLAGE SQUELETTE SUPPORT
! RESULT  /O/ : NOM K8 DU MODE_MECA QUE L'ON VEUT PROJETER
!
!-----------------------------------------------------------------------
!
!
    character(len=19) :: pfchn1, pfchn2
    character(len=19) :: chexin, chexou, chamno
    integer(kind=8) :: i, iadnew, iadold, ieq, igd
    integer(kind=8) :: iold, iord, iret, j, k, ldeeq
    integer(kind=8) :: lnunew, lprnew, lprold
    integer(kind=8) :: lvnew, nbord, ncmp, nddl, ndeeq
    integer(kind=8) :: ndi, nec, nnodes
    character(len=24) :: nequ
    character(len=24), pointer :: refe(:) => null()
    integer(kind=8), pointer :: corres(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: p_nequ(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    call rsexch('F', result, 'DEPL', 1, chamno, &
                iret)
    call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=pfchn1)
    call nueq_chck(pfchn1)
!
!
    call copisd('NUME_EQUA', 'G', pfchn1, pfchn2)
    call copisd('RESULTAT', 'G', result, nomres)
!
!
    call jeveuo(mailsk//'.CORRES', 'L', vi=corres)
    call jelira(mailsk//'.CORRES', 'LONUTI', nnodes)
    call jeveuo(jexnum(pfchn1//'.PRNO', 1), 'L', lprold)
    call jeveuo(pfchn1//'.NUEQ', 'L', vi=nueq)
    call jelira(nomres//'           .ORDR', 'LONUTI', nbord)
    call jeveuo(nomres//'           .ORDR', 'L', vi=ordr)
!
    call dismoi('NUM_GD', chamno, 'CHAM_NO', repi=igd)
!
    nec = nbec(igd)
    ndi = nec+2
!
! --- CALCUL DU NOMBRE DE DDL ---
    nddl = 0
    do i = 1, nnodes
        iold = corres(i)
        nddl = nddl+zi(lprold+(iold-1)*ndi+1)
    end do
!
! - Create nume_equa
!
    call profchno_crsd(pfchn2, 'G', nddl, prno_lengthz=nnodes*ndi)
    call jeveuo(pfchn2//'.DEEQ', 'E', ldeeq)
    call jeveuo(pfchn2//'.NUEQ', 'E', lnunew)
    call jeveuo(jexnum(pfchn2//'.PRNO', 1), 'E', lprnew)

!
! - THis is a NUME_EQUA object
!
    nequ = pfchn2(1:19)//'.NEQU'
    call jeveuo(nequ, 'E', vi=p_nequ)
    p_nequ(1) = nddl

!
    do iord = 1, nbord
        call rsexch('F', result, 'DEPL', ordr(iord), chexin, &
                    iret)
        call rsexch('F', nomres, 'DEPL', ordr(iord), chexou, &
                    iret)
!
!     --- MISE A JOUR DU .REFE
        call jeveuo(chexou//'.REFE', 'E', vk24=refe)
        refe(1) = mailsk
        call detrsd('NUME_EQUA', refe(2))
        refe(2) = nomres//'.PROFC.NUME'
!
        ieq = 1
        do i = 1, nnodes
            iold = corres(i)
            ncmp = zi(lprold+(iold-1)*ndi+1)
            zi(lprnew+(i-1)*ndi) = ieq
            zi(lprnew+(i-1)*ndi+1) = ncmp
            zi(lprnew+(i-1)*ndi+2) = zi(lprold+(iold-1)*ndi+2)
            do k = 1, ncmp
                zi(lnunew-1+ieq) = ieq
                ieq = ieq+1
            end do
        end do
!
!     --- MISE A JOUR DU .VALE (DEPL_R) ---
        call jeexin(chexou//'.VALE', iret)
        if (iret .ne. 0) then
            call jedetr(chexou//'.VALE')
            call wkvect(chexou//'.VALE', 'G V R', nddl, lvnew)
            call jeveuo(chexin//'.VALE', 'L', vr=vale)
            do i = 1, nnodes
                iold = corres(i)
                ncmp = zi(lprold+(iold-1)*ndi+1)
                iadold = nueq(zi(lprold+(iold-1)*ndi))
                iadnew = zi(lnunew-1+zi(lprnew+(i-1)*ndi))
                do j = 1, ncmp
                    zr(lvnew-1+iadnew+j-1) = vale(iadold+j-1)
                end do
            end do
        end if
    end do
!

    ndeeq = 0
    do i = 1, nnodes
        ncmp = zi(lprnew-1+(i-1)*ndi+2)
        do j = 1, ncmp
            ndeeq = ndeeq+1
            zi(ldeeq-1+ndeeq) = i
            ndeeq = ndeeq+1
            zi(ldeeq-1+ndeeq) = j
        end do
    end do
!
    call jedema()
end subroutine
