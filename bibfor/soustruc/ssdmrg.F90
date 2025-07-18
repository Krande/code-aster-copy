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
subroutine ssdmrg(mag)
    implicit none
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ssdmu1.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: mag
! ----------------------------------------------------------------------
!     BUT:
!        - TRAITER LE MOTS CLEF "RECO_GLOBAL"
!          DE LA COMMANDE DEFI_MAILLAGE.
!
!     IN:
!        MAG : NOM DU MAILLAGE QUE L'ON DEFINIT.
!     VAR:
!        --MODIFICATION DE L'OBJET .NOEUD_CONF CREE DANS SSDMRC
!
    character(len=8) :: kbid, crit
    real(kind=8) :: prec, di, dj
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacoo2
    integer(kind=8) :: iasupi, iasupj, iconf, ii, inoi, inoj
    integer(kind=8) :: iocc, isma, j, jj, jsma, n1, nbnoi
    integer(kind=8) :: nbnoj, nbsma, nbsmar, nnnoe, nocc
    integer(kind=8), pointer :: liis(:) => null()
    character(len=8), pointer :: lik8(:) => null()
    real(kind=8), pointer :: para_r(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: dime_2(:) => null()
    integer(kind=8), pointer :: noeud_conf(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getfac('RECO_GLOBAL', nocc)
    if (nocc .eq. 0) goto 999
!
!     -- ON RECUPERE CERTAINES DIMENSIONS:
!     ------------------------------------
    call jeveuo(mag//'.DIME', 'L', vi=dime)
    nbsma = dime(4)
    nnnoe = dime(1)
!
    call jeveuo(mag//'.NOEUD_CONF', 'E', vi=noeud_conf)
!
    call jeveuo(mag//'.COORDO_2', 'L', iacoo2)
    call jeveuo(mag//'.DIME_2', 'L', vi=dime_2)
    call jeveuo(mag//'.PARA_R', 'L', vr=para_r)
    AS_ALLOCATE(vk8=lik8, size=nbsma)
    AS_ALLOCATE(vi=liis, size=nbsma)
!
!
!     -- BOUCLE SUR LES OCCURENCES DU MOT-CLEF:
!     -----------------------------------------
    do iocc = 1, nocc
!
!     -- ON RECUPERE LA LISTE DES MAILLES A TRAITER :
!     -----------------------------------------------
        call getvtx('RECO_GLOBAL', 'TOUT', iocc=iocc, scal=kbid, nbret=n1)
        if (n1 .eq. 1) then
            nbsmar = nbsma
            do i = 1, nbsmar
                liis(i) = i
            end do
        else
            call getvem(mag, 'MAILLE', 'RECO_GLOBAL', 'SUPER_MAILLE', iocc, &
                        nbsma, lik8, n1)
            if (n1 .lt. 0) then
                call utmess('F', 'SOUSTRUC_63')
            end if
            nbsmar = n1
            do i = 1, nbsmar
                call jenonu(jexnom(mag//'.SUPMAIL', lik8(i)), isma)
                liis(i) = isma
            end do
        end if
!
        call getvr8('RECO_GLOBAL', 'PRECISION', iocc=iocc, scal=prec, nbret=n1)
        call getvtx('RECO_GLOBAL', 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
!
        do i = 1, nbsmar
            isma = liis(i)
            call jeveuo(jexnum(mag//'.SUPMAIL', isma), 'L', iasupi)
            nbnoi = dime_2(4*(isma-1)+1)+dime_2(4*(isma-1)+2)
            di = para_r(14*(isma-1)+13)
            do j = i+1, nbsmar
                jsma = liis(j)
                call jeveuo(jexnum(mag//'.SUPMAIL', jsma), 'L', iasupj)
                nbnoj = dime_2(4*(jsma-1)+1)+dime_2(4*(jsma-1) &
                                                    +2)
                dj = para_r(14*(jsma-1)+13)
                dj = min(di, dj)
                do ii = 1, nbnoi
                    inoi = zi(iasupi-1+ii)
!               -- SI C'EST UN NOEUD DE LAGRANGE, ON SAUTE :
                    if (inoi .gt. nnnoe) goto 7
                    do jj = 1, nbnoj
                        inoj = zi(iasupj-1+jj)
                        if (inoj .gt. nnnoe) goto 8
                        call ssdmu1(dj, crit, prec, zr(iacoo2+3*(inoi-1)), zr(iacoo2+3*(inoj-1)), &
                                    iconf)
                        if (iconf .eq. 0) then
                            if (inoi .lt. inoj) then
                                noeud_conf(inoj) = inoi
                            else
                                noeud_conf(inoi) = inoj
                            end if
                        end if
8                       continue
                    end do
7                   continue
                end do
            end do
        end do
!
    end do
!
!
999 continue
! --- MENAGE
    AS_DEALLOCATE(vk8=lik8)
    AS_DEALLOCATE(vi=liis)
!
    call jedema()
end subroutine
