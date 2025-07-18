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

subroutine ssdmrm(mag)
    implicit none
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/ssdmu1.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: mag
! ----------------------------------------------------------------------
!     BUT:
!        - TRAITER LE MOT CLEF "RECO_SUPER_MAILLE"
!          DE LA COMMANDE DEFI_MAILLAGE.
!
!     IN:
!        MAG : NOM DU MAILLAGE QUE L'ON DEFINIT.
!     VAR:
!     -- MODIFICATION DE L'OBJET .NOEUD_CONF CREE DANS SSDMRC
!
! ----------------------------------------------------------------------
! INSPI SSDMRM  SSDMRG
    character(len=8) :: crit, nomacr, mal, nosma
    real(kind=8) :: prec, di, dj
    character(len=16) :: option
    character(len=24) :: valk(2), nognoi, nognoj
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacoo2, iagno, ialiii, ialiij
    integer(kind=8) :: ibid(1), ico, iconf, ii, inoi
    integer(kind=8) :: inoii, inoj, inojj, iocc, ismai, ismaj, j
    integer(kind=8) :: jj, kk, longi, longj, n1, n2, n3
    integer(kind=8) :: nbexti, nbextj, nbid, nbngno, nbnore, nbnori, nbnorj
    integer(kind=8) :: nbsma, nbsmar, nocc
    character(len=24), pointer :: likg(:) => null()
    character(len=8), pointer :: likm(:) => null()
    integer(kind=8), pointer :: noeud_conf(:) => null()
    real(kind=8), pointer :: para_r(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: dime_2(:) => null()
    integer(kind=8), pointer :: lini(:) => null()
    integer(kind=8), pointer :: linj(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getfac('RECO_SUPER_MAILLE', nocc)
    if (nocc .eq. 0) goto 999
!
!     -- ON RECUPERE CERTAINES DIMENSIONS:
!     ------------------------------------
    call jeveuo(mag//'.DIME', 'L', vi=dime)
    nbsma = dime(4)
!
    call jeveuo(mag//'.NOMACR', 'L', vk8=vnomacr)
    call jeveuo(mag//'.NOEUD_CONF', 'E', vi=noeud_conf)
!
    call jeveuo(mag//'.COORDO_2', 'L', iacoo2)
    call jeveuo(mag//'.DIME_2', 'L', vi=dime_2)
    call jeveuo(mag//'.PARA_R', 'L', vr=para_r)
    AS_ALLOCATE(vk8=likm, size=nbsma)
    AS_ALLOCATE(vk24=likg, size=nbsma)
!
!
!     -- BOUCLE SUR LES OCCURENCES DU MOT-CLEF:
!     -----------------------------------------
    longi = 0
    longj = 0
    do iocc = 1, nocc
!
!     -- ON RECUPERE LA LISTE DES MAILLES ET LA LISTE DES GROUP_NO:
!     -------------------------------------------------------------
        call getvtx('RECO_SUPER_MAILLE', 'SUPER_MAILLE', iocc=iocc, nbval=nbsma, &
                    vect=likm, nbret=n1)
        call getvtx('RECO_SUPER_MAILLE', 'GROUP_NO', iocc=iocc, nbval=nbsma, vect=likg, &
                    nbret=n2)
        if (n1 .lt. 0) then
            call utmess('F', 'SOUSTRUC_64')
        end if
        if (n1 .ne. n2) then
            call utmess('F', 'SOUSTRUC_65')
        end if
        if (n1 .lt. 2) then
            call utmess('F', 'SOUSTRUC_66')
        end if
!
        nbsmar = n1
!
        call getvtx('RECO_SUPER_MAILLE', 'OPTION', iocc=iocc, scal=option, nbret=n1)
        if (option(1:11) .eq. 'GEOMETRIQUE') then
            call getvr8('RECO_SUPER_MAILLE', 'PRECISION', iocc=iocc, scal=prec, nbret=n1)
            call getvtx('RECO_SUPER_MAILLE', 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
        end if
!
        do i = 1, nbsmar
            nosma = likm(i)
            nognoi = likg(i)
            call jenonu(jexnom(mag//'.SUPMAIL', nosma), ismai)
            di = para_r(14*(ismai-1)+13)
            nomacr = vnomacr(ismai)
            call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
            call jeveuo(nomacr//'.LINO', 'L', vi=lini)
            call jelira(nomacr//'.LINO', 'LONUTI', nbexti)
            call jeveuo(jexnom(mal//'.GROUPENO', nognoi), 'L', iagno)
            call jelira(jexnom(mal//'.GROUPENO', nognoi), 'LONUTI', nbngno)
            call utlisi('INTER', zi(iagno), nbngno, lini, nbexti, &
                        ibid, 0, n3)
            nbnori = -n3
            if (nbnori .gt. longi) then
                if (longi .ne. 0) call jedetr('&&SSDMRM.LIII')
!           POUR NE PAS LE DETRUIRE A CHAQUE FOIS, ON ALLOUE PLUS GRAND
                longi = (nbnori+10)*2
                call wkvect('&&SSDMRM.LIII', 'V V I', longi, ialiii)
            end if
            call utlisi('INTER', zi(iagno), nbngno, lini, nbexti, &
                        zi(ialiii), nbnori, nbid)
            do j = i+1, nbsmar
                nosma = likm(j)
                nognoj = likg(j)
                call jenonu(jexnom(mag//'.SUPMAIL', nosma), ismaj)
                dj = para_r(14*(ismaj-1)+13)
                nomacr = vnomacr(ismaj)
                call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
                call jeveuo(nomacr//'.LINO', 'L', vi=linj)
                call jelira(nomacr//'.LINO', 'LONUTI', nbextj)
                call jeveuo(jexnom(mal//'.GROUPENO', nognoj), 'L', iagno)
                call jelira(jexnom(mal//'.GROUPENO', nognoj), 'LONUTI', nbngno)
                call utlisi('INTER', zi(iagno), nbngno, linj, nbextj, &
                            ibid, 0, n3)
                nbnorj = -n3
                if (nbnorj .gt. longj) then
                    if (longj .ne. 0) call jedetr('&&SSDMRM.LIIJ')
                    longj = (nbnorj+10)*2
                    call wkvect('&&SSDMRM.LIIJ', 'V V I', longj, ialiij)
                end if
                call utlisi('INTER', zi(iagno), nbngno, linj, nbextj, &
                            zi(ialiij), nbnorj, nbid)
!
                dj = min(di, dj)
!
                if (nbnori .ne. nbnorj) then
                    valk(1) = nognoi
                    valk(2) = nognoj
                    call utmess('A', 'SOUSTRUC_67', nk=2, valk=valk)
                end if
                nbnore = min(nbnori, nbnorj)
!
!
                if (option(1:13) .eq. 'NOEUD_A_NOEUD') then
!           ---------------------------------
                    do ii = 1, nbnore
                        inoii = indiis(lini, zi(ialiii-1+ii), 1, &
                                       nbexti)
                        inojj = indiis(linj, zi(ialiij-1+ii), 1, &
                                       nbextj)
                        inoi = dime_2(4*(ismai-1)+3)+inoii
                        inoj = dime_2(4*(ismaj-1)+3)+inojj
                        if (inoi .lt. inoj) then
                            noeud_conf(inoj) = inoi
                        else
                            noeud_conf(inoi) = inoj
                        end if
                    end do
!
!
                else if (option(1:7) .eq. 'INVERSE') then
!           -----------------------------------
                    do ii = 1, nbnore
                        if (i .eq. 1) then
                            kk = nbnore+1-ii
                            inoii = indiis(lini, zi(ialiii-1+kk), 1, &
                                           nbexti)
                        else
                            inoii = indiis(lini, zi(ialiii-1+ii), 1, &
                                           nbexti)
                        end if
                        inojj = indiis(linj, zi(ialiij-1+ii), 1, &
                                       nbextj)
                        inoi = dime_2(4*(ismai-1)+3)+inoii
                        inoj = dime_2(4*(ismaj-1)+3)+inojj
                        if (inoi .lt. inoj) then
                            noeud_conf(inoj) = inoi
                        else
                            noeud_conf(inoi) = inoj
                        end if
                    end do
!
!
                else if (option(1:11) .eq. 'GEOMETRIQUE') then
!           --------------------------------------
                    ico = 0
                    do ii = 1, nbnore
                        inoii = indiis(lini, zi(ialiii-1+ii), 1, &
                                       nbexti)
                        inoi = dime_2(4*(ismai-1)+3)+inoii
                        do jj = 1, nbnore
                            inojj = indiis(linj, zi(ialiij-1+jj), 1, &
                                           nbextj)
                            inoj = dime_2(4*(ismaj-1)+3)+inojj
                            call ssdmu1(dj, crit, prec, zr(iacoo2+3*(inoi-1)), &
                                        zr(iacoo2+3*(inoj-1)), iconf)
                            if (iconf .eq. 0) then
                                if (inoi .lt. inoj) then
                                    noeud_conf(inoj) = inoi
                                else
                                    noeud_conf(inoi) = inoj
                                end if
                                ico = ico+1
                                goto 63
                            end if
                        end do
63                      continue
                    end do
!
                    if (nbnore .ne. ico) then
                        valk(1) = nognoi
                        valk(2) = nognoj
                        call utmess('A', 'SOUSTRUC_68', nk=2, valk=valk)
                    end if
!
                end if
!
            end do
        end do
!
    end do
!
!
999 continue
!
! --- MENAGE
!
    AS_DEALLOCATE(vk8=likm)
    AS_DEALLOCATE(vk24=likg)
    call jedetr('&&SSDMRM.LIII')
    call jedetr('&&SSDMRM.LIIJ')
!
    call jedema()
end subroutine
