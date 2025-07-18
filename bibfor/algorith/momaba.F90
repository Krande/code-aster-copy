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

subroutine momaba(mailla)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/barhex.h"
#include "asterfort/barpen.h"
#include "asterfort/barpyr.h"
#include "asterfort/barqua.h"
#include "asterfort/bartet.h"
#include "asterfort/bartri.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: mailla
!
    integer(kind=8) :: jtyma, nbmc, nbma, jnuma, i, ityp, n1, n2, i1, nbno
    integer(kind=8) :: nbmat, jpoin, ifm, niv, jcon, ndim, nn, jnbma, nbnmf
    integer(kind=8) :: ilcnx2, nbm1, kk, numa, i1sauv, i2sauv, nbnsf
    aster_logical :: lnmf, lmodi, lconx
    parameter(nbmc=2)
    character(len=8) :: k8b, type
    character(len=16) :: tymocl(nbmc), motcle(nbmc)
    character(len=24) :: connex, nomjv
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: coninv(:) => null()
    real(kind=8), pointer :: conm(:) => null()
    real(kind=8), pointer :: coor(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
    connex = mailla//'.CONNEX'
    call jeveuo(mailla//'.TYPMAIL        ', 'L', jtyma)
    call jeveuo(mailla//'.COORDO    .VALE', 'E', vr=coor)
    call jeveuo(mailla//'.DIME           ', 'L', vi=dime)
    call dismoi('NB_MA_MAILLA', mailla, 'MAILLAGE', repi=nbmat)
    lmodi = .false.
    lconx = .false.
!     ------------------------------------------------------------------
!
! --- LECTURE DE LA LISTE DE MAILLES
!
    motcle(1) = 'GROUP_MA_FOND'
    tymocl(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE_FOND'
    tymocl(2) = 'MAILLE'
    nomjv = '&&MOMABA.LISTE_MA'
    call reliem(' ', mailla, 'NU_MAILLE', 'MODI_MAILLE', 1, &
                nbmc, motcle, tymocl, nomjv, nbma)
!
    if (nbma .eq. 0) goto 150
!
    call jeveuo(nomjv, 'L', jnuma)
!
! --- TRAITEMENT DES MAILLES
!
!     ON INTERDIT 'GROUP_MA_FOND' ET 'MAILLE_FOND'
!     SI LE MAILLAGE EST DE DIMENSION 2.
    if (dime(6) .eq. 2) then
        call utmess('F', 'ALGORITH6_16')
    end if
!
! --- CREATION DE LA CONNECTIVITE INVERSE
!
    call cncinv(mailla, [0], 0, 'V', '&&MOMABA.CONINV')
    call jeveuo('&&MOMABA.CONINV', 'L', vi=coninv)
    call jeveuo(jexatr('&&MOMABA.CONINV', 'LONCUM'), 'L', ilcnx2)
    lconx = .true.
!
    do i = 1, nbma
        ityp = jtyma-1+zi(jnuma+i-1)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
!
        if (niv .eq. 2) then
            k8b = int_to_char8(zi(jnuma+i-1))
            write (ifm, *) 'TRAITEMENT DE LA MAILLE ', k8b
        end if
        if (type(1:4) .eq. 'SEG3') then
            call jeveuo(jexnum(connex, zi(jnuma+i-1)), 'L', jpoin)
            n1 = zi(jpoin)
            n2 = zi(jpoin+1)
        else if (type(1:4) .eq. 'POI1') then
            call jeveuo(jexnum(connex, zi(jnuma+i-1)), 'L', jpoin)
            n1 = zi(jpoin)
            n2 = 0
        else
            call utmess('F', 'ALGORITH6_17')
        end if
!
        nbm1 = zi(ilcnx2+n1)-zi(ilcnx2-1+n1)
!
! ----- BOUCLE SUR LES MAILLES AUXQUELLES SONT LIEES CE NOEUD
        do kk = 1, nbm1
            numa = coninv(1+zi(ilcnx2-1+n1)-1+kk-1)
            call jeveuo(jexnum(connex, numa), 'L', jpoin)
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtyma-1+numa)), type)
!
            nbno = 0
! --------- TRIA6, TRIA7
            if (type .eq. 'TRIA6' .or. type .eq. 'TRIA7') then
                nbno = 3
! --------- QUAD8, QUAD9
            else if (type .eq. 'QUAD8' .or. type .eq. 'QUAD9') then
                nbno = 4
! --------- TETRA10
            else if (type .eq. 'TETRA10') then
                nbno = 4
! --------- PENTA15, PENTA18
            else if (type .eq. 'PENTA15' .or. type .eq. 'PENTA18') then
                nbno = 6
! --------- PYRAM13
            else if (type .eq. 'PYRAM13') then
                nbno = 5
! --------- HEXA20, HEXA27
            else if (type .eq. 'HEXA20' .or. type .eq. 'HEXA27') then
                nbno = 8
            end if
!
            i1sauv = 0
            i2sauv = 0
            do i1 = 1, nbno
                if (zi(jpoin+i1-1) .eq. n1) i1sauv = i1
                if (zi(jpoin+i1-1) .eq. n2) i2sauv = i1
            end do
!
            if ((i1sauv .ne. 0 .and. i2sauv .ne. 0) .or. &
                (i1sauv .ne. 0 .and. n2 .eq. 0)) then
! ------------- TRIA6, TRIA7
                if (type .eq. 'TRIA6' .or. type .eq. 'TRIA7') then
                    lmodi = .true.
                    call bartri(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- QUAD8, QUAD9
                else if (type .eq. 'QUAD8' .or. type .eq. 'QUAD9') then
                    lmodi = .true.
                    call barqua(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- TETRA10
                else if (type .eq. 'TETRA10') then
                    lmodi = .true.
                    call bartet(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- PENTA15, PENTA18
                else if (type .eq. 'PENTA15' .or. type .eq. 'PENTA18') then
                    lmodi = .true.
                    call barpen(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- PYRAM13
                else if (type .eq. 'PYRAM13') then
                    lmodi = .true.
                    call barpyr(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- HEXA20, HEXA27
                else if (type .eq. 'HEXA20' .or. type .eq. 'HEXA27') then
                    lmodi = .true.
                    call barhex(i1sauv, i2sauv, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
                end if
            end if
        end do
!
    end do
!
    call jedetr(nomjv)
    call jedetr('&&MOMABA_MAILLE')
!
150 continue
!     ------------------------------------------------------------------
!
! --- LECTURE DE LA LISTE DE NOEUDS
!
    call jeveuo(mailla//'.COORDO    .VALE', 'L', vr=conm)
    call jelira(mailla//'.COORDO    .VALE', 'LONMAX', ndim)
!
!     ON STOCKE LES COORDONNEES DES NOEUDS DU FOND DE FISSURE AVANT
!     LEURS MODIFICATIONS
    call wkvect('&&COORD_NOEUDS', 'V V R', ndim, jcon)
    do i = 1, ndim
        zr(jcon+i-1) = conm(i)
    end do
!
    motcle(1) = 'GROUP_NO_FOND'
    tymocl(1) = 'GROUP_NO'
    motcle(2) = 'NOEUD_FOND'
    tymocl(2) = 'NOEUD'
    nomjv = '&&MOMABA.LISTE_NO'
    call reliem(' ', mailla, 'NU_NOEUD', 'MODI_MAILLE', 1, &
                nbmc, motcle, tymocl, nomjv, nbma)
    if (nbma .eq. 0) goto 260
!
    if (.not. lconx) then
!
! ------- CREATION DE LA CONNECTIVITE INVERSE
!
        call cncinv(mailla, [0], 0, 'V', '&&MOMABA.CONINV')
        call jeveuo('&&MOMABA.CONINV', 'L', vi=coninv)
        call jeveuo(jexatr('&&MOMABA.CONINV', 'LONCUM'), 'L', ilcnx2)
        lconx = .true.
    end if
!
!     ON VERIFIE L'UNICITE DU NOEUD DU FOND DE FISSURE POUR UN
!     MAILLAGE DE DIMENSION 2
    if (dime(6) .eq. 2 .and. nbma .gt. 1) then
        call utmess('F', 'ALGORITH6_18')
    end if
!
!     ON VERIFIE QU'IL Y A AU MOINS 2 NOEUDS POUR UN
!     MAILLAGE DE DIMENSION 3
    if (dime(6) .eq. 3 .and. nbma .lt. 2) then
        call utmess('F', 'ALGORITH16_74')
    end if
!
    call jeveuo(nomjv, 'L', jnuma)
    call wkvect('&&NOEU_MIL_FISS', 'V V I', nbma, jnbma)
!
! --- TRAITEMENT DES NOEUDS
!
!   nombre de noeuds sommet du fond de fissure
    nbnsf = 0
!   nombre de noeuds milieu du fond de fissure
    nbnmf = 0
    do i = 1, nbma
        n1 = zi(jnuma+i-1)
        n2 = 0
        if (niv .eq. 2) then
            k8b = int_to_char8(zi(jnuma+i-1))
            write (ifm, *) 'TRAITEMENT DU NOEUD ', k8b
        end if
        lnmf = .true.
!
        nbm1 = zi(ilcnx2+n1)-zi(ilcnx2-1+n1)
!
! ----- BOUCLE SUR LES MAILLES AUXQUELLES SONT LIEES CE NOEUD
        do kk = 1, nbm1
            numa = coninv(1+zi(ilcnx2-1+n1)-1+kk-1)
            call jeveuo(jexnum(connex, numa), 'L', jpoin)
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtyma-1+numa)), type)
!
            nbno = 0
! --------- TRIA6, TRIA7
            if (type .eq. 'TRIA6' .or. type .eq. 'TRIA7') then
                nbno = 3
! --------- QUAD8, QUAD9
            else if (type .eq. 'QUAD8' .or. type .eq. 'QUAD9') then
                nbno = 4
! --------- TETRA10
            else if (type .eq. 'TETRA10') then
                nbno = 4
! --------- PENTA15, PENTA18
            else if (type .eq. 'PENTA15' .or. type .eq. 'PENTA18') then
                nbno = 6
! --------- PYRAM13
            else if (type .eq. 'PYRAM13') then
                nbno = 5
! --------- HEXA20, HEXA27
            else if (type .eq. 'HEXA20' .or. type .eq. 'HEXA27') then
                nbno = 8
            end if
!
            i1sauv = 0
            do i1 = 1, nbno
                if (zi(jpoin+i1-1) .eq. n1) i1sauv = i1
            end do
!
            if (i1sauv .ne. 0) then
! ------------- TRIA6, TRIA7
                if (type .eq. 'TRIA6' .or. type .eq. 'TRIA7') then
                    lnmf = .false.
                    lmodi = .true.
                    call bartri(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- QUAD8, QUAD9
                else if (type .eq. 'QUAD8' .or. type .eq. 'QUAD9') then
                    lnmf = .false.
                    lmodi = .true.
                    call barqua(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- TETRA10
                else if (type .eq. 'TETRA10') then
                    lnmf = .false.
                    lmodi = .true.
                    call bartet(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- PENTA15
                else if (type .eq. 'PENTA15') then
                    lnmf = .false.
                    lmodi = .true.
                    call barpen(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- PYRAM13
                else if (type .eq. 'PYRAM13') then
                    lnmf = .false.
                    lmodi = .true.
                    call barpyr(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
! ------------- HEXA20, HEXA27
                else if (type .eq. 'HEXA20' .or. type .eq. 'HEXA27') then
                    lnmf = .false.
                    lmodi = .true.
                    call barhex(i1sauv, n2, coor, zi(jpoin))
                    if (niv .eq. 2) then
                        k8b = int_to_char8(numa)
                        write (ifm, *) '   MAILLE MODIFIEE ', k8b
                    end if
                end if
            end if
        end do
!
        if (lnmf) then
!         ON STOCKE LES NOEUDS MILIEU DU FOND DE FISSURE
            zi(jnbma+nbnmf) = n1
            nbnmf = nbnmf+1
        else
            nbnsf = nbnsf+1
        end if
!
    end do
!
    do i = 1, nbnmf
!       ON REAJUSTE LES COORDONNEES DES NOEUDS MILIEU
!       DU FOND DE FISSURE
        nn = 3*(zi(jnbma+i-1)-1)
        coor(nn+1) = zr(jcon+nn)
        coor(1+nn+1) = zr(jcon+nn+1)
        coor(1+nn+2) = zr(jcon+nn+2)
    end do
!
    if (.not. lmodi) then
        call utmess('F', 'ALGORITH16_72')
    elseif (nbma .gt. 1 .and. nbnmf .ne. nbnsf-1) then
        call utmess('F', 'ALGORITH16_73', si=nbnsf-1-nbnmf)
    end if
!
    call jedetr('&&NOEU_MIL_FISS')
    call jedetr('&&COORD_NOEUDS')
    call jedetr(nomjv)
    call jedetr('&&MOMABA_MAILLE')
    call jedetr('&&MOMABA.CONINV')
!
260 continue
!
    call jedema()
end subroutine
