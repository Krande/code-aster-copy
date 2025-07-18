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
subroutine xstam1(noma, nbma, nmafis, mafis, stano, &
                  mafon, maen1, maen2, maen3, nmafon, &
                  nmaen1, nmaen2, nmaen3, typdis, cnslt)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/panbno.h"
!
    integer(kind=8) :: nmafis, nmafon, nmaen1, nmaen2, nmaen3, nbma
    integer(kind=8) :: stano(*), mafis(nmafis)
    integer(kind=8) :: mafon(nmafis), maen1(nbma), maen2(nbma), maen3(nbma)
    character(len=8) :: noma
    character(len=16) :: typdis
    character(len=19) :: cnslt
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
! CALCUL DU STATUT DES MAILLES SUIVANT LE STATUT DES NOEUDS EN COURS
!
! ----------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NBMA   : NOMBRE DE MAILLES DU MAILLAGE
! IN  NMAFIS : NOMBRE DE MAILLES DE LA ZONE FISSURE
! IN  MAFIS  : VECTUER DES MAILLES DE LA ZONE FISSURE
! IN  STANO  : VECTEUR STATUT DES NOEUDS
!
! OUT  NMAFON : NOMBRE DE MAILLES CONTENANT LE FOND DE FISSURE
! OUT  NMAEN1 : NOMBRE DE MAILLES 'HEAVISIDE'
! OUT  NMAEN2 : NOMBRE DE MAILLES 'CRACKTIP'
! OUT  NMAEN3 : NOMBRE DE MAILLES 'HEAVISIDE-CRACKTIP'
! OUT  MAFON  : VECTEUR DES MAILLES 'CONTENANT LE FOND DE FISSURE
! OUT  MAEN1  : VECTEUR DES MAILLES 'HEAVISIDE'
! OUT  MAEN2  : VECTEUR DES MAILLES 'CRACKTIP'
! OUT  MAEN3  : VECTEUR DES MAILLES 'HEAVISIDE-CRACKTIP'
!
!
!
!
    integer(kind=8) :: jma, nuno, jconx2, jltsv
    integer(kind=8) :: i, im1, im2, im3, ima, itypma, in, imae, nunop
    integer(kind=8) :: em, em1, em2, nmaabs, nbnott(3), nno, en
    integer(kind=8) :: ndim, dime_topo
    character(len=8) :: typma
    character(len=19) :: mai
    aster_logical :: lstch
    integer(kind=8), pointer :: connex(:) => null()
!
    call jemarq()
!
    i = 0
    im1 = 0
    im2 = 0
    im3 = 0
!
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    if (cnslt(3:8) .eq. 'OP0010') call jeveuo(cnslt//'.CNSV', 'L', jltsv)
!
!   dimension du maillage
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
!
!     BOUCLE SUR LES MAILLES DU MAILLAGE
    do ima = 1, nbma
!
        nmaabs = ima
!
        itypma = zi(jma-1+nmaabs)
        call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
        if (typma(1:3) .eq. 'POI') goto 310
!
        em = 0
        em1 = 0
        em2 = 0
        call panbno(itypma, nbnott)
        nno = nbnott(1)+nbnott(2)+nbnott(3)
!
!       BOUCLE SUR LES NOEUDS DE LA MAILLE
        lstch = .false.
        do in = 1, nno
            if (in .gt. 1) nunop = nuno
            nuno = connex(zi(jconx2+nmaabs-1)+in-1)
            en = abs(stano(nuno))
            if (en .eq. 1 .or. en .eq. 3) em1 = em1+1
            if (en .eq. 2 .or. en .eq. 3) em2 = em2+1
            if (in .gt. 1 .and. cnslt(3:8) .eq. 'OP0010') lstch = lstch .or. &
                                                                  ( &
                                                                  ( &
                                                                  zr(jltsv-1+nuno)*zr(jltsv-1+nun&
                                                                  &op) &
                                                                  ) &
                                                                  .le. 0.d0 &
                                                                  )
        end do
        if (em1 .ge. 1) em = 1
        if (em2 .ge. 1) em = 2
        if (em1 .ge. 1 .and. em2 .ge. 1) em = 3
!
!         MAILLE RETENUE POUR MAFOND (TS LS NOEUDS SONT 'CARRÉS')
!         SOUS RÉSERVE QUE CE SOIT UNE MAILLE PRINCIPALE DE MAFIS
        if (typdis .ne. 'COHESIF') then
            if (em2 .eq. nno) then
!
                call dismoi('DIM_TOPO', typma, 'TYPE_MAILLE', repi=dime_topo)
                if (dime_topo .eq. ndim) then
!
                    do imae = 1, nmafis
                        if (nmaabs .eq. mafis(imae)) then
                            i = i+1
                            ASSERT(i .le. nmafis)
                            mafon(i) = nmaabs
!                       ON SORT DE LA BOUCLE 312
                            goto 313
                        end if
                    end do
313                 continue
                end if
            end if
        else if (typdis .eq. 'COHESIF') then
!
            if (lstch .or. (em1 .lt. nno .and. em1 .gt. 0) .and. typma(1:3) .ne. 'SEG') then
                do imae = 1, nmafis
                    if (nmaabs .eq. mafis(imae)) then
                        i = i+1
                        ASSERT(i .le. nmafis)
                        mafon(i) = nmaabs
                        goto 413
                    end if
                end do
413             continue
            end if
        end if
!
!       ON RÉCUPÈRE LES NUMEROS DES MAILLES ENRICHIES
        if (em .eq. 1) then
            im1 = im1+1
            ASSERT(im1 .le. nbma)
            maen1(im1) = nmaabs
        else if (em .eq. 2) then
            im2 = im2+1
            ASSERT(im2 .le. nbma)
            maen2(im2) = nmaabs
        else if (em .eq. 3) then
            im3 = im3+1
            ASSERT(im3 .le. nbma)
            maen3(im3) = nmaabs
        end if
310     continue
    end do
!
    nmafon = i
    nmaen1 = im1
    nmaen2 = im2
    nmaen3 = im3
!
    call jedema()
end subroutine
