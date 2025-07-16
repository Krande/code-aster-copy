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
subroutine ornofd(mafour, nomail, nbma, noeord, ndorig, &
                  ndextr, base, vecori)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/i2extf.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: nbma
    character(len=24) :: mafour
    character(len=8) :: nomail, ndorig, ndextr
    character(len=24) :: noeord
    character(len=1) :: base
    real(kind=8) :: vecori(3)
! FONCTION REALISEE:
!
!       ORNOFD -- ORDONNANCEMENT D'UNE LISTE DE NOEUD
!                 A PARTIR D'UN NOEUD ORIGINE
!                 UTILISE DANS DEFI_GROUP ET DEFI_FOND_FISS
!
!     ENTREES:
!        MAFOUR     : LISTE DES MAILLES SEG
!        NOMAIL     : NOM DU MAILLAGE
!        NBMA       : NOMBRE DE MAILLES TRAITEES
!        NOEORD     : NOM DE L'OBJET
!        NDORIG     : NOM DU NOEUD ORIGINE
!        NDEXTR     : NOM DU NOEUD EXTREMITE
!        BASE       : TYPE DE BASE DE SAUVEGARDE
!
!-----------------------------------------------------------------------
!
    real(kind=8) :: vecta(3), ps1, ps2
!
    integer(kind=8) :: iatyma, jtypm, jmail
    integer(kind=8) :: im, nid, nig, njonc, n, i, k, nbno
    integer(kind=8) :: jrdm, jnoe, ntemp
    character(len=8) :: typm
    character(len=8) :: noeud
    character(len=24) :: conec, typp
    integer(kind=8), pointer :: mailles_triee(:) => null()
    integer(kind=8), pointer :: noeuds_extrem(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
! DEB-------------------------------------------------------------------
    call jemarq()
!
    conec = nomail//'.CONNEX'
    typp = nomail//'.TYPMAIL'
!
!     RECUPERATION DES NOEUDS DESORDONNES
    call jeveuo(mafour, 'L', jmail)
!
!     ------------------------------------------------------------------
!     RECUPERATION DU TYPE DE MAILLE
!     ------------------------------------------------------------------
    call jeveuo(typp, 'L', iatyma)
    jtypm = iatyma-1+zi(jmail-1+1)
    call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypm)), typm)
!
!     ------------------------------------------------------------------
!     CONSTRUCTION D'UN VECTEUR DE TRAVAIL LOCAL CONTENANT
!     LES NOEUDS EXTREMITES  DE CHAQUE MAILLE
!     ------------------------------------------------------------------
    AS_ALLOCATE(vi=noeuds_extrem, size=2*nbma)
    do im = 1, nbma
        call i2extf(zi(jmail-1+im), 1, conec(1:15), typp(1:16), nig, &
                    nid)
        noeuds_extrem(im) = nig
        noeuds_extrem(nbma+im) = nid
    end do
!
!
!     ------------------------------------------------------------------
!     --- ORDONNANCEMENT DES MAILLES EN PARTANT DU NOEUD ORIGINE
!     ------------------------------------------------------------------
    njonc = char8_to_int(ndorig)
    n = 1
!     ------------------------------------------------------------------
!     CONSTRUCTION D'UN VECTEUR DE TRAVAIL LOCAL POUR
!     TRIER LES NOEUDS ET CONTENANT
!     LES MAILLES, LES NOEUDS SOMMET 1 ET LES NOEUDS SOMMET 2
!     ------------------------------------------------------------------
    AS_ALLOCATE(vi=mailles_triee, size=3*nbma)
!     EQUIVALENT D'UNE BOUCLE WHILE
550 continue
    do i = n, nbma
        if (noeuds_extrem(i) .eq. njonc) then
            mailles_triee(n) = zi(jmail-1+i)
            mailles_triee(nbma+n) = noeuds_extrem(i)
            mailles_triee(2*nbma+n) = noeuds_extrem(nbma+i)
            njonc = noeuds_extrem(nbma+i)
            goto 555
        end if
        if (noeuds_extrem(nbma+i) .eq. njonc) then
            mailles_triee(n) = zi(jmail-1+i)
            mailles_triee(nbma+n) = noeuds_extrem(nbma+i)
            mailles_triee(2*nbma+n) = noeuds_extrem(i)
            njonc = noeuds_extrem(i)
            goto 555
        end if
    end do
!
555 continue
    do k = n, i-1
        mailles_triee(1+k) = zi(jmail-1+k)
        mailles_triee(nbma+1+k) = noeuds_extrem(k)
        mailles_triee(2*nbma+1+k) = noeuds_extrem(nbma+k)
    end do
    do k = i+1, nbma
        mailles_triee(k) = zi(jmail-1+k)
        mailles_triee(nbma+k) = noeuds_extrem(k)
        mailles_triee(2*nbma+k) = noeuds_extrem(nbma+k)
    end do
    do k = n, nbma
        zi(jmail-1+k) = mailles_triee(k)
        noeuds_extrem(k) = mailles_triee(nbma+k)
        noeuds_extrem(nbma+k) = mailles_triee(2*nbma+k)
    end do
    n = n+1
    if (n .gt. nbma) goto 560
    goto 550
!
!
560 continue
!
!
!     ------------------------------------------------------------------
!     --- SAUVEGARDE DES NOEUDS ORDONNES DANS LA STRUCTURE DE DONNEES
!     --- AVEC RAJOUT DES NOEUDS MILIEUX SI SEG3
!     ------------------------------------------------------------------
    if (typm(1:4) .eq. 'SEG2') then
!
        nbno = nbma+1
        call wkvect(noeord, base//' V I', nbno, jnoe)
        do i = 1, nbma
            zi(jnoe-1+i) = noeuds_extrem(i)
        end do
        zi(jnoe-1+nbma+1) = noeuds_extrem(2*nbma)
!
    else if (typm(1:4) .eq. 'SEG3') then
!
        nbno = 2*nbma+1
        call wkvect(noeord, base//' V I', nbno, jnoe)
        do i = 1, nbma
            zi(jnoe-1+2*i-1) = noeuds_extrem(i)
            call jeveuo(jexnum(conec, zi(jmail-1+i)), 'L', jrdm)
            zi(jnoe-1+2*i) = zi(jrdm-1+3)
        end do
        zi(jnoe-1+2*nbma+1) = noeuds_extrem(2*nbma)
!
    else if (typm(1:4) .eq. 'SEG4') then
!
        nbno = 3*nbma+1
        call wkvect(noeord, base//' V I', nbno, jnoe)
        do i = 1, nbma
            zi(jnoe-1+3*i-2) = noeuds_extrem(i)
            call jeveuo(jexnum(conec, zi(jmail-1+i)), 'L', jrdm)
            ASSERT((zi(jrdm-1+1) .eq. noeuds_extrem(i)) .or. (zi(jrdm-1+2) .eq. noeuds_extrem(i)))
            if (zi(jrdm-1+1) .eq. noeuds_extrem(i)) then
                zi(jnoe-1+3*i-1) = zi(jrdm-1+3)
                zi(jnoe-1+3*i) = zi(jrdm-1+4)
            else if (zi(jrdm-1+2) .eq. noeuds_extrem(i)) then
                zi(jnoe-1+3*i-1) = zi(jrdm-1+4)
                zi(jnoe-1+3*i) = zi(jrdm-1+3)
            end if
        end do
        zi(jnoe-1+3*nbma+1) = noeuds_extrem(2*nbma)
!
    end if
!
!
!     ------------------------------------------------------------------
!     --- VERIFICATION DU NOEUD EXTREMITE LORSQU'IL EST DONNE
!     --- DANS LE CAS D UNE COURBE NON FERMEE
!     ------------------------------------------------------------------
    if (ndextr .ne. ' ') then
        noeud = int_to_char8(zi(jnoe-1+nbno))
        if (noeud .ne. ndextr) then
            call utmess('F', 'ELEMENTS_77', sk=ndextr)
        end if
    end if
!
!
!     -- SI VECORI EST RENSEIGNE (I.E. != 0),
!        IL FAUT EVENTUELLEMENT RETOURNER LA LISTE
!     ------------------------------------------------------------------
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ps1 = ddot(b_n, vecori, b_incx, vecori, b_incy)
    if (ps1 .gt. 0.d0) then
        ASSERT(nbno .ge. 3)
        ASSERT(zi(jnoe-1+1) .eq. zi(jnoe-1+nbno))
        call jeveuo(nomail//'.COORDO    .VALE', 'L', vr=vale)
!
!       PS1 : DDOT(VECORI,(1,2))/NORME((1,2))
        do k = 1, 3
            vecta(k) = vale(3*(zi(jnoe-1+2)-1)+k)
            vecta(k) = vecta(k)-vale(3*(zi(jnoe-1+1)-1)+k)
        end do
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ps1 = ddot(b_n, vecta, b_incx, vecori, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ps1 = ps1/sqrt(ddot(b_n, vecta, b_incx, vecta, b_incy))
!
!       PS2 : DDOT(VECORI,(N,N-1))/NORME((N,N-1))
        do k = 1, 3
            vecta(k) = vale(3*(zi(jnoe-1+nbno-1)-1)+k)
            vecta(k) = vecta(k)-vale(3*(zi(jnoe-1+nbno)-1)+k)
        end do
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ps2 = ddot(b_n, vecta, b_incx, vecori, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        ps2 = ps2/sqrt(ddot(b_n, vecta, b_incx, vecta, b_incy))
!
!       -- SI PS2 > PS1 : ON RETOURNE LA LISTE :
        if (ps2 .gt. ps1) then
            do k = 1, nbno/2
                ntemp = zi(jnoe-1+k)
                zi(jnoe-1+k) = zi(jnoe-1+nbno+1-k)
                zi(jnoe-1+nbno+1-k) = ntemp
            end do
        end if
    end if
!
!
!     -- MENAGE :
    AS_DEALLOCATE(vi=mailles_triee)
    AS_DEALLOCATE(vi=noeuds_extrem)
!
    call jedema()
end subroutine
