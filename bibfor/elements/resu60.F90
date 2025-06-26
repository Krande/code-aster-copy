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
subroutine resu60(resu1, resu2)
    implicit none
!
!     CETTE ROUTINE PERMET LA CONCATENATION DE DEUX CONCEPTS DYNA_GENE
!     DE TYPE HARMONIQUE CALCULES PAR DEUX COMMANDE DYNA_VIBRA//HARM/GENE
!     RESU1 ET RESU2 SONT COPIES DANS RESU1
!
!     LA ROUTINE RESU74 FAIT LA MEME CHOSE MAIS POUR DES CALCULS TRANS.
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/refdcp.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
    character(len=8) :: resu1, resu2
!
! IN  : RESU1 : PREMIER CONCEPT DYNA_GENE HARMONIQUE
! IN  : RESU2 : SECOND CONCEPT DYNA_GENE HARMONIQUE
!
    integer(kind=8) :: nbsto1, nbsau1, nbsto2, nbsau2
    integer(kind=8) :: nbstoc, nbsauv
    character(len=8) :: resu
    integer(kind=8) :: i
    integer(kind=8) :: flagd1, flagv1, flaga1, flagd2, flagv2, flaga2
    aster_logical :: flagd, flagv, flaga
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: jacce, jvite
    integer(kind=8) :: jdepl, jdesc
    integer(kind=8) :: jfreq, jordr
    complex(kind=8), pointer :: acce1(:) => null()
    complex(kind=8), pointer :: acce2(:) => null()
    complex(kind=8), pointer :: depl1(:) => null()
    complex(kind=8), pointer :: depl2(:) => null()
    complex(kind=8), pointer :: vite1(:) => null()
    complex(kind=8), pointer :: vite2(:) => null()
    real(kind=8), pointer :: freq1(:) => null()
    real(kind=8), pointer :: freq2(:) => null()
    integer(kind=8), pointer :: ordr1(:) => null()
    integer(kind=8), pointer :: ordr2(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    resu = '88888'
!
!     --- VERIFICATION DE LA COMPATIBILITE DES RESULTATS A MERGER
    call jeexin(resu1//'           .DEPL', flagd1)
    call jeexin(resu1//'           .VITE', flagv1)
    call jeexin(resu1//'           .ACCE', flaga1)
    call jeexin(resu2//'           .DEPL', flagd2)
    call jeexin(resu2//'           .VITE', flagv2)
    call jeexin(resu2//'           .ACCE', flaga2)
!
    flagd = ((flagd1 .eq. 0) .and. (flagd2 .eq. 0)) .or. ((flagd1 .ne. 0) .and. (flagd2 .ne. 0))
    flagv = ((flagv1 .eq. 0) .and. (flagv2 .eq. 0)) .or. ((flagv1 .ne. 0) .and. (flagv2 .ne. 0))
    flaga = ((flaga1 .eq. 0) .and. (flaga2 .eq. 0)) .or. ((flaga1 .ne. 0) .and. (flaga2 .ne. 0))
!
!     CONDITION POUR SAVOIR SI LES FLAGS SONT BIEN TOUS LES 2 ZEROS
!     OU BIEN DIFFERENTS DE ZERO = COMPATIBILITE DES RESUS
    if (.not. (flagd .and. flagv .and. flaga)) then
        call utmess('F', 'ALGORITH17_25')
    end if
!
    call jeveuo(resu1//'           .DESC', 'E', jdesc)
!
!     --- RECHERCHE DU NUMERO D'ORDRE DE LA FREQUENCE DE REPRISE
    if (flagd1 .gt. 0) then
        call jelira(resu1//'           .DEPL', 'LONUTI', nbsto1)
        call jelira(resu2//'           .DEPL', 'LONUTI', nbsto2)
    else if (flagv1 .gt. 0) then
        call jelira(resu1//'           .VITE', 'LONUTI', nbsto1)
        call jelira(resu2//'           .VITE', 'LONUTI', nbsto2)
    else
        call jelira(resu1//'           .ACCE', 'LONUTI', nbsto1)
        call jelira(resu2//'           .ACCE', 'LONUTI', nbsto2)
    end if
!
    nbstoc = nbsto1+nbsto2
!
!     --- RECUPERATION DES CHAMPS DEPL VITE ET ACCE ---
    if (flagd1 .ne. 0) then
        call jeveuo(resu1//'           .DEPL', 'E', vc=depl1)
        call jeveuo(resu2//'           .DEPL', 'E', vc=depl2)
        call wkvect(resu//'           .DEPL', 'G V C', nbstoc, jdepl)
        b_n = to_blas_int(nbsto1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, depl1, b_incx, zc(jdepl), b_incy)
        b_n = to_blas_int(nbsto2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, depl2, b_incx, zc(jdepl+nbsto1), b_incy)
    end if
!                         ^
!                         |______________
!       --- VALEURS COMPLEXES = COPIER |2| FOIS PLUS DE REELS
!
    if (flagv1 .ne. 0) then
        call jeveuo(resu1//'           .VITE', 'E', vc=vite1)
        call jeveuo(resu2//'           .VITE', 'E', vc=vite2)
        call wkvect(resu//'           .VITE', 'G V C', nbstoc, jvite)
        b_n = to_blas_int(nbsto1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, vite1, b_incx, zc(jvite), b_incy)
        b_n = to_blas_int(nbsto2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, vite2, b_incx, zc(jvite+nbsto1), b_incy)
    end if
!
    if (flaga1 .ne. 0) then
        call jeveuo(resu1//'           .ACCE', 'E', vc=acce1)
        call jeveuo(resu2//'           .ACCE', 'E', vc=acce2)
        call wkvect(resu//'           .ACCE', 'G V C', nbstoc, jacce)
        b_n = to_blas_int(nbsto1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, acce1, b_incx, zc(jacce), b_incy)
        b_n = to_blas_int(nbsto2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, acce2, b_incx, zc(jacce+nbsto1), b_incy)
    end if
!
!     --- RECUPERATION DES CHAMPS ORDR
!
    call jeveuo(resu1//'           .ORDR', 'E', vi=ordr1)
    call jelira(resu1//'           .ORDR', 'LONUTI', nbsau1)
!
    call jeveuo(resu2//'           .ORDR', 'E', vi=ordr2)
    call jelira(resu2//'           .ORDR', 'LONUTI', nbsau2)
!     --- CUMULER LES NUMEROS D'ORDRE POUR CONSERVER LA MONOTONIE
    do i = 0, nbsau2-1
        ordr2(1+i) = ordr2(1+i)+ordr1(nbsau1)+1
    end do
!
    nbsauv = nbsau1+nbsau2
!
    call wkvect(resu//'           .ORDR', 'G V I', nbsauv, jordr)
    call copvis(nbsau1, ordr1, zi(jordr))
    call copvis(nbsau2, ordr2, zi(jordr+nbsau1))
!
    call jeveuo(resu1//'           .DISC', 'E', vr=freq1)
    call jeveuo(resu2//'           .DISC', 'E', vr=freq2)
    call wkvect(resu//'           .DISC', 'G V R', nbsauv, jfreq)
    b_n = to_blas_int(nbsau1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, freq1, b_incx, zr(jfreq), b_incy)
    b_n = to_blas_int(nbsau2)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, freq2, b_incx, zr(jfreq+nbsau1), b_incy)
!
!     --- DUPLICATION ---
!
    if (flagd1 .ne. 0) call jedupo(resu//'           .DEPL', 'G', resu1//'           .DEPL', &
                                   .false._1)
    if (flagv1 .ne. 0) call jedupo(resu//'           .VITE', 'G', resu1//'           .VITE', &
                                   .false._1)
    if (flaga1 .ne. 0) call jedupo(resu//'           .ACCE', 'G', resu1//'           .ACCE', &
                                   .false._1)
    call jedupo(resu//'           .ORDR', 'G', resu1//'           .ORDR', .false._1)
    call jedupo(resu//'           .DISC', 'G', resu1//'           .DISC', .false._1)
!
!
!     --- COPIE DU NOUVEAU .REFD DANS LA SD FINALE ---
!
    call refdcp(resu2, resu1)
!
!     --- DESTRUCTION DES OBJETS PROVISOIRES
!
    call jedetc('G', resu//'           ', 1)
    call jedetc('G', resu2//'           ', 1)
!
    call jedema()
!
end subroutine
