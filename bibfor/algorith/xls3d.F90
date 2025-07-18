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
subroutine xls3d(callst, grille, jltsv, jltsl, jlnsv, &
                 jlnsl, nbno, jcoor, jcoorg, nbmaf, &
                 jdlima, nbsef, jdlise, jconx1, jconx2, &
                 noma)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/normev.h"
#include "asterfort/padist.h"
#include "asterfort/panbno.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/xorima.h"
#include "asterfort/int_to_char8.h"
#include "blas/ddot.h"
!
    character(len=8) :: noma
    integer(kind=8) :: jltsv, jltsl, jlnsv, jlnsl, nbno, jcoor, jcoorg
    aster_logical :: callst, grille
!
! person_in_charge: samuel.geniaut at edf.fr
!
!
    real(kind=8) :: dmin, eps, eps1, eps2, eps3
    integer(kind=8) :: imafis, inoma, inose, isefis, itri, jconx1, jconx2, jma
    integer(kind=8) :: jdlima, jdlise, n1, n2, nbnoma, nbsef, nmaabs
    integer(kind=8) :: nseabs, ntri, num, nunoc, itypma, jcrd
    real(kind=8) :: xln, xlt
    integer(kind=8) :: ino, nbmaf, nuno(4), nunose(2), i, nbnott(3)
    real(kind=8) :: ab(3), ac(3), ap(3), vn(3), vnt(3), bc(3)
    real(kind=8) :: a(3), p(3), b(3), c(3), m(3), pm(3)
    real(kind=8) :: norme, ps, ps1, ps2, d
    aster_logical :: ma2ff
    character(len=19) :: mai, sens
    character(len=8) :: nomail
    real(kind=8) :: mprim(3), pmprim(3), cos, sin, vect(3), nove, pronor, angle
    real(kind=8) :: anglem
!
!-----------------------------------------------------------------------
    integer(kind=8) :: jsens
    integer(kind=8), pointer :: nbno_ma_fondfiss(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
!
    if (grille) then
        jcrd = jcoorg
    else
        jcrd = jcoor
    end if
!
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
!
!     TABLEAU POUR STOCKER LE NOMBRE DE NOEUDS SOMMETS
!     DES MAILLES DE FISSURE
    AS_ALLOCATE(vi=nbno_ma_fondfiss, size=nbmaf)
    do imafis = 1, nbmaf
        nmaabs = zi(jdlima+imafis-1)
        itypma = zi(jma-1+nmaabs)
        call panbno(itypma, nbnott)
        nbno_ma_fondfiss(imafis) = nbnott(1)
    end do
!
!     VERIFICATION DE L'ORIENTATION DES MAILLES DE LA FISSURES
    sens = '&&XLS3D.ORI_MAFIS'
    call xorima(noma, nbmaf, jdlima, jconx1, jconx2, &
                jcoor, sens)
    call jeveuo(sens, 'L', jsens)
!
!
!     BOUCLE SUR TOUS LES NOEUDS P DU MAILLAGE
    do ino = 1, nbno
!
        p(1) = zr(jcrd-1+3*(ino-1)+1)
        p(2) = zr(jcrd-1+3*(ino-1)+2)
        p(3) = zr(jcrd-1+3*(ino-1)+3)
!
!       CALCUL DE LSN
!       -------------
        dmin = r8maem()
!       RECHERCHE DE LA MAILLE LA PLUS PROCHE :
!       BOUCLE SUR NOEUDS DE MAFIS
        do imafis = 1, nbmaf
            nmaabs = zi(jdlima-1+(imafis-1)+1)
            nbnoma = zi(jconx2+nmaabs)-zi(jconx2+nmaabs-1)
            if ((nbnoma .eq. 4) .or. (nbnoma .eq. 8)) ntri = 4
            if ((nbnoma .eq. 3) .or. (nbnoma .eq. 6)) ntri = 1
!
!         BOUCLE SUR LE NOMBRE DE TRIANGLES DE LA MAILLE
            do itri = 1, ntri
!
                inoma = 1
                if (itri .eq. 4) inoma = 4
                nuno(inoma) = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
                a(1) = zr(jcoor-1+3*(nuno(inoma)-1)+1)
                a(2) = zr(jcoor-1+3*(nuno(inoma)-1)+2)
                a(3) = zr(jcoor-1+3*(nuno(inoma)-1)+3)
!
                inoma = 2
                if (itri .eq. 2) inoma = 3
                nuno(inoma) = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
                b(1) = zr(jcoor-1+3*(nuno(inoma)-1)+1)
                b(2) = zr(jcoor-1+3*(nuno(inoma)-1)+2)
                b(3) = zr(jcoor-1+3*(nuno(inoma)-1)+3)
!
                inoma = 3
                if (itri .eq. 2 .or. itri .eq. 3) inoma = 4
                nuno(inoma) = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
                c(1) = zr(jcoor-1+3*(nuno(inoma)-1)+1)
                c(2) = zr(jcoor-1+3*(nuno(inoma)-1)+2)
                c(3) = zr(jcoor-1+3*(nuno(inoma)-1)+3)
!
                do i = 1, 3
                    ab(i) = b(i)-a(i)
                    bc(i) = c(i)-b(i)
                    ap(i) = p(i)-a(i)
                    ac(i) = c(i)-a(i)
                end do
!
!           CALCUL DE LA NORMALE A LA MAILLE TRIA3
!           PROJECTION DE P SUR LA MAILLE VOIR R5.03.50-B
                call provec(ab, ac, vn)
                call normev(vn, norme)
                call provec(ap, vn, vnt)
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                ps = ddot(b_n, vnt, b_incx, ac, b_incy)
                eps1 = -1*ps/norme
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                ps = ddot(b_n, vnt, b_incx, ab, b_incy)
                eps2 = ps/norme
                eps3 = 1-eps1-eps2
!
!           SI M EST DS LE SECTEUR 1
                if (eps1 .lt. 0.d0) then
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps = ddot(b_n, ac, b_incx, ac, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps1 = ddot(b_n, ab, b_incx, ac, b_incy)
                    eps2 = eps2+eps1*ps1/ps
                    eps1 = 0.d0
                end if
!           SI M EST DS LE SECTEUR 2
                if (eps2 .lt. 0.d0) then
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps = ddot(b_n, ab, b_incx, ab, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps1 = ddot(b_n, ab, b_incx, ac, b_incy)
                    eps1 = eps1+eps2*ps1/ps
                    eps2 = 0.d0
                end if
!           SI M EST DS LE SECTEUR 3
                if (eps3 .lt. 0.d0) then
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps = ddot(b_n, bc, b_incx, bc, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps1 = ddot(b_n, ab, b_incx, bc, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps2 = ddot(b_n, ac, b_incx, bc, b_incy)
                    eps1 = (-1.d0*eps1*ps1+(1.d0-eps2)*ps2)/ps
                    eps2 = 1.d0-eps1
                end if
!
!           ON FINIT DE RAMENER LES POINTS ENCORE DEHORS
                if (eps1 .lt. 0.d0) eps1 = 0.d0
                if (eps2 .lt. 0.d0) eps2 = 0.d0
                if (eps1 .gt. 1.d0) eps1 = 1.d0
                if (eps2 .gt. 1.d0) eps2 = 1.d0
!
                do i = 1, 3
                    m(i) = a(i)+eps1*ab(i)+eps2*ac(i)
                    pm(i) = m(i)-p(i)
                end do
!
!           CALCUL DE LA DISTANCE PM
                d = padist(3, p, m)
!
!           MISE EN MEMOIRE DE LSN POUR LA MAILLE LA PLUS PROCHE
                if ((dmin-d) .gt. r8prem()*1.d04) then
                    dmin = d
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    xln = zi(jsens-1+imafis)*ddot(b_n, vn, b_incx, pm, b_incy)
                end if
!
            end do
!
        end do
!
        zr(jlnsv-1+(ino-1)+1) = xln
        zl(jlnsl-1+(ino-1)+1) = .true.
!
!
!       CALCUL DE LST
!       -------------
!
        if (.not. callst) then
            xlt = -1.d0
            goto 888
        end if
!
        dmin = r8maem()
        anglem = r8maem()
!
!       RECHERCHE DU SEGMENT LE PLUS PROCHE : BOUCLE SUR SEG DE FONFIS
        do isefis = 1, nbsef
!
            nseabs = zi(jdlise-1+(isefis-1)+1)
            nomail = int_to_char8(nseabs)
!
            inose = 1
            nunose(inose) = zi(jconx1-1+zi(jconx2+nseabs-1)+inose-1)
            inose = 2
            nunose(inose) = zi(jconx1-1+zi(jconx2+nseabs-1)+inose-1)
!
!         BOUCLE SUR LES MAILLES DE MAFIS POUR TROUVER LA BONNE MAILLE
            ma2ff = .false.
            do imafis = 1, nbmaf
!
                nmaabs = zi(jdlima-1+(imafis-1)+1)
!           ON RECUPERE LES NUMEROS DS NOEUDS DE LA MAILLE ET ON TESTE
                n1 = 0
                n2 = 0
!
                do inoma = 1, nbno_ma_fondfiss(imafis)
                    num = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
                    if (nunose(1) .eq. num) n1 = 1
                    if (nunose(2) .eq. num) n2 = 1
!             POUR RECUPERER UN 3EME POINT (SOMMET) DE LA MAILLE
!             QUI NE SOIT PAS SUR LE FOND
                    if ((nunose(1) .ne. num) .and. (nunose(2) .ne. num)) nunoc = num
                end do
!
                if ((n1*n2) .eq. 1) then
!
                    ma2ff = .true.
                    do i = 1, 3
                        a(i) = zr(jcoor-1+3*(nunose(1)-1)+i)
                        b(i) = zr(jcoor-1+3*(nunose(2)-1)+i)
                        c(i) = zr(jcoor-1+3*(nunoc-1)+i)
                        ab(i) = b(i)-a(i)
                        ap(i) = p(i)-a(i)
                        ac(i) = c(i)-a(i)
                    end do
!
!             CALCUL DE LA NORMALE A LA MAILLE
                    call provec(ab, ac, vn)
                    call normev(vn, norme)
!
!             CALCUL DE LA NORMALE INTERIEURE AU SEGMENT
                    call provec(ab, vn, vnt)
                    call normev(vnt, norme)
                    vn(1) = -1.d0*vnt(1)
                    vn(2) = -1.d0*vnt(2)
                    vn(3) = -1.d0*vnt(3)
!
!             PROJECTION SUR LE SEGMENT
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps = ddot(b_n, ap, b_incx, ab, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps1 = ddot(b_n, ab, b_incx, ab, b_incy)
                    eps = ps/ps1
!
                    do i = 1, 3
                        mprim(i) = a(i)+eps*ab(i)
                    end do
!
!             ON RAMENE M SUR LES BORDS S'IL LE FAUT
                    if (eps .gt. 1.d0) eps = 1.d0
                    if (eps .lt. 0.d0) eps = 0.d0
!
                    do i = 1, 3
                        m(i) = a(i)+eps*ab(i)
                        pm(i) = m(i)-p(i)
                        pmprim(i) = mprim(i)-p(i)
                    end do
!
!              CALCUL DE L'ANGLE (PM,PM')
!                  OU M EST LE PROJETE RAMENE
!                  ET M' LE PROJETE AVANT RAMENAGE
!              COS A = <U,V> / (||U|| * ||V||)
!              SIN A = ||U^V|| / (||U|| * ||V||)
                    call provec(pm, pmprim, vect)
                    call normev(vect, nove)
                    pronor = sqrt( &
                             pm(1)**2+pm(2)**2+pm(3)**2+pmprim(1)**2+pmprim(2)**2+pmprim(3)**2)
                    if (pronor .ne. 0.d0) then
                        cos = (pm(1)*pmprim(1)+pm(2)*pmprim(2)+pm(3)*pmprim(3))/pronor
                        sin = nove/pronor
                        angle = atan2(sin, cos)
                    else
                        cos = 0.d0
                        sin = 0.d0
                        angle = 0.d0
                    end if
!
!             CALCUL DE LA DISTANCE PM
                    d = padist(3, p, m)
!
!             MISE EN MEMOIRE DE LSN=PM.N POUR LE SEG LE PLUS PROCHE
                    if ((dmin-d) .gt. r8prem()*1.d04 .or. &
                        (abs(dmin-d) .le. r8prem()*1.d04 .and. angle .lt. anglem)) then
                        dmin = d
                        anglem = angle
                        b_n = to_blas_int(3)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        xlt = ddot(b_n, vn, b_incx, pm, b_incy)
                    end if
!
                end if
!
            end do
!
            if (.not. ma2ff) then
                call utmess('F', 'XFEM2_17')
            end if
        end do
!
888     continue
        zr(jltsv-1+(ino-1)+1) = xlt
        zl(jltsl-1+(ino-1)+1) = .true.
!
    end do
!
    AS_DEALLOCATE(vi=nbno_ma_fondfiss)
    call jedetr('&&XLS3D.ORI_MAFIS')
!
    call jedema()
end subroutine
