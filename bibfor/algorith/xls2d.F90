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
subroutine xls2d(callst, grille, jltsv, jltsl, jlnsv, &
                 jlnsl, nbno, jcoor, jcoorg, nbmaf, &
                 jdlima, nbsef, jdlise, jconx1, jconx2)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/padist.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    integer(kind=8) :: nbno, jcoor, jcoorg, nbmaf, nbsef, jdlima, jdlise
    integer(kind=8) :: jlnsv, jlnsl, jltsv, jltsl, jconx1, jconx2
    aster_logical :: callst, grille
!
! person_in_charge: samuel.geniaut at edf.fr
    integer(kind=8) :: ino, imafis, nmaabs, inoma, nuno(2), jcrd
    real(kind=8) :: p(2), dmin, a(2), b(2), m(2), ap(2), ab(2), norcab, ps, eps
    real(kind=8) :: d, oriabp, xln, ps1, xlt, tole
    integer(kind=8) :: isefis, nseabs, inose, nunose, n1, nbnoma, num, nunoc, i
    aster_logical :: ma2ff
    integer(kind=8) :: ir, ir2, ir3, jmafit, jmafif, jmaori, nuno1, nuno2, nunoi, ori
    aster_logical :: finfis
    aster_logical, pointer :: is_pt_fond(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    parameter(tole=1.d-12)
!
    call jemarq()
!
    if (grille) then
        jcrd = jcoorg
    else
        jcrd = jcoor
    end if
!
!     REORGANISATION DES MAILLES DE LA FISSURE
!     VECTEURS INTERMEDIAIRE ET ORIENTATION DES MAILLES
    call wkvect('&&XINILS.LIMFISO', 'V V I', nbmaf, jmafit)
    call wkvect('&&XINILS.ORIENT', 'V V I', nbmaf, jmaori)
    AS_ALLOCATE(vl=is_pt_fond, size=nbno)
    is_pt_fond(1:nbno) = .false.
    do isefis = 1, nbsef
        nseabs = zi(jdlise-1+(isefis-1)+1)
        is_pt_fond(zi(jconx1-1+zi(jconx2+nseabs-1))) = .true.
    end do
!     INITIALISATION PREMIERE MAILLE
    ori = 1
    zi(jmaori) = 1
    zi(jmafit) = zi(jdlima)
!     PARCOURS DES MAILLES CONTIGUES DANS 1 SENS
    do ir = 2, nbmaf
        nunoi = zi(jconx1-1+zi(jconx2+zi(jmafit+ir-1-1)-1)+1+ori-1)
        finfis = .true.
        do ir2 = 1, nbmaf
            nuno1 = zi(jconx1-1+zi(jconx2+zi(jdlima+ir2-1)-1)+1-1)
            nuno2 = zi(jconx1-1+zi(jconx2+zi(jdlima+ir2-1)-1)+2-1)
            if (zi(jdlima+ir2-1) .ne. zi(jmafit+ir-1-1) .and. &
                (nunoi .eq. nuno1 .or. nunoi .eq. nuno2)) then
                if (nunoi .eq. nuno1) ori = 1
                if (nunoi .eq. nuno2) ori = 0
                zi(jmaori+ir-1) = ori
                zi(jmafit+ir-1) = zi(jdlima+ir2-1)
                finfis = .false.
                goto 1135
            end if
        end do
1135    continue
        if (finfis) goto 1145
    end do
1145 continue
!
!     VECTEUR FINAL REORGANISE
    call wkvect('&&XINILS.LIMFISOF', 'V V I', nbmaf, jmafif)
!     DECALAGE DES MAILLES TROUVEES
    do ir3 = 1, ir-1
        zi(jmafif+nbmaf-ir+ir3) = zi(jmafit-1+ir3)
    end do
!     PARCOURS DANS L'AUTRE SENS A PARTIR DE LA PREMIERE MAILLE
    ori = 0
    do ir2 = ir, nbmaf+1
        nunoi = zi(jconx1-1+zi(jconx2+zi(jmafif+nbmaf-ir2+1)-1)+1+ori-1)
        finfis = .true.
        do ir3 = 1, nbmaf
            nuno1 = zi(jconx1-1+zi(jconx2+zi(jdlima+ir3-1)-1)+1-1)
            nuno2 = zi(jconx1-1+zi(jconx2+zi(jdlima+ir3-1)-1)+2-1)
            if (zi(jdlima+ir3-1) .ne. zi(jmafif+nbmaf-ir2+1) .and. &
                (nunoi .eq. nuno1 .or. nunoi .eq. nuno2)) then
                if (nunoi .eq. nuno1) ori = 1
                if (nunoi .eq. nuno2) ori = 0
                zi(jmaori+nbmaf-ir2) = ori
                zi(jmafif+nbmaf-ir2) = zi(jdlima+ir3-1)
                finfis = .false.
                goto 1165
            end if
        end do
1165    continue
        if (finfis) goto 1175
    end do
1175 continue
!
!      DO 118 IR3=1,NBMAF
!      WRITE (6,*) IR3, ZI(JCONX1-1+ZI(JCONX2+ZI(JDLIMA+IR3-1)-1)+1-1)
!     & , ZI(JCONX1-1+ZI(JCONX2+ZI(JDLIMA+IR3-1)-1)+2-1)
! 118  CONTINUE
!      DO 1182 IR3=1,NBMAF
!      WRITE (6,*) IR3, ZI(JCONX1-1+ZI(JCONX2+ZI(JMAFIF+IR3-1)-1)+1-1)
!     & , ZI(JCONX1-1+ZI(JCONX2+ZI(JMAFIF+IR3-1)-1)+2-1)
!     & , ZI(JMAORI-1+IR3)
! 1182 CONTINUE
!
    if (ir2-1 .ne. nbmaf) then
        write (6, *) ir, ir2, nbmaf
!     MAILLES MANQUANTES
        ASSERT(ir2-1 .eq. nbmaf)
    end if
!
!     BOUCLE SUR LES NOEUDS P DU MAILLAGE
    do ino = 1, nbno
        p(1) = zr(jcrd-1+3*(ino-1)+1)
        p(2) = zr(jcrd-1+3*(ino-1)+2)
!
!     CALCUL DE LSN
!     -------------
        dmin = r8maem()
        xln = r8maem()
!         RECHERCHE DE LA MAILLE LA PLUS PROCHE :
!         BOUCLE SUR NOEUDS DE MAFIS
        do imafis = 1, nbmaf
            nmaabs = zi(jmafif-1+(imafis-1)+1)
            inoma = 1
            nuno(inoma) = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
            a(1) = zr(jcoor-1+3*(nuno(inoma)-1)+1)
            a(2) = zr(jcoor-1+3*(nuno(inoma)-1)+2)
!
            inoma = 2
            nuno(inoma) = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
            b(1) = zr(jcoor-1+3*(nuno(inoma)-1)+1)
            b(2) = zr(jcoor-1+3*(nuno(inoma)-1)+2)
!
            do i = 1, 2
                ab(i) = b(i)-a(i)
                ap(i) = p(i)-a(i)
            end do
!
!           CALCUL DE EPS TEL QUE AM=EPS*AB
            norcab = ab(1)*ab(1)+ab(2)*ab(2)
            b_n = to_blas_int(2)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            ps = ddot(b_n, ap, b_incx, ab, b_incy)
            eps = ps/norcab
!
!           ON RAMENE LES POINTS EN DEHORS DU SEGMENT
!             > SI LE POINT N EST PAS SUR LA DROITE DIRECTRICE AU SEGMENT
            b_n = to_blas_int(2)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            if (abs(ddot(b_n, ap, b_incx, [-ab(2), ab(1)], b_incy)) .gt. norcab*tole) then
                if (eps .lt. -tole .and. .not. is_pt_fond(nuno(1))) eps = 0.d0
                if (eps .gt. (1.d0+tole) .and. .not. is_pt_fond(nuno(2))) eps = 1.d0
            end if
!
            do i = 1, 2
                m(i) = a(i)+eps*ab(i)
            end do
!
!           CALCUL DE LA DISTANCE PM
            d = padist(2, p, m)
!
!           MISE EN MEMOIRE DE LSN POUR LA MAILLE LA PLUS PROCHE
!           EN VERIFIANT DE QUEL COTE DE LA FISSURE SE TROUVE P
            if ((dmin-d) .gt. (r8prem()*1.d02)) then
                dmin = d
                oriabp = ab(1)*ap(2)-ab(2)*ap(1)
                do i = 1, 2
                    m(i) = a(i)+ps/norcab*ab(i)
                end do
                d = padist(2, p, m)
                if (oriabp .gt. 0.d0) then
                    xln = d
                else
                    xln = -1.d0*d
                end if
!                mp(1:2)=p(1:2)-m(1:2)
!                xln=(ab(1)*mp(2)-ab(2)*mp(1))/sqrt(norcab)
                if (zi(jmaori-1+imafis) .eq. 0) then
                    xln = -1.d0*xln
                end if
            end if
!
        end do
!
        zr(jlnsv-1+(ino-1)+1) = xln
        zl(jlnsl-1+(ino-1)+1) = .true.
!
!        CALCUL DE LST
!        -------------
!
        if (.not. callst) then
            xlt = -1.d0
            goto 888
        end if
!
        dmin = r8maem()
        xlt = r8maem()
!
!         RECHERCHE DU POINT LE PLUS PROCHE : BOUCLE SUR POINT DE FONFIS
        do isefis = 1, nbsef
!
            nseabs = zi(jdlise-1+(isefis-1)+1)
            inose = 1
            nunose = zi(jconx1-1+zi(jconx2+nseabs-1)+inose-1)
!
!           BOUCLE SUR LES MAILLES DE MAFIS POUR TROUVER LA BONNE MAILLE
            ma2ff = .false.
            do imafis = 1, nbmaf
!
                nmaabs = zi(jmafif-1+(imafis-1)+1)
                nbnoma = zi(jconx2+nmaabs)-zi(jconx2+nmaabs-1)
!             ON RECUPERE LES NUMEROS DS NOEUDS DE LA MAILLE ET ON TESTE
                n1 = 0
!
                do inoma = 1, nbnoma
                    num = zi(jconx1-1+zi(jconx2+nmaabs-1)+inoma-1)
                    if (nunose .eq. num) n1 = 1
!               POUR RECUPERER UN 2EME POINT DE LA MAILLE QUI NE SOIT
!               PAS SUR LE FOND
                    if ((nunose .ne. num)) nunoc = num
                end do
!
                if (n1 .eq. 1) then
!
                    ma2ff = .true.
                    do i = 1, 2
                        a(i) = zr(jcoor-1+3*(nunose-1)+i)
                        b(i) = zr(jcoor-1+3*(nunoc-1)+i)
                        ab(i) = b(i)-a(i)
                        ap(i) = p(i)-a(i)
                    end do
!
!               PROJECTION SUR LE SEGMENT
                    b_n = to_blas_int(2)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps = ddot(b_n, ap, b_incx, ab, b_incy)
                    b_n = to_blas_int(2)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    ps1 = ddot(b_n, ab, b_incx, ab, b_incy)
                    eps = ps/ps1
!
!               CALCUL DE LA DISTANCE PA
                    d = padist(2, p, a)
!               MISE EN MEMOIRE DE LSN=PA.N POUR LE SEG LE PLUS PROCHE
                    if ((dmin-d) .gt. (r8prem()*1.d02)) then
                        dmin = d
                        xlt = -1.d0*eps*sqrt(ab(1)*ab(1)+ab(2)*ab(2))
                    end if
!
                end if
!
            end do
!
            if (.not. ma2ff) then
                call utmess('F', 'XFEM2_15')
            end if
        end do
!
888     continue
        zr(jltsv-1+(ino-1)+1) = xlt
        zl(jltsl-1+(ino-1)+1) = .true.
!
    end do
!
    AS_DEALLOCATE(vl=is_pt_fond)
    call jedema()
!
end subroutine
