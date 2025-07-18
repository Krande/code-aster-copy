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
subroutine op0031()
    implicit none
!
!     COMBINAISON LINEAIRE DE MATRICE
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/amogen.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/gettco.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/utmess.h"
#include "asterfort/vrrefe.h"
#include "asterfort/wkvect.h"
    character(len=1) :: typres, base
    character(len=6) :: combrc
    character(len=8) :: matres, matri1, partie
    character(len=19) :: matr19, nomi
    character(len=8) :: nomddl
    character(len=16) :: concep, nomcmd, typrep
    character(len=24) :: cnom, ccoef, ctypec, valk(2)
    real(kind=8) :: r8val(2)
    aster_logical :: lcoefc, lreent
    complex(kind=8) :: cval
    integer(kind=8) :: nbocag, nboccr, nboccc, nbocc, ldesc, l, lnom, iocc, i
    integer(kind=8) :: ibid, lcoef, ltypec, nbcst, lr, lc, iret, ides1
    integer(kind=8) :: jrefe, jpomr
! ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call getres(matres, concep, nomcmd)
    call gettco(matres, typrep)
    matr19 = matres
!
!     -- LA MATRICE RESULTAT  EST-ELLE REENTRANTE ?
    call jeexin(matr19//'.REFA', iret)
    lreent = (iret .ne. 0)
    base = 'G'
    if (lreent) then
        matr19 = '&&OP0031.MATRES'
        base = 'V'
    end if
!
    call getfac('CALC_AMOR_GENE', nbocag)
    if (nbocag .ne. 0) then
        call amogen(matr19)
        goto 999
    end if
!
    call getfac('COMB_R', nboccr)
    call getfac('COMB_C', nboccc)
    if (nboccr .ne. 0) then
        nbocc = nboccr
        typres = 'R'
        combrc = 'COMB_R'
    else
        nbocc = nboccc
        typres = 'C'
        combrc = 'COMB_C'
    end if
!
!
!
!   -- Creation du .DESC pour la MATR_ASSE_GENE resultat :
!   -----------------------------------------------------
    if (typrep(1:14) .eq. 'MATR_ASSE_GENE') then
        call wkvect(matr19//'.DESC', base//' V I', 3, ldesc)
        do iocc = 1, nbocc
            call getvid(combrc, 'MATR_ASSE', iocc=iocc, scal=matri1, nbret=l)
            call jeveuo(matri1//'           .DESC', 'L', ides1)
            if (iocc .eq. 1) then
                do i = 1, 3
                    zi(ldesc+i-1) = zi(ides1+i-1)
                end do
            else
                ASSERT(zi(ldesc) .eq. zi(ides1))
                if (zi(ldesc-1+2) .ne. zi(ides1-1+2)) call utmess('F', 'ALGELINE2_29')
! si l'une des matr_asse_gene est pleine, le resultat le sera aussi :
                if (zi(ides1-1+3) .eq. 2) zi(ldesc-1+3) = 2
                if (zi(ides1-1+3) .eq. 3) then
                    if (zi(ldesc-1+3) .eq. 1) zi(ldesc-1+3) = 3
                end if
            end if
        end do
    end if
!
!
    cnom = '&&OP0031.LISTE_MATRICE'
    call wkvect(cnom, 'V V K8', nbocc, lnom)
    jpomr = 0
    do iocc = 0, nbocc-1
        call getvid(combrc, 'MATR_ASSE', iocc=iocc+1, scal=zk8(lnom+iocc), nbret=l)
!       -- on recherche une eventuelle matrice non symetrique
        nomi = zk8(lnom+iocc)
        call jeveuo(nomi//'.REFA', 'L', jrefe)
        if (zk24(jrefe-1+9) .eq. 'MR') then
            jpomr = iocc
        end if
    end do
!
!
    nomddl = ' '
    call getvtx(' ', 'SANS_CMP', scal=nomddl, nbret=ibid)
!
!
!   --- recuperation des coefficients :
!   ------------------------------------
!   remarque : pour partie='imag', on fait coef=-j*coef
!
    ccoef = '&&OP0031.COEF_VALEURS'
    call wkvect(ccoef, 'V V R', 2*nbocc, lcoef)
    ctypec = '&&OP0031.COEF_TYPE'
    call wkvect(ctypec, 'V V K8', nbocc, ltypec)
!
    nbcst = 0
    do iocc = 0, nbocc-1
        call getvr8(combrc, 'COEF_R', iocc=iocc+1, scal=r8val(1), nbret=lr)
        if (lr .eq. 1) then
            lcoefc = .false.
            if (combrc .eq. 'COMB_R') then
                partie = ' '
                call getvtx(combrc, 'PARTIE', iocc=iocc+1, scal=partie, nbret=ibid)
                if (partie .eq. 'IMAG') lcoefc = .true.
            end if
!
            if (.not. lcoefc) then
                zr(lcoef+nbcst) = r8val(1)
                nbcst = nbcst+1
                zk8(ltypec+iocc) = 'R'
            else
                zr(lcoef+nbcst) = 0.d0
                zr(lcoef+nbcst+1) = -1.d0*r8val(1)
                nbcst = nbcst+2
                zk8(ltypec+iocc) = 'C'
            end if
        else
            call getvc8(combrc, 'COEF_C', iocc=iocc+1, scal=cval, nbret=lc)
            ASSERT(lc .eq. 1)
            zr(lcoef+nbcst) = dble(cval)
            zr(lcoef+nbcst+1) = dimag(cval)
            nbcst = nbcst+2
            zk8(ltypec+iocc) = 'C'
        end if
    end do
!
!
!   --- controle des references :
!   --------------------------------
    do iocc = 0, nbocc-2
        call vrrefe(zk8(lnom+iocc), zk8(lnom+iocc+1), iret)
        if (iret .ne. 0) then
            valk(1) = zk8(lnom+iocc)
            valk(2) = zk8(lnom+iocc+1)
            call utmess('F', 'ALGELINE2_28', nk=2, valk=valk)
        end if
    end do
!
!
!
!   -- combinaison des matrices :
!   ------------------------------------------------------------------
! initialisation de la matrice resultat :
    call mtdefs(matr19, zk8(lnom+jpomr), base, typres)
    call mtcmbl(nbocc, zk8(ltypec), zr(lcoef), zk8(lnom), matr19, &
                nomddl, ' ', 'ELIM=')
!
!
!   -- si la matrice est reentrante, on la detruit et on recopie
!      la matrice intermediaire :
!   -------------------------------------------------------------
    if (lreent) then
        ASSERT(matr19(1:8) .eq. '&&OP0031')
        call detrsd('MATR_ASSE', matres)
        call copisd('MATR_ASSE', 'G', matr19, matres)
        call detrsd('MATR_ASSE', matr19)
    end if
!
999 continue
    call jedema()
end subroutine
