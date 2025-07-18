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
subroutine wpnorm(norm, para, lmatr, neq, nbmode, &
                  ddlexc, vecpro, resufr, coef)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcmult.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: norm, para
    integer(kind=8) :: nbmode, neq, lmatr(*), ddlexc(*)
    complex(kind=8) :: vecpro(neq, *)
    real(kind=8) :: resufr(nbmode, *), coef(*)
!     NORMALISATION DE VECTEURS COMPLEXES ET DE GRANDEURS MODALES
!     ------------------------------------------------------------------
! IN  NORM   : TYPE DE NORMALISATION
!          = 'AVEC_CMP'
!          = 'MASS_GENE'
!          = 'RIGI_GENE'
!          = 'EUCL'
! IN  PARA   : ON REPERCUTE LA NORMALISATION SUR LES PARAMETRES MODAUX
!          = 'OUI' DANS CE CAS ILS DOIVENT DEJA AVOIR ETE CALCULES
!          = 'NON' ON NE NORMALISE QUE LES VECTEURS PROPRES
! IN  LMTR   : DESCRIPTEUR D'UNE MATRICE
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  NBMODE : NOMBRE DE MODES
! IN  DDLEXC : TABLEAU DES DDL EXCLUS
!          = 0 SI EXCLUS
!          = 1 SI NON EXCLUS
! VAR VECPRO : TABLEAU DES VECTEURS PROPRES
! VAR RESUFR : TABLEAU DES GRANDEURS MODALES RANGEES SELON
!        'FREQ'            , 'OMEGA2'          , 'AMOR_REDUIT'     ,
!        'MASS_GENE'       , 'RIGI_GENE'       , 'AMOR_GENE'       ,
!        'MASS_EFFE_DX'    , 'MASS_EFFE_DY'    , 'MASS_EFFE_DZ'    ,
!        'FACT_PARTICI_DX' , 'FACT_PARTICI_DY' , 'FACT_PARTICI_DZ' ,
!        'MASS_EFFE_UN_DX' , 'MASS_EFFE_UN_DY' , 'MASS_EFFE_UN_DZ'
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: im, ieq, ldynam
    character(len=1) :: typcst(2)
    character(len=19) :: matmod
    real(kind=8) :: rnorm, rx1, rx2, constr(4), fr, am, zero
    complex(kind=8) :: xx1, cmpl, xnorm, dconjg, czero
    character(len=24) :: nmatr(2), ndynam
    character(len=24) :: valk
    complex(kind=8), pointer :: xxxx_gene_2(:) => null()
!     ------------------------------------------------------------------
    data typcst/'C', 'C'/
!     ------------------------------------------------------------------
!
    call jemarq()
    zero = 0.d0
    czero = dcmplx(zero, zero)
!
    if (norm .eq. 'AVEC_CMP' .or. norm .eq. 'EUCL') then
!
!        --- NORMALISATION SUR LES DDL NON EXCLUS
        do im = 1, nbmode
            rnorm = 0.0d0
            if (norm .eq. 'EUCL') then
                do ieq = 1, neq
                    xx1 = vecpro(ieq, im)*ddlexc(ieq)
                    rnorm = rnorm+dble(xx1*dconjg(xx1))
                end do
                rnorm = sqrt(rnorm)
            else
                do ieq = 1, neq
                    rx1 = abs(vecpro(ieq, im)*ddlexc(ieq))
                    rnorm = max(rx1, rnorm)
                end do
            end if
            rx1 = 1.0d0/rnorm
            coef(im) = rx1
            do ieq = 1, neq
                vecpro(ieq, im) = vecpro(ieq, im)*rx1
            end do
            if (para .eq. 'OUI') then
                rx2 = rx1*rx1
                resufr(im, 4) = resufr(im, 4)*rx2
                resufr(im, 5) = resufr(im, 5)*rx2
                resufr(im, 6) = resufr(im, 6)*rx2
            end if
        end do
!
    else if (norm .eq. 'MASS_GENE' .or. norm .eq. 'RIGI_GENE') then
!
!        --- ON NORMALISE LA MASSE OU LA RAIDEUR GENERALISEE A 1 ---
!        --- DU PROBLEME GENERALISE ASSOCIE AU PROBLEME QUADRATIQUE ---
        matmod = zk24(zi(lmatr(1)+1))
        nmatr(1) = zk24(zi(lmatr(1)+1))
        nmatr(2) = zk24(zi(lmatr(2)+1))
        AS_ALLOCATE(vc=xxxx_gene_2, size=neq)
        call mtdefs('&&WPNORM.MATR.DYNAMIC', matmod, 'V', 'C')
        call mtdscr('&&WPNORM.MATR.DYNAM')
        ndynam = '&&WPNORM.MATR.DYNAM'//'.&INT'
        call jeveuo(ndynam, 'E', ldynam)
        if (norm .eq. 'MASS_GENE') then
            constr(3) = 1.d0
            constr(4) = 0.d0
        else if (norm .eq. 'RIGI_GENE') then
            constr(1) = -1.d0
            constr(2) = 0.d0
        end if
        do im = 1, nbmode
            fr = sqrt(resufr(im, 2))
            am = resufr(im, 3)
            am = -abs(am*fr)/sqrt(1.0d0-am*am)
            if (norm .eq. 'MASS_GENE') then
                constr(1) = 2.d0*am
                constr(2) = 2.d0*fr
            else if (norm .eq. 'RIGI_GENE') then
                cmpl = dcmplx(am, fr)
                cmpl = cmpl*cmpl
                constr(3) = dble(cmpl)
                constr(4) = dimag(cmpl)
            end if
            call mtcmbl(2, typcst, constr, nmatr, ndynam, &
                        ' ', ' ', 'ELIM=')
            call mcmult('ZERO', ldynam, vecpro(1, im), xxxx_gene_2, 1, &
                        .true._1)
            xnorm = czero
            do ieq = 1, neq
                xnorm = xnorm+vecpro(ieq, im)*xxxx_gene_2(ieq)
            end do
            xnorm = 1.d0/sqrt(xnorm)
            coef(im) = dble(xnorm)
            do ieq = 1, neq
                vecpro(ieq, im) = vecpro(ieq, im)*xnorm
            end do
            if (para .eq. 'OUI') then
                xnorm = xnorm*xnorm
                resufr(im, 4) = resufr(im, 4)*dble(xnorm)
                resufr(im, 5) = resufr(im, 5)*dble(xnorm)
                resufr(im, 6) = resufr(im, 6)*dble(xnorm)
            end if
        end do
! --- MENAGE
        AS_DEALLOCATE(vc=xxxx_gene_2)
        call detrsd('MATR_ASSE', '&&WPNORM.MATR.DYNAMIC')
        call detrsd('MATR_ASSE', '&&WPNORM.MATR.DYNAM')
        call jedetr('&&WPNORM.MATR.DYNAM'//'.&INT')
!
    else
!
        valk = norm
        call utmess('F', 'ALGELINE4_77', sk=valk)
!
    end if
!
    call jedema()
end subroutine
