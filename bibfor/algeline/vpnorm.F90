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
subroutine vpnorm(norm, para, lmatr, neq, nbmode, &
                  ddlexc, vecpro, resufr, xmastr, isign, &
                  numddl, coef)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: norm, para
    integer(kind=8) :: nbmode, neq, lmatr, ddlexc(*)
    real(kind=8) :: vecpro(neq, *), resufr(nbmode, *), xmastr(3), coef(*)
!     NORMALISATION DE VECTEURS ET DE GRANDEURS MODALES
!     ------------------------------------------------------------------
! IN  NORM   : TYPE DE NORMALISATION
!          = 'AVEC_CMP'
!          = 'MASS_GENE'
!          = 'RIGI_GENE'
!          = 'EUCL', 'EUCL_TRAN', ...
! IN  PARA   : ON REPERCUTE LA NORMALISATION SUR LES PARAMETRES MODAUX
!          = 'OUI' DANS CE CAS ILS DOIVENT DEJA AVOIR ETE CALCULES
!          = 'NON' ON NE NORMALISE QUE LES VECTEURS PROPRES
! IN  LMATR   : DESCRIPTEUR D'UNE MATRICE
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
! IN  XMASTR : MASSE DE LA STRUCTURE
! OUT COEF   : COEFFICIENTS
!     ------------------------------------------------------------------
!
    character(len=24) :: valk
!
    real(kind=8) :: xmn, xx1, xx2, xx3, xnorm, epsi
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ie, im, indg, isign, numddl
    real(kind=8), pointer :: poi1(:) => null()
    real(kind=8), pointer :: poi2(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    epsi = r8prem()
    if (norm .eq. 'AVEC_CMP' .or. norm(1:4) .eq. 'EUCL') then
!
!     --- NORMALISATION SUR LES DDL NON EXCLUS
!
        do im = 1, nbmode
            xnorm = 0.0d0
            if (norm(1:4) .eq. 'EUCL') then
                do ie = 1, neq
                    xx1 = vecpro(ie, im)*ddlexc(ie)
                    xnorm = xnorm+xx1*xx1
                end do
                xnorm = sqrt(xnorm)
            else
                do ie = 1, neq
                    xx1 = vecpro(ie, im)*ddlexc(ie)
                    if (abs(xnorm) .lt. abs(xx1)) then
                        xnorm = xx1
                    end if
                end do
            end if
            if (abs(xnorm) .lt. epsi) then
                call utmess('F', 'MODAL_23', si=im, sr=xnorm)
            end if
            xx1 = 1.0d0/xnorm
            coef(im) = xx1
            do ie = 1, neq
                vecpro(ie, im) = vecpro(ie, im)*xx1
            end do
            if (para .eq. 'OUI') then
                xx2 = xx1*xx1
                resufr(im, 4) = resufr(im, 4)*xx2
                resufr(im, 5) = resufr(im, 5)*xx2
!-PROV        RESUFR(IM,6)  = RESUFR(IM,6)  * XX2
                resufr(im, 10) = resufr(im, 10)*xnorm
                resufr(im, 11) = resufr(im, 11)*xnorm
                resufr(im, 12) = resufr(im, 12)*xnorm
            else
                xx2 = xx1*xx1
                resufr(im, 4) = resufr(im, 4)*xx2
                resufr(im, 5) = resufr(im, 5)*xx2
!-PROV        RESUFR(IM,6)  = RESUFR(IM,6)  * XX2
            end if
        end do
!
    else if (norm .eq. 'MASS_GENE' .or. norm .eq. 'RIGI_GENE') then
!
!     --- ON NORMALISE LA MASSE OU LA RAIDEUR GENERALISEE A 1 ---
!
        indg = 4
        if (norm .eq. 'RIGI_GENE') indg = 5
        if (para .eq. 'OUI') then
            do im = 1, nbmode
                xmn = resufr(im, indg)
                xx1 = 1.0d0/xmn
                xx2 = sqrt(xmn)
                xx3 = 1.0d0/xx2
                resufr(im, 4) = resufr(im, 4)*xx1
                resufr(im, 5) = resufr(im, 5)*xx1
                resufr(im, 10) = resufr(im, 10)*xx2
                resufr(im, 11) = resufr(im, 11)*xx2
                resufr(im, 12) = resufr(im, 12)*xx2
                coef(im) = xx3
                do ie = 1, neq
                    vecpro(ie, im) = vecpro(ie, im)*xx3
                end do
            end do
        else
            AS_ALLOCATE(vr=poi1, size=neq)
            AS_ALLOCATE(vr=poi2, size=neq)
            do im = 1, nbmode
                do ie = 1, neq
                    poi1(ie) = vecpro(ie, im)
                end do
                call mrmult('ZERO', lmatr, poi1, poi2, 1, &
                            .true._1)
                xmn = 0.0d0
                do ie = 1, neq
                    xmn = xmn+(poi1(ie)*poi2(ie))
                end do
                xx1 = 1.0d0/sqrt(xmn)
                coef(im) = xx1
                do ie = 1, neq
                    vecpro(ie, im) = vecpro(ie, im)*xx1
                end do
            end do
            AS_DEALLOCATE(vr=poi2)
            AS_DEALLOCATE(vr=poi1)
        end if
!
    else
!
        valk = norm
        call utmess('F', 'ALGELINE4_77', sk=valk)
!
    end if
!
    if (para .eq. 'OUI') then
        do im = 1, nbmode
            resufr(im, 13) = 0.d0
            resufr(im, 14) = 0.d0
            resufr(im, 15) = 0.d0
            if (xmastr(1) .gt. epsi) resufr(im, 13) = resufr(im, 7)/xmastr(1)
            if (xmastr(2) .gt. epsi) resufr(im, 14) = resufr(im, 8)/xmastr(2)
            if (xmastr(3) .gt. epsi) resufr(im, 15) = resufr(im, 9)/xmastr(3)
        end do
    end if
!
    if (isign .eq. 0) then
    else if (isign .eq. 1) then
        do im = 1, nbmode
            xx1 = vecpro(numddl, im)
            if (xx1 .lt. 0.0d0) then
                coef(im) = -coef(im)
                do ie = 1, neq
                    vecpro(ie, im) = -vecpro(ie, im)
                end do
                if (para .eq. 'OUI') then
                    resufr(im, 10) = -resufr(im, 10)
                    resufr(im, 11) = -resufr(im, 11)
                    resufr(im, 12) = -resufr(im, 12)
                end if
            end if
        end do
    else if (isign .eq. -1) then
        do im = 1, nbmode
            xx1 = vecpro(numddl, im)
            if (xx1 .gt. 0.0d0) then
                coef(im) = -coef(im)
                do ie = 1, neq
                    vecpro(ie, im) = -vecpro(ie, im)
                end do
                if (para .eq. 'OUI') then
                    resufr(im, 10) = -resufr(im, 10)
                    resufr(im, 11) = -resufr(im, 11)
                    resufr(im, 12) = -resufr(im, 12)
                end if
            end if
        end do
    end if
!
    call jedema()
end subroutine
