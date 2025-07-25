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
subroutine ascomb(lischa, vecelz, typres, nompar, valpar, &
                  cnchar)
!
!
    implicit none
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/corich.h"
#include "asterfort/fointc.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liscpp.h"
#include "asterfort/lislnf.h"
#include "asterfort/lisltf.h"
#include "asterfort/lisnnb.h"
#include "asterfort/vtcmbl.h"
#include "asterfort/wkvect.h"
    character(len=1) :: typres
    character(len=8) :: nompar
    real(kind=8) :: valpar, tval(1)
    character(len=19) :: lischa
    character(len=*) :: vecelz
    character(len=19) :: cnchar
!
! ----------------------------------------------------------------------
!
! COMBINER LES CHAM_NO
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : SD LISTE DES CHARGES
! IN  VECELE : NOM DU VECT_ELEM
! IN  TYPRES : TYPE DES VECTEURS ET DU CHAM_NO RESULTANT 'R' OU 'C'
! IN  NOMPAR : NOM DU PARAMETRE
! IN  VALPAR : VALEUR DU PARAMETRE
! OUT CNCHAR : CHAM_NO RESULTAT
!
!
!
!
    integer(kind=8) :: nbchar
    integer(kind=8) :: iret, ichar
    integer(kind=8) :: jcoef, jtype
    character(len=24) :: vachar
    integer(kind=8) :: ivec, ivecc, nbvec, jvacha
    character(len=8) :: nomfct
    character(len=24) :: chamno
    real(kind=8) :: valres, valre, valim
    complex(kind=8) :: calpha
    real(kind=8) :: phase, omega, dgrd
    integer(kind=8) :: npuis
    character(len=19) :: vecele
    character(len=16) :: typfct
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- VERIFICATION DU VACHAR
!
    vecele = vecelz
    vachar = vecele(1:19)//'.CHNO'
    call jeexin(vachar, iret)
    ASSERT(iret .ne. 0)
    call jelira(vachar, 'LONMAX', nbvec)
    ASSERT(nbvec .ne. 0)
    call jeveuo(vachar, 'L', jvacha)
    ASSERT(typres .eq. 'R' .or. typres .eq. 'C')
    call lisnnb(lischa, nbchar)
    dgrd = r8dgrd()
!
! --- CALCUL DES COEFFICIENTS - CAS REEL
!
    if (typres .eq. 'R') then
        call wkvect('&&ASCOMB.COEF', 'V V R8', nbvec, jcoef)
        call wkvect('&&ASCOMB.TYPE', 'V V K8', nbvec, jtype)
        do ivec = 1, nbvec
!
! ------- NOM DU CHAMNO
!
            chamno = zk24(jvacha+ivec-1)
!
! ------- NUMERO DE LA CHARGE
!
            call corich('L', chamno, ichout_=ichar)
            ASSERT((ichar .ne. 0) .and. (ichar .ge. -2))
!
! ------- FONCTION MULTIPLICATRICE
!
            if (ichar .gt. 0) then
                call lislnf(lischa, ichar, nomfct)
                call lisltf(lischa, ichar, typfct)
            end if
!
! ------- VALEUR DU COEFFICIENT
!
            if (ichar .eq. -1) then
                valres = 1.d0
            else if (ichar .eq. -2) then
                valres = 0.d0
            else if (ichar .gt. 0) then
                valres = 1.d0
                if (nomfct .ne. ' ') then
                    ASSERT(typfct(7:10) .eq. 'REEL')
                    tval(1) = valpar
                    call fointe('F', nomfct, 1, nompar, tval, &
                                valres, iret)
                end if
            else
                ASSERT(.false.)
            end if
!
            zr(jcoef+ivec-1) = valres
            zk8(jtype+ivec-1) = 'R'
        end do
    end if
!
! --- CALCUL DES COEFFICIENTS - CAS COMPLEXE
!
    if (typres .eq. 'C') then
        omega = r8depi()*valpar
        call wkvect('&&ASCOMB.COEF', 'V V R8', 2*nbvec, jcoef)
        call wkvect('&&ASCOMB.TYPE', 'V V K8', nbvec, jtype)
        ivecc = 0
        do ivec = 1, nbvec
!
! ------- NOM DU CHAMNO
!
            chamno = zk24(jvacha+ivec-1)
!
! ------- NUMERO DE LA CHARGE
!
            call corich('L', chamno, ichout_=ichar)
            ASSERT((ichar .ne. 0) .and. (ichar .ge. -2))
!
! ------- FONCTION MULTIPLICATRICE
!
            if (ichar .gt. 0) then
                call lislnf(lischa, ichar, nomfct)
                call lisltf(lischa, ichar, typfct)
            end if
!
! ------- MULTIPLICATEUR COMPLEXE
!
            if (ichar .gt. 0) then
                call liscpp(lischa, ichar, phase, npuis)
            end if
!
! ------- VALEUR DU COEFFICIENT
!
            if (ichar .eq. -1) then
                valre = 1.d0
                valim = 0.d0
                calpha = 1.d0
            else if (ichar .eq. -2) then
                valre = 0.d0
                valim = 0.d0
                calpha = 1.d0
            else if (ichar .gt. 0) then
                valre = 1.d0
                valim = 0.d0
                calpha = exp(dcmplx(0.d0, phase*dgrd))
                if (npuis .ne. 0) calpha = calpha*omega**npuis
                if (nomfct .ne. ' ') then
                    tval(1) = valpar
                    if (typfct(7:10) .eq. 'REEL') then
                        call fointe('F', nomfct, 1, nompar, tval, &
                                    valre, iret)
                        valim = 0.d0
                    else if (typfct(7:10) .eq. 'COMP') then
                        call fointc('F', nomfct, 1, nompar, tval, &
                                    valre, valim, iret)
                    else
                        ASSERT(.false.)
                    end if
                end if
            else
                ASSERT(.false.)
            end if
!
            zk8(jtype+ivec-1) = 'C'
            ivecc = ivecc+1
            zr(jcoef+ivecc-1) = valre*dble(calpha)-valim*dimag(calpha)
            ivecc = ivecc+1
            zr(jcoef+ivecc-1) = valim*dble(calpha)+valre*dimag(calpha)
!
        end do
    end if
!
! --- COMBINAISON LINEAIRE DES CHAM_NO
!
    call vtcmbl(nbvec, zk8(jtype), zr(jcoef), zk8(jtype), zk24(jvacha), &
                zk8(jtype), cnchar)
!
    call jedetr('&&ASCOMB.COEF')
    call jedetr('&&ASCOMB.TYPE')
!
    call jedema()
end subroutine
