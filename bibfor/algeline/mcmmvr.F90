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
subroutine mcmmvr(cumul, lmat, smdi, smhc, neq, &
                  vect, xsol, nbvect, vectmp, prepos)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mcconl.h"
!
    character(len=*) :: cumul
    integer(kind=4) :: smhc(*)
    integer(kind=8) :: smdi(*), neq, nbvect, lmat
    complex(kind=8) :: vect(neq, nbvect), xsol(neq, nbvect), vectmp(neq)
    aster_logical :: prepos
!                   MULTIPLICATION MATRICE PAR N VECTEURS
!         XSOL(1..NEQ,1..NBVECT) = MATRICE  * VECT(1..NEQ,1..NBVECT)
!     ------------------------------------------------------------------
!     VERSION : LA MATRICE EST REELLE SYMETRIQUE OU NON (MORSE)
!             : LES VECTEURS SONT COMPLEXES
!     ------------------------------------------------------------------
!
!
!
    integer(kind=8) :: jmat1, jmat2, nbloc, jcol, i, j, kdeb, kfin, ki, jvec, k
    character(len=19) :: nom19
    character(len=24) :: valm
    complex(kind=8) :: czero
    aster_logical :: nonsym
    integer(kind=8) :: keta, iexi, ieq
    integer(kind=8), pointer :: ccid(:) => null()
!     ------------------------------------------------------------------
!
!
!
    call jemarq()
    nom19 = zk24(zi(lmat+1))
    valm = nom19//'.VALM'
    call jelira(valm, 'NMAXOC', nbloc)
    ASSERT(nbloc .eq. 1 .or. nbloc .eq. 2)
    nonsym = (nbloc .eq. 2)
    czero = dcmplx(0.d0, 0.d0)
    if (cumul .eq. 'ZERO') then
        do i = 1, nbvect
            do j = 1, neq
                xsol(j, i) = czero
            end do
        end do
    end if
!     -- VALM(1) : AU DESSUS DE LA DIAGONALE
    call jeveuo(jexnum(valm, 1), 'L', jmat1)
    if (nonsym) then
!        -- VALM(2) : AU DESSOUS DE LA DIAGONALE
        call jeveuo(jexnum(valm, 2), 'L', jmat2)
    else
        jmat2 = jmat1
    end if
!
!
    do jvec = 1, nbvect
        do k = 1, neq
            vectmp(k) = vect(k, jvec)
        end do
!        -- LES LAGRANGE DOIVENT ETRE MIS A L'ECHELLE AVANT LA
!           MULTIPLICATION :
        if (prepos) call mcconl('DIVI', lmat, 0, 'C', vectmp, &
                                1)
!
!        -- PREMIERE LIGNE
        xsol(1, jvec) = xsol(1, jvec)+zr(jmat1-1+1)*vectmp(1)
!
!        -- LIGNES SUIVANTES
        do i = 2, neq
            kdeb = smdi(i-1)+1
            kfin = smdi(i)-1
            do ki = kdeb, kfin
                jcol = smhc(ki)
                xsol(jcol, jvec) = xsol(jcol, jvec)+zr(jmat1-1+ki)*vectmp( &
                                   i)
                xsol(i, jvec) = xsol(i, jvec)+zr(jmat2-1+ki)*vectmp(jcol)
            end do
            xsol(i, jvec) = xsol(i, jvec)+zr(jmat1+kfin)*vectmp(i)
        end do
        if (prepos) call mcconl('DIVI', lmat, 0, 'C', xsol(1, jvec), &
                                1)
    end do
!
!
!     -- POUR LES DDLS ELIMINES PAR AFFE_CHAR_CINE, ON NE PEUT PAS
!        CALCULER F=K*U. CES DDLS SONT MIS A ZERO.
!     -------------------------------------------------------------
    call jeexin(nom19//'.CCID', iexi)
    if (iexi .ne. 0) then
        call jeveuo(nom19//'.CCID', 'L', vi=ccid)
        do jvec = 1, nbvect
            do ieq = 1, neq
                keta = ccid(ieq)
                ASSERT(keta .eq. 1 .or. keta .eq. 0)
                if (keta .eq. 1) xsol(ieq, jvec) = dcmplx(0.d0, 0.d0)
            end do
        end do
    end if
!
!
    call jedema()
end subroutine
