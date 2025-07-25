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
subroutine mrmtmv(cumul, lmat, smdi, smhc, lmatd, &
                  neq, neql, vect, xsol, nbvect, &
                  vectmp, prepos)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mrconl.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: cumul
    integer(kind=4) :: smhc(*)
    integer(kind=8) :: smdi(*), neq, nbvect, neql, lmat
    real(kind=8) :: vect(neq, nbvect), xsol(neq, nbvect), vectmp(neq)
    aster_logical :: lmatd, prepos
!                   MULTIPLICATION MATRICE PAR N VECTEURS
!         XSOL(1..NEQ,1..NBVECT) = MATRICE ^T * VECT(1..NEQ,1..NBVECT)
!     ------------------------------------------------------------------
! IN  CUMUL  : K4 :
!              / 'ZERO' : XSOL =        MAT ^T *VECT
!              / 'CUMU' : XSOL = XSOL + MAT ^T *VECT
!     ------------------------------------------------------------------
!
!
!
    real(kind=8) :: zero
    character(len=14) :: numddl
    character(len=19) :: nom19
    character(len=24) :: valm, refa, kxfem
    integer(kind=8) :: kfin, jvalms, jvalmi, jvec, ki, kdeb, nbloc
    integer(kind=8) :: ilig, jcol, jrefa, jnulg, iligg, jcolg, numglo, k
    integer(kind=8) :: keta, iexi, jccid, ieq
    aster_logical :: nonsym
!     ------------------------------------------------------------------
!
!
!
    call jemarq()
    nom19 = zk24(zi(lmat+1)) (1:19)
!
    call dismoi('XFEM', nom19, 'MATR_ASSE', repk=kxfem)
    if (kxfem .eq. 'XFEM_PRECOND') call utmess('A', 'XFEMPRECOND_4', nk=1, valk=nom19)
!
    valm = nom19//'.VALM'
    call jelira(valm, 'NMAXOC', nbloc)
    ASSERT(nbloc .eq. 1 .or. nbloc .eq. 2)
    nonsym = (nbloc .eq. 2)
!
    zero = 0.d0
    ASSERT(cumul .eq. 'ZERO' .or. cumul .eq. 'CUMU')
    if (cumul .eq. 'ZERO') then
        do jvec = 1, nbvect
            do ilig = 1, neq
                xsol(ilig, jvec) = zero
            end do
        end do
    end if
!
!
!     -- VALM(1) : AU DESSUS DE LA DIAGONALE
    call jeveuo(jexnum(valm, 1), 'L', jvalms)
    if (nonsym) then
!        -- VALM(2) : AU DESSOUS DE LA DIAGONALE
        call jeveuo(jexnum(valm, 2), 'L', jvalmi)
    else
        jvalmi = jvalms
    end if
!
!
!     -- CAS D'UNE MATRICE NON DISTRIBUEE :
!     ----------------------------------------
    if (.not. lmatd) then
        do jvec = 1, nbvect
            do k = 1, neq
                vectmp(k) = vect(k, jvec)
            end do
!         -- LES LAGRANGE DOIVENT ETRE MIS A L'ECHELLE AVANT LA
!            MULTIPLICATION :
            if (prepos) call mrconl('DIVI', lmat, 0, 'R', vectmp, &
                                    1)
            xsol(1, jvec) = xsol(1, jvec)+zr(jvalms-1+1)*vectmp(1)
            do ilig = 2, neq
                kdeb = smdi(ilig-1)+1
                kfin = smdi(ilig)-1
                do ki = kdeb, kfin
                    jcol = smhc(ki)
                    xsol(jcol, jvec) = xsol(jcol, jvec)+zr(jvalmi-1+ki)* &
                                       vectmp(ilig)
                    xsol(ilig, jvec) = xsol(ilig, jvec)+zr(jvalms-1+ki)* &
                                       vectmp(jcol)
                end do
                xsol(ilig, jvec) = xsol(ilig, jvec)+zr(jvalms+kfin)* &
                                   vectmp(ilig)
            end do
            if (prepos) call mrconl('DIVI', lmat, 0, 'R', xsol(1, jvec), &
                                    1)
        end do
!
!
!     -- CAS D'UNE MATRICE DISTRIBUEE :
!     ----------------------------------------
    else
        refa = nom19//'.REFA'
        call jeveuo(refa, 'L', jrefa)
        numddl = zk24(jrefa+2-1) (1:14)
        call jeveuo(numddl//'.NUML.NULG', 'L', jnulg)
        do jvec = 1, nbvect
            do k = 1, neq
                vectmp(k) = vect(k, jvec)
            end do
            if (prepos) call mrconl('DIVI', lmat, 0, 'R', vectmp, &
                                    1)
            numglo = zi(jnulg+1-1)
            xsol(numglo, jvec) = xsol(numglo, jvec)+zr(jvalms-1+1)* &
                                 vectmp(numglo)
            do ilig = 2, neql
                iligg = zi(jnulg+ilig-1)
                kdeb = smdi(ilig-1)+1
                kfin = smdi(ilig)-1
                do ki = kdeb, kfin
                    jcol = smhc(ki)
                    jcolg = zi(jnulg+jcol-1)
                    xsol(jcolg, jvec) = xsol(jcolg, jvec)+zr(jvalmi-1+ki) &
                                        *vectmp(iligg)
                    xsol(iligg, jvec) = xsol(iligg, jvec)+zr(jvalms-1+ki) &
                                        *vectmp(jcolg)
                end do
                xsol(iligg, jvec) = xsol(iligg, jvec)+zr(jvalms+kfin)* &
                                    vectmp(iligg)
            end do
            if (prepos) call mrconl('DIVI', lmat, 0, 'R', xsol(1, jvec), &
                                    1)
        end do
    end if
!
!
!     -- POUR LES DDLS ELIMINES PAR AFFE_CHAR_CINE, ON NE PEUT PAS
!        CALCULER F=K*U. CES DDLS SONT MIS A ZERO.
!     -------------------------------------------------------------
    call jeexin(nom19//'.CCID', iexi)
    if (iexi .ne. 0) then
        call jeveuo(nom19//'.CCID', 'L', jccid)
        do jvec = 1, nbvect
            do ieq = 1, neql
                if (lmatd) then
                    keta = zi(jccid-1+zi(jnulg+ieq-1))
                else
                    keta = zi(jccid-1+ieq)
                end if
                ASSERT(keta .eq. 1 .or. keta .eq. 0)
                if (keta .eq. 1) xsol(ieq, jvec) = 0.d0
            end do
        end do
    end if
!
!
!
    call jedema()
end subroutine
