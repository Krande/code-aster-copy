! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine modsta(motcle, matfac, matpre, solveu, lmatm, &
                  nume, iddl, coef, neq, nbmode, &
                  zrmod)
    implicit none
#include "jeveux.h"
#include "asterfort/ddllag.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/pteddl.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer :: lmatm, iddl(*), neq, nbmode
    real(kind=8) :: coef(*), zrmod(neq, *)
    character(len=*) :: motcle, nume, matfac, matpre, solveu
    complex(kind=8) :: cbid
!
!     CALCUL DE MODES STATIQUES
!
!     SI MOTCLE = 'DEPL' : CALCUL DE MODES CONTRAINTS
!                                    ( DEPLACEMENT UNITAIRE )
!                          LE TABLEAU IDDL EST CELUI DES NOEUDS BLOQUES
!                          ON APPLIQUE UNE FORCE UNITAIRE SUR LES LAGR
!     SI MOTCLE = 'FORC' : CALCUL DE MODES D'ATTACHE
!                                    ( FORCE UNITAIRE )
!                          LE TABLEAU IDDL EST CELUI DES NOEUDS ACTIFS
!     SI MOTCLE = 'ACCE' : CALCUL DE DEFORMEES STATIQUES
!                                    ( ACCELERATION VECTEUR UNITAIRE )
!     SI MOTCLE = 'ACCD' : CALCUL DE DEFORMEES STATIQUES
!                                    ( ACCELERATION DDL UNITAIRE )
!-----------------------------------------------------------------------
!  IN  : MOTCLE : CALCUL DE MODES CONTRAINTS OU D'ATTACHE
!  IN  : MATFAC : MATRICE DE RAIDEUR FACTORISEE
!  IN  : MATPRE : MATRICE DE PRECONDIONNEMENT POUR LA RAIDEUR (GCPC)
!  IN  : LMATM  : POINTEUR SUR LE DESCRIPTEUR DE LA MATRICE DE MASSE
!  IN  : NUME   : NOM DU NUME_DDL
!  IN  : IDDL   : TABLEAU DES DDL
!                 IDDL(I) = 0  PAS DE CALCUL DU MODE
!                 IDDL(I) = 1  CALCUL DU MODE
!  IN  : COEF   : COEFFICIENTS A APPLIQUER
!  IN  : NEQ    : NOMBRE D'EQUATIONS DU NUME
!  IN  : NBMODE : NOMBRE DE MODES STATIQUES
!  OUT : ZRMOD  : TABLEAU DES MODES STATIQUES CALCULES
!-----------------------------------------------------------------------
!
!   Right-hand side are passed by batches of ICMPL(27) to MUMPS
!   Whose value is currently 8
    integer, parameter :: icmpl27 = 8
    real(kind=8) :: zrbuff(neq, icmpl27)
    real(kind=8) :: un
    character(len=8) :: nomcmp(3)
    character(len=19) :: numeq
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer :: batch_size
    integer :: ib, ic, ie, ila1, ila2, im, imod, in
    integer :: in2, ind, jddr
    integer :: iret
    integer :: nbatch, n_last_batch
    integer, pointer :: position_ddl(:) => null()
    integer, pointer :: deeq(:) => null()
!
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    nomcmp = (/'DX', 'DY', 'DZ'/)
!     ------------------------------------------------------------------
    call jemarq()
    un = 1.d0
    imod = 0
!
    numeq = nume(1:14)//'.NUME'
    call jeveuo(numeq//'.DEEQ', 'L', vi=deeq)
!
    if (motcle(1:4) .eq. 'ACCE') then
        AS_ALLOCATE(vi=position_ddl, size=3*neq)
        call pteddl('NUME_DDL', nume, 3, nomcmp, neq, &
                    tabl_equa=position_ddl)
        do im = 1, nbmode
            imod = imod+1
            in2 = 3*(im-1)
            call wkvect('&&MODSTA.POSITION_DDR', 'V V R', neq, jddr)
            do ic = 1, 3
                ind = neq*(ic-1)
                do in = 0, neq-1
                    zr(jddr+in) = zr(jddr+in)+position_ddl(1+ind+in)*coef(in2+ic)
                end do
            end do
            call mrmult('ZERO', lmatm, zr(jddr), zrmod(1, imod), 1, &
                        .true._1)
            call jedetr('&&MODSTA.POSITION_DDR')
!
        end do
        AS_DEALLOCATE(vi=position_ddl)
    else
        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
                if (motcle(1:4) .eq. 'DEPL') then
                    call ddllag(nume, ie, neq, ila1, ila2)
                    if (ila1 .eq. 0 .or. ila2 .eq. 0) then
                        call utmess('F', 'ALGELINE2_4')
                    end if
                    zrmod(ila1, imod) = un
                    zrmod(ila2, imod) = un
                else if (motcle(1:4) .eq. 'FORC') then
                    zrmod(ie, imod) = un
                else
                    call wkvect('&&MODSTA.POSITION_DDR', 'V V R', neq, jddr)
                    call ddllag(nume, ie, neq, ila1, ila2)
                    if (ila1 .eq. 0 .or. ila2 .eq. 0) then
                        call utmess('F', 'ALGELINE2_4')
                    end if
                    zr(jddr+ila1-1) = un
                    zr(jddr+ila2-1) = un
                    call resoud(matfac, matpre, solveu, ' ', 1, &
                                ' ', ' ', ' ', zr(jddr), [cbid], &
                                ' ', .true._1, 0, iret)
                    call mrmult('ZERO', lmatm, zr(jddr), zrmod(1, imod), 1, &
                                .true._1)
                    call jedetr('&&MODSTA.POSITION_DDR')
                end if
            end if
        end do
    end if
!
!     --- RESOLUTION, BY BATCH OF ICMPL(27), SEE ISSUE32438 ---
    if (imod .gt. 0) then
        n_last_batch = mod(imod, icmpl27)
        if (n_last_batch == 0) then
            nbatch = imod/icmpl27
            n_last_batch = icmpl27
        else
            nbatch = imod/icmpl27+1
        end if
        do ib = 1, nbatch
            if (ib == nbatch) batch_size = n_last_batch
            if (ib /= nbatch) batch_size = icmpl27
            do in = 1, batch_size
                do ie = 1, neq
                    zrbuff(ie, in) = zrmod(ie, icmpl27*(ib-1)+in)
                end do
            end do
            call resoud(matfac, matpre, solveu, ' ', batch_size, &
                        ' ', ' ', ' ', zrbuff, [cbid], &
                        ' ', .true._1, 0, iret)
            do in = 1, batch_size
                ! Setting Lagrange dofs to 0
                call zerlag(neq, deeq, vectr=zrbuff(1, in))
                do ie = 1, neq
                    zrmod(ie, icmpl27*(ib-1)+in) = zrbuff(ie, in)
                end do
            end do
        end do
    end if
    call jedema()
end subroutine
