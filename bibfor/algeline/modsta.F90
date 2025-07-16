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
subroutine modsta(motcle, matfac, matpre, solveu, lmatm, &
                  nume, iddl, coef, neq, nbmode, &
                  zrmod)
    implicit none
#include "jeveux.h"
#include "asterfort/ddllag.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/pteddl.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "asterfort/vecini.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: lmatm, iddl(*), neq, nbmode
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
!   Right-hand side are passed by batches of ICNTL(27) to MUMPS
    integer(kind=8) :: icmpl27
    real(kind=8) :: un, zero
    character(len=8) :: nomcmp(3)
    character(len=19) :: numeq
    character(len=24) :: kmumps_sparse(3)
!
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ib, ic, ie, ila1, ila2, im, imod, in, ila1o, nblag, jlag
    integer(kind=8) :: in2, ind, jddr, iret, nbatch, n_last_batch, batch_size, iaux
    integer(kind=8) :: coef_sparse, nzrhs_max, js1, js2, js3
    integer(kind=8), pointer :: position_ddl(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    aster_logical :: lverbose, lrhs_sparse
!
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    nomcmp = (/'DX', 'DY', 'DZ'/)
!     ------------------------------------------------------------------
    call jemarq()
! pour forcer le mode VERBOSE:
    lverbose = .true.
    lverbose = .false.
    un = 1.d0
    zero = 0.d0
    imod = 0
    ila1o = -9999
    numeq = nume(1:14)//'.NUME'
    call jeveuo(numeq//'.DEEQ', 'L', vi=deeq)
    call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
    lrhs_sparse = .false.
! taille de blocs de rhs traites en meme temps: valeur par defaut
    icmpl27 = min(32, nbmode)
    if (slvk(1) (1:5) .eq. 'MUMPS') then
! si MUMPS: taille de blocs de rhs traites en meme temps, valeur du .comm
        call jeveuo(solveu//'.SLVI', 'L', vi=slvi)
        if ((slvi(8) .lt. 0) .and. &
            ((motcle(1:4) .eq. 'DEPL') .or. (motcle(1:4) .eq. 'FORC'))) then
            lrhs_sparse = .true.
        end if
! au cas ou slvi(8) serait mal renseigne
        icmpl27 = min(max(abs(slvi(8)), 1), nbmode)
    end if
!
! Pour appel a MUMPS solve en mode sparse
    if (lrhs_sparse) then
        coef_sparse = 0
        if (motcle(1:4) .eq. 'DEPL') then
            coef_sparse = 2
        else if (motcle(1:4) .eq. 'FORC') then
            coef_sparse = 1
        end if
        nzrhs_max = nbmode*coef_sparse
        kmumps_sparse(1) = '&&MUMPS.SPARSE1'
        kmumps_sparse(2) = '&&MUMPS.SPARSE2'
        kmumps_sparse(3) = '&&MUMPS.SPARSE3'
! objet MUMPS RHS_SPARSE
        call jeexin(kmumps_sparse(1), iret)
        if (iret .ne. 0) call jedetr(kmumps_sparse(1))
        call wkvect(kmumps_sparse(1), 'V V R', nzrhs_max, js1)
        call vecini(nzrhs_max, un, zr(js1))
! objet MUMPS IRHS_SPARSE
        call jeexin(kmumps_sparse(2), iret)
        if (iret .ne. 0) call jedetr(kmumps_sparse(2))
        call wkvect(kmumps_sparse(2), 'V V I', nzrhs_max, js2)
! objet indice premier mode du patch + MUMPS IRHS_PTR
        call jeexin(kmumps_sparse(3), iret)
        if (iret .ne. 0) call jedetr(kmumps_sparse(3))
        call wkvect(kmumps_sparse(3), 'V V I', nbmode+2, js3)
        zi(js3+1) = 1
    end if
    if (lverbose) write (6, *) '<modsta entree> motcle/nbmode=', motcle, nbmode, lrhs_sparse
!
!***********************************
! Traitements locaux a chaque option
!***********************************
!
! Option ACCE*****************************************************************************
    if (motcle(1:4) .eq. 'ACCE') then
        AS_ALLOCATE(vi=position_ddl, size=3*neq)
        call pteddl('NUME_DDL', nume, 3, nomcmp, neq, &
                    tabl_equa=position_ddl)
        call wkvect('&&MODSTA.POSITION_DDR', 'V V R', neq, jddr)
        do im = 1, nbmode
            imod = imod+1
            in2 = 3*(im-1)
            call vecini(neq, zero, zr(jddr))
            do ic = 1, 3
                ind = neq*(ic-1)
                do in = 0, neq-1
                    zr(jddr+in) = zr(jddr+in)+position_ddl(1+ind+in)*coef(in2+ic)
                end do
            end do
            call mrmult('ZERO', lmatm, zr(jddr), zrmod(1, imod), 1, &
                        .true._1)
        end do
        call jedetr('&&MODSTA.POSITION_DDR')
        AS_DEALLOCATE(vi=position_ddl)
!
! Option DEPL*****************************************************************************
    else if (motcle(1:4) .eq. 'DEPL') then
        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
                call ddllag(nume, ie, neq, ila1, ila2, &
                            ila1o)
                ila1o = ila1
                if (ila1 .eq. 0 .or. ila2 .eq. 0) call utmess('F', 'ALGELINE2_4')
                zrmod(ila1, imod) = un
                zrmod(ila2, imod) = un
                if (lrhs_sparse) then
! Pour appel a MUMPS solve en mode sparse
                    zi(js2+2*(imod-1)) = ila1
                    zi(js2+2*(imod-1)+1) = ila2
                    zi(js3+1+imod) = 2*imod+1
                end if
            end if
        end do
!
! Option FORC*****************************************************************************
    else if (motcle(1:4) .eq. 'FORC') then
        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
                zrmod(ie, imod) = un
                if (lrhs_sparse) then
! Pour appel a MUMPS solve en mode sparse
                    zi(js2+imod-1) = ie
                    zi(js3+1+imod) = imod+1
                end if
            end if
        end do
!
    else
! Option ACCD*****************************************************************************

        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
            end if
        end do
        call wkvect('&&MODSTA.POSITION_DDR', 'V V R', neq*imod, jddr)
        imod = 0
        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
                call ddllag(nume, ie, neq, ila1, ila2, &
                            ila1o)
                ila1o = ila1
                if (ila1 .eq. 0 .or. ila2 .eq. 0) call utmess('F', 'ALGELINE2_4')
                !call vecini(neq, zero, zr(jddr))
                zr(jddr+(imod-1)*neq+ila1-1) = un
                zr(jddr+(imod-1)*neq+ila2-1) = un
            end if
        end do
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
            call resoud(matfac, matpre, solveu, ' ', batch_size, &
                        ' ', ' ', ' ', zr(jddr+neq*icmpl27*(ib-1)), [cbid], &
                        ' ', .true._1, 0, iret)
        end do
        imod = 0
        do ie = 1, neq
            if (iddl(ie) .eq. 1) then
                imod = imod+1
                call mrmult('ZERO', lmatm, zr(jddr+neq*(imod-1)), zrmod(1, imod), 1, &
                            .true._1)
            end if
        end do
        call jedetr('&&MODSTA.POSITION_DDR')
    end if
!
!
!*****************************************************************************************
! Traitement global: resolution par blocs (en mode dense ou sparse avec blocking si MUMPS)
!*****************************************************************************************
!
! Vecteur contenant les indices des Lagranges (au lieu de faire zerlag vecteur par vecteur)
    nblag = 0
    do ib = 1, neq
        if (deeq(2*ib) .le. 0) nblag = nblag+1
    end do
    if (nblag .ge. 1) then
        call wkvect('&&MODSTA.LAG', 'V V I', nblag, jlag)
        ie = 0
        do ib = 1, neq
            if (deeq(2*ib) .le. 0) then
                zi(jlag+ie) = ib
                ie = ie+1
            end if
        end do
    end if
!
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
!
            if (lrhs_sparse) then
! Pour appel a MUMPS solve en mode sparse: indice du premier rhs du patch
                zi(js3) = icmpl27*(ib-1)+1
                if (lverbose) then
                    write (6, *) '<modsta> ibatch/batch_size=', icmpl27*(ib-1)+1, batch_size
                    do im = 1, nzrhs_max
                        write (6, *) im, zr(js1-1+im), zi(js2-1+im)
                    end do
                    do im = 1, nbmode+2
                        write (6, *) im, zi(js3-1+im)
                    end do
                end if
            end if
!
            call resoud(matfac, matpre, solveu, ' ', batch_size, &
                        ' ', ' ', ' ', zrmod(1, icmpl27*(ib-1)+1), [cbid], &
                        ' ', .true._1, 0, iret)
!
        end do
    end if
!
!*****************************************************************************************
! Post-traitements
!*****************************************************************************************
! On annule les composantes de Lagranges des vecteurs solutions: vecteur par vecteur pour
! fluidifier les acces memoire.
    if (nblag .ge. 1) then
        do ib = 1, imod
            do ie = 1, nblag
                iaux = zi(jlag+ie-1)
                zrmod(iaux, ib) = zero
            end do
        end do
        call jedetr('&&MODSTA.LAG')
    end if
    call jedema()
end subroutine
