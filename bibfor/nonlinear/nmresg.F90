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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmresg(numedd, sddyna, instap, cndonn, accsol)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdgeph.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmdebg.h"
#include "asterfort/vtzero.h"
#include "blas/ddot.h"
!
    real(kind=8) :: instap
    character(len=19) :: cndonn, sddyna
    character(len=24) :: numedd
    character(len=19) :: accsol
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! RESOLUTION SUR DDLS GENERALISES
!
! --------------------------------------------------------------------------------------------------
!
! IN  INSTAP : INSTANT COURANT
! IN  NUMEDD : NUME_DDL
! IN  CNDONN : CHAM_NO POUR LE SECOND MEMBRE
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE
! I/O ACCSOL : ACCELERATION CALCULEE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ier
    integer(kind=8) :: ifonc, imode, imode2
    integer(kind=8) :: neq, nbgene, nbmodp
    integer(kind=8) :: j2memb, jaccp, jaccg
    aster_logical :: lexge, lacce
    character(len=19) :: fmodal, valfon
    integer(kind=8) :: jfmoda, jvalfo
    character(len=19) :: depgep, vitgep, accgep
    integer(kind=8) :: jdepgp, jvitgp, jaccgp
    character(len=19) :: basmod, masgen, amogen, riggen
    integer(kind=8) :: jbasmo, jmasge, jamoge, jrigge
    character(len=19) :: fongen, forgen
    integer(kind=8) :: jfonge, jforge
    character(len=19) :: accgcn
    integer(kind=8) :: jacccn
    integer(kind=8) :: ifm, niv
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE><RESO> RESOLUTION SUR BASE MODALE'
    end if
!
! --- INITIALISATIONS
!
    call vtzero(accsol)
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- FONCTIONNALITES ACTIVEES
!
    lexge = ndynlo(sddyna, 'EXPL_GENE')
    lacce = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 3
    if (.not. lacce) then
        ASSERT(.false.)
    end if
!
! --- OBJETS PROJECTION MODALE
!
    call ndynkk(sddyna, 'PRMO_DEPGEP', depgep)
    call ndynkk(sddyna, 'PRMO_VITGEP', vitgep)
    call ndynkk(sddyna, 'PRMO_ACCGEP', accgep)
    call ndynkk(sddyna, 'PRMO_BASMOD', basmod)
    call ndynkk(sddyna, 'PRMO_MASGEN', masgen)
    call ndynkk(sddyna, 'PRMO_AMOGEN', amogen)
    call ndynkk(sddyna, 'PRMO_RIGGEN', riggen)
    call ndynkk(sddyna, 'PRMO_FONGEN', fongen)
    call ndynkk(sddyna, 'PRMO_FORGEN', forgen)
    call ndynkk(sddyna, 'PRMO_ACCGCN', accgcn)
    call ndynkk(sddyna, 'PRMO_VALFON', valfon)
    call ndynkk(sddyna, 'PRMO_FMODAL', fmodal)
    call jeveuo(masgen, 'L', jmasge)
    call jeveuo(basmod, 'L', jbasmo)
    call jeveuo(fmodal, 'E', jfmoda)
!
! --- NOMBRE DE MODES
!
    nbmodp = ndynin(sddyna, 'NBRE_MODE_PROJ')
!
! --- VECTEUR SECOND MEMBRE
!
    call jeveuo(cndonn(1:19)//'.VALE', 'E', j2memb)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE><RESO> -> SECOND MEMBRE DONNE'
        call nmdebg('VECT', cndonn, 6)
    end if
!
! --- RECUPERATION VECTEUR DES FORCES GENERALISEES
!
    if (lexge) then
        call jeveuo(riggen, 'L', jrigge)
        call jeveuo(amogen, 'L', jamoge)
        call jeveuo(vitgep, 'L', jvitgp)
        call jeveuo(depgep, 'L', jdepgp)
        call jeveuo(accgep, 'E', jaccgp)
!
! --- FORCES GENERALISEES ?
!
        nbgene = ndynin(sddyna, 'NBRE_EXCIT_GENE')
!
! --- EVALUATION DES FONCTIONS MULTIPLICATRICES
!
        if (nbgene .gt. 0) then
            call jeveuo(fongen, 'L', jfonge)
            call jeveuo(forgen, 'L', jforge)
            call jeveuo(valfon, 'E', jvalfo)
            do ifonc = 1, nbgene
                call fointe('F ', zk24(jfonge+ifonc-1) (1:8), 1, ['INST'], [instap], &
                            zr(jvalfo+ifonc-1), ier)
            end do
        end if
!
! --- CALCUL DES FORCES MODALES
!
        do imode = 1, nbmodp
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(jfmoda+imode-1) = ddot(b_n, zr(jbasmo+(imode-1)*neq), b_incx, zr(j2memb), b_incy)
            do imode2 = 1, nbmodp
                zr(jfmoda+imode-1) = zr(jfmoda+imode-1)-zr(jrigge+(imode2-1)*nbmodp+imode-1)*zr(j&
                                     &depgp+imode2-1)-zr(jamoge+(imode2-1)*nbmodp+imode-1)*zr(jvi&
                                     &tgp+imode2-1)
            end do
            do ifonc = 1, nbgene
                zr(jfmoda+imode-1) = zr(jfmoda+imode-1)+zr(jforge+(ifonc-1)*nbmodp+imode-1)*zr(jv&
                                     &alfo+ifonc-1)
            end do
        end do
!
! --- CALCUL DES ACCELERATIONS GENERALISEES
!
        do imode = 1, nbmodp
            zr(jaccgp+imode-1) = zr(jfmoda+imode-1)/zr(jmasge+imode-1)
        end do
!
        jaccg = jaccgp
!
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE><RESO> -> SOLUTION (MODALE):'
            call nmdebg(' ', accgep, 6)
        end if
    else
        call jeveuo(accgcn, 'E', jacccn)
!
! --- CALCUL DES FORCES GENERALISEES
!
        do imode = 1, nbmodp
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(jfmoda+imode-1) = ddot(b_n, zr(jbasmo+(imode-1)*neq), b_incx, zr(j2memb), b_incy)
        end do
!
! --- CALCUL DES ACCELERATIONS GENERALISEES
!
        do imode = 1, nbmodp
            zr(jacccn+imode-1) = zr(jfmoda+imode-1)/zr(jmasge+imode-1)
        end do
!
        jaccg = jacccn
!
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE><RESO> -> SOLUTION (MODALE):'
            call nmdebg(' ', accgcn, 6)
        end if
    end if
!
! --- CALCUL DES ACCELERATIONS PHYSIQUES
!
    call jeveuo(accsol(1:19)//'.VALE', 'E', jaccp)
    call mdgeph(neq, nbmodp, zr(jbasmo), zr(jaccg), zr(jaccp))
!
! --- AFFICHAGE DES SOLUTIONS
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE><RESO> -> SOLUTION (PHYSIQUE):'
        call nmdebg('VECT', accsol, 6)
    end if
!
    call jedema()
end subroutine
