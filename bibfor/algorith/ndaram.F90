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
subroutine ndaram(result, sddyna, numarc)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/rsadpa.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    character(len=8) :: result
    integer(kind=8) :: numarc
    character(len=19) :: sddyna
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE)
!
! ARCHIVAGE DEPL/VITE/ACCE GENERALISES DANS SD_RESULT
!
! ----------------------------------------------------------------------
!
!
! IN  RESULT : NOM UTILISATEUR DU CONCEPT RESULTAT
! IN  NUMARC : NUMERO DE L'ARCHIVAGE
! IN  SDDYNA : SD DYNAMIQUE
!
!
!
!
    integer(kind=8) :: jpara
    character(len=8) :: k8bid
    character(len=6) :: chford
    character(len=19) :: dgen, vgen, agen
    integer(kind=8) :: jrestd, jrestv, jresta
    integer(kind=8) :: nbmodp
    character(len=24) :: trgene
    character(len=19) :: depgep, vitgep, accgep
    integer(kind=8) :: jdepgp, jvitgp, jaccgp
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOM TRAN_GENE_NOLI
!
    call codent(numarc, 'D0', chford)
    trgene = result(1:8)//'.TRG'//chford(1:6)
!
! --- ARCHIVAGE NOM TRAN_GENE_NOLI
!
    call rsadpa(result, 'E', 1, 'TRAN_GENE_NOLI', numarc, &
                0, sjv=jpara, styp=k8bid)
    zk24(jpara) = trgene
    nbmodp = ndynin(sddyna, 'NBRE_MODE_PROJ')
!
! --- NOM VECTEURS DEPL/VITE/ACCE GENERALISES
!
    dgen = trgene(1:18)//'D'
    vgen = trgene(1:18)//'V'
    agen = trgene(1:18)//'A'
!
! --- CREATION DEPL/VITE/ACCE GENERALISES
!
    call wkvect(dgen, 'G V R', nbmodp, jrestd)
    call wkvect(vgen, 'G V R', nbmodp, jrestv)
    call wkvect(agen, 'G V R', nbmodp, jresta)
!
! --- ARCHIVAGE
!
    call ndynkk(sddyna, 'PRMO_DEPGEP', depgep)
    call ndynkk(sddyna, 'PRMO_VITGEP', vitgep)
    call ndynkk(sddyna, 'PRMO_ACCGEP', accgep)
    nbmodp = ndynin(sddyna, 'NBRE_MODE_PROJ')
    call jeveuo(accgep, 'L', jaccgp)
    call jeveuo(vitgep, 'L', jvitgp)
    call jeveuo(depgep, 'L', jdepgp)
    b_n = to_blas_int(nbmodp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdepgp), b_incx, zr(jrestd), b_incy)
    b_n = to_blas_int(nbmodp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jvitgp), b_incx, zr(jrestv), b_incy)
    b_n = to_blas_int(nbmodp)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jaccgp), b_incx, zr(jresta), b_incy)
!
    call jedema()
!
end subroutine
