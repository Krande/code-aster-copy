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
subroutine mxmoam(sddyna, nbmodp)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mginfo.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
!
    character(len=19) :: sddyna
    integer(kind=8) :: nbmodp
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EXPLICITE)
!
! PROJECTION MODALE EN EXPLICITE
!
! ----------------------------------------------------------------------
!
!
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE (CF NDLECT)
! OUT NBMODP : NOMBRE DE MODES DE PROJECTION
!
!
!
!
    integer(kind=8) :: nbmd, nbmg, neq, nbmax, nbrg, nbag
    integer(kind=8) :: nbgene
    integer(kind=8) :: iddeeq, jval
    integer(kind=8) :: ldblo, ldblo1, ldblo2
    integer(kind=8) :: imode, ifonc, imode2
    integer(kind=8) :: iret, ibid, nf, lpar, vali(3)
    character(len=8) :: k8bid
    character(len=8) :: modmec, magene, amgene, rigene
    character(len=14) :: numddl
    character(len=19) :: fmodal, valfon
    integer(kind=8) :: jfmoda, jvalfo
    character(len=19) :: depgem, vitgem, accgem
    integer(kind=8) :: jdepgm, jvitgm, jaccgm
    character(len=19) :: depgep, vitgep, accgep
    integer(kind=8) :: jdepgp, jvitgp, jaccgp
    character(len=19) :: basmod
    integer(kind=8) :: jbasmo
    character(len=19) :: riggen, masgen, amogen, fongen, forgen
    integer(kind=8) :: jrigge, jmasge, jamoge, jfonge, jforge
    character(len=19) :: accgcn
    integer(kind=8) :: jacccn
    character(len=24) :: deeq
    character(len=24) :: nomcha
    character(len=24), pointer :: lifoge(:) => null()
    real(kind=8), pointer :: fge(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOM DES OBJETS POUR PROJECTION MODALE
!
    call ndynkk(sddyna, 'PRMO_DEPGEM', depgem)
    call ndynkk(sddyna, 'PRMO_VITGEM', vitgem)
    call ndynkk(sddyna, 'PRMO_ACCGEM', accgem)
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
!
! --- MATRICE DES MODES MECA
!
    call getvid('PROJ_MODAL', 'MODE_MECA', iocc=1, scal=modmec, nbret=nbmd)
    if (nbmd .eq. 0) then
        ASSERT(.false.)
    end if
!
! --- MASSE, RIGIDITE ET AMORTISSEMENT GENERALISES
!
    call getvid('PROJ_MODAL', 'MASS_GENE', iocc=1, scal=magene, nbret=nbmg)
    call getvid('PROJ_MODAL', 'RIGI_GENE', iocc=1, scal=rigene, nbret=nbrg)
    call getvid('PROJ_MODAL', 'AMOR_GENE', iocc=1, scal=amgene, nbret=nbag)
!
! --- IL FAUT MASS_GENE _ET_ RIGI_GENE (VOIR CAPY)
!
    if ((nbmg .gt. 0) .and. (nbrg .eq. 0)) then
        ASSERT(.false.)
    end if
!
! --- INFORMATIONS SUR MATRICE DES MODES MECANIQUES
!
    call mginfo(modmec, numddl, nbmodp, neq)
    deeq = numddl//'.NUME.DEEQ'
    call jeveuo(deeq, 'L', iddeeq)
!
! --- NOMBRE DE MODES
!
    call getvis('PROJ_MODAL', 'NB_MODE', iocc=1, scal=nbmax, nbret=ibid)
!
    if (ibid .eq. 0) then
        nbmax = nbmodp
    end if
!
    if (nbmax .ne. nbmodp) then
        vali(1) = nbmodp
        vali(2) = nbmax
        vali(3) = min(nbmodp, nbmax)
        call utmess('I', 'MECANONLINE5_29', ni=3, vali=vali)
        nbmodp = min(nbmodp, nbmax)
    end if
!
! --- CREATION VECTEUR DES FORCES MODALES
!
    call wkvect(fmodal, 'V V R', nbmodp, jfmoda)
!
! --- CREATION MASSES GENERALISEES
!
    call wkvect(masgen, 'V V R', nbmodp, jmasge)
!
! --- CREATION BASE MODALE
!
    call wkvect(basmod, 'V V R', nbmodp*neq, jbasmo)
!
! --- SI MASS_GENE NON DONNE
!
    if (nbmg .eq. 0) then
!
! ---   ON RECUPERE MODES DANS MODE_MECA
!
        do imode = 1, nbmodp
            call rsexch('F', modmec, 'DEPL', imode, nomcha, &
                        iret)
            call jeveuo(nomcha(1:19)//'.VALE', 'L', jval)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jval), b_incx, zr(jbasmo+(imode-1)*neq), b_incy)
            call zerlag(neq, zi(iddeeq), vectr=zr(jbasmo+(imode-1)*neq))
        end do
!
! ---   ON RECUPERE MASSES GENERALISEES DANS MODE_MECA
!
        do imode = 1, nbmodp
            call rsadpa(modmec, 'L', 1, 'MASS_GENE', imode, &
                        0, sjv=lpar, styp=k8bid)
            zr(jmasge+imode-1) = zr(lpar)
        end do
!
! --- CREATION ACCELERATION DE REFERENCE
!
        call wkvect(accgcn, 'V V R', nbmodp, jacccn)
    else
!
! --- CREATION DEPL/VITE/ACCE GENERALISES T- ET T+
!
        call wkvect(accgem, 'V V R', nbmodp, jaccgm)
        call wkvect(accgep, 'V V R', nbmodp, jaccgp)
        call wkvect(vitgem, 'V V R', nbmodp, jvitgm)
        call wkvect(vitgep, 'V V R', nbmodp, jvitgp)
        call wkvect(depgem, 'V V R', nbmodp, jdepgm)
        call wkvect(depgep, 'V V R', nbmodp, jdepgp)
!
! --- ON RECUPERE MODES DANS MODE_MECA
!
        do imode = 1, nbmodp
            call rsexch('F', modmec, 'DEPL', imode, nomcha, &
                        iret)
            call jeveuo(nomcha(1:19)//'.VALE', 'L', jval)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jval), b_incx, zr(jbasmo+(imode-1)*neq), b_incy)
            call zerlag(neq, zi(iddeeq), vectr=zr(jbasmo+(imode-1)*neq))
        end do
!
! --- CREATION RIGIDITES GENERALISEES
!
        call wkvect(riggen, 'V V R', nbmodp*nbmodp, jrigge)
!
! --- CREATION AMORTISSEMENTS GENERALISES
!
        call wkvect(amogen, 'V V R', nbmodp*nbmodp, jamoge)
!
! --- CREATION FORCES/FONC_MULT GENERALISEES
!
        nbgene = ndynin(sddyna, 'NBRE_EXCIT_GENE')
        if (nbgene .ne. 0) then
            call wkvect(fongen, 'V V K24', nbgene, jfonge)
            AS_ALLOCATE(vk24=lifoge, size=nbgene)
            call wkvect(forgen, 'V V R', nbgene*nbmodp, jforge)
        end if
!
! --- CREATION VECTEUR DE RESOLUTION FORCES
!
        if (nbgene .ne. 0) then
            call wkvect(valfon, 'V V R', nbgene, jvalfo)
        end if
!
! --- RECUPERATION MASSES GENERALISEES
!
        call jeveuo(jexnum(magene//'           .VALM', 1), 'L', ldblo)
        do imode = 1, nbmodp
            zr(jmasge+imode-1) = zr(ldblo+imode-1)
        end do
!
! --- RECUPERATION RIGIDITES GENERALISEES
!
        call jeveuo(jexnum(rigene//'           .VALM', 1), 'L', ldblo1)
        do imode = 1, nbmodp
            do imode2 = 1, imode
                zr(jrigge+(imode-1)*nbmodp+imode2-1) = zr(ldblo1+imode*(imode-1)/2+imode2-1)
                zr(jrigge+(imode2-1)*nbmodp+imode-1) = zr(ldblo1+imode*(imode-1)/2+imode2-1)
            end do
        end do
!
! --- RECUPERATION AMORTISSEMENTS GENERALISES
!
        call jeveuo(jexnum(amgene//'           .VALM', 1), 'L', ldblo2)
        do imode = 1, nbmodp
            do imode2 = 1, imode
                if (nbag .ne. 0) then
                    zr(jamoge+(imode-1)*nbmodp+imode2-1) = zr(ldblo2+imode*(imode-1)/2+imode2-1)
                    zr(jamoge+(imode2-1)*nbmodp+imode-1) = zr(ldblo2+imode*(imode-1)/2+imode2-1)
                end if
            end do
        end do
!
! --- RECUPERATION FORCES/FONC_MULT GENERALISEES
!
        if (nbgene .ne. 0) then
            do ifonc = 1, nbgene
                call getvid('EXCIT_GENE', 'FONC_MULT', iocc=ifonc, scal=zk24(jfonge+ifonc-1), &
                            nbret=nf)
                call getvid('EXCIT_GENE', 'VECT_GENE', iocc=ifonc, scal=lifoge(ifonc), nbret=nf)
                call jeveuo(lifoge(ifonc) (1:19)//'.VALE', 'L', vr=fge)
                do imode = 1, nbmodp
                    zr(jforge+(ifonc-1)*nbmodp+imode-1) = fge(imode)
                end do
            end do
        end if
    end if
!
! --- MENAGE
!
    AS_DEALLOCATE(vk24=lifoge)
    call jedema()
end subroutine
