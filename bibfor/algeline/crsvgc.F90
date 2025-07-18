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

subroutine crsvgc(motfac, solveu, kellag)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/gcncon.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"

    character(len=3) :: kellag
    character(len=16) :: motfac
    character(len=19) :: solveu
!  BUT : REMPLISSAGE SD_SOLVEUR GCPC
!
! IN K19 SOLVEU  : NOM DU SOLVEUR DONNE EN ENTREE
! OUT    SOLVEU  : LE SOLVEUR EST CREE ET INSTANCIE
! IN  K3 KELLAG  : ELIM_LAGR
! ----------------------------------------------------------
!
!
!
!
    integer(kind=8) :: ibid, nmaxit, niremp, reacpr, pcpiv, redmpi
    real(kind=8) :: resire, blreps
    character(len=8) :: precon
    character(len=19) :: solvbd
    character(len=24) :: usersm
    character(len=8) :: renum
    real(kind=8), pointer :: slvr(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
!
!------------------------------------------------------------------
    call jemarq()
!
! --- LECTURES PARAMETRES DEDIES AU SOLVEUR
    call getvtx(motfac, 'PRE_COND', iocc=1, scal=precon, nbret=ibid)
    ASSERT(ibid .eq. 1)
    call getvtx(motfac, 'RENUM', iocc=1, scal=renum, nbret=ibid)
    ASSERT(ibid .eq. 1)
    call getvis(motfac, 'NMAX_ITER', iocc=1, scal=nmaxit, nbret=ibid)
    ASSERT(ibid .eq. 1)
    call getvr8(motfac, 'RESI_RELA', iocc=1, scal=resire, nbret=ibid)
    ASSERT(ibid .eq. 1)
!
!
! --- LECTURES PARAMETRES LIES A LDLT_INC
!
!     -- INITIALISATION
    niremp = -9999
    reacpr = -9999
    pcpiv = -9999
    blreps = 0.d0
    solvbd = 'XXXXXXXXXXXXXXXXXXX'
!
!     -- LECTURE
    if (precon .eq. 'LDLT_INC') then
        call getvis(motfac, 'NIVE_REMPLISSAGE', iocc=1, scal=niremp, nbret=ibid)
        ASSERT(ibid .eq. 1)
    else if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
        call getvis(motfac, 'REAC_PRECOND', iocc=1, scal=reacpr, nbret=ibid)
        ASSERT(ibid .eq. 1)
        call getvis(motfac, 'PCENT_PIVOT', iocc=1, scal=pcpiv, nbret=ibid)
        ASSERT(ibid .eq. 1)
        call getvtx(motfac, 'GESTION_MEMOIRE', iocc=1, scal=usersm, nbret=ibid)
        ASSERT(ibid .eq. 1)
        call getvr8(motfac, 'LOW_RANK_SEUIL', iocc=1, scal=blreps, nbret=ibid)
        ASSERT(ibid .eq. 1)
        redmpi = -9999
! a finir de tester avant restitution
!        call getvis(motfac, 'REDUCTION_MPI', iocc=1, scal=redmpi, nbret=ibid)
!
!       NOM DE SD SOLVEUR BIDON QUI SERA PASSEE A MUMPS
!       POUR LE PRECONDITIONNEMENT
        call gcncon('.', solvbd)
!
    else
        ASSERT(.false.)
    end if
!
! --- ON REMPLIT LA SD_SOLVEUR
    call jeveuo(solveu//'.SLVK', 'E', vk24=slvk)
    call jeveuo(solveu//'.SLVR', 'E', vr=slvr)
    call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
!
    slvk(1) = 'GCPC'
    slvk(2) = precon
    slvk(3) = solvbd
    slvk(4) = renum
    slvk(5) = 'XXXX'
    slvk(6) = 'XXXX'
    slvk(7) = 'XXXX'
    slvk(8) = 'XXXX'
    slvk(9) = usersm
    slvk(10) = 'XXXX'
    slvk(11) = 'XXXX'
    slvk(12) = 'XXXX'
    slvk(13) = kellag
    slvk(14) = 'XXXX'
!
!     POUR NEWTON_KRYLOV LE RESI_RELA VARIE A CHAQUE
!     ITERATION DE NEWTON, CEPENDANT LE RESI_RELA DONNE
!     PAR L'UTILISATEUR TOUT DE MEME NECESSAIRE
!     C'EST POURQUOI ON EN FAIT UNE COPIE EN POSITION 1
    slvr(1) = resire
    slvr(2) = resire
    slvr(3) = 0.d0
    slvr(4) = blreps
    slvr(5) = 0.d0
!
    slvi(1) = redmpi
    slvi(2) = nmaxit
    slvi(3) = -9999
    slvi(4) = niremp
    slvi(5) = 0
    slvi(6) = reacpr
    slvi(7) = pcpiv
    slvi(8) = 0
!
!
    call jedema()
end subroutine
