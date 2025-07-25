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
subroutine ndmuap(numins, numedd, sddyna, sddisc)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/nmdebg.h"
#include "asterfort/r8inir.h"
!
    integer(kind=8) :: numins
    character(len=24) :: numedd
    character(len=19) :: sddyna, sddisc
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE)
!
! INITIALISATION DES CHAMPS D'ENTRAINEMENT EN MULTI-APPUI
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  NUMINS : NUMERO INSTANT COURANT
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  SDDYNA : SD DYNAMIQUE
!
!
!
!
    real(kind=8), parameter :: zero = 0.d0
!
    integer(kind=8) :: neq, ie, iex
    real(kind=8) :: instap
    real(kind=8) :: coef1, coef2, coef3
    character(len=19) :: depent, vitent, accent
    character(len=19) :: mafdep, mafvit, mafacc, mamula, mapsid
    integer(kind=8) :: jnodep, jnovit, jnoacc, jmltap, jpsdel
    integer(kind=8) :: nbexci
    integer(kind=8) :: ifm, niv
    real(kind=8), pointer :: accen(:) => null()
    real(kind=8), pointer :: depen(:) => null()
    real(kind=8), pointer :: viten(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... INIT. MULTI-APPUI'
    end if
!
! --- INITIALISATIONS
!
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    nbexci = ndynin(sddyna, 'NBRE_EXCIT')
    instap = diinst(sddisc, numins)
!
! --- ACCES SD MULTI-APPUI
!
    call ndynkk(sddyna, 'MUAP_MAFDEP', mafdep)
    call ndynkk(sddyna, 'MUAP_MAFVIT', mafvit)
    call ndynkk(sddyna, 'MUAP_MAFACC', mafacc)
    call ndynkk(sddyna, 'MUAP_MAMULA', mamula)
    call ndynkk(sddyna, 'MUAP_MAPSID', mapsid)
!
    call jeveuo(mafdep, 'L', jnodep)
    call jeveuo(mafvit, 'L', jnovit)
    call jeveuo(mafacc, 'L', jnoacc)
    call jeveuo(mamula, 'L', jmltap)
    call jeveuo(mapsid, 'L', jpsdel)
!
! --- ACCES DEPL/VITE/ACCE ENTRAINEMENT
!
    call ndynkk(sddyna, 'DEPENT', depent)
    call ndynkk(sddyna, 'VITENT', vitent)
    call ndynkk(sddyna, 'ACCENT', accent)
    call jeveuo(depent(1:19)//'.VALE', 'E', vr=depen)
    call jeveuo(vitent(1:19)//'.VALE', 'E', vr=viten)
    call jeveuo(accent(1:19)//'.VALE', 'E', vr=accen)
!
! --- EVALUATION DEPL/VITE/ACCE ENTRAINEMENT
!
    call r8inir(neq, zero, depen, 1)
    call r8inir(neq, zero, viten, 1)
    call r8inir(neq, zero, accen, 1)
!
    do iex = 1, nbexci
        if (zi(jmltap+iex-1) .eq. 1) then
            call fointe('F ', zk8(jnodep+iex-1), 1, ['INST'], [instap], &
                        coef1, ie)
            call fointe('F ', zk8(jnovit+iex-1), 1, ['INST'], [instap], &
                        coef2, ie)
            call fointe('F ', zk8(jnoacc+iex-1), 1, ['INST'], [instap], &
                        coef3, ie)
        else
            coef1 = zero
            coef2 = zero
            coef3 = zero
        end if
        do ie = 1, neq
            depen(ie) = depen(ie)+zr(jpsdel+(iex-1)*neq+ie-1)*coef1
            viten(ie) = viten(ie)+zr(jpsdel+(iex-1)*neq+ie-1)*coef2
            accen(ie) = accen(ie)+zr(jpsdel+(iex-1)*neq+ie-1)*coef3
        end do
    end do
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ...... DEPL. ENTRAINEMENT'
        call nmdebg('VECT', depent, ifm)
        write (ifm, *) '<MECANONLINE> ...... VITE. ENTRAINEMENT'
        call nmdebg('VECT', vitent, ifm)
        write (ifm, *) '<MECANONLINE> ...... ACCE. ENTRAINEMENT'
        call nmdebg('VECT', accent, ifm)
    end if
!
    call jedema()
!
end subroutine
