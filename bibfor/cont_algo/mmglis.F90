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

subroutine mmglis(ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES CONTINUES - ALGORITHME)
!
! GESTION DE LA GLISSIERE
!
! ----------------------------------------------------------------------
!
!     ON MET LE POINT EN GLISSIERE SI LGLISS=.TRUE. ET
!     SI LA CONVERGENCE EN CONTRAINTE ACTIVE EST ATTEINTE
! In  ds_contact       : datastructure for contact management
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ztabf
    character(len=24) :: tabfin
    integer(kind=8) :: jtabf
    integer(kind=8) :: nzoco, nbmae, nptm
    aster_logical :: lveri, lgliss
    integer(kind=8) :: izone, imae, iptc, iptm
    integer(kind=8) :: xs
    integer(kind=8) :: posmae, jdecme
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... GESTION GLISSIERE'
    end if
!
! --- ACCES SD CONTACT
!
    tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    call jeveuo(tabfin, 'L', jtabf)
    ztabf = cfmmvd('ZTABF')
!
! --- INITIALISATIONS
!
    nzoco = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    iptc = 1
!
! --- BOUCLE SUR LES ZONES
!
    do izone = 1, nzoco
!
! ----- MODE VERIF: ON SAUTE LES POINTS
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        if (lveri) then
            goto 25
        end if
!
! --- OPTIONS SUR LA ZONE DE CONTACT
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        nbmae = mminfi(ds_contact%sdcont_defi, 'NBMAE', izone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', izone)
        lgliss = mminfl(ds_contact%sdcont_defi, 'GLISSIERE_ZONE', izone)
!
! ----- BOUCLE SUR LES MAILLES ESCLAVES
!
        do imae = 1, nbmae
!
! ------- NUMERO ABSOLU DE LA MAILLE ESCLAVE
!
            posmae = jdecme+imae
!
! ------- NOMBRE DE POINTS SUR LA MAILLE ESCLAVE
!
            call mminfm(posmae, ds_contact%sdcont_defi, 'NPTM', nptm)
!
! ------- BOUCLE SUR LES POINTS
!
            if (lgliss) then
                do iptm = 1, nptm
                    xs = nint(zr(jtabf+ztabf*(iptc-1)+22))
                    if (xs .eq. 1) then
                        zr(jtabf+ztabf*(iptc-1)+17) = 1.d0
                    end if
                    iptc = iptc+1
                end do
            else
                iptc = iptc+nptm
            end if
        end do
25      continue
    end do
!
    call jedema()
end subroutine
