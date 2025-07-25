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

subroutine mmreas(mesh, ds_contact, hval_incr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnumm.h"
#include "asterfort/mmfield_prep.h"
#include "asterfort/detrsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mmelty.h"
#include "asterfort/mmextm.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
#include "asterfort/mmvalp_scal.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: hval_incr(*)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Update triggers for friction
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  hval_incr        : hat-variable for incremental values fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ztabf
    integer(kind=8) :: ibid
    integer(kind=8) :: posmae, jdecme, nummae
    integer(kind=8) :: iptc
    integer(kind=8) :: izone, imae, iptm
    integer(kind=8) :: nne, nbmae, nptm
    integer(kind=8) :: ndimg, nzoco
    aster_logical :: lveri
    real(kind=8) :: lambdc, ksipc1, ksipc2
    real(kind=8) :: mlagc(9)
    character(len=8) :: aliase
    character(len=19) :: cnslbd, depplu
    character(len=24) :: tabfin
    integer(kind=8) :: jtabf
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... MISE A JOUR DES SEUILS DE FROTTEMENT'
    end if
!
! --- INITIALISATIONS
!
    ndimg = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nzoco = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    ibid = 0
!
! --- RECUPERATION DES QCQS DONNEES
!
    tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    call jeveuo(tabfin, 'E', jtabf)
    ztabf = cfmmvd('ZTABF')
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', depplu)
!
! --- TRANSFORMATION DEPPLU EN CHAM_NO_S ET REDUCTION SUR LES LAGRANGES
!
    cnslbd = '&&REACLM.CNSLBD'
    !print*, 'dans mmreas'
    !call nmdebg('VECT', depplu, ifm)
    call mmfield_prep(depplu, cnslbd, &
                      l_sort_=.true._1, nb_cmp_=1, list_cmp_=['LAGS_C  '])
!
! --- BOUCLE SUR LES ZONES
!
    iptc = 1
    do izone = 1, nzoco
!
! --- OPTIONS SUR LA ZONE DE CONTACT
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        nbmae = mminfi(ds_contact%sdcont_defi, 'NBMAE', izone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', izone)
!
! ----- MODE VERIF: ON SAUTE LES POINTS
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        if (lveri) then
            goto 25
        end if
!
! ----- BOUCLE SUR LES MAILLES ESCLAVES
!
        do imae = 1, nbmae
!
! ------- NUMERO ABSOLU DE LA MAILLE ESCLAVE
!
            posmae = jdecme+imae
            call cfnumm(ds_contact%sdcont_defi, posmae, nummae)
!
! ------- INFOS SUR LA MAILLE
!
            call mmelty(mesh, nummae, aliase, nne)
!
! ------- MULTIPLICATEURS DE CONTACT SUR LES NOEUDS DE LA MAILLE ESCLAVE
!
            call mmextm(ds_contact%sdcont_defi, cnslbd, posmae, mlagc)
!
! ------- NOMBRE DE POINTS SUR LA MAILLE ESCLAVE
!
            call mminfm(posmae, ds_contact%sdcont_defi, 'NPTM', nptm)
!
! ------- BOUCLE SUR LES POINTS
!
            do iptm = 1, nptm
!
! --------- COORDONNEES ACTUALISEES DU POINT DE CONTACT
!
                ksipc1 = zr(jtabf+ztabf*(iptc-1)+3)
                ksipc2 = zr(jtabf+ztabf*(iptc-1)+4)
!
! --------- MULTIPLICATEUR DE LAGRANGE DE CONTACT DU POINT
!
                call mmvalp_scal(ndimg, aliase, nne, ksipc1, &
                                 ksipc2, mlagc, lambdc)
!
! --------- SAUVEGARDE
!
                zr(jtabf+ztabf*(iptc-1)+16) = lambdc
!
! --------- LIAISON DE CONTACT SUIVANTE
!
                iptc = iptc+1
            end do
        end do
25      continue
    end do
!
    call detrsd('CHAM_NO_S', cnslbd)
    call jedema()
end subroutine
