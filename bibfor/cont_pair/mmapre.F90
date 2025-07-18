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
subroutine mmapre(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/apinfi.h"
#include "asterfort/apinfr.h"
#include "asterfort/apvect.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mmapma.h"
#include "asterfort/mmelty.h"
#include "asterfort/mmexcl.h"
#include "asterfort/mmimp1.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Save pairing in contact datastructures
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: sdappa
    integer(kind=8) :: nzoco, ndimg
    character(len=24) :: tabfin
    integer(kind=8) :: jtabf
    integer(kind=8) :: izone, ip, imae, iptm, ztabf
    integer(kind=8) :: iptc
    integer(kind=8) :: ntpc, nbpt, nbmae, nptm, nnomae
    real(kind=8) :: tau1m(3), tau2m(3), norm(3)
    real(kind=8) :: ksipr1, ksipr2
    character(len=8) :: aliase, nommam
    aster_logical :: lveri
    integer(kind=8) :: jdecme
    integer(kind=8) :: typint, typapp, entapp
    integer(kind=8) :: posmae, nummae, posmam, nummam
    aster_logical :: lappar, l_excl_frot, l_node_excl
    integer(kind=8) :: ndexfr
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... Save pairing in contact datastructures'
    end if
!
! - Pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! --- INITIALISATIONS
!
    iptc = 1
    ntpc = 0
!
! --- NOUVEL APPARIEMENT
!
    lappar = .true.
!
! --- PARAMETRES
!
    nzoco = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    ndimg = cfdisi(ds_contact%sdcont_defi, 'NDIM')
!
! --- ACCES SD CONTACT
!
    tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    call jeveuo(tabfin, 'E', jtabf)
    ztabf = cfmmvd('ZTABF')
!
! --- BOUCLE SUR LES ZONES
!
    ip = 1
    do izone = 1, nzoco
!
! ----- INFORMATION SUR LA ZONE
!
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', izone)
        nbmae = mminfi(ds_contact%sdcont_defi, 'NBMAE', izone)
        typint = mminfi(ds_contact%sdcont_defi, 'INTEGRATION', izone)
!
! ----- MODE VERIF: ON SAUTE LES POINTS
!
        lveri = mminfl(ds_contact%sdcont_defi, 'VERIF', izone)
        if (lveri) then
            nbpt = mminfi(ds_contact%sdcont_defi, 'NBPT', izone)
            ip = ip+nbpt
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
! ------- NOMBRE DE POINTS SUR LA MAILLE ESCLAVE
!
            call mminfm(posmae, ds_contact%sdcont_defi, 'NPTM', nptm)
!
! ------- INFOS SUR LA MAILLE ESCLAVE
!
            call mmelty(mesh, nummae, aliase, nnomae)
!
! ------- NOEUDS EXCLUS PAR SANS_GROUP_NO_FR OU SANS_NOEUD_FR
!
            call mminfm(posmae, ds_contact%sdcont_defi, 'NDEXFR', ndexfr)
!
! ------- BOUCLE SUR LES POINTS
!
            do iptm = 1, nptm
!
! --------- INFOS APPARIEMENT
!
                call apinfi(sdappa, 'APPARI_TYPE', ip, typapp)
                call apinfi(sdappa, 'APPARI_ENTITE', ip, entapp)
                call apinfr(sdappa, 'APPARI_PROJ_KSI1', ip, ksipr1)
                call apinfr(sdappa, 'APPARI_PROJ_KSI2', ip, ksipr2)
                call apvect(sdappa, 'APPARI_TAU1', ip, tau1m)
                call apvect(sdappa, 'APPARI_TAU2', ip, tau2m)
!
! --------- APPARIEMENT NODAL INTERDIT !
!
                if (typapp .eq. 1) then
                    ASSERT(.false.)
                end if
!
! ------------- Excluded nodes
!
                call mmexcl(typint, typapp, iptm, ndexfr, l_node_excl, &
                            l_excl_frot)
                zr(jtabf+ztabf*(iptc-1)+18) = 0.d0
                zr(jtabf+ztabf*(iptc-1)+19) = 0.d0
                if (l_node_excl) then
                    zr(jtabf+ztabf*(iptc-1)+18) = 1.d0
                end if
                if (l_excl_frot) then
                    zr(jtabf+ztabf*(iptc-1)+19) = ndexfr
                end if
!
! ------------- Excluded nodes => no contact !
!
                if (l_node_excl) then
                    zr(jtabf+ztabf*(iptc-1)+22) = 0.d0
                end if
!
! --------- NUMEROS DE LA MAILLE MAITRE
!
                posmam = entapp
                call cfnumm(ds_contact%sdcont_defi, posmam, nummam)
!
! --------- SAUVEGARDE APPARIEMENT
!
                call mmapma(mesh, ds_contact, ndimg, izone, l_excl_frot, &
                            typint, aliase, posmae, nummae, nnomae, &
                            posmam, nummam, ksipr1, ksipr2, tau1m, &
                            tau2m, iptm, iptc, norm, nommam)
!
! --------- LIAISON DE CONTACT EFFECTIVE
!
                iptc = iptc+1
                ntpc = ntpc+1
!
! --------- POINT SUIVANT
!
                ip = ip+1
!
            end do
        end do
25      continue
    end do
!
! --- NOMBRE TOTAL DE NOEUDS EN CONTACT
!
    zr(jtabf-1+1) = ntpc
    ASSERT(ntpc .eq. cfdisi(ds_contact%sdcont_defi, 'NTPC'))
!
! - Flag for (re) numbering
!
    ds_contact%l_renumber = lappar
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        call mmimp1(ifm, mesh, ds_contact)
    end if
!
    call jedema()
!
end subroutine
