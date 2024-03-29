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

subroutine xreacl(mesh, model, hval_incr, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/copisd.h"
#include "asterfort/dbgcal.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/xmchex.h"
!
! person_in_charge: samuel.geniaut at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: hval_incr(*)
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - ALGORITHME)
!
! MISE À JOUR DU SEUIL DE FROTTEMENT
!
! ----------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  ds_contact       : datastructure for contact management
! In  hval_incr        : hat-variable for incremental values fields
!
!
!
    integer, parameter :: nbout = 1
    integer, parameter :: nbin = 10
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=19) :: ligrmo, xdonco, xseuco, cseuil
    character(len=19) :: lnno, ltno
    character(len=16) :: option
    character(len=24) :: ainter, cface, faclon, pinter, chgeom, baseco
    character(len=19) :: depplu
    aster_logical :: debug, lcontx
    integer :: ifm, niv, ifmdbg, nivdbg
    integer, pointer :: xfem_cont(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! --- INITIALISATIONS
!
    ligrmo = model(1:8)//'.MODELE'
    cseuil = '&&XREACL.SEUIL'
    xdonco = ds_contact%sdcont_solv(1:14)//'.XFDO'
    xseuco = ds_contact%sdcont_solv(1:14)//'.XFSE'
    option = 'XREACL'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', depplu)
!
! --- SI PAS DE CONTACT ALORS ON ZAPPE LA VÉRIFICATION
!
    call jeveuo(model(1:8)//'.XFEM_CONT', 'L', vi=xfem_cont)
    lcontx = xfem_cont(1) .ge. 1
    if (.not. lcontx) then
        goto 999
    end if
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! --- RECUPERATION DES DONNEES XFEM
!
    lnno = model(1:8)//'.LNNO'
    ltno = model(1:8)//'.LTNO'
    ainter = model(1:8)//'.TOPOFAC.AI'
    cface = model(1:8)//'.TOPOFAC.CF'
    faclon = model(1:8)//'.TOPOFAC.LO'
    pinter = model(1:8)//'.TOPOFAC.OE'
    baseco = model(1:8)//'.TOPOFAC.BA'
!
! --- CREATION DU CHAM_ELEM_S VIERGE
!
    call xmchex(mesh, xseuco, cseuil)
!
! --- RECUPERATION DES COORDONNEES DES NOEUDS
!
    chgeom = mesh(1:8)//'.COORDO'
!
! --- CREATION DES LISTES DES CHAMPS IN
!
    lpain(1) = 'PDEPL_P'
    lchin(1) = depplu(1:19)
    lpain(2) = 'PAINTER'
    lchin(2) = ainter(1:19)
    lpain(3) = 'PCFACE'
    lchin(3) = cface(1:19)
    lpain(4) = 'PLONGCO'
    lchin(4) = faclon(1:19)
    lpain(5) = 'PDONCO'
    lchin(5) = xdonco(1:19)
    lpain(6) = 'PPINTER'
    lchin(6) = pinter(1:19)
    lpain(7) = 'PGEOMER'
    lchin(7) = chgeom(1:19)
    lpain(8) = 'PLSN'
    lchin(8) = lnno(1:19)
    lpain(9) = 'PLST'
    lchin(9) = ltno(1:19)
    lpain(10) = 'PBASECO'
    lchin(10) = baseco(1:19)
!
!
! --- CREATION DES LISTES DES CHAMPS OUT
!
    lpaout(1) = 'PSEUIL'
    lchout(1) = cseuil(1:19)
!
! --- APPEL A CALCUL
!
    call calcul('S', option, ligrmo, nbin, lchin, &
                lpain, nbout, lchout, lpaout, 'V', &
                'OUI')
!
    if (debug) then
        call dbgcal(option, ifmdbg, nbin, lpain, lchin, &
                    nbout, lpaout, lchout)
    end if
!
! --- ON COPIE CSEUIL DANS RESOCO.SE
!
    call copisd('CHAMP_GD', 'V', lchout(1), xseuco)
!
999 continue
!
    call jedema()
end subroutine
