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

subroutine nmcerr(sddisc, iter_glob_maxi, iter_glob_elas, pas_mini_elas, resi_glob_maxi, &
                  resi_glob_rela, inikry, ds_contact_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/ceil.h"
#include "asterfort/cfdisi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmlerr.h"
#include "asterfort/utdidt.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iter_glob_maxi
    integer(kind=8), intent(in) :: iter_glob_elas
    real(kind=8), intent(in) :: pas_mini_elas
    real(kind=8), intent(in) :: inikry
    real(kind=8), intent(in) :: resi_glob_maxi
    real(kind=8), intent(in) :: resi_glob_rela
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (STRUCTURES DE DONNES - SD DISCRETISATION)
!
! CREATION SD STOCKAGE DES INFOS EN COURS DE CALCUL
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! IN  iter_glob_maxi   : ITER_GLOB_MAXI
! IN  iter_glob_rela   : ITER_GLOB_ELAS
! IN  pas_mini_elas    : PAS_MINI_ELAS
! IN  INIKRY           : PRECISION INITIALE POUR NEWTON-KRYLOV
! IN  resi_glob_maxi   : RESI_GLOB_MAXI
! IN  resi_glob_rela   : RESI_GLOB_RELA
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: r8bid
    real(kind=8) :: pcplus
    integer(kind=8) :: typres, ibid, nbiter, iter_maxi, iter_mini
    integer(kind=8) :: maxgeo, nbreag
    integer(kind=8) :: nmax, nplus
    integer(kind=8) :: i_echec, nb_echec, itesup, nbitct
    character(len=24) :: infocv, infore
    integer(kind=8) :: jifcv, jifre
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    iter_maxi = max(iter_glob_maxi, iter_glob_elas)
    iter_mini = min(iter_glob_maxi, iter_glob_elas)
    itesup = 0
    nmax = 0
    call utdidt('L', sddisc, 'LIST', 'NECHEC', &
                vali_=nb_echec)
!
! --- NOMBRE D'ITERATIONS AUTORISEES EN PLUS
!
    do i_echec = 1, nb_echec
        call utdidt('L', sddisc, 'ECHE', 'PCENT_ITER_PLUS', index_=i_echec, &
                    valr_=pcplus)
        nplus = nint(pcplus)
        nmax = max(nmax, nplus)
    end do
!
! --- NOMBRE MAXIMUM D'ITERATIONS
!
    nbiter = ceil(iter_maxi*(1.d0+nplus/100.0d0))
!
! --- CREATION DU VECTEUR D'INFORMATIONS SUR LA CONVERGENCE
!
    infocv = sddisc(1:19)//'.IFCV'
    call wkvect(infocv, 'V V R', 10, jifcv)
!
! --- MXITER      : MAX( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
! --- MNITER      : MIN( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
! --- NBITER      : MXITER + ITERATION EN PLUS
! --- ITESUP      :  =1 ON AUTORISE DES ITERATIONS EN PLUS
! ---                <>1 ON N'AUTORISE PAS D'ITERATIONS EN PLUS
! --- RESIDUS     : RELA ET MAXI
!
    call nmlerr(sddisc, 'E', 'MXITER', r8bid, iter_maxi)
    call nmlerr(sddisc, 'E', 'MNITER', r8bid, iter_mini)
    call nmlerr(sddisc, 'E', 'NBITER', r8bid, nbiter)
    call nmlerr(sddisc, 'E', 'PAS_MINI_ELAS', pas_mini_elas, ibid)
    call nmlerr(sddisc, 'E', 'ITERSUP', r8bid, itesup)
    call nmlerr(sddisc, 'E', 'RESI_GLOB_RELA', resi_glob_rela, ibid)
    call nmlerr(sddisc, 'E', 'RESI_GLOB_MAXI', resi_glob_maxi, ibid)
!
! --- TYPE_RESI   :  =1 ON A DONNE RESI_GLOB_RELA
! ---                =2 ON A DONNE RESI_GLOB_MAXI
! ---                =3 ON A DONNE RESI_GLOB_RELA ET RESI_GLOB_MAXI
! ---                =0 ON A RIEN DONNE ==> =1
!
    typres = 0
    if (resi_glob_rela .ne. r8vide()) then
        typres = typres+1
    end if
    if (resi_glob_maxi .ne. r8vide()) then
        typres = typres+2
    end if
    if (typres .eq. 0) then
        typres = 1
    end if
    call nmlerr(sddisc, 'E', 'TYPE_RESI', r8bid, typres)
!
! --- RESIDU INITIAL POUR NEWTON-KRYLOV
!
    call nmlerr(sddisc, 'E', 'INIT_NEWTON_KRYLOV', inikry, ibid)
!
! --- RESIDU COURANT POUR NEWTON-KRYLOV
!
    call nmlerr(sddisc, 'E', 'ITER_NEWTON_KRYLOV', inikry, ibid)
!
! --- RECUPERATION NOMBRE DE REAC_GEOM EN CONTACT DISCRET
!
    nbitct = 0
    if (present(ds_contact_)) then
        if (ds_contact_%l_meca_cont) then
            maxgeo = cfdisi(ds_contact_%sdcont_defi, 'ITER_GEOM_MAXI')
            nbreag = cfdisi(ds_contact_%sdcont_defi, 'NB_ITER_GEOM')
            nbitct = max(maxgeo, nbreag)
        end if
    end if
    nbiter = nbiter+nbitct
!
! --- CREATION DU VECTEUR DE STOCKAGE DES RESIDUS
!
    nbiter = nbiter+1
    infore = sddisc(1:19)//'.IFRE'
    call wkvect(infore, 'V V R', 3*nbiter, jifre)
!
    call jedema()
end subroutine
