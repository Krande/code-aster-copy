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

subroutine apimpr_c(ifm, mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/apcopt.h"
#include "asterfort/apinfi.h"
#include "asterfort/apinfr.h"
#include "asterfort/apnomp.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfnumn.h"
#include "asterfort/cfdisi.h"
#include "asterfort/apvect.h"
#include "asterfort/mminfi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: ifm
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Node to segment - Debug print
!
! --------------------------------------------------------------------------------------------------
!
! In  ifm              : unit for message
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: sdappa
    integer(kind=8) :: nb_cont_zone, nt_poin, nb_poin_zone
    integer(kind=8) :: pair_type, pair_enti
    real(kind=8) :: poin_coor(3)
    real(kind=8) :: dist, ksi1, ksi2, tau1(3), tau2(3)
    character(len=16) :: poin_name
    integer(kind=8) :: i_zone, i_poin, k, i_poin_zone
    integer(kind=8) :: node_mast_nume(1), elem_mast_nume
    integer(kind=8) :: node_mast_indx(1), elem_mast_indx
    character(len=8) :: node_mast_name, elem_mast_name
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    i_poin = 1
!
! - Pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - Get parameters
!
    nt_poin = cfdisi(ds_contact%sdcont_defi, 'NTPT')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
!
! ----------------------------------------------------------------------
! --- INFOS SUR LES ZONES
! ----------------------------------------------------------------------
!
    write (ifm, *) '<APPARIEMENT> ------ ZONES ------ '
!
    write (ifm, 100) nb_cont_zone
    write (ifm, 101) nt_poin
!
100 format(' <APPARIEMENT> NOMBRE DE ZONES                   : ', i6)
101 format(' <APPARIEMENT> NOMBRE MAX. DE POINTS A APPARIER  : ', i6)
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Get parameters on current contact zone
!
        nb_poin_zone = mminfi(ds_contact%sdcont_defi, 'NBPT', i_zone)
!
! ----- Loop on points
!
        do i_poin_zone = 1, nb_poin_zone
!
! --------- Get parameters of current point
!
            call apnomp(sdappa, i_poin, poin_name)
            call apcopt(sdappa, i_poin, poin_coor)
!
! --------- Get parameters from pairing
!
            call apinfi(sdappa, 'APPARI_TYPE', i_poin, pair_type)
            call apinfi(sdappa, 'APPARI_ENTITE', i_poin, pair_enti)
            call apinfr(sdappa, 'APPARI_PROJ_KSI1', i_poin, ksi1)
            call apinfr(sdappa, 'APPARI_PROJ_KSI2', i_poin, ksi2)
            call apinfr(sdappa, 'APPARI_DIST', i_poin, dist)
            call apvect(sdappa, 'APPARI_TAU1', i_poin, tau1)
            call apvect(sdappa, 'APPARI_TAU2', i_poin, tau2)
!
! --------- Print parameters of current point
!
            write (ifm, 400) i_poin, poin_name
            if (pair_type .eq. -1) then
                write (ifm, 501)
            else if (pair_type .eq. -2) then
                write (ifm, 502)
            else if (pair_type .eq. -3) then
                write (ifm, 503)
            else if (pair_type .eq. 0) then
                write (ifm, 504)
            else if (pair_type .eq. 1) then
                write (ifm, 401) poin_coor(1), poin_coor(2), poin_coor(3)
                node_mast_indx = pair_enti
                call cfnumn(ds_contact%sdcont_defi, 1, node_mast_indx(1), node_mast_nume(1))
                node_mast_name = int_to_char8(node_mast_nume(1))
                write (ifm, 601) node_mast_name
                write (ifm, 801) dist
            else if (pair_type .eq. 2) then
                write (ifm, 401) poin_coor(1), poin_coor(2), poin_coor(3)
                elem_mast_indx = pair_enti
                call cfnumm(ds_contact%sdcont_defi, elem_mast_indx, elem_mast_nume)
                elem_mast_name = int_to_char8(elem_mast_nume)
                write (ifm, 602) elem_mast_name
                write (ifm, 701) ksi1, ksi2
                write (ifm, 801) dist
                write (ifm, 901) (tau1(k), k=1, 3)
                write (ifm, 902) (tau2(k), k=1, 3)
            else
                write (ifm, 504)
            end if
!
! --------- Next point
!
            i_poin = i_poin+1
        end do
    end do
!
400 format(' <APPARIEMENT> POINT            ', i6, ' (', a16, ')')
401 format(' <APPARIEMENT> ** DE COORDONNEES ', 1pe15.8, 1pe15.8, 1pe15.8)
!
501 format(' <APPARIEMENT> -> EXCLU - PAR SANS_NOEUD')
502 format(' <APPARIEMENT> -> EXCLU - PAR TOLE_APPA')
503 format(' <APPARIEMENT> -> EXCLU - PAR TOLE_PROJ_EXT')
504 format(' <APPARIEMENT> -> NON APPARIE (ERREUR)')
!
!
601 format(' <APPARIEMENT> -> APPARIEMENT AVEC NOEUD  ', a8)
602 format(' <APPARIEMENT> -> APPARIEMENT AVEC MAILLE ', a8)
!
701 format(' <APPARIEMENT>      SUR POINT KSI1,KSI2: ', 1pe15.8, 1pe15.8)
801 format(' <APPARIEMENT>      DISTANCE: ', 1pe15.8)
901 format(' <APPARIEMENT>      TANGENTE BRUTE  DIR. 1   : ', 3(1pe15.8, 2x))
902 format(' <APPARIEMENT>      TANGENTE BRUTE  DIR. 2   : ', 3(1pe15.8, 2x))
!
    call jedema()
!
end subroutine
