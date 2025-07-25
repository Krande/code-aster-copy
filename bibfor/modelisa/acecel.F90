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

subroutine acecel(noma, nomo, nbocc, ele_sup_num, ele_sup_typ, nb_ty_el, zjdlm, ier)
!
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!     COMPTEUR D'ELEMENTS
!
!       IN
!           ele_sup_num : numéro des éléments dans les catalogues &CATA.TE.NOMTE
!           ele_sup_typ : numéro interne à ACE du type de l'élément
!                           cf cara_elem_parameter_module
!       OUT
!           nb_ty_el    : nombre du type d'élément  nb_ty_el ( [ele_sup_typ] )
!           zjdlm       : numéro de l'élément porté par la maille
!
! --------------------------------------------------------------------------------------------------
!
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : NOMO   : NOM DU MODELE
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use cara_elem_parameter_module
    implicit none
    character(len=8) :: noma, nomo
    integer(kind=8) :: nbocc(*), ele_sup_num(*), ele_sup_typ(*), nb_ty_el(*), zjdlm(*), ier
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ii, ifm, ixma, tt, jdme, nbmail, nummai, nutyel, nb_elem_p
!
    character(len=24) :: mlgnma, modmai
    aster_logical :: l_pmesh
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    modmai = nomo//'.MAILLE'
    mlgnma = noma//'.TYPMAIL'
    call jeexin(modmai, ixma)
    call jelira(mlgnma, 'LONMAX', nbmail)
    ASSERT(ixma .ne. 0)
    call jeveuo(modmai, 'L', jdme)
!
    l_pmesh = isParallelMesh(noma)
!
    ifm = iunifi('MESSAGE')
!
    nb_ty_el(1:ACE_NB_ELEMENT) = 0
!
    do nummai = 1, nbmail
        nutyel = zi(jdme+nummai-1)
        zjdlm(nummai) = nutyel
!
        do ii = 1, ACE_NB_TYPE_ELEM
            if (nutyel .eq. ele_sup_num(ii)) then
                tt = ele_sup_typ(ii)
                nb_ty_el(tt) = nb_ty_el(tt)+1
                exit
            end if
        end do
    end do
!
    write (ifm, 100) nomo
    do ii = 1, ACE_NB_ELEMENT
        if (nb_ty_el(ii) .gt. 0) then
            write (ifm, 110) nb_ty_el(ii), ACE_NM_ELEMENT(ii)
        end if
    end do
100 format(/, 5x, 'LE MODELE ', a8, ' CONTIENT : ')
110 format(35x, i6, ' ELEMENT(S) ', A16)
!
!   Vérification de la cohérence des affectations
    if (nbocc(ACE_POUTRE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_POUTRE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_29', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_COQUE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_COQUE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_30', sk=nomo)
            ier = ier+1
        end if
    end if
    if ((nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_DISCRET)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_31', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_ORIENTATION) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_POUTRE)+nb_ty_el(ACE_NU_DISCRET)+nb_ty_el(ACE_NU_BARRE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_32', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_CABLE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_CABLE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_33', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_BARRE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_BARRE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_34', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_MASSIF) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_MASSIF)+nb_ty_el(ACE_NU_THHMM)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_35', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_GRILLE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_GRILLE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_36', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_MEMBRANE) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_MEMBRANE)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_55', sk=nomo)
            ier = ier+1
        end if
    end if
    if (nbocc(ACE_MASS_REP) .ne. 0) then
        nb_elem_p = nb_ty_el(ACE_NU_DISCRET)
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nb_elem_p)
        end if
        if (nb_elem_p .eq. 0) then
            call utmess('E', 'MODELISA_31', sk=nomo)
            ier = ier+1
        end if
    end if
!
    call jedema()
end subroutine
