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

subroutine ace_verif_affe(infoconcept, nbocc, nb_ty_el, zjdlm)
!
    use cara_elem_parameter_module
    use cara_elem_info_type
!
    implicit none
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
!
    type(cara_elem_info), intent(in) :: infoconcept
    integer(kind=8), intent(in)      :: nbocc(*)
    integer(kind=8), intent(out)     :: nb_ty_el(*), zjdlm(*)
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!     COMPTEUR D'ELEMENTS
!
! --------------------------------------------------------------------------------------------------
!

    integer(kind=8) :: ii, ifm, tt, jdme, nbmail, nummai, nutyel, nbelem, ier
!
    logical :: TestOcc
!
    character(len=8)    :: noma, nomo
    character(len=32)   :: vmessk(3)
!
    aster_logical :: l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    nomo = infoconcept%modele
    noma = infoconcept%maillage
    nbmail = infoconcept%nbmail
    jdme = infoconcept%jmodmail
!
    l_pmesh = infoconcept%IsParaMesh
    ifm = iunifi('MESSAGE')
!
    nb_ty_el(1:ACE_NB_ELEMENT) = 0
!
    do nummai = 1, nbmail
!       Le type de l'élément venant de affe_modele
        nutyel = zi(jdme+nummai-1)
        zjdlm(nummai) = nutyel
!       Comptage des éléments en fonction de leur classification dans cara_elem_parameter_module
        doii: do ii = 1, ACE_NB_TYPE_ELEM
            if (nutyel .eq. elem_supp%catanum(ii)) then
                tt = elem_supp%acenum(ii)
                nb_ty_el(tt) = nb_ty_el(tt)+1
                exit doii
            end if
        end do doii
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
!   Si pas de vérification
    if (.not. infoconcept%VerifMaille) then
        goto 999
    end if
!
!   Vérification de la cohérence des affectations
    ier = 0
    vmessk(1) = nomo
!
!   Les éléments : POUTRE COQUE CABLE BARRE GRILLE MEMBRANE
!       ne peuvent affectés que par leur mot clef facteur
!       Si MOT_F si pas d'élément ==> E
!       Sinon    si des éléments  ==> A
    nbelem = nb_ty_el(ACE_NU_POUTRE)
    if (nbocc(ACE_POUTRE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'POUTRE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='POUTRE')
        ier = ier+1
    end if
!
    nbelem = nb_ty_el(ACE_NU_COQUE)
    if (nbocc(ACE_COQUE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'COQUE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='COQUE')
        ier = ier+1
    end if
!
    nbelem = nb_ty_el(ACE_NU_CABLE)
    if (nbocc(ACE_CABLE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'CABLE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='CABLE')
        ier = ier+1
    end if
!
    nbelem = nb_ty_el(ACE_NU_BARRE)
    if (nbocc(ACE_BARRE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'BARRE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='BARRE')
        ier = ier+1
    end if
!
    nbelem = nb_ty_el(ACE_NU_GRILLE)
    if (nbocc(ACE_GRILLE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'GRILLE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='GRILLE')
        ier = ier+1
    end if
!
    nbelem = nb_ty_el(ACE_NU_MEMBRANE)
    if (nbocc(ACE_MEMBRANE) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'MEMBRANE'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    else if (nbelem .ne. 0) then
        call utmess('A', 'AFFECARAELEM_32', sk='MEMBRANE')
        ier = ier+1
    end if
!
!   Les éléments concernés par ORIENTATION
!       L'orientation est facultative, on ne teste donc pas (nbelem .ne. 0)
    nbelem = nb_ty_el(ACE_NU_POUTRE)+nb_ty_el(ACE_NU_DISCRET)+nb_ty_el(ACE_NU_BARRE)
    if (nbocc(ACE_ORIENTATION) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'ORIENTATION'; vmessk(3) = 'POUTRE | DISCRET | BARRE'
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    end if
!
!   Les éléments DISCRET peuvent être affectés de plusieurs façons
!       DISCRET DISCRET_2D MASS_REP RIGI_PARASOL RIGI_MISS_3D
    TestOcc = .false.
    nbelem = nb_ty_el(ACE_NU_DISCRET)
    if ((nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'DISCRET'; vmessk(3) = vmessk(2)
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
        TestOcc = TestOcc .or. .true.
    end if
    if (nbocc(ACE_MASS_REP) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'MASS_REP'; vmessk(3) = 'DISCRET'
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
        TestOcc = TestOcc .or. .true.
    end if
    if (nbocc(ACE_RIGI_PARASOL) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'RIGI_PARASOL'; vmessk(3) = 'DISCRET'
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
        TestOcc = TestOcc .or. .true.
    end if
    if (nbocc(ACE_RIGI_MISS_3D) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'RIGI_MISS_3D'; vmessk(3) = 'DISCRET'
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
        TestOcc = TestOcc .or. .true.
    end if
!   Si on à des DISCRET (nbelem <>0) et TestOcc = False
    if ((nbelem .ne. 0) .and. (.not. TestOcc)) then
        call utmess('A', 'AFFECARAELEM_32', sk='DISCRET')
        ier = ier+1
    end if
!
!   L'affectation des éléments MASSIF est facultative, on ne teste donc pas (nbelem .ne. 0)
    nbelem = nb_ty_el(ACE_NU_MASSIF)+nb_ty_el(ACE_NU_THHMM)+nb_ty_el(ACE_NB_HHO)
    if (nbocc(ACE_MASSIF) .ne. 0) then
        if (l_pmesh) then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=nbelem)
        end if
        if (nbelem .eq. 0) then
            vmessk(2) = 'MASSIF'; vmessk(3) = 'THERMIQUE | MECANIQUE | HHO'
            call utmess('E', 'AFFECARAELEM_31', nk=3, valk=vmessk)
            ier = ier+1
        end if
    end if
!
    if (ier .ne. 0) call utmess('F', 'AFFECARAELEM_33')
999 continue
    call jedema()
end subroutine
