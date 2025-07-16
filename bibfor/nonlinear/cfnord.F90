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
subroutine cfnord(noma, typent, nument, itype, vector, &
                  tau1, tau2, lnfixe)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
#include "blas/dcopy.h"
!
    character(len=8) :: noma
    character(len=4) :: typent
    integer(kind=8) :: nument
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: vector(3)
    integer(kind=8) :: itype
    aster_logical :: lnfixe
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - APPARIEMENT)
!
! MODIFIE LES VECTEURS TANGENTS LOCAUX QUAND NORMALE DONNEE PAR
! UTILISATEUR
!
! ----------------------------------------------------------------------
!
!  NB: LE REPERE EST ORTHORNORME ET TEL QUE LA NORMALE POINTE VERS
!  L'EXTERIEUR DE LA MAILLE
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  TYPENT : TYPE DE L'ENTITE
!               'MAIL' UNE MAILLE
!               'NOEU' UN NOEUD
! IN  NUMENT : NUMERO ABSOLU DE L'ENTITE DANS LE MAILLAGE
! IN  ITYPE  : TYPE DE NORMALE
!                0 AUTO
!                1 FIXE   (DONNE PAR VECTOR)
!                2 VECT_Y (DONNE PAR VECTOR)
! IN  VECTOR : VALEUR DE LA NORMALE FIXE OU VECT_Y
! I/O TAU1   : PREMIER VECTEUR TANGENT LOCAL
! I/O TAU2   : SECOND VECTEUR TANGENT LOCAL
! OUT LNFIXE : VAUT .TRUE. SI NORMALE='FIXE' OU 'VECT_Y'
!                   .FALSE. SI NORMALE='AUTO'
!
!
!
!
    character(len=8) :: noment
    real(kind=8) :: norm(3), noor, noor2
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    lnfixe = .false.
!
! --- NOM DE L'ENTITE (NOEUD OU MAILLE)
!
    if (typent .eq. 'MAIL') then
        noment = int_to_char8(nument)
    else if (typent .eq. 'NOEU') then
        noment = int_to_char8(nument)
    else
        ASSERT(.false.)
    end if
!
! --- NORMALE AUTOMATIQUE: ON SORT
!
    if (itype .eq. 0) then
        lnfixe = .false.
        goto 999
    else
        call normev(vector, noor)
        if (noor .le. r8prem()) then
            ASSERT(.false.)
        end if
        lnfixe = .true.
    end if
!
! --- REDEFINITION SI VECT_ == 'FIXE' (ON GARDE T1 COMME REFERENCE)
!
    if (itype .eq. 1) then
        call provec(vector, tau1, tau2)
        call normev(tau2, noor)
        if (noor .le. r8prem()) then
            if (typent .eq. 'MAIL') then
                call utmess('F', 'CONTACT_14', sk=noment)
            else if (typent .eq. 'NOEU') then
                call utmess('F', 'CONTACT_13', sk=noment)
            else
                ASSERT(.false.)
            end if
        end if
    end if
!
! --- REDEFINITION SI VECT_ == 'VECT_Y'
!
    if (itype .eq. 2) then
!
! --- VECTEUR TAU2 NUL OU POUTRE !
!
        call normev(tau2, noor)
        if (noor .le. r8prem()) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vector, b_incx, tau2, b_incy)
            call provec(tau1, tau2, norm)
            call normev(norm, noor2)
            if (noor2 .le. r8prem()) then
                if (typent .eq. 'MAIL') then
                    call utmess('F', 'CONTACT3_27', sk=noment)
                else if (typent .eq. 'NOEU') then
                    call utmess('F', 'CONTACT3_26', sk=noment)
                else
                    ASSERT(.false.)
                end if
            end if
        else
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vector, b_incx, tau2, b_incy)
            call provec(tau1, tau2, norm)
            call normev(norm, noor2)
            if (noor2 .le. r8prem()) then
                if (typent .eq. 'MAIL') then
                    call utmess('F', 'CONTACT3_27', sk=noment)
                else if (typent .eq. 'NOEU') then
                    call utmess('F', 'CONTACT3_26', sk=noment)
                else
                    ASSERT(.false.)
                end if
            end if
            call provec(tau2, norm, tau1)
        end if
    end if
!
    if (itype .ge. 3) then
        ASSERT(.false.)
    end if
!
999 continue
!
    call jedema()
!
end subroutine
