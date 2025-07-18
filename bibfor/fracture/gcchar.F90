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

subroutine gcchar(ichar, iprec, time, carteo, lfchar, &
                  lpchar, lformu, lfmult, lccomb, cartei, &
                  nomfct, newfct, oldfon)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcharf.h"
#include "asterfort/gcharm.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
!
!
    aster_logical :: lfchar
    aster_logical :: lpchar
    aster_logical :: lformu
    aster_logical :: lfmult
    aster_logical :: lccomb
    character(len=24) :: oldfon
    character(len=8) :: nomfct
    character(len=8) :: newfct
    integer(kind=8) :: ichar
    integer(kind=8) :: iprec
    real(kind=8) :: time
    character(len=19) :: cartei
    character(len=19) :: carteo
!
! ----------------------------------------------------------------------
!
! ROUTINE CALC_G
!
! CONSTRUIT LA CARTE A PARTIR DU CHARGEMENT
!
! ----------------------------------------------------------------------
!
!
! IN  LFCHAR : .TRUE.  SI LE CHARGEMENT EST 'FONCTION'
! IN  LFORMU : .TRUE.  SI LE CHARGEMENT 'FONCTION' UTILISE UNE FORMULE
! IN  LFMULT : .TRUE.  S'IL Y A UNE FONCTION MULTIPLICATRICE
! IN  LCCOMB : .TRUE.  SI LE CHARGEMENT EST COMBINABLE
! IN  LPCHAR : .TRUE.  SI C'EST LA PREMIERE FOIS QU'ON A UNE CHARGE DU STYLE COURANT
! IN  ICHAR  : INDICE DU CHARGEMENT
! I/O OLDFON : LISTE DES TYPES DE CHARGEMENTS
! IN  NOMFCT : NOM DE LA FONCTION MULTIPLICATRICE
! IN  TIME   : INSTANT
! I/O IPREC  : INDICE DU CHARGEMENT PRECEDENT DU MEME TYPE
! I/O NEWFCT : FONCTION MULTIPLICATRICE MODIFIEE DANS LA CARTE DE SORTIE
!              PRODUIT DE LA FONC_MULT ET DE LA DEPENDANCE EVENTUELLE
!              VENUE D'AFFE_CHAR_MECA_F
! IN  CARTEI : CARTE DU CHARGEMENT AVANT LA PRISE EN COMPTE
!              DE LA FONCTION MULTIPLICATRICE
! OUT CARTEO : CARTE DU CHARGEMENT APRES LA PRISE EN COMPTE
!              DE LA FONCTION MULTIPLICATRICE
!
! ----------------------------------------------------------------------
!
    character(len=19) :: chtmp1, chtmp2
    aster_logical :: fonc1, fonc2
    integer(kind=8) :: jfonci
!
! ----------------------------------------------------------------------
!
    chtmp1 = '&&GCCHAR_INTERM1'
    chtmp2 = '&&GCCHAR_INTERM2'
    call jeveuo(oldfon, 'L', jfonci)
!
    if (lpchar) then
        call copisd('CHAMP_GD', 'V', cartei, carteo)
        if (lfmult .and. (.not. lformu)) then
            call gcharm(lfchar, cartei, nomfct, newfct, time, &
                        carteo)
        end if
    else
        if (.not. lccomb) then
            call utmess('F', 'RUPTURE2_3')
        end if
        if (lformu) then
            call utmess('F', 'RUPTURE2_2')
        end if
        call copisd('CHAMP_GD', 'V', carteo, chtmp1)
        call copisd('CHAMP_GD', 'V', cartei, chtmp2)
        call detrsd('CHAMP_GD', carteo)
        if (lfmult) then
            call gcharm(lfchar, cartei, nomfct, newfct, time, &
                        chtmp2)
        end if
        fonc1 = zl(jfonci+iprec-1)
        fonc2 = zl(jfonci+ichar-1)
        call tecart(chtmp1)
        call tecart(chtmp2)
!
! ----- EFFECTUE LA FUSION DE 2 CHARGES DE MEME TYPE
!
        call gcharf(ichar, fonc1, chtmp1, fonc2, chtmp2, &
                    carteo, oldfon)
    end if
    iprec = ichar
    call detrsd('CHAMP_GD', chtmp1)
    call detrsd('CHAMP_GD', chtmp2)
!
end subroutine
