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

subroutine asret2(lmasym, jtmp2, lgtmp2, nbterm, jsmhc, &
                  jsmdi, nbi1, ti1, ti2)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jeveut.h"
#include "asterfort/juveca.h"
    aster_logical :: lmasym
    integer(kind=8) :: jtmp2, lgtmp2, nbterm, jsmhc, jsmdi
    integer(kind=8) :: ideb, ifin, imil
    integer(kind=8) :: nbi1, ti1(nbi1), ti2(nbi1)
!     ROUTINE SERVANT A RETENIR OU S'ACCUMULENT LES TERMES ELEMENTAIRES:
!     DANS LE CAS D'UN STOCKAGE MORSE SYMETRIQUE
! -----------------------------------------------------------------
! IN/OUT I JTMP2   : ADRESSE JEVEUX DE L'OBJET ".TMP2"
! IN    I4 JSMHC   : ADRESSE DE ".SMHC".
! IN     I JSMDI   : ADRESSE DE ".SMDI".
! IN     I NBI1    : NOMBRES DE TERMES DANS TI1 (ET TI2)
! IN     I TI1(*),TI2(*)   : NUMEROS GLOBAUX (LIGNE ET COLONNE)
! IN/OUT I NBTERM   : NOMBRE DE TERMES (R/C) A RECOPIER
!                     (ISSU DE LA MATRICE ELEMENTAIRE)
! -----------------------------------------------------------------
    integer(kind=8) :: ili, jco, icoefc, icoefl, i, ncoefc, nubloc, k, i1, i2
! -----------------------------------------------------------------
    do k = 1, nbi1
        i1 = ti1(k)
        i2 = ti2(k)
!
        if (i1 .le. i2) then
            ili = i1
            jco = i2
            nubloc = 1
        else
            ili = i2
            jco = i1
            nubloc = 2
        end if
        if (lmasym) nubloc = 1
!
        if (jco .eq. 1) then
            icoefc = 0
        else
            icoefc = zi(jsmdi+jco-2)
        end if
        ncoefc = zi(jsmdi+jco-1)-icoefc
!
!
!     -- CALCUL DE ICOEFL :
!     ------------------------------------------
        icoefl = 0
        if (.false._1) then
!         -- RECHERCHE BESTIALE :
            do i = 1, ncoefc
                if (zi4(jsmhc-1+icoefc+i) .eq. ili) then
                    icoefl = i
                    goto 40
!
                end if
            end do
!
        else
!          -- RECHERCHE PAR DICHOTOMIE :
            ideb = 1
            ifin = ncoefc
20          continue
            if (ifin-ideb .lt. 5) then
                do i = ideb, ifin
                    if (zi4(jsmhc-1+icoefc+i) .eq. ili) then
                        icoefl = i
                        goto 40
!
                    end if
                end do
            end if
            imil = (ideb+ifin)/2
            if (zi4(jsmhc-1+icoefc+imil) .gt. ili) then
                ifin = imil
            else
                ideb = imil
            end if
            goto 20
!
        end if
!       IF (ICOEFL.EQ.0) CALL UTMESS('F','MODELISA_67')
!
!
40      continue
!
!     -- NBTERM COMPTE LES REELS TRAITES:
        nbterm = nbterm+1
        if (2*nbterm .gt. lgtmp2) then
            lgtmp2 = 2*lgtmp2
            call juveca('&&ASSMAM.TMP2', lgtmp2)
!         -- IL NE FAUT PAS QUE .TMP2 SOIT LIBERE :
            call jeveut('&&ASSMAM.TMP2', 'E', jtmp2)
        end if
        zi(jtmp2-1+(nbterm-1)*2+1) = nubloc
        zi(jtmp2-1+(nbterm-1)*2+2) = icoefc+icoefl
    end do
end subroutine
