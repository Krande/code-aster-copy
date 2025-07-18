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

subroutine vpstur(lmatk, valshi, lmatm, lmatsh, mantis, &
                  expo, pivot, ier, solveu, caldet, &
                  calfac)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/freqom.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdete.h"
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
#include "asterfort/vpshif.h"
!
    real(kind=8) :: valshi, mantis
    integer(kind=8) :: lmatk, lmatm, lmatsh, expo, pivot, ier
    character(len=19) :: solveu
!     EFFECTUE L'OPERATION DE STURM
!        1) COMBINAISON LINEAIRE MSH =  K - W * M    (W ETANT LE SHIFT)
!        2) DECOMPOSITION DE  MSH
!        3) NOMBRE DE PIVOT NEGATIF
!        4) CALCUL DU DETERMINANT
!     ------------------------------------------------------------------
! IN  VALSHI : R8 : VALEUR DU DECALAGE
! IN  LMATK  : IS : ADRESSE ATTRIBUT MATRICE K
! IN  LMATM  : IS : ADRESSE ATTRIBUT MATRICE M
! IN  LMATSH : IS : ADRESSE ATTRIBUT MATRICE SHIFTEE
! OUT MANTIS : R8 : MANTISSE DU DETERMINANT
! OUT EXPO   : IS : EXPOSANT DU DETERMINANT
! OUT PIVOT  : IS : NOMBRE DE TERMES DIAGONAUX NEGATIFS
! OUT IER    : IS : CODE RETOUR  /= 0 ==> LE SHIFT EST UNE VALEUR PROPRE
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
! IN  CALDET : LOG : SI TRUE ON CALCULE LE DETERMINANT, SI FALSE ON NE
!                     LE CALCULE PAS (GAIN TEMPS)
! IN  CALFAC : LOG : SI MUMPS: SI TRUE, ON GARDE LES TERMES DE LA FACTO
!                    RISEE, SI FALSE, ON NE LES GARDE PAS (GAIN ESPACE)
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: iret, npvneg, iold, iold2
    real(kind=8) :: valr
    complex(kind=8) :: cbid
    character(len=19) :: matpre, matass
    character(len=24) :: metres
    aster_logical :: caldet, calfac
    integer(kind=8), pointer :: slvi(:) => null()
    character(len=24), pointer :: slvk(:) => null()
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call jemarq()
!
!     --- INITIALISATION ---
    call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
    metres = slvk(1)
    call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
!
!
! ---- OPTIMISATION MUMPS VIA CALFAC/CALDET: PART 1/2
! ---- SI CALFAC=.FALSE., ON NE STOCKE PAS LES FACTEURS, SEUL LE
! ---- CARACTERE SINGULIER ET LE NBRE DE TERMES <0 DE LA DIAGONALE NOUS
! ---- INTERESSENT
! ---- SI CALDET=.TRUE.: A PARTIR DE MUMPS.4.10.0 ON CALCULE LE DET
    iold = -9999
    iold2 = -9999
    if (metres(1:5) .eq. 'MUMPS') then
        if (.not. calfac) then
            iold = slvi(4)
            slvi(4) = 1
        end if
        if (caldet) then
            iold2 = slvi(5)
            slvi(5) = 1
        end if
    end if
!
!     --- DECALAGE SPECTRAL  K - W * M    (W ETANT LE SHIFT) ---
    call vpshif(lmatk, valshi, lmatm, lmatsh)
!
!     --- FACTORISATION LDLT DE LA MATRICE SHIFTEE---
    ier = 0
!
    matpre = ' '
    matass = zk24(zi(lmatsh+1))
    call preres(solveu, 'V', iret, matpre, matass, &
                npvneg, 2)
!
    if (iret .ge. 1) ier = 1
    if (iret .gt. 1) then
        valr = freqom(valshi)
        call utmess('A', 'ALGELINE5_27', sr=valr)
    end if
    pivot = -npvneg
!
! ---  CALCUL OPTIONNEL DU DETERMINANT
! ---- OPTIMISATION VIA CALDET (SI MF OU LDLT OU MUMPS V4.10.0 ET PLUS)
    if (caldet) then
        if ((metres(1:10) .ne. 'MULT_FRONT') .and. (metres(1:4) .ne. 'LDLT') .and. &
            (metres(1:5) .ne. 'MUMPS')) then
            call utmess('F', 'ALGELINE5_73')
        else
            call mtdete(1, metres, lmatsh, mantis, expo, &
                        cbid)
        end if
    end if
!
! ---- OPTIMISATION MUMPSVIA CALFAC/CALDET: PART 2/2
! ---- ON REMET DANS LA SD_SOLVEUR.SLVI(4) L'ANCIENNE VALEUR
! ---- (NORMALEMENT LA VALEUR INITIALISEE PAR DEFAUT -9999).
    if (metres(1:5) .eq. 'MUMPS') then
        if (.not. calfac) slvi(4) = iold
        if (caldet) slvi(5) = iold2
    end if
!
    call jedema()
end subroutine
