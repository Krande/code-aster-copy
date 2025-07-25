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

subroutine uttcpl(dim, nbmesu, nomc, noml, prpal)
    implicit none
#include "asterfort/assert.h"
    integer(kind=8) :: dim, nbmesu
    character(len=1) :: prpal(dim)
    character(len=24) :: nomc(dim)
    character(len=80) :: noml(dim)
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!  RETOURNE LA LISTE DES DIFFERENTES MESURES DE TEMPS GENERALES
!
!  IN DIM : DIMENSION DES TABLEAUX NOMC, NOML ET PRPAL
!  OUT NBMESU : NOMBRE DE MESURES
!  OUT NOMC(IMES) : NOM COURT DE LA MESURE
!  OUT NOML(IMES) : NOM LONG DE LA MESURE
!  OUT PRPAL(IMES) :  'P' : MESURE PRICIPALE
!                     'S' : MESURE SECONDAIRE
! ----------------------------------------------------------------------
    integer(kind=8) :: nbmax, k, i1, i2
    parameter(nbmax=30)
    character(len=80) :: d1(nbmax)
!
!     -- COMMONS POUR MESURE DE TEMPS :
    integer(kind=8) :: mtpniv, mtpsta, indmax
    parameter(indmax=5)
    character(len=80) :: snolon(indmax)
    real(kind=8) :: valmes(indmax*7), valmei(indmax*7)
    common/mestp1/mtpniv, mtpsta
    common/mestp2/snolon
    common/mestp3/valmes, valmei
!     ------------------------------------------------------------------
!
!     -- SI L'UTILISATEUR NE VEUT PAS DE MESURE, NBMESU=0
    if (mtpniv .eq. 0) then
        nbmesu = 0
        goto 999
    end if
!
!     ON ECRIT LES DONNEES DANS LE TABLEAU D1 :
!     -----------------------------------------------------
    d1(1) = 'CPU.RESO.1|P|1#Resolution des systemes lineaires'
    d1(2) = 'CPU.RESO.2|S|1.1#Numerotation, connectivite de la matrice'
    d1(3) = 'CPU.RESO.3|S|1.2#Factorisation symbolique'
    d1(4) = 'CPU.RESO.4|S|1.3#Factorisation numerique (ou precond.)'
    d1(5) = 'CPU.RESO.5|S|1.4#Resolution'
!
    d1(6) = 'CPU.CALC.1|P|2#Calculs elementaires et assemblages'
    d1(7) = 'CPU.CALC.2|S|2.1#Routine calcul'
    d1(8) = 'CPU.CALC.3|S|2.1.1#Routines te00ij'
    d1(9) = 'CPU.ASSE.1|S|2.2#Assemblages'
    d1(10) = 'CPU.ASSE.2|S|2.2.1#Assemblage matrices'
    d1(11) = 'CPU.ASSE.3|S|2.2.2#Assemblage seconds membres'
!
!     -- LES 2 MESURES SUIVANTES SONT SPECIALES : ON S'INTERDIT DE
!        FAIRE APPEL A JEVEUX. VOIR UTTCPU.F, UTTCPR.F
!        C'EST LEUR NOM QUI EST CONNU DE TOUS : 'CPU.MEMD.1/2'
    d1(12) = 'CPU.MEMD.1|P|3#Dechargement de la memoire sur disque'
    d1(13) = 'CPU.MEMD.2|P|3#??? libre pour mesure interne jeveux'
!
    d1(14) = 'CPU.CMPI.1|P|4#Communications MPI'
!
    nbmesu = 14
    ASSERT(nbmax .ge. nbmesu)
    ASSERT(dim .ge. nbmesu)
!
!
!
!     ON "SPLITE" D1 DANS NOMC, NOML ET PRPAL :
!     -----------------------------------------------------
    do k = 1, nbmesu
        i1 = index(d1(k), '|')
        nomc(k) = d1(k) (1:i1-1)
        i2 = index(d1(k) (i1+1:), '|')
        ASSERT(i2 .eq. 2)
        prpal(k) = d1(k) (i1+1:i1+2)
        ASSERT(prpal(k) .eq. 'P' .or. prpal(k) .eq. 'S')
        noml(k) = d1(k) (i1+i2+1:)
!
    end do
!
999 continue
end subroutine
