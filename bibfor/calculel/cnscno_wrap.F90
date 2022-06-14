! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine cnscno_wrap(cnsz, prchnz, prol0, basez, cnoz, kstop, iret)
! A_UTIL
  implicit none
#include "asterf_types.h"
#include "asterfort/cnscno.h"
!
character(len=*) :: cnsz, cnoz, basez, prchnz, prol0
character(len=1) :: kstop
integer :: iret
!
! ------------------------------------------------------------------
! BUT : TRANSFORMER UN CHAM_NO_S (CNSZ) EN CHAM_NO (CNOZ)
! ------------------------------------------------------------------
!     ARGUMENTS:
! CNSZ    IN/JXIN  K19 : SD CHAM_NO_S A TRANSFORMER
! PRCHNZ  IN/JXVAR K19 : SD PROF_CHNO  (OU ' ')
!          SI PRCHNZ EXISTE ON CREE CNOZ CONFORMEMENT A PRCHNZ :
!             => SI CNSZ CONTIENT DES VALEURS QUE L'ON NE SAIT PAS
!                STOCKER DANS PRCHNZ, ON LES "OUBLIE"
!             => SI PRCHNZ EXIGE DES VALEURS QUE L'ON NE TROUVE PAS
!                DANS CNSZ :
!                  - SI PROL0='OUI' : ON PRENDS LA VALEUR "ZERO"
!                  - SI PROL0='NON' : ERREUR <F>
!
!          SI PRCHNZ N'EXISTE PAS ON CREE CNOZ EN FONCTION
!             DU CONTENU DE CNSZ
!             SI PRCHNZ  = ' ' ON CREE UN PROF_CHNO "SOUS-TERRAIN"
!             SI PRCHNZ /= ' ' ON CREE UN PROF_CHNO DE NOM PRCHNZ
! PROL0   IN   K3  :  POUR PROLONGER (OU NON) LE CHAMP PAR "ZERO"
!        /OUI /NON  ( CET ARGUMENT N'EST UTILISE QUE SI PRCHNZ /= ' ')
!        "ZERO" : / 0       POUR LES CHAMPS NUMERIQUES (R/C/I)
!                 / ' '     POUR LES CHAMPS "KN"
!                 / .FALSE. POUR LES CHAMPS DE "L"
!
! BASEZ   IN       K1  : BASE DE CREATION POUR CNOZ : G/V/L
! CNOZ    IN/JXOUT K19 : SD CHAM_NO A CREER
! KSTOP   IN       K1  : COMPORTEMENT EN CAS DE PROBLEME :
!              / 'A' : ON EMET UNE ALARME ET ON REND IRET > 0
!              / 'F' : ON EMET UNE ERREUR FATALE
!              / ' ' : ON N'EMET PAS DE MESSAGE
! IRET    OUT       I  : CODE DE RETOUR :
!              / 0 : OK
!              / 1 : LE CHAM_NO N'A PAS PU ETRE CREE
!-----------------------------------------------------------------------
!
    call cnscno(cnsz, prchnz, prol0, basez, cnoz, kstop, iret)
end subroutine
