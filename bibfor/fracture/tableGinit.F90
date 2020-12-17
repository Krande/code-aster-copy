! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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

subroutine tableGinit(table, option, ndim, nxpara,&
                  lmoda, nbpara, linopa, litypa)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/cgajpa.h"
#include "asterfort/assert.h"
!
    integer :: nbpara, nxpara, ndim
    aster_logical :: lmoda
    character(len=*) :: litypa(nxpara), linopa(nxpara)
    character(len=8) :: table, option
!
! person_in_charge: samuel.geniaut at edf.fr
!
!     SOUS-ROUTINE DE L'OPERATEUR CALC_G
!
!     BUT : CREATION DE LA TABLE ISSUE DE CALC_G
!           ET AFFECTATION DES PARAMETRES
!
! ----------------------------------------------
!  IN :
!     TABLE : NOM DE LA TABLE
!     OPTION : OPTION DE CALCUL
!     NDIM    : DIMENSION DU CALCUL : 2 OU 3
!     TYPFIS : TYPE D'OBJET POUR DECRIRE LE FOND DE FISSURE
!              'FONDFISS' OU 'FISSURE' OU 'THETA'
!     NXPARA : NOMBRE MAXI DE PARAMETRES DE LA TABLE
!     LMODA  : .TRUE.  SI TYPE SD RESULTAT = MODE_MECA
!              .FALSE. SINON
!
!  OUT :
!     NBPARA : NOMBRE DE PARAMETRES
!     linopa : NOMS DES PARAMETRES
!     litypa : TYPES DES PARAMETRES
! ----------------------------------------------
!
    integer :: i
    
    aster_logical :: debug
    nbpara = 0
    debug = .false.
!
!   --------------------    
!   1.2 TEMPOREL/CHARGEMENT
!   ---------------------           
    if (lmoda) then
        call cgajpa('NUME_MODE', 'I', nbpara, linopa, litypa, nxpara)   
    else
        call cgajpa('NUME_ORDRE', 'I', nbpara, linopa, litypa, nxpara)
        call cgajpa('INST', 'R', nbpara, linopa, litypa, nxpara)    
    endif
!
    call cgajpa('TEMP', 'R', nbpara, linopa, litypa, nxpara)
    call cgajpa('COMPORTEMENT', 'K8', nbpara, linopa, litypa, nxpara)
!
!   --------------------    
!   1.3 POINT DU FOND DE FISSURE
!   ---------------------    

    call cgajpa('COOR_X', 'R', nbpara, linopa, litypa, nxpara)
    call cgajpa('COOR_Y', 'R', nbpara, linopa, litypa, nxpara)  
    if (ndim.eq.3) then
        call cgajpa('COOR_Z', 'R', nbpara, linopa, litypa, nxpara)       
        call cgajpa('ABSC_CURV_NORM', 'R', nbpara, linopa, litypa, nxpara)
    endif
!
!   --------------------
!   2. OPTIONS DE CALCUL
!   --------------------
!   --------------------
!   2.1 G COMMUN A TOUTES LES OPTIONS
!   ---------------------
    call cgajpa('G', 'R', nbpara, linopa, litypa, nxpara)
!   --------------------
!   2.2 CALC_K_G
!   ---------------------
    if (option.eq.'K') then
        call cgajpa('K1', 'R', nbpara, linopa, litypa, nxpara)
        call cgajpa('K2', 'R', nbpara, linopa, litypa, nxpara)
!            call cgajpa('G_IRWIN', 'R', nbpara, linopa, litypa, nxpara)
        if (ndim.eq.3) then
            call cgajpa('K3', 'R', nbpara, linopa, litypa, nxpara)
        endif
    endif

!   --------------------
!   3. CREATION DE LA TABLE
!   --------------------
!    call tbcrsd(table, 'G')
    call tbajpa(table, nbpara, linopa, litypa)

!   --------
!   4. DEBUG
!   --------
    if (debug) then
        write(6,*)'OPTION = ', option
        write(6,*)'NOMBRE DE PARAMETRES DE LA TABLE = ', nbpara
        write(6,*)'NO_PARA, NOM_PARA, TYP_PARA'
        do i = 1, nbpara
            write(6,*)i, linopa(i), litypa(i)
        enddo
    endif

end subroutine
