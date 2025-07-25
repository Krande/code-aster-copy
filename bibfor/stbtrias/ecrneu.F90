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
subroutine ecrneu(imod, nbnode, ama, bma, cma, &
                  ami, bmi, cmi, min, man, &
                  ites)
    implicit none
!     ==============================================================
!A PRESUPER
!
!   ================================================================
!   !                                                              !
!   !  FONCTION: ECRITURE DES COORDONNEES DES NOEUDS SUR LE FICHIER!
!   !            MODELE A PARTIR DU FICHIER BUFFER (NOEUDS-BUFFER) !
!   !                                                              !
!   ================================================================
!   !                                                              !
!   !  ROUTINES APPELES: CODENT                                    !
!   !                         : IUNIFI (FONCTION)                  !
!   !                         : JJMMAA                             !
!   !                                                              !
!   !                                                              !
!   !  ROUTINE APPELANTE : PRESUP                                  !
!   !                                                              !
!   ================================================================
!   !                                                              !
!   !                 ***************                              !
!   !                 *  ARGUMENTS  *                              !
!   !                 ***************                              !
!   !                                                              !
!   !  **********************************************************  !
!   !  *   NOM    * TYPE  *  MODE  *ALTERE *       ROLE         *  !
!   !  **********************************************************  !
!   !  *          *       *        *       *                    *  !
!   !  * NBNODE   *INTEGER*ENTREE  * NON   *NBRE DE NOEUDS TOTAL*  !
!   !  *          *       *        *       *                    *  !
!   !  * AMA      *D.PRECI*ENTREE  * NON   * X(MAXIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * BMA      *D.PRECI*ENTREE  * NON   * Y(MAXIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * CMA      *D.PRECI*ENTREE  * NON   * Z(MAXIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * AMI      *D.PRECI*ENTREE  * NON   * X(MINIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * BMI      *D.PRECI*ENTREE  * NON   * Y(MINIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * CMI      *D.PRECI*ENTREE  * NON   * Z(MINIMUM)         *  !
!   !  *          *       *        *       *                    *  !
!   !  * MAN      *INTEGER*ENTREE  * NON   * N DE NOEUD MAXI    *  !
!   !  *          *       *        *       *                    *  !
!   !  * MIN      *INTEGER*ENTREE  * NON   * N DE NOEUD MINI    *  !
!   !  *          *       *        *       *                    *  !
!   !  * ITES     *INTEGER*ENTREE  * NON   * INDIQUE L'EXISTENCE*  !
!   !  *          *       *        *       * D'AU MOINS DEUX    *  !
!   !  *          *       *        *       * SYSTEMES DE COORD. *  !
!   !  *          *       *        *       *                    *  !
!   !  **********************************************************  !
!   !                                                              !
!   ================================================================
!
!
!
!
! ---> DECLARATION DES VARIABLES POUR LE TYPE D'ECRITURE
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/codnop.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jjmmaa.h"
    character(len=1) :: prfnoe, prfnsy
!
!
!  --> DECLARATION DES ARGUMENTS
!
    integer(kind=8) :: nbnode, min, man, ites
    real(kind=8) :: ama, bma, cma, ami, bmi, cmi
!
!  --> DECLARATION DES VARIABLES LOCALES
!
    character(len=4) :: ct(3)
    character(len=8) :: chnode, chsc
    character(len=12) :: chenti, chnomi, chnoma, aut
    character(len=13) :: chlign, chlige
    character(len=80) :: chfone(4)
    real(kind=8) :: x, y, z
    integer(kind=8) :: nblit, nblie, nblif
    integer(kind=8) :: node, isc
!
!  --------- FIN DECLARATION ---------
!
!  --> N  D'UNITE LOGIQUE ASSOCIE AU FICHIER
!-----------------------------------------------------------------------
    integer(kind=8) :: i, imod
    real(kind=8), pointer :: coor(:) => null()
    integer(kind=8), pointer :: info(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    prfnoe = 'N'
    prfnsy = ' '
    chfone(1) = '%FORMAT=(1*NOM_DE_NOEUD,3*COORD)'
    chfone(2) = '%FORMAT=(1*NOM_DE_NOEUD,3*COORD,1*NOM_DE_SYSTEME)'
    chfone(3) = '%FORMAT=(1*NOM_DE_NOEUD,2*COORD)'
    chfone(4) = '%FORMAT=(1*NOM_DE_NOEUD,2*COORD,1*NOM_DE_SYSTEME)'
    chenti = 'NBOBJ=      '
    chlign = 'NBLIGT=      '
    chnode = '        '
    chlige = 'NBLIGE=      '
    chnomi = 'NUMIN=      '
    chnoma = 'NUMAX=      '
    chsc = '        '
!
    nblif = 1
!
! --- RECUPERATION DES VECTEURS DE TRAVAIL :
!     ------------------------------------
    call jeveuo('&&PRESUP.INFO.NOEUDS', 'L', vi=info)
    call jeveuo('&&PRESUP.COOR.NOEUDS', 'L', vr=coor)
!
    call codent(nbnode, 'G', chenti(7:12))
    call codent(min, 'G', chnomi(7:12))
    call codent(man, 'G', chnoma(7:12))
!
    if (ites .eq. 0) then
        nblie = 5
        nblit = nbnode+nblie+nblif+1
        call codent(nblie, 'G', chlige(8:13))
        call codent(nblit, 'G', chlign(8:13))
    else
        nblie = 3
        nblit = nbnode+nblie+nblif+1
        call codent(nblie, 'G', chlige(8:13))
        call codent(nblit, 'G', chlign(8:13))
    end if
!
    call jjmmaa(ct, aut)
!
    write (imod, '(A,4X,A,4X,A,3X,A,3X,A)') 'COOR_3D', 'NOM=INDEFINI',&
     &         chenti, chlige, chlign
!
    write (imod, '(11X,A,19X,A)') chnomi, chnoma
    write (imod, '(11X,2A,12X,A,A2,A,A2,A,A4)') 'AUTEUR=', aut, 'DATE=',&
     &          ct(1) (1:2), '/', ct(2) (1:2), '/', ct(3)
    if (ites .eq. 0) then
        write (imod, '(A,7X,3(A,E15.8,1X))') '%', 'XMAX=', ama, 'YMAX=', &
            bma, 'ZMAX=', cma
        write (imod, '(A,7X,3(A,E15.8,1X))') '%', 'XMIN=', ami, 'YMIN=', &
            bmi, 'ZMIN=', cmi
    end if
!
    if (ites .eq. 0) then
        write (imod, '(A)') chfone(1)
    else
        write (imod, '(A)') chfone(2)
    end if
!
    do i = 1, nbnode
        node = info((i-1)*3+1)
        isc = info((i-1)*3+2)
        x = coor((i-1)*3+1)
        y = coor((i-1)*3+2)
        z = coor((i-1)*3+3)
!
        call codnop(chnode, prfnoe, 1, 1)
        call codent(node, 'G', chnode(2:8))
!
        if (ites .eq. 0) then
            write (imod, '(2X,A,2X,3(1PE21.14,1X))') chnode, x, y, z
        else
            call codnop(chsc, prfnsy, 1, 1)
!
! ---> RENUMEROTATION DES SYSTEMES DE COORDONNEES C.A.D. ON
!      INCREMENTE DE 1 LES NUMEROS DONNEES PAR SUPERTAB
!
            isc = isc+1
!
            call codent(isc, 'G', chsc(2:8))
            write (imod, '(2X,A,2X,3(1PE21.14,1X),A)') chnode, x, y, z, &
                chsc
        end if
    end do
    write (imod, '(A)') 'FINSF'
    write (imod, '(A)') '%'
    call jedema()
end subroutine
