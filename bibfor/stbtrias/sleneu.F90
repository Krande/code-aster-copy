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
subroutine sleneu(iunv, nbnode, ama, bma, cma, &
                  ami, bmi, cmi, mix, man, &
                  ites, datset)
! aslint: disable=
    implicit none
!     ==============================================================
!A PRESUPER
!
!     ==============================================================
!     !                                                            !
!     !  FONCTION:LECTURE SUR LE FICHIER UNIVERSEL ISSU DE SUPER-  !
!     !           TAB I-DEAS 4.0, 6.0 OU 7.0   DES COORDONNEES DES !
!     !           DES NOEUDS ET STOCKAGE DANS OBJETS JEVEUX        !
!     !                                                            !
!     ==============================================================
!     !                                                            !
!     !  ROUTINES APPELES : CODENT                                 !
!     !                          : IUNIFI (FONCTION)               !
!     !                                                            !
!     !  ROUTINE APPELANTE : PRESUP                                !
!     !                                                            !
!     ==============================================================
!     !                                                            !
!     !                  **************                            !
!     !                  *  ARGUMENT  *                            !
!     !                  **************                            !
!     !                                                            !
!     !  ********************************************************  !
!     !  *   NOM    *  TYPE * MODE *ALTERE *      ROLE          *  !
!     !  ********************************************************  !
!     !  *          *       *      *       *                    *  !
!     !  * NBNODE   *INTEGER*SORTIE* NON   *NBRE TOTAL DE NOEUDS*  !
!     !  *          *       *      *       *                    *  !
!     !  * AMA      *D.PRECI*SORTIE* NON   * X(MAXIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * BMA      *D.PRECI*SORTIE* NON   * Y(MAXIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * CMA      *D.PRECI*SORTIE* NON   * Z(MAXIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * AMI      *D.PRECI*SORTIE* NON   * X(MINIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * BMI      *D.PRECI*SORTIE* NON   * Y(MINIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * CMI      *D.PRECI*SORTIE* NON   * Z(MINIMUM)         *  !
!     !  *          *       *      *       *                    *  !
!     !  * MAN      *INTEGER*SORTIE* NON   * N DE NOEUD MAXI    *  !
!     !  *          *       *      *       *                    *  !
!     !  * MIX      *INTEGER*SORTIE* NON   * N DE NOEUD MINI    *  !
!     !  *          *       *      *       *                    *  !
!     !  * ITES     *INTEGER*SORTIE* NON   * INDIQUE S'IL EXISTE*  !
!     !  *          *       *      *       * AU MOINS DE SYST.  *  !
!     !  *          *       *      *       * DE COORDONNEES     *  !
!     !  *          *       *      *       *                    *  !
!     !  * DATSET   *INTEGER*ENTREE* NON   * NUMERO DU DATASET  *  !
!     !  *          *       *      *       * TRAITE             *  !
!     !  *          *       *      *       *                    *  !
!     !  ********************************************************  !
!     !                                                            !
!     ==============================================================
!  --> DECLARATION DES ARGUMENTS
#include "jeveux.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbnode, mix, man, ites, datset
    real(kind=8) :: ama, bma, cma, ami, bmi, cmi
!  --> DECLARATION DES VARIABLES LOCALES
    character(len=80) :: cbuf
    integer(kind=8) :: node, i, icnode, ind, j, itmp, inus
    real(kind=8) :: x, y, z
!
!  ---------- FIN DECLARATION -----------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: imes, ire1, ire2, iret, isyst, iter
    integer(kind=8) :: iunv, jcoor, jinfo, ndeca, niter, ntail
    integer(kind=8), pointer :: syst(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    nbnode = 0
    ites = 0
!
!  --> N  D'UNITE LOGIQUE ASSOCIE AUX FICHIERS
    imes = iunifi('MESSAGE')
!
    ntail = 1000
    niter = 1000
    ndeca = 0
    itmp = -6
    inus = 10
!
    call jeexin('&&PRESUP.INFO.NOEUDS', ire1)
    if (ire1 .ne. 0) call jedetr('&&PRESUP.INFO.NOEUDS')
    call wkvect('&&PRESUP.INFO.NOEUDS', 'V V I', 3*ntail, jinfo)
    call jeexin('&&PRESUP.COOR.NOEUDS', ire2)
    if (ire2 .ne. 0) call jedetr('&&PRESUP.COOR.NOEUDS')
    call wkvect('&&PRESUP.COOR.NOEUDS', 'V V R', 3*ntail, jcoor)
!
! -->  GESTION DES SYSTEMES DE COORDONNEES - PREMIERE PARTIE
! -->  ON TESTE L'EXISTENCE D'UN SYSTEME DE COORDONNEES (DATASET 2420)
!
    call jeexin('&&IDEAS.SYST', iret)
    if (iret .eq. 0) then
        call utmess('I', 'STBTRIAS_9')
!     Il n'y a pas de sys de coord defini dans le fichier, pour ne pas
!     planter on en cree un bidon ici qu'on declare comme cartesien
        isyst = 0
    end if
!
1   continue
    do iter = 1, niter
        read (iunv, '(A)') cbuf
        read (unit=cbuf, fmt='(4X,I2)') ind
        if (ind .eq. -1) goto 99
!
! --> LES COORDONNEES DES NOEUDS SONT EN SIMPLE PRECISION (SUPERTAB 4
!     OU 6) OU EN REAL*8 (SUPERTAB 6 OU 7)
!
        if (datset .eq. 15) then
            read (cbuf, '(4I10,3E13.6)') node, i, j, icnode, x, y, z
        else if (datset .eq. 781 .or. datset .eq. 2411) then
            read (cbuf, '(4I10)') node, i, j, icnode
            read (iunv, '(3E25.16)') x, y, z
        end if
!
!
! -->  GESTION DES SYSTEMES DE COORDONNEES - BIS
! -->  SI UN SYSTEME EST DEFINI, ON LE RECUPERE
!
        call jeexin('&&IDEAS.SYST', iret)
        if (iret .ne. 0) then
            call jeveuo('&&IDEAS.SYST', 'L', vi=syst)
            isyst = syst(i)
        end if
!
!        On ne teste ici que si le systeme de coordonnne est cartesien,
!        cylindrique ou autre
        if (isyst .ne. 0) then
            call utmess('F', 'STBTRIAS_10')
        end if
!        On ne teste ici que si les noeuds font reference a plusieurs
!        systeme de coordonnne
!        On ne teste pas si ces systemes sont identiques juste si leur
!        label est different
        if ((i .ne. itmp) .and. (itmp .ne. -6)) then
            call utmess('A', 'STBTRIAS_11')
        end if
!
!
!  --> INITIALISATION POUR LA RECHERCHE DES MINI ET MAXI
        if (nbnode .eq. 0) then
            ama = x
            bma = y
            cma = z
            ami = x
            bmi = y
            cmi = z
        else
            ama = max(ama, x)
            bma = max(bma, y)
            cma = max(cma, z)
            ami = min(ami, x)
            bmi = min(bmi, y)
            cmi = min(cmi, z)
        end if
!
        if (nbnode .eq. 0) then
            mix = node
        else
            man = max(mix, node)
        end if
!
        nbnode = nbnode+1
!
        zi(jinfo-1+ndeca+(iter-1)*3+1) = node
        zi(jinfo-1+ndeca+(iter-1)*3+2) = i
        zi(jinfo-1+ndeca+(iter-1)*3+3) = icnode
        zr(jcoor-1+ndeca+(iter-1)*3+1) = x
        zr(jcoor-1+ndeca+(iter-1)*3+2) = y
        zr(jcoor-1+ndeca+(iter-1)*3+3) = z
    end do
    ntail = ntail+niter
    ndeca = ndeca+3000
    call juveca('&&PRESUP.INFO.NOEUDS', 3*ntail)
    call jeveuo('&&PRESUP.INFO.NOEUDS', 'E', jinfo)
    call juveca('&&PRESUP.COOR.NOEUDS', 3*ntail)
    call jeveuo('&&PRESUP.COOR.NOEUDS', 'E', jcoor)
    goto 1
99  continue
!
    write (imes, *) 'NOMBRE DE NOEUDS :', nbnode
    call jedema()
end subroutine
