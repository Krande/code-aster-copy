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
subroutine ecrelt(imod, maxnod, nbtyma, nomail, nbmail, &
                  mint, mant, limail, nbmtot)
    implicit none
!     ===============================================================
!A PRESUPER
!
!  ==================================================================
!  !                                                                !
!  ! FONCTION: ECRITURE DES MAILLES DANS LE FICHIER MODELE          !
!  !                                                                !
!  ==================================================================
!  !                                                                !
!  !  ROUTINE APPELES : CODENT                                      !
!  !                         : IUNIFI (FONCTION)                    !
!  !                         : JJMMAA                               !
!  !                                                                !
!  !  ROUTINE APPELANTE : PRESUP                                    !
!  !                                                                !
!  ==================================================================
!  !                                                                !
!  !                 ***************                                !
!  !                 *  ARGUMENTS  *                                !
!  !                 ***************                                !
!  !                                                                !
!  !  ************************************************************* !
!  !  *   NOM   *  TYPE * MODE *ALTERE *        ROLE              * !
!  !  ************************************************************* !
!  !  *         *       *      *       *                          * !
!  !  * MAXNOD  *INTEGER*ENTREE* NON   * NBRE MAXI DE NOEUDS POUR * !
!  !  *         *       *      *       *  UNE MAILLE SUPERTAB     * !
!  !  *         *       *      *       *                          * !
!  !  * NBTYMA  *INTEGER*ENTREE* NON   * NBRE DE TYPES DE MAILLES * !
!  !  *         *       *      *       *     SUPERTAB             * !
!  !  *         *       *      *       *                          * !
!  !  * NOMAIL  *CHARACT*ENTREE* NON   * NOMS DES DIFFERENTS TYPES* !
!  !  *         *       *      *       * DE MAILLES LUS/LE FICHIER* !
!  !  *         *       *      *       *                          * !
!  !  * NBMAIL  *INTEGER*ENTREE* NON   * NBRE DE MAILLES DE CHAQUE* !
!  !  *         *       *      *       *      TYPE                * !
!  !  *         *       *      *       *                          * !
!  !  * MINT    *INTEGER*ENTREE* NON   * TABLEAU CONTENANT LE N DE* !
!  !  *         *       *      *       * MAILLE MIN POUR CHAQUE   * !
!  !  *         *       *      *       * TYPE DE MAILLE           * !
!  !  *         *       *      *       *                          * !
!  !  * MANT    *INTEGER*ENTREE* NON   * TABLEAU CONTENANT LE N DE* !
!  !  *         *       *      *       * MAILLE MAX POUR CHAQUE   * !
!  !  *         *       *      *       * TYPE DE MAILLE           * !
!  !  *         *       *      *       *                          * !
!  !  * LIMAIL  *INTEGER*ENTREE* NON   * TABLEAU CONTENANT LE NBRE* !
!  !  *         *       *      *       * DE LIGNE POUR L'ECRITURE * !
!  !  *         *       *      *       * DE CHAQUE TYPE DE MAILLE * !
!  !  *         *       *      *       *                          * !
!  !  ************************************************************* !
!  !                                                                !
!  ==================================================================
!
!
!  --> DECLARATION DES ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/codnop.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jjmmaa.h"
    integer(kind=8) :: maxnod, nbtyma
    character(len=1) :: prfnoe, prfmai
    character(len=8) :: nomail(nbtyma)
    integer(kind=8) :: nbmail(nbtyma)
    integer(kind=8) :: mint(nbtyma), mant(nbtyma), limail(nbtyma)
!
!  --> DECLARATION DES VARIABLES LOCALES
!
    character(len=4) :: ct(3)
    character(len=8) :: chmail, chtab(32)
    character(len=12) :: chnomi, chnoma, chenti, aut
    character(len=13) :: chlign, chlige
    character(len=80) :: chfoma(35)
    integer(kind=8) :: irest, idiv, nblit, nblie, nblif
    integer(kind=8) :: neu2(32), inum, ityp, neul
    integer(kind=8) :: l
!
!  --> DECLARATION DES INDICES DE BOUCLES
!
    integer(kind=8) :: in, nte
!
!  ------------- FIN DECLARATION ------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifo, ima, imod, ipos
    integer(kind=8) :: k, nbmtot
    integer(kind=8), pointer :: conn(:) => null()
    integer(kind=8), pointer :: info(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    prfnoe = 'N'
    prfmai = 'M'
    do ifo = 1, 35
        chfoma(ifo) = ' '
    end do
!
    chfoma(1) = '%FORMAT=(1*NOM_DE_MAILLE,2*NOM_DE_NOEUD)'
    chfoma(2) = '%FORMAT=(1*NOM_DE_MAILLE,3*NOM_DE_NOEUD)'
    chfoma(3) = '%FORMAT=(1*NOM_DE_MAILLE,6*NOM_DE_NOEUD)'
    chfoma(4) = '%FORMAT=(1*NOM_DE_MAIILE,9*NOM_DE_NOEUD)'
    chfoma(5) = '%FORMAT=(1*NOM_DE_MAILLE,4*NOM_DE_NOEUD)'
    chfoma(6) = '%FORMAT=(1*NOM_DE_MAILLE,8*NOM_DE_NOEUD)'
    chfoma(7) = '%FORMAT=(1*NOM_DE_MAILLE,12*NOM_DE_NOEUD)'
    chfoma(8) = '%FORMAT=(1*NOM_DE_MAILLE,6*NOM_DE_NOEUD)'
    chfoma(9) = '%FORMAT=(1*NOM_DE_MAILLE,12*NOM_DE_NOEUD)'
    chfoma(10) = '%FORMAT=(1*NOM_DE_MAILLE,18*NOM_DE_NOEUD)'
    chfoma(11) = '%FORMAT=(1*NOM_DE_MAILLE,8*NOM_DE_NOEUD)'
    chfoma(12) = '%FORMAT=(1*NOM_DE_MAILLE,16*NOM_DE_NOEUD)'
    chfoma(13) = '%FORMAT=(1*NOM_DE_MAILLE,24*NOM_DE_NOEUD)'
    chfoma(14) = '%FORMAT=(1*NOM_DE_MAILLE,4*NOM_DE_NOEUD)'
    chfoma(15) = '%FORMAT=(1*NOM_DE_MAILLE,10*NOM_DE_NOEUD)'
    chfoma(16) = '%FORMAT=(1*NOM_DE_MAILLE,6*NOM_DE_NOEUD)'
    chfoma(17) = '%FORMAT=(1*NOM_DE_MAILLE,15*NOM_DE_NOEUD)'
    chfoma(18) = '%FORMAT=(1*NOM_DE_MAILLE,24*NOM_DE_NOEUD)'
    chfoma(19) = '%FORMAT=(1*NOM_DE_MAILLE,8*NOM_DE_NOEUD)'
    chfoma(20) = '%FORMAT=(1*NOM_DE_MAILLE,20*NOM_DE_NOEUD)'
    chfoma(21) = '%FORMAT=(1*NOM_DE_MAILLE,32*NOM_DE_NOEUD)'
    chfoma(29) = '%FORMAT=(1*NOM_DE_MAILLE,2*NOM_DE_NOEUD)'
    chfoma(30) = '%FORMAT=(1*NOM_DE_MAILLE,1*NOM_DE_NOEUD)'
    chfoma(31) = '%FORMAT=(1*NOM_DE_MAILLE,2*NOM_DE_NOEUD)'
    chfoma(32) = '%FORMAT=(1*NOM_DE_MAILLE,1*NOM_DE_NOEUD)'
    chfoma(33) = '%FORMAT=(1*NOM_DE_MAILLE,1*NOM_DE_NOEUD)'
    chfoma(34) = '%FORMAT=(1*NOM_DE_MAILLE,2*NOM_DE_NOEUD)'
    chfoma(35) = '%FORMAT=(1*NOM_DE_MAILLE,3*NOM_DE_NOEUD)'
!
!     ECRITURE DES MAILLES
!     ====================
!
    chmail = '        '
    chenti = 'NBOBJ=      '
    chlign = 'NBLIGT=      '
    chlige = 'NBLIGE=      '
    chnomi = 'NUMIN=      '
    chnoma = 'NUMAX=      '
!
!  --> N  D'UNITES LOGIQUES DES FICHIERS
!
    call jeveuo('&&PRESUP.INFO.MAILLE', 'L', vi=info)
    call jeveuo('&&PRESUP.CONN.MAILLE', 'L', vi=conn)
    do in = 1, maxnod
        chtab(in) = '        '
    end do
!
    do nte = 1, nbtyma
        ipos = 0
        if (nbmail(nte) .eq. 0) goto 10
        nblie = 3
        nblif = 1
        nblit = nbmail(nte)*limail(nte)+nblie+nblif+1
!
        call codent(nbmail(nte), 'G', chenti(7:12))
        call codent(nblit, 'G', chlign(8:13))
        call codent(nblie, 'G', chlige(8:13))
        call codent(mint(nte), 'G', chnomi(7:12))
        call codent(mant(nte), 'G', chnoma(7:12))
!
!
!  --> ECRITURE DE LA DATE (IBM&CRAY)
        call jjmmaa(ct, aut)
!
!  --> ECRITURE DE L'ENTETE DANS LE FICHIER NEUTRE
        write (unit=imod, fmt='(A,3X,A,3X,A,3X,A,3X,A,3X,A)') nomail( &
            nte), 'NOM=INDEFINI', chenti, chlige, chlign
        write (imod, '(11X,A,18X,A)') chnomi, chnoma
        write (imod, '(11X,2A,11X,A,A2,A,A2,A,A4)') 'AUTEUR=', aut, &
            'DATE=', ct(1) (1:2), '/', ct(2) (1:2), '/', ct(3)
!
!  --> ECRITURE DU FORMAT DANS LE FICHIER NEUTRE
        write (imod, '(A)') chfoma(nte)
!
        do ima = 1, nbmtot
            ityp = info((ima-1)*4+2)
            neul = info((ima-1)*4+3)
            if (ityp .eq. nte) then
!
                inum = info((ima-1)*4+1)
                call codnop(chmail, prfmai, 1, 1)
                call codent(inum, 'G', chmail(2:8))
!
                do in = 1, neul
                    neu2(in) = conn(ipos+in)
                    call codnop(chtab(in), prfnoe, 1, 1)
                    call codent(neu2(in), 'G', chtab(in) (2:8))
                end do
!
                idiv = int(neul/8)
                irest = mod(neul, 8)
                if (irest .ne. 0) then
                    write (imod, 202) chmail, (chtab(i), i=1, neul)
                else
                    do k = 1, idiv
                        l = 8*(k-1)
                        if (idiv .eq. 1) then
                            write (imod, '(A,8(1X,A))') chmail, (chtab( &
                                                                 i), i=1, 8)
                        else
                            write (imod, '(8X,8(1X,A))') (chtab(i), i=1+ &
                                                          l, 8+l)
                        end if
                    end do
                end if
            end if
            ipos = ipos+neul
        end do
        write (imod, '(A)') 'FINSF'
        write (imod, '(A)') '%'
10      continue
    end do
202 format(a, 8(1x, a), /, (8x, 8(1x, a)))
    call jedema()
end subroutine
