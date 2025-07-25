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
subroutine erglth(champ, inst, niveau, iordr, resuco)
! ERREUR GLOBALE AU MAILLAGE - THERMIQUE
! **     **                    **
! =====================================================================
!     BUT :  EN THERMIQUE, CALCULER LES ESTIMATEURS GLOBAUX
!            A PARTIR DES ESTIMATEURS LOCAUX CONTENUS DANS CHAMP
!
! IN  CHAMP    :  NOM DU CHAM_ELEM_ERREUR
! IN  INSTANT  :  INSTANT DE CALCUL
! IN  NIVEAU   :  NIVEAU DE L'ESTIMATEUR
! IN  IORDR    :  NUMERO D'ORDRE
! IN  RESUCO   :  SD RESULTAT.
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES:
!       MESSAGE : INFNIV.
!       JEVEUX  : JEMARQ,JELIRA,JEVEUO,JEDEMA.
!       ASTER   : IUNIFI,CELVER.
!       ENVIMA  : R8MIEM.
!
!     FONCTIONS INTRINSEQUES:
!       ABS,SQRT.
!  -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       05/07/01 (OB) : CREATION EN S'INSPIRANT DE ERGLOB.F.
!       12/09/02 (OB) : MODIF. MSG D'ALARME DE LA DIVISION PAR ZERO.
! --------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/celver.h"
#include "asterfort/digdel.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/utmess.h"
    real(kind=8) :: inst
    integer(kind=8) :: niveau, iordr
    character(len=8) :: resuco
    character(len=*) :: champ
!
! ------------------------------------------------------------------
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: ifi, longt, long2, mode, j, ibid, nbgr, icoef, nel
    integer(kind=8) :: idecgr, k, iad, iavale, nbel
    real(kind=8) :: termvo, termsa, termfl, termec, ovfl, terms1, termf1, terme1
    real(kind=8) :: termv1, termv2, terms2, termf2, terme2, err0, nors, nu0
    character(len=4) :: docu
    character(len=19) :: champ2, ligrel
    aster_logical :: first
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
!
! INIT.
    call jemarq()
    ovfl = r8miem()
    ifi = iunifi('RESULTAT')
    champ2 = champ
!
! ON RETROUVE LE NOM DU LIGREL:
!     -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
    call celver(champ2, 'NBVARI_CST', 'STOP', ibid)
    call celver(champ2, 'NBSPT_1', 'STOP', ibid)
    call jelira(champ2//'.CELD', 'DOCU', cval=docu)
    if (docu .ne. 'CHML') then
        call utmess('F', 'CALCULEL5_44')
    end if
    call jeveuo(champ2//'.CELK', 'L', vk24=celk)
    ligrel = celk(1) (1:19)
!
    call jeveuo(champ2//'.CELD', 'L', vi=celd)
!
!     -- ON VERIFIE LES LONGUEURS:
    first = .true.
    nbgr = nbgrel(ligrel)
    do j = 1, nbgr
        mode = celd(celd(4+j)+2)
        if (mode .eq. 0) goto 1
        long2 = digdel(mode)
        icoef = max(1, celd(4))
        long2 = long2*icoef
        if (first) then
            longt = long2
        else
            if (longt .ne. long2) then
                call utmess('F', 'CALCULEL5_45')
            end if
        end if
        first = .false.
1       continue
    end do
!
!        -- ON CUMULE :
    call jeveuo(champ2//'.CELV', 'E', iavale)
!
    err0 = 0.d0
    nors = 0.d0
    nbel = 0
    if (niveau .eq. 2) then
        termvo = 0.d0
        termsa = 0.d0
        termfl = 0.d0
        termec = 0.d0
        termv1 = 0.d0
        terms1 = 0.d0
        termf1 = 0.d0
        terme1 = 0.d0
        termv2 = 0.d0
        terms2 = 0.d0
        termf2 = 0.d0
        terme2 = 0.d0
    end if
!
    do j = 1, nbgr
        mode = celd(celd(4+j)+2)
        if (mode .eq. 0) goto 2
        nel = nbelem(ligrel, j)
        idecgr = celd(celd(4+j)+8)
        do k = 1, nel
            iad = iavale-1+idecgr+(k-1)*longt
            err0 = err0+zr(iad)**2
            nors = nors+zr(iad+2)**2
            if (niveau .eq. 2) then
                termvo = termvo+zr(iad+3)**2
                termv1 = termv1+zr(iad+5)**2
                termsa = termsa+zr(iad+6)**2
                terms1 = terms1+zr(iad+8)**2
                termfl = termfl+zr(iad+9)**2
                termf1 = termf1+zr(iad+11)**2
                termec = termec+zr(iad+12)**2
                terme1 = terme1+zr(iad+14)**2
            end if
            nbel = nbel+1
        end do
2       continue
    end do
    err0 = sqrt(err0)
    nors = sqrt(nors)
    if (niveau .eq. 2) then
! ERREURS PARTIELLES ABSOLUES
        termvo = sqrt(termvo)
        termsa = sqrt(termsa)
        termfl = sqrt(termfl)
        termec = sqrt(termec)
        termv1 = sqrt(termv1)
        terms1 = sqrt(terms1)
        termf1 = sqrt(termf1)
        terme1 = sqrt(terme1)
! ERREURS PARTIELLES RELATIVES
        if (termv1 .gt. ovfl) termv2 = 100.d0*(termvo/termv1)
        if (terms1 .gt. ovfl) terms2 = 100.d0*(termsa/terms1)
        if (termf1 .gt. ovfl) termf2 = 100.d0*(termfl/termf1)
        if (terme1 .gt. ovfl) terme2 = 100.d0*(termec/terme1)
    end if
    if (nors .gt. ovfl) then
        nu0 = 100.d0*err0/nors
    else
        call utmess('I', 'CALCULEL5_46')
        nu0 = 0.d0
    end if
    write (ifi, *) ' '
    write (ifi, *) '**********************************************'
    write (ifi, *) ' THERMIQUE: ESTIMATEUR D''ERREUR EN RESIDU '
    write (ifi, *) '**********************************************'
    write (ifi, *)
    write (ifi, *) '   IMPRESSION DES NORMES GLOBALES :'
    write (ifi, *)
!
! ESTIMATEURS D'ERREURS EN THERMIQUE LINEAIRE
    write (ifi, 111) ' SD EVOL_THER    ', resuco
    write (ifi, 110) ' NUMERO D''ORDRE  ', iordr
    write (ifi, 109) ' INSTANT         ', inst
    write (ifi, *) 'ERREUR             ABSOLUE   /  RELATIVE '//&
     &      '/ NORMALISATION'
    write (ifi, 108) ' TOTAL           ', err0, nu0, '%', nors
    if (niveau .eq. 2) then
        write (ifi, 108) ' TERME VOLUMIQUE ', termvo, termv2, '%', termv1
        write (ifi, 108) ' TERME SAUT      ', termsa, terms2, '%', terms1
        write (ifi, 108) ' TERME FLUX      ', termfl, termf2, '%', termf1
        write (ifi, 108) ' TERME ECHANGE   ', termec, terme2, '%', terme1
    end if
108 format(a17, d16.8, 1x, d16.8, a2, 1x, d16.8)
109 format(a17, d16.8)
110 format(a17, i5)
111 format(a17, a8)
    write (ifi, *)
    write (ifi, *) '**********************************************'
    call jedema()
end subroutine
