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

subroutine crelrl(typcoz, typvaz, basez, lisrez)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
    character(len=1) :: base
    character(len=4) :: typcoe, typval
    character(len=19) :: lisrel
    character(len=*) :: typcoz, typvaz, basez, lisrez
! ------------------------------------------------------------
!     CREATION DE L'OBJET DE TYPE LISTE_DE_RELATIONS
!     DE NOM LISREL
!        LE NOM LISREL EST FOURNI EN ARGUMENT
!     SI L'OBJET LISREL EXISTE DEJA, ON LE DETRUIT
!                               PUIS ON LE RECREE
!
!     LA TAILLE DES VECTEURS CONTENANT LES COMPOSANTES
!     DES RELATIONS EST DIMENSIONNEE A LVECRL = 10000
!
!     LA TAILLE DES VECTEURS QUI CORRESPOND AUX NOMBRES
!     DE RELATIONS EST DIMENSIONNEE A NBRELA = 1000
!
!
!-------------------------------------------------------------
! TYPCOZ        - IN - K4  - : TYPE DES COEFFICIENTS DE LA RELATION
!               -    -     -   = 'REEL' OU 'COMP'
!-------------------------------------------------------------
!  TYPVAZ       - IN - K4  - : INDICATEUR DU TYPE DU SECOND MEMBRE
!               -    -     -     = 'REEL' OU 'COMP' OU 'FONC'
!-------------------------------------------------------------
!  BASEZ        - IN - K1  - : INDICATEUR DE LA BASE SUR LAQUELLE
!                              DOIT ETRE CREE L'OBJET DE TYPE
!                              LISTE_DE_RELATIONS
!-------------------------------------------------------------
!  LISREZ     - IN    - K19 - : NOM DE LA LISTE_RELA
!             - JXOUT -     -
!-------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: idbeta, idcoef, iddl, idnbre, idnoeu, idpoin, idsurc
    integer(kind=8) :: idterm, idtyco, idtyva, iret, lvecrl, nbrela
!
!-----------------------------------------------------------------------
    parameter(lvecrl=10000)
    parameter(nbrela=1000)
!
    call jemarq()
!
! --- SI L'OBJET LISREL EXISTE , ON LE DETRUIT ---
!
    typcoe = typcoz
    typval = typvaz
    base = basez
    lisrel = lisrez
!
    call jeexin(lisrel//'.RLCO', iret)
    if (iret .ne. 0) then
        call jedetr(lisrel//'.RLCO')
        call jedetr(lisrel//'.RLDD')
        call jedetr(lisrel//'.RLNO')
        call jedetr(lisrel//'.RLBE')
        call jedetr(lisrel//'.RLNT')
        call jedetr(lisrel//'.RLPO')
        call jedetr(lisrel//'.RLNR')
        call jedetr(lisrel//'.RLSU')
        call jedetr(lisrel//'.RLTC')
        call jedetr(lisrel//'.RLTV')
        call jedetr(lisrel//'.RLBE')
    end if
!
! ---  VECTEUR DES COEFFICIENTS DES TERMES DES RELATIONS
!
    if (typcoe .eq. 'COMP') then
        call wkvect(lisrel//'.RLCO', base//' V C', lvecrl, idcoef)
    else
        call wkvect(lisrel//'.RLCO', base//' V R', lvecrl, idcoef)
    end if
!
! ---  VECTEUR DES NOMS DES DDLS IMPLIQUES DANS LES RELATIONS
!
    call wkvect(lisrel//'.RLDD', base//' V K8', lvecrl, iddl)
!
! ---  VECTEUR DES NOMS DES NOEUDS IMPLIQUES DANS LES RELATIONS
!
    call wkvect(lisrel//'.RLNO', base//' V K8', lvecrl, idnoeu)
!
! ---  VECTEUR DES VALEURS DES SECONDS MEMBRES DE CHAQUE RELATION
!
    if (typval .eq. 'REEL') then
        call wkvect(lisrel//'.RLBE', base//' V R', nbrela, idbeta)
    else if (typval .eq. 'COMP') then
        call wkvect(lisrel//'.RLBE', base//' V C', nbrela, idbeta)
    else if (typval .eq. 'FONC') then
        call wkvect(lisrel//'.RLBE', base//' V K24', nbrela, idbeta)
    end if
!
! ---  VECTEUR DES NOMBRES DE TERMES DE CHAQUE RELATION
!
    call wkvect(lisrel//'.RLNT', base//' V I', nbrela, idterm)
!
! ---  VECTEUR DES NOMBRES DE TERMES CUMULES DES RELATIONS
! ---  ( SERT DE POINTEUR DANS LES TABLEAUX RELATIFS AUX TERMES
! ---    DES RELATIONS)
!
    call wkvect(lisrel//'.RLPO', base//' V I', nbrela, idpoin)
!
! ---  VECTEUR D'INDICATEURS DE PRISE EN COMPTE DES RELATIONS
! ---  ( SERT A APPLIQUER LA REGLE DE SURCHARGE)
!
    call wkvect(lisrel//'.RLSU', base//' V I', nbrela, idsurc)
!
! ---  NOMBRE DE RELATIONS
!
    call wkvect(lisrel//'.RLNR', base//' V I', 1, idnbre)
    zi(idnbre) = 0
!
! ---  TYPE DES COEFFICIENTS DES RELATIONS (TYPCOE)
!
    call wkvect(lisrel//'.RLTC', base//' V K8', 1, idtyco)
    zk8(idtyco) = typcoe(1:4)//'    '
!
! ---  TYPE DES VALEURS DES RELATIONS (TYPVAL)
!
    call wkvect(lisrel//'.RLTV', base//' V K8', 1, idtyva)
    zk8(idtyva) = typval(1:4)//'    '
!
    call jedema()
end subroutine
