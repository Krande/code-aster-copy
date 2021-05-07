! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine te0047(optioz, nomtez)
!
! --------------------------------------------------------------------------------------------------
!
!         COMPORTEMENT LINÉAIRE ET NON-LINÉAIRE POUR LES DISCRETS
!
! --------------------------------------------------------------------------------------------------
!
! Éléments concernés
!     MECA_DIS_TR_L     : sur une maille a 2 noeuds
!     MECA_DIS_T_L      : sur une maille a 2 noeuds
!     MECA_DIS_TR_N     : sur une maille a 1 noeud
!     MECA_DIS_T_N      : sur une maille a 1 noeud
!     MECA_2D_DIS_TR_L  : sur une maille a 2 noeuds
!     MECA_2D_DIS_T_L   : sur une maille a 2 noeuds
!     MECA_2D_DIS_TR_N  : sur une maille a 1 noeud
!     MECA_2D_DIS_T_N   : sur une maille a 1 noeud
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       optioz : nom de l'option a calculer
!       nomtez : nom du type_element
!
! --------------------------------------------------------------------------------------------------
!
use Behaviour_module, only : behaviourOption
use te0047_type
!
implicit none
character(len=*) :: optioz, nomtez
!
#include "jeveux.h"
#include "asterfort/diarm0.h"
#include "asterfort/dibili.h"
#include "asterfort/dis_contact_frot.h"
#include "asterfort/dicho0.h"
#include "asterfort/dicora.h"
#include "asterfort/didashpot.h"
#include "asterfort/diecci.h"
#include "asterfort/dielas.h"
#include "asterfort/digou2.h"
#include "asterfort/digric.h"
#include "asterfort/dichoc_endo_pena.h"
#include "asterfort/diisotrope.h"
#include "asterfort/difondabb.h"
#include "asterfort/dichoc_endo_ldc.h"
#include "asterfort/dizeng.h"
#include "asterfort/disjvp.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dcopy.h"
!
! --------------------------------------------------------------------------------------------------
!
integer :: neq, ideplm, ideplp, icompo
integer :: ii, jcret, itype, lorien
integer :: iadzi, iazk24, ibid, infodi, codret
!
real(kind=8) :: ang(3), ugm(12), dug(12), r8bid
!
character(len=8)  :: k8bid
character(len=24) :: messak(5)
!
type(te0047_dscr) :: for_discret
!
! --------------------------------------------------------------------------------------------------
!
    for_discret%option = optioz
    for_discret%nomte  = nomtez
    ! on vérifie que les caractéristiques ont été affectées
    ! le code du discret
    call infdis('CODE', ibid, r8bid, for_discret%nomte)
    ! le code stocké dans la carte
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=for_discret%nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    endif
    ! discret de type raideur
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=for_discret%nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    endif
    ! les discrets sont obligatoirement symétriques
    call infdis('SYMK', infodi, r8bid, k8bid)
    if (infodi .ne. 1) then
        call utmess('F', 'DISCRETS_40')
    endif
    ! informations sur les discrets :
    !     nbt   = nombre de coefficients dans k
    !     nno   = nombre de noeuds
    !     nc    = nombre de composante par noeud
    !     ndim  = dimension de l'élément
    !     itype = type de l'élément
    call infted(for_discret%nomte, infodi, &
                for_discret%nbt, for_discret%nno, for_discret%nc, for_discret%ndim, itype)
    neq = for_discret%nno*for_discret%nc
    ! récupération des adresses jeveux
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    !
    ! Récupération des infos concernant les comportements :
    !     rela_comp   zk16(icompo-1+RELA_NAME)   NOM_DU_COMPORTEMENT
    !     nbvar       zk16(icompo-1+NVAR)        nbvar = read (zk16(icompo-1+NVAR),'(i16)')
    !     defo_comp   zk16(icompo-1+DEFO)        PETIT   PETIT_REAC  GROT_GDEP
    !     type_comp   zk16(icompo-1+INCRELAS)    COMP_ELAS   COMP_INCR
    !
    ! Properties of behaviour
    call jevech('PCOMPOR', 'L', icompo)
    for_discret%rela_comp = zk16(icompo-1+RELA_NAME)
    for_discret%defo_comp = zk16(icompo-1+DEFO)
    for_discret%type_comp = zk16(icompo-1+INCRELAS)
    !
    ! Select objects to construct from option name
    call behaviourOption(for_discret%option, zk16(icompo), &
                         for_discret%lMatr, for_discret%lVect, &
                         for_discret%lVari, for_discret%lSigm , codret)
    for_discret%lMatrPred = for_discret%option .eq. 'RIGI_MECA_TANG'
    !
    if (for_discret%defo_comp .ne. 'PETIT') then
        call utmess('A', 'DISCRETS_18')
    endif
    ! si COMP_ELAS alors comportement elas
    if ((for_discret%type_comp.eq.'COMP_ELAS') .and. (for_discret%rela_comp.ne.'ELAS')) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_8', nk=5, valk=messak)
    endif
    ! dans les cas *_ELAS, les comportements qui ont une matrice de
    ! décharge sont : elas DIS_GRICRA. pour tous les autres cas : <f>
    if ((for_discret%option(10:14).eq.'_ELAS') .and. (for_discret%rela_comp.ne.'ELAS') .and.&
        (for_discret%rela_comp.ne.'DIS_GRICRA')) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_10', nk=5, valk=messak)
    endif
    !
    ! récupération des orientations (angles nautiques -> vecteur ang)
    ! orientation de l'élément et déplacements dans les repères g et l
    call tecach('ONO', 'PCAORIE', 'L', codret, iad=lorien)
    if (codret .ne. 0) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_6', nk=5, valk=messak)
    endif
    call dcopy(3, zr(lorien), 1, ang, 1)
    ! déplacements dans le repère global :
    !     ugm = déplacement précédent
    !     dug = incrément de déplacement
    do ii = 1, neq
        ugm(ii) = zr(ideplm+ii-1)
        dug(ii) = zr(ideplp+ii-1)
    enddo
    !
    ! matrice pgl de passage repère global -> repère local
    call matrot(ang, for_discret%pgl)
    ! déplacements dans le repère local :
    !     ulm = déplacement précédent    = pgl * ugm
    !     dul = incrément de déplacement = pgl * dug
    if (for_discret%ndim .eq. 3) then
        call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugm, for_discret%ulm)
        call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, dug, for_discret%dul)
    else if (for_discret%ndim.eq.2) then
        call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugm, for_discret%ulm)
        call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, dug, for_discret%dul)
    endif
    !
    codret = 0
    if (for_discret%rela_comp .eq. 'ELAS') then
        ! comportement élastique
        call dielas(for_discret, codret)
    else if (for_discret%rela_comp .eq. 'DASHPOT') then
        ! comportement dashpot
        call didashpot(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_VISC') then
        ! comportement de type zener (raideurs série et parallèle + amortisseur)
        call dizeng(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_ECRO_TRAC') then
        ! comportement avec ecrouissage isotrope à partir de la courbe force-déplacement
        call diisotrope(for_discret, codret)
    else if (for_discret%rela_comp(1:10).eq.'DIS_GOUJ2E') then
        ! comportement pour les goujons :  DIS_GOUJ2E_PLAS DIS_GOUJ2E_ELAS
        call digou2(for_discret, codret)
    else if (for_discret%rela_comp.eq.'ARME') then
        ! comportement pour les armements
        call diarm0(for_discret, codret)
    else if (for_discret%rela_comp.eq.'ASSE_CORN') then
        ! comportement pour les assemblages boulonés de cornière
        call dicora(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_GRICRA') then
        ! comportement des liaisons grille-crayon combustible
        call digric(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_CHOC') then
        ! comportement choc sans frottement de coulomb et sans amortissement
        call dicho0(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_CONTACT') then
        ! comportement choc avec frottement de coulomb avec amortissement
        call dis_contact_frot(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_ECRO_CINE') then
        ! comportement avec écrouissage cinématique
        call diecci(for_discret, codret)
    else if (for_discret%rela_comp.eq.'DIS_BILI_ELAS') then
        ! comportement élastique bi-linéaire
        call dibili(for_discret, codret)
    else if (for_discret%rela_comp.eq.'FONDATION') then
        ! comportement fondation : superficielle
        call difondabb(for_discret, codret)
    else if (for_discret%rela_comp.eq.'CHOC_ENDO') then
        ! comportement de choc avec déformation résiduelle
        call dichoc_endo_ldc(for_discret, codret)
    else if (for_discret%rela_comp.eq.'CHOC_ENDO_PENA') then
        ! comportement de choc avec déformation résiduelle par pénalisation
        call dichoc_endo_pena(for_discret, codret)
    else if (for_discret%rela_comp.eq.'JONC_ENDO_PLAS') then
        ! comportement élasto-plastique endommageable : jonction voile-plancher
        call disjvp(for_discret, codret)
    else
        ! si on passe par ici c'est qu'aucun comportement n'est valide
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_7', nk=5, valk=messak)
    endif
    ! les comportements valides passent par ici
    if (for_discret%lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = codret
    endif
!
end subroutine
