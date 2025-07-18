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

subroutine lissav(lischa, ichar, charge, typech, genrec, &
                  motclc, prefob, typapp, nomfct, typfct, &
                  phase, npuis)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=19) :: lischa
    character(len=13) :: prefob
    integer(kind=8) :: ichar
    character(len=8) :: charge
    character(len=16) :: typapp, typfct
    integer(kind=8) :: genrec, motclc(2)
    character(len=8) :: typech, nomfct
    real(kind=8) :: phase
    integer(kind=8) :: npuis
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! SAUVEGARDE DANS LA SD LISTE_CHARGES
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : NOM DE LA SD LISTE_CHARGES
! IN  ICHAR  : INDICE DE LA CHARGE
! IN  CHARGE : NOM DE LA CHARGE (AFFE_CHAR_*)
! IN  TYPECH : TYPE DE LA CHARGE
!               'REEL'    - CHARGE CONSTANTE REELLE
!               'COMP'    - CHARGE CONSTANTE COMPLEXE
!               'FONC_F0' - CHARGE FONCTION QUELCONQUE
!               'FONC_FT' - CHARGE FONCTION DU TEMPS
! IN  GENREC : CODE (ENTIER CODE) CONTENANT LES GENRES
! IN  MOTCLC : CODE (ENTIER CODE) CONTENANT LES MOTS-CLEFS
! IN  PREFOB : PREFIXE DE L'OBJET DE LA CHARGE
! IN  TYPAPP : TYPE D'APPLICATION DE LA CHARGE
!              FIXE_CSTE
!              FIXE_PILO
!              SUIV
!              DIDI
! IN  NOMFCT : NOM DE LA FONCTION MULTIPLICATRICE
! IN  TYPFCT : TYPE DE LA FONCTION MULTIPLICATRICE
!              'FONCT_REEL' FONCTION MULTIPLICATRICE REELLE
!              'FONCT_COMP' FONCTION MULTIPLICATRICE COMPLEXE
!              'CONST_REEL' FONCTION MULTIPLICATRICE CONSTANTE REELLE
!              'CONST_COMP' FONCTION MULTIPLICATRICE CONSTANTE COMPLEXE
! IN  PHASE  : PHASE POUR LES FONCTIONS MULTIPLICATRICES COMPLEXES
! IN  NPUIS  : PUISSANCE POUR LES FONCTIONS MULTIPLICATRICES COMPLEXES
!
! ----------------------------------------------------------------------
!
    character(len=24) :: nomcha, genrch, typcha, typeap, precha, mocfch
    integer(kind=8) :: jncha, jgenc, jtypc, jtypa, jprec, jmcfc
    character(len=24) :: nomfon, typfon, valfon
    integer(kind=8) :: jnfon, jtfon, jvfon
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nomcha = lischa(1:19)//'.NCHA'
    genrch = lischa(1:19)//'.GENC'
    mocfch = lischa(1:19)//'.MCFC'
    typcha = lischa(1:19)//'.TYPC'
    typeap = lischa(1:19)//'.TYPA'
    precha = lischa(1:19)//'.PREO'
    nomfon = lischa(1:19)//'.NFON'
    typfon = lischa(1:19)//'.TFON'
    valfon = lischa(1:19)//'.VFON'
!
    call jeveuo(nomcha, 'E', jncha)
    call jeveuo(genrch, 'E', jgenc)
    call jeveuo(mocfch, 'E', jmcfc)
    call jeveuo(typcha, 'E', jtypc)
    call jeveuo(typeap, 'E', jtypa)
    call jeveuo(precha, 'E', jprec)
    call jeveuo(nomfon, 'E', jnfon)
    call jeveuo(typfon, 'E', jtfon)
    call jeveuo(valfon, 'E', jvfon)
!
    zk8(jncha-1+ichar) = charge
    zi(jgenc-1+ichar) = genrec
    zi(jmcfc-1+2*(ichar-1)+1) = motclc(1)
    zi(jmcfc-1+2*(ichar-1)+2) = motclc(2)
    zk8(jtypc-1+ichar) = typech
    zk16(jtypa-1+ichar) = typapp
    zk24(jprec-1+ichar) = prefob
    zk8(jnfon-1+ichar) = nomfct
    zk16(jtfon-1+ichar) = typfct
    zr(jvfon-1+2*(ichar-1)+1) = phase
    zr(jvfon-1+2*(ichar-1)+2) = npuis
!
    call jedema()
end subroutine
