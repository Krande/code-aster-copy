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

subroutine w190af(modele, chmar1)
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: modele
    character(len=19) :: chmar1
!
! BUT : CREER LE CHAMP DE DONNEES POUR CALC_FERRAILLAGE
!
!-------------------------------------------------------------------------------------------------
    integer(kind=8) :: gd, nocc, ncmpmx, nbtou
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15
    integer(kind=8) :: n16, n17, n18, n20
    integer(kind=8) ::   jmail, iocc, nbmail
    real(kind=8) :: valrcb, valrco, valruc
    character(len=8) :: k8b, typmcl(2), noma, typcb, clacier
    character(len=8) :: typdiag, typstru, unitc
    character(len=16) :: motcls(2), typco
    character(len=24) :: mesmai
    character(len=8), pointer :: ncmp(:) => null()
    real(kind=8), pointer :: valv(:) => null()
!     ---------------------------------------------------------------------------------------------
    call jemarq()
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    ASSERT(noma .ne. ' ')
!
    call getfac('AFFE', nocc)
!
    mesmai = '&&W190AF.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
!     1- ALLOCATION DU CHAMP CHMAR1 (CARTE)
!     --------------------------------------------
    call alcart('V', chmar1, noma, 'VFER1_R')
    call jeveuo(chmar1//'.NCMP', 'E', vk8=ncmp)
    call jeveuo(chmar1//'.VALV', 'E', vr=valv)
!
    call jenonu(jexnom('&CATA.GD.NOMGD', 'VFER1_R'), gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    ASSERT(ncmpmx .eq. 20)

    ncmp(1) = 'TYPCOMB'
    ncmp(2) = 'CODIF'
    ncmp(3) = 'THITER'
    ncmp(4) = 'TYPSTRU'
    ncmp(5) = 'ENROBI'
    ncmp(6) = 'ENROBS'
    ncmp(7) = 'SIGS'
    ncmp(8) = 'SIGCI'
    ncmp(9) = 'SIGCS'
    ncmp(10) = 'ALPHACC'
    ncmp(11) = 'GAMMAS'
    ncmp(12) = 'GAMMAC'
    ncmp(13) = 'FACIER'
    ncmp(14) = 'EYS'
    ncmp(15) = 'TYPDIAG'
    ncmp(16) = 'FBETON'
    ncmp(17) = 'CLACIER'
    ncmp(18) = 'UC'
    ncmp(19) = 'UM'
    ncmp(20) = 'CEQUI'

!
!     2. MOTS CLES GLOBAUX :
!     ----------------------
!     2.1 TYPE_COMB :
    call getvtx(' ', 'TYPE_COMB', scal=typcb, nbret=n1)
    if (typcb .eq. 'ELU') valrcb = 0.d0
    if (typcb .eq. 'ELS') valrcb = 1.d0
    valv(1) = valrcb
!
!     2.2 CODIFICATION :
    call getvtx(' ', 'CODIFICATION', scal=typco, nbret=n2)
    if (typco .eq. 'BAEL91') valrco = 1.d0
    if (typco .eq. 'EC2') valrco = 2.d0
    valv(2) = valrco
!
!
!     2.3 CHOIX CRITERES DE PRECISION :
    call getvr8(' ', 'PAS_THETA', scal=valv(3), nbret=n3)
!
!     2.5 UNITES :
    call getvtx(' ', 'UNITE_CONTRAINTE', scal=unitc, nbret=n18)
!        UC = 0 CONTRAINTES EN Pa
!        UC = 1 CONTRAINTES EN MPa
    if (unitc .eq. 'Pa') valruc = 0.d0
    if (unitc .eq. 'MPa') valruc = 1.d0
    valv(18) = valruc
!   Pa avec m
!   MPa avec mm
    valv(19) = valruc

!     3- BOUCLE SUR LES OCCURENCES DU MOT CLE AFFE
!     --------------------------------------------
    do iocc = 1, nocc
!
        if (typco .eq. 'BAEL91') then
!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'BAEL91'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n4)
            if (typstru .eq. '2D') valv(4) = 0.d0
            call getvr8('AFFE', 'N', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(5), nbret=n5)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(6), nbret=n6)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(7), nbret=n7)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(8), nbret=n8)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(9), nbret=n9)
            call getvr8('AFFE', 'FE', iocc=iocc, scal=valv(13), nbret=n13)
            call getvr8('AFFE', 'FCJ', iocc=iocc, scal=valv(16), nbret=n16)
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(10), nbret=n10)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(11), nbret=n11)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(12), nbret=n12)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(14), nbret=n14)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n15)
            if (typdiag .eq. 'B1') valv(15) = 1.d0
            if (typdiag .eq. 'B2') valv(15) = 2.d0

        else if (typco .eq. 'EC2') then

!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'EC2'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n4)
            if (typstru .eq. '2D') valv(4) = 0.d0
            call getvr8('AFFE', 'ALPHA_E', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(5), nbret=n5)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(6), nbret=n6)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(7), nbret=n7)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(8), nbret=n8)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(9), nbret=n9)
            call getvr8('AFFE', 'FYK', iocc=iocc, scal=valv(13), nbret=n13)
            call getvr8('AFFE', 'FCK', iocc=iocc, scal=valv(16), nbret=n16)
            call getvtx('AFFE', 'CLASSE_ACIER', iocc=iocc, scal=clacier, nbret=n17)
            if (clacier .eq. 'A') valv(17) = 0.d0
            if (clacier .eq. 'B') valv(17) = 1.d0
            if (clacier .eq. 'C') valv(17) = 2.d0
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(10), nbret=n10)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(11), nbret=n11)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(12), nbret=n12)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(14), nbret=n14)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n15)
            if (typdiag .eq. 'B1') valv(15) = 1.d0
            if (typdiag .eq. 'B2') valv(15) = 2.d0
        end if

!
!       VERIFICATION DE LA COHERENCE DES PARAMETRES

!
        if (valv(1) .eq. 0.d0) then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELU
            if (n11 .eq. 0 .or. n12 .eq. 0 .or. n13 .eq. 0 &
                & .or. n14 .eq. 0 .or. n15 .eq. 0 .or. n16 .eq. 0) then
                call utmess('F', 'VERIFERRAILLAGE_8')
            end if

        else if (valv(1) .eq. 1.d0) then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELS
            if (n20 .eq. 0 .or. n7 .eq. 0 .or. n14 .eq. 0) then
                call utmess('F', 'VERIFERRAILLAGE_9')
            end if
            if (valv(4) .eq. 0.d0) then
                if (n7 .eq. 0 .or. n8 .eq. 0 .or. n9 .eq. 0) then
                    call utmess('F', 'VERIFERRAILLAGE_9')
                end if
            end if

        end if
!
        call getvtx('AFFE', 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(chmar1, 1, ncmpmx)
!
        else
            call reliem(' ', noma, 'NU_MAILLE', 'AFFE', iocc, &
                        2, motcls, typmcl, mesmai, nbmail)
            call jeveuo(mesmai, 'L', jmail)
            call nocart(chmar1, 3, ncmpmx, mode='NUM', nma=nbmail, &
                        limanu=zi(jmail))
            call jedetr(mesmai)
        end if
    end do
!
    call jedetr(chmar1//'.NCMP')
    call jedetr(chmar1//'.VALV')
!
    call jedema()
end subroutine
