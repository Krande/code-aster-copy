! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine w175af(modele, chfer1)
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
    character(len=19) :: chfer1
!
! BUT : CREER LE CHAMP DE DONNEES POUR CALC_FERRAILLAGE
!
!-------------------------------------------------------------------------------------------------
    integer :: gd, nocc, ncmpmx, nbtou
    integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15
    integer :: n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32
    integer :: n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49
    integer :: n50, n51, n52, n53, n54, n55, n56
    integer ::   jmail, iocc, nbmail
    real(kind=8) :: valrcb, valrco
    character(len=8) :: k8b, typmcl(2), noma, typcb, clacier, uc, compress
    character(len=8) :: epucisa, ferrcomp, ferrsyme, typdiag, typstru
    character(len=16) :: motcls(2), typco, ferrmin
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
    mesmai = '&&W175AF.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
!     1- ALLOCATION DU CHAMP CHFER1 (CARTE)
!     --------------------------------------------
    call alcart('V', chfer1, noma, 'FER1_R')
    call jeveuo(chfer1//'.NCMP', 'E', vk8=ncmp)
    call jeveuo(chfer1//'.VALV', 'E', vr=valv)
!
    call jenonu(jexnom('&CATA.GD.NOMGD', 'FER1_R'), gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
!
    ASSERT(ncmpmx .eq. 56)
    ncmp(1) = 'TYPCOMB'
    ncmp(2) = 'CODIF'
    ncmp(3) = 'TYPSTRU'
    ncmp(4) = 'FERRSYME'
    ncmp(5) = 'SLSYME'
    ncmp(6) = 'FERRCOMP'
    ncmp(7) = 'EPUCISA'
    ncmp(8) = 'FERRMIN'
    ncmp(9) = 'RHOLMIN'
    ncmp(10) = 'RHOTMIN'
    ncmp(11) = 'COMPRESS'
    ncmp(12) = 'CEQUI'
    ncmp(13) = 'ENROBI'
    ncmp(14) = 'ENROBS'
    ncmp(15) = 'ENROBYI'
    ncmp(16) = 'ENROBYS'
    ncmp(17) = 'ENROBZI'
    ncmp(18) = 'ENROBZS'
    ncmp(19) = 'SIGS'
    ncmp(20) = 'SIGCI'
    ncmp(21) = 'SIGCS'
    ncmp(22) = 'SIGCYI'
    ncmp(23) = 'SIGCYS'
    ncmp(24) = 'SIGCZI'
    ncmp(25) = 'SIGCZS'
    ncmp(26) = 'ALPHACC'
    ncmp(27) = 'GAMMAS'
    ncmp(28) = 'GAMMAC'
    ncmp(29) = 'FACIER'
    ncmp(30) = 'EYS'
    ncmp(31) = 'TYPDIAG'
    ncmp(32) = 'FBETON'
    ncmp(33) = 'CLACIER'
    ncmp(34) = 'UC'
    ncmp(35) = 'UM'
    ncmp(36) = 'RHOACIER'
    ncmp(37) = 'AREINF'
    ncmp(38) = 'ASHEAR'
    ncmp(39) = 'ASTIRR'
    ncmp(40) = 'RHOCRIT'
    ncmp(41) = 'DATCRIT'
    ncmp(42) = 'LCRIT'
    ncmp(43) = 'WMAXI'
    ncmp(44) = 'WMAXS'
    ncmp(45) = 'WMAXYI'
    ncmp(46) = 'WMAXYS'
    ncmp(47) = 'WMAXZI'
    ncmp(48) = 'WMAXZS'
    ncmp(49) = 'SIGELSQP'
    ncmp(50) = 'KT'
    ncmp(51) = 'PHIXI'
    ncmp(52) = 'PHIXS'
    ncmp(53) = 'PHIYI'
    ncmp(54) = 'PHIYS'
    ncmp(55) = 'PHIZI'
    ncmp(56) = 'PHIZS'
!
!     2. MOTS CLES GLOBAUX :
!     ----------------------
!     2.1 TYPE_COMB :
    call getvtx(' ', 'TYPE_COMB', scal=typcb, nbret=n1)
    if (typcb .eq. 'ELU') valrcb = 0.d0
    if (typcb .eq. 'ELS') valrcb = 1.d0
    if (typcb .eq. 'ELS_QP') valrcb = 2.d0
    valv(1) = valrcb
!
!     2.2 CODIFICATION :
    call getvtx(' ', 'CODIFICATION', scal=typco, nbret=n2)
    if (typco .eq. 'BAEL91') valrco = 1.d0
    if (typco .eq. 'EC2') valrco = 2.d0
    valv(2) = valrco
!
!
!     3- BOUCLE SUR LES OCCURENCES DU MOT CLE AFFE
!     --------------------------------------------
    do iocc = 1, nocc
!
        if (typco .eq. 'BAEL91') then
!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'BAEL91'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n3)
            if (typstru .eq. '2D') valv(3) = 0.d0
            if (typstru .eq. '1D') valv(3) = 1.d0
            call getvtx('AFFE', 'FERR_SYME', iocc=iocc, scal=ferrsyme, nbret=n4)
            if (ferrsyme .eq. 'NON') valv(4) = 0.d0
            if (ferrsyme .eq. 'OUI') valv(4) = 1.d0
            call getvr8('AFFE', 'SEUIL_SYME', iocc=iocc, scal=valv(5), nbret=n5)
            call getvtx('AFFE', 'FERR_COMP', iocc=iocc, scal=ferrcomp, nbret=n6)
            if (ferrcomp .eq. 'NON') valv(6) = 0.d0
            if (ferrcomp .eq. 'OUI') valv(6) = 1.d0
            call getvtx('AFFE', 'EPURE_CISA', iocc=iocc, scal=epucisa, nbret=n7)
            if (epucisa .eq. 'NON') valv(7) = 0.d0
            if (epucisa .eq. 'OUI') valv(7) = 1.d0
            call getvtx('AFFE', 'FERR_MIN', iocc=iocc, scal=ferrmin, nbret=n8)
            if (ferrmin .eq. 'NON') valv(8) = 0.d0
            if (ferrmin .eq. 'OUI') valv(8) = 1.d0
            if (ferrmin .eq. 'CODE') valv(8) = 2.d0
            call getvr8('AFFE', 'RHO_LONGI_MIN', iocc=iocc, scal=valv(9), nbret=n9)
            call getvr8('AFFE', 'RHO_TRNSV_MIN', iocc=iocc, scal=valv(10), nbret=n10)
            call getvr8('AFFE', 'N', iocc=iocc, scal=valv(12), nbret=n12)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(13), nbret=n13)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(14), nbret=n14)
            call getvr8('AFFE', 'C_INF_Y', iocc=iocc, scal=valv(15), nbret=n15)
            call getvr8('AFFE', 'C_SUP_Y', iocc=iocc, scal=valv(16), nbret=n16)
            call getvr8('AFFE', 'C_INF_Z', iocc=iocc, scal=valv(17), nbret=n17)
            call getvr8('AFFE', 'C_SUP_Z', iocc=iocc, scal=valv(18), nbret=n18)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(19), nbret=n19)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(21), nbret=n21)
            call getvr8('AFFE', 'SIGC_INF_Y_ELS', iocc=iocc, scal=valv(22), nbret=n22)
            call getvr8('AFFE', 'SIGC_SUP_Y_ELS', iocc=iocc, scal=valv(23), nbret=n23)
            call getvr8('AFFE', 'SIGC_INF_Z_ELS', iocc=iocc, scal=valv(24), nbret=n24)
            call getvr8('AFFE', 'SIGC_SUP_Z_ELS', iocc=iocc, scal=valv(25), nbret=n25)
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(26), nbret=n26)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(27), nbret=n27)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(28), nbret=n28)
            call getvr8('AFFE', 'FE', iocc=iocc, scal=valv(29), nbret=n29)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(30), nbret=n30)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n31)
            if (uc .eq. 'B1') valv(31) = 1.d0
            if (uc .eq. 'B2') valv(31) = 2.d0
            call getvr8('AFFE', 'FCJ', iocc=iocc, scal=valv(32), nbret=n32)
            call getvtx('AFFE', 'UNITE_CONTRAINTE', iocc=iocc, scal=uc, nbret=n34)
            if (uc .eq. 'Pa') valv(34) = 0.d0
            if (uc .eq. 'MPa') valv(34) = 1.d0
            call getvtx('AFFE', 'UNITE_DIMENSION', iocc=iocc, scal=uc, nbret=n35)
            if (uc .eq. 'm') valv(35) = 0.d0
            if (uc .eq. 'mm') valv(35) = 1.d0
            call getvr8('AFFE', 'RHO_ACIER', iocc=iocc, scal=valv(36), nbret=n36)
            call getvr8('AFFE', 'ALPHA_REINF', iocc=iocc, scal=valv(37), nbret=n37)
            call getvr8('AFFE', 'ALPHA_SHEAR', iocc=iocc, scal=valv(38), nbret=n38)
            call getvr8('AFFE', 'ALPHA_STIRRUPS', iocc=iocc, scal=valv(39), nbret=n39)
            call getvr8('AFFE', 'RHO_CRIT', iocc=iocc, scal=valv(40), nbret=n40)
            call getvr8('AFFE', 'DNSTRA_CRIT', iocc=iocc, scal=valv(41), nbret=n41)
            call getvr8('AFFE', 'L_CRIT', iocc=iocc, scal=valv(42), nbret=n42)
            call getvr8('AFFE', 'WMAX_INF', iocc=iocc, scal=valv(43), nbret=n43)
            call getvr8('AFFE', 'WMAX_SUP', iocc=iocc, scal=valv(44), nbret=n44)
            call getvr8('AFFE', 'WMAX_INF_Y', iocc=iocc, scal=valv(45), nbret=n45)
            call getvr8('AFFE', 'WMAX_SUP_Y', iocc=iocc, scal=valv(46), nbret=n46)
            call getvr8('AFFE', 'WMAX_INF_Z', iocc=iocc, scal=valv(47), nbret=n47)
            call getvr8('AFFE', 'WMAX_SUP_Z', iocc=iocc, scal=valv(48), nbret=n48)
            call getvr8('AFFE', 'SIGC_ELS_QP', iocc=iocc, scal=valv(49), nbret=n49)
            call getvr8('AFFE', 'KT', iocc=iocc, scal=valv(50), nbret=n50)
            call getvr8('AFFE', 'PHI_INF_X', iocc=iocc, scal=valv(51), nbret=n51)
            call getvr8('AFFE', 'PHI_SUP_X', iocc=iocc, scal=valv(52), nbret=n52)
            call getvr8('AFFE', 'PHI_INF_Y', iocc=iocc, scal=valv(53), nbret=n53)
            call getvr8('AFFE', 'PHI_SUP_Y', iocc=iocc, scal=valv(54), nbret=n54)
            call getvr8('AFFE', 'PHI_INF_Z', iocc=iocc, scal=valv(55), nbret=n55)
            call getvr8('AFFE', 'PHI_SUP_Z', iocc=iocc, scal=valv(56), nbret=n56)
        else if (typco .eq. 'EC2') then
!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'EC2'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n3)
            if (typstru .eq. '2D') valv(3) = 0.d0
            if (typstru .eq. '1D') valv(3) = 1.d0
            call getvtx('AFFE', 'FERR_SYME', iocc=iocc, scal=ferrsyme, nbret=n4)
            if (ferrsyme .eq. 'NON') valv(4) = 0.d0
            if (ferrsyme .eq. 'OUI') valv(4) = 1.d0
            call getvr8('AFFE', 'SEUIL_SYME', iocc=iocc, scal=valv(5), nbret=n5)
            call getvtx('AFFE', 'FERR_COMP', iocc=iocc, scal=ferrcomp, nbret=n6)
            if (ferrcomp .eq. 'NON') valv(6) = 0.d0
            if (ferrcomp .eq. 'OUI') valv(6) = 1.d0
            call getvtx('AFFE', 'EPURE_CISA', iocc=iocc, scal=epucisa, nbret=n7)
            if (epucisa .eq. 'NON') valv(7) = 0.d0
            if (epucisa .eq. 'OUI') valv(7) = 1.d0
            call getvtx('AFFE', 'FERR_MIN', iocc=iocc, scal=ferrmin, nbret=n8)
            if (ferrmin .eq. 'NON') valv(8) = 0.d0
            if (ferrmin .eq. 'OUI') valv(8) = 1.d0
            if (ferrmin .eq. 'CODE') valv(8) = 2.d0
            call getvr8('AFFE', 'RHO_LONGI_MIN', iocc=iocc, scal=valv(9), nbret=n9)
            call getvr8('AFFE', 'RHO_TRNSV_MIN', iocc=iocc, scal=valv(10), nbret=n10)
            call getvtx('AFFE', 'UTIL_COMPR', iocc=iocc, scal=compress, nbret=n11)
            if (compress .eq. 'NON') valv(11) = 0.d0
            if (compress .eq. 'OUI') valv(11) = 1.d0
            call getvr8('AFFE', 'ALPHA_E', iocc=iocc, scal=valv(12), nbret=n12)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(13), nbret=n13)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(14), nbret=n14)
            call getvr8('AFFE', 'C_INF_Y', iocc=iocc, scal=valv(15), nbret=n15)
            call getvr8('AFFE', 'C_SUP_Y', iocc=iocc, scal=valv(16), nbret=n16)
            call getvr8('AFFE', 'C_INF_Z', iocc=iocc, scal=valv(17), nbret=n17)
            call getvr8('AFFE', 'C_SUP_Z', iocc=iocc, scal=valv(18), nbret=n18)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(19), nbret=n19)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(21), nbret=n21)
            call getvr8('AFFE', 'SIGC_INF_Y_ELS', iocc=iocc, scal=valv(22), nbret=n22)
            call getvr8('AFFE', 'SIGC_SUP_Y_ELS', iocc=iocc, scal=valv(23), nbret=n23)
            call getvr8('AFFE', 'SIGC_INF_Z_ELS', iocc=iocc, scal=valv(24), nbret=n24)
            call getvr8('AFFE', 'SIGC_SUP_Z_ELS', iocc=iocc, scal=valv(25), nbret=n25)
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(26), nbret=n26)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(27), nbret=n27)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(28), nbret=n28)
            call getvr8('AFFE', 'FYK', iocc=iocc, scal=valv(29), nbret=n29)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(30), nbret=n30)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n31)
            if (typdiag .eq. 'B1') valv(31) = 1.d0
            if (typdiag .eq. 'B2') valv(31) = 2.d0
            call getvr8('AFFE', 'FCK', iocc=iocc, scal=valv(32), nbret=n32)
            call getvtx('AFFE', 'CLASSE_ACIER', iocc=iocc, scal=clacier, nbret=n33)
            if (clacier .eq. 'A') valv(33) = 0.d0
            if (clacier .eq. 'B') valv(33) = 1.d0
            if (clacier .eq. 'C') valv(33) = 2.d0
            call getvtx('AFFE', 'UNITE_CONTRAINTE', iocc=iocc, scal=uc, nbret=n34)
            if (uc .eq. 'Pa') valv(34) = 0.d0
            if (uc .eq. 'MPa') valv(34) = 1.d0
            call getvtx('AFFE', 'UNITE_DIMENSION', iocc=iocc, scal=uc, nbret=n35)
            if (uc .eq. 'm') valv(35) = 0.d0
            if (uc .eq. 'mm') valv(35) = 1.d0
            call getvr8('AFFE', 'RHO_ACIER', iocc=iocc, scal=valv(36), nbret=n36)
            call getvr8('AFFE', 'ALPHA_REINF', iocc=iocc, scal=valv(37), nbret=n37)
            call getvr8('AFFE', 'ALPHA_SHEAR', iocc=iocc, scal=valv(38), nbret=n38)
            call getvr8('AFFE', 'ALPHA_STIRRUPS', iocc=iocc, scal=valv(39), nbret=n39)
            call getvr8('AFFE', 'RHO_CRIT', iocc=iocc, scal=valv(40), nbret=n40)
            call getvr8('AFFE', 'DNSTRA_CRIT', iocc=iocc, scal=valv(41), nbret=n41)
            call getvr8('AFFE', 'L_CRIT', iocc=iocc, scal=valv(42), nbret=n42)
            call getvr8('AFFE', 'WMAX_INF', iocc=iocc, scal=valv(43), nbret=n43)
            call getvr8('AFFE', 'WMAX_SUP', iocc=iocc, scal=valv(44), nbret=n44)
            call getvr8('AFFE', 'WMAX_INF_Y', iocc=iocc, scal=valv(45), nbret=n45)
            call getvr8('AFFE', 'WMAX_SUP_Y', iocc=iocc, scal=valv(46), nbret=n46)
            call getvr8('AFFE', 'WMAX_INF_Z', iocc=iocc, scal=valv(47), nbret=n47)
            call getvr8('AFFE', 'WMAX_SUP_Z', iocc=iocc, scal=valv(48), nbret=n48)
            call getvr8('AFFE', 'SIGC_ELS_QP', iocc=iocc, scal=valv(49), nbret=n49)
            call getvr8('AFFE', 'KT', iocc=iocc, scal=valv(50), nbret=n50)
            call getvr8('AFFE', 'PHI_INF_X', iocc=iocc, scal=valv(51), nbret=n51)
            call getvr8('AFFE', 'PHI_SUP_X', iocc=iocc, scal=valv(52), nbret=n52)
            call getvr8('AFFE', 'PHI_INF_Y', iocc=iocc, scal=valv(53), nbret=n53)
            call getvr8('AFFE', 'PHI_SUP_Y', iocc=iocc, scal=valv(54), nbret=n54)
            call getvr8('AFFE', 'PHI_INF_Z', iocc=iocc, scal=valv(55), nbret=n55)
            call getvr8('AFFE', 'PHI_SUP_Z', iocc=iocc, scal=valv(56), nbret=n56)
        end if

!
!       VERIFICATION DE LA COHERENCE DES PARAMETRES

        if (typstru .eq. '1D') then
!           VERIFICATION DES ENROBAGES 1D
            if (n15 .eq. 0 .or. n16 .eq. 0 .or. n17 .eq. 0 &
                & .or. n18 .eq. 0) then
                call utmess('F', 'CALCULEL7_18')
            end if
        elseif (typstru .eq. '2D') then
!           VERIFICATION DES ENROBAGES 2D
            if (n13 .eq. 0 .or. n14 .eq. 0) then
                call utmess('F', 'CALCULEL7_19')
            end if
        end if

        if (ferrsyme .eq. 'OUI') then
!           VERIFICATION DES SYMETRIES
            if (n5 .eq. 0) then
                call utmess('F', 'CALCULEL7_30')
            end if
            if (typstru .eq. '1D') then
                if (valv(15) .ne. valv(16) &
                    & .or. valv(17) .ne. valv(18)) then
                    call utmess('F', 'CALCULEL7_20')
                end if
                if (typcb .eq. 'ELS') then
                    if (valv(22) .ne. valv(23) &
                        & .or. valv(24) .ne. valv(25)) then
                        call utmess('F', 'CALCULEL7_21')
                    end if
                elseif (typcb .eq. 'ELS_QP') then
                    if (valv(45) .ne. valv(46) &
                        & .or. valv(47) .ne. valv(48)) then
                        call utmess('F', 'CALCULEL7_22')
                    end if
                    if (valv(53) .ne. valv(54) &
                        & .or. valv(55) .ne. valv(56)) then
                        call utmess('F', 'CALCULEL7_23')
                    end if
                end if
            elseif (typstru .eq. '2D') then
                if (valv(13) .ne. valv(14)) then
                    call utmess('F', 'CALCULEL7_24')
                end if
                if (typcb .eq. 'ELS') then
                    if (valv(20) .ne. valv(21)) then
                        call utmess('F', 'CALCULEL7_25')
                    end if
                elseif (typcb .eq. 'ELS_QP') then
                    if (valv(43) .ne. valv(44)) then
                        call utmess('F', 'CALCULEL7_26')
                    end if
                    if (valv(51) .ne. valv(52) &
                        & .or. valv(53) .ne. valv(54)) then
                        call utmess('F', 'CALCULEL7_27')
                    end if
                end if
            end if
        end if

        if (valv(8) .eq. (1.d0)) then
            if ((n9 .eq. 0) .or. (n10 .eq. 0)) then
                call utmess('F', 'CALCULEL7_11')
            end if
        end if

        if (valv(36) .lt. 0.d0) then
!           MASSE VOLUMIQUE NEGATIVE
            call utmess('I', 'CALCULEL_89')
        end if
!
        if (typcb .eq. 'ELU') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELU
            if (n27 .eq. 0 .or. n28 .eq. 0 .or. n29 .eq. 0 &
                & .or. n30 .eq. 0 .or. n31 .eq. 0 .or. n32 .eq. 0) then
                call utmess('F', 'CALCULEL_74')
            end if

        else if (typcb .eq. 'ELS') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELS
            if (n12 .eq. 0 .or. n19 .eq. 0) then
                call utmess('F', 'CALCULEL_82')
            end if
            if (typstru .eq. '2D') then
                if (n20 .eq. 0 .or. n21 .eq. 0) then
                    call utmess('F', 'CALCULEL_82')
                end if
            elseif (typstru .eq. '1D') then
                if (n22 .eq. 0 .or. n23 .eq. 0 .or. n24 .eq. 0 .or. n25 .eq. 0) then
                    call utmess('F', 'CALCULEL_82')
                end if
            end if
            if (typco .eq. 'EC2' .and. n32 .eq. 0) then
                call utmess('F', 'CALCULEL_82')
            end if
            if (typco .eq. 'BAEL91') then
!               MESSAGE D'INFORMATION : PAS DE CALCUL DES ACIERS
!               TRANSVERSAUX POUR LA CODIFICATION BAEL91
                call utmess('I', 'CALCULEL_80')
            end if

        else if (typcb .eq. 'ELS_QP') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELS QP
            if (n12 .eq. 0 .or. n29 .eq. 0 .or. n30 .eq. 0 .or. n32 .eq. 0 &
                & .or. n49 .eq. 0 .or. n50 .eq. 0) then
                call utmess('F', 'CALCULEL7_8')
            end if
            if (typstru .eq. '2D') then
                if (n43 .eq. 0 .or. n44 .eq. 0 .or. n51 .eq. 0 &
                    & .or. n52 .eq. 0 .or. n53 .eq. 0 .or. n54 .eq. 0) then
                    call utmess('F', 'CALCULEL7_8')
                end if
            elseif (typstru .eq. '1D') then
                if (n45 .eq. 0 .or. n46 .eq. 0 .or. n47 .eq. 0 &
                    & .or. n48 .eq. 0 .or. n53 .eq. 0 .or. n54 .eq. 0 &
                    & .or. n55 .eq. 0 .or. n56 .eq. 0) then
                    call utmess('F', 'CALCULEL7_8')
                end if
            end if
            if (typco .eq. 'BAEL91') then
!               MESSAGE D'INFORMATION : CALCUL RÉALISÉ SUR LA
!               THEORIE DE EC2
                call utmess('I', 'CALCULEL7_9')
            end if
        end if
!
        call getvtx('AFFE', 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(chfer1, 1, ncmpmx)
!
        else
            call reliem(' ', noma, 'NU_MAILLE', 'AFFE', iocc, &
                        2, motcls, typmcl, mesmai, nbmail)
            call jeveuo(mesmai, 'L', jmail)
            call nocart(chfer1, 3, ncmpmx, mode='NUM', nma=nbmail, &
                        limanu=zi(jmail))
            call jedetr(mesmai)
        end if
    end do
!
    call jedetr(chfer1//'.NCMP')
    call jedetr(chfer1//'.VALV')
!
    call jedema()
end subroutine
