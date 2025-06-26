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
    integer(kind=8) :: gd, nocc, ncmpmx, nbtou
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15
    integer(kind=8) :: n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32
   integer(kind=8) :: n33, n34, n35, n36, n37, n38, n39, n41, n42, n43, n44, n45, n46, n47, n48, n49
    integer(kind=8) :: n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61
    integer(kind=8) ::   jmail, iocc, nbmail
    real(kind=8) :: valrcb, valrco, valrmt, valrcd, valruc
    character(len=8) :: k8b, typmcl(2), noma, typcb, clacier, compress
    character(len=8) :: epucisa, ferrcomp, ferrsyme, typdiag, typstru, cond109, unitc
    character(len=16) :: meth2D
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
    ASSERT(ncmpmx .eq. 61)
    ncmp(1) = 'TYPCOMB'
    ncmp(2) = 'CODIF'
    ncmp(3) = 'METH2D'
    ncmp(4) = 'THITER'
    ncmp(5) = 'EPITER'
    ncmp(6) = 'APHITER'
    ncmp(7) = 'COND109'
    ncmp(8) = 'TYPSTRU'
    ncmp(9) = 'FERRSYME'
    ncmp(10) = 'SLSYME'
    ncmp(11) = 'FERRCOMP'
    ncmp(12) = 'EPUCISA'
    ncmp(13) = 'FERRMIN'
    ncmp(14) = 'RHOLMIN'
    ncmp(15) = 'RHOTMIN'
    ncmp(16) = 'COMPRESS'
    ncmp(17) = 'CEQUI'
    ncmp(18) = 'ENROBI'
    ncmp(19) = 'ENROBS'
    ncmp(20) = 'ENROBYI'
    ncmp(21) = 'ENROBYS'
    ncmp(22) = 'ENROBZI'
    ncmp(23) = 'ENROBZS'
    ncmp(24) = 'SIGS'
    ncmp(25) = 'SIGCI'
    ncmp(26) = 'SIGCS'
    ncmp(27) = 'SIGCYI'
    ncmp(28) = 'SIGCYS'
    ncmp(29) = 'SIGCZI'
    ncmp(30) = 'SIGCZS'
    ncmp(31) = 'ALPHACC'
    ncmp(32) = 'GAMMAS'
    ncmp(33) = 'GAMMAC'
    ncmp(34) = 'FACIER'
    ncmp(35) = 'EYS'
    ncmp(36) = 'TYPDIAG'
    ncmp(37) = 'FBETON'
    ncmp(38) = 'CLACIER'
    ncmp(39) = 'UC'
    ncmp(40) = 'UM'
    ncmp(41) = 'RHOACIER'
    ncmp(42) = 'AREINF'
    ncmp(43) = 'ASHEAR'
    ncmp(44) = 'ASTIRR'
    ncmp(45) = 'RHOCRIT'
    ncmp(46) = 'DATCRIT'
    ncmp(47) = 'LCRIT'
    ncmp(48) = 'WMAXI'
    ncmp(49) = 'WMAXS'
    ncmp(50) = 'WMAXYI'
    ncmp(51) = 'WMAXYS'
    ncmp(52) = 'WMAXZI'
    ncmp(53) = 'WMAXZS'
    ncmp(54) = 'SIGELSQP'
    ncmp(55) = 'KT'
    ncmp(56) = 'PHIXI'
    ncmp(57) = 'PHIXS'
    ncmp(58) = 'PHIYI'
    ncmp(59) = 'PHIYS'
    ncmp(60) = 'PHIZI'
    ncmp(61) = 'PHIZS'

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
!     2.3 CHOIX ALGO CALCUL 2D :
    call getvtx(' ', 'METHODE_2D', scal=meth2D, nbret=n3)
    if (meth2D .eq. 'CAPRA-MAURY') valrmt = 1.d0
    if (meth2D .eq. 'SANDWICH') valrmt = 2.d0
    valv(3) = valrmt
!
!     2.4 CHOIX CRITERES DE PRECISION :
    call getvr8(' ', 'PAS_THETA', scal=valv(4), nbret=n4)
    call getvr8(' ', 'PAS_EPAI', scal=valv(5), nbret=n5)
    call getvr8(' ', 'PAS_SIGM', scal=valv(6), nbret=n6)
    call getvtx(' ', 'COND_109', scal=cond109, nbret=n7)
!        COND109 = 0 (NON)
!        COND109 = 1 (OUI)
    if (cond109 .eq. 'OUI') valrcd = 1.d0
    if (cond109 .eq. 'NON') valrcd = 0.d0
    valv(7) = valrcd
!
!     2.5 UNITES :
    call getvtx(' ', 'UNITE_CONTRAINTE', scal=unitc, nbret=n39)
!        UC = 0 CONTRAINTES EN Pa
!        UC = 1 CONTRAINTES EN MPa
    if (unitc .eq. 'Pa') valruc = 0.d0
    if (unitc .eq. 'MPa') valruc = 1.d0
    valv(39) = valruc
!   Pa avec m
!   MPa avec mm
    valv(40) = valruc

!     3- BOUCLE SUR LES OCCURENCES DU MOT CLE AFFE
!     --------------------------------------------
    do iocc = 1, nocc
!
        if (typco .eq. 'BAEL91') then
!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'BAEL91'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n8)
            if (typstru .eq. '2D') valv(8) = 0.d0
            if (typstru .eq. '1D') valv(8) = 1.d0
            call getvtx('AFFE', 'FERR_SYME', iocc=iocc, scal=ferrsyme, nbret=n9)
            if (ferrsyme .eq. 'NON') valv(9) = 0.d0
            if (ferrsyme .eq. 'OUI') valv(9) = 1.d0
            call getvr8('AFFE', 'SEUIL_SYME', iocc=iocc, scal=valv(10), nbret=n10)
            call getvtx('AFFE', 'FERR_COMP', iocc=iocc, scal=ferrcomp, nbret=n11)
            if (ferrcomp .eq. 'NON') valv(11) = 0.d0
            if (ferrcomp .eq. 'OUI') valv(11) = 1.d0
            call getvtx('AFFE', 'EPURE_CISA', iocc=iocc, scal=epucisa, nbret=n12)
            if (epucisa .eq. 'NON') valv(12) = 0.d0
            if (epucisa .eq. 'OUI') valv(12) = 1.d0
            call getvtx('AFFE', 'FERR_MIN', iocc=iocc, scal=ferrmin, nbret=n13)
            if (ferrmin .eq. 'NON') valv(13) = 0.d0
            if (ferrmin .eq. 'OUI') valv(13) = 1.d0
            if (ferrmin .eq. 'CODE') valv(13) = 2.d0
            call getvr8('AFFE', 'RHO_LONGI_MIN', iocc=iocc, scal=valv(14), nbret=n14)
            call getvr8('AFFE', 'RHO_TRNSV_MIN', iocc=iocc, scal=valv(15), nbret=n15)
            call getvr8('AFFE', 'N', iocc=iocc, scal=valv(17), nbret=n17)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(18), nbret=n18)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(19), nbret=n19)
            call getvr8('AFFE', 'C_INF_Y', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'C_SUP_Y', iocc=iocc, scal=valv(21), nbret=n21)
            call getvr8('AFFE', 'C_INF_Z', iocc=iocc, scal=valv(22), nbret=n22)
            call getvr8('AFFE', 'C_SUP_Z', iocc=iocc, scal=valv(23), nbret=n23)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(24), nbret=n24)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(25), nbret=n25)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(26), nbret=n26)
            call getvr8('AFFE', 'SIGC_INF_Y_ELS', iocc=iocc, scal=valv(27), nbret=n27)
            call getvr8('AFFE', 'SIGC_SUP_Y_ELS', iocc=iocc, scal=valv(28), nbret=n28)
            call getvr8('AFFE', 'SIGC_INF_Z_ELS', iocc=iocc, scal=valv(29), nbret=n29)
            call getvr8('AFFE', 'SIGC_SUP_Z_ELS', iocc=iocc, scal=valv(30), nbret=n30)
            call getvr8('AFFE', 'FE', iocc=iocc, scal=valv(34), nbret=n34)
            call getvr8('AFFE', 'FCJ', iocc=iocc, scal=valv(37), nbret=n37)
            call getvr8('AFFE', 'WMAX_INF', iocc=iocc, scal=valv(48), nbret=n48)
            call getvr8('AFFE', 'WMAX_SUP', iocc=iocc, scal=valv(49), nbret=n49)
            call getvr8('AFFE', 'WMAX_INF_Y', iocc=iocc, scal=valv(50), nbret=n50)
            call getvr8('AFFE', 'WMAX_SUP_Y', iocc=iocc, scal=valv(51), nbret=n51)
            call getvr8('AFFE', 'WMAX_INF_Z', iocc=iocc, scal=valv(52), nbret=n52)
            call getvr8('AFFE', 'WMAX_SUP_Z', iocc=iocc, scal=valv(53), nbret=n53)
            call getvr8('AFFE', 'SIGC_ELS_QP', iocc=iocc, scal=valv(54), nbret=n54)
            call getvr8('AFFE', 'KT', iocc=iocc, scal=valv(55), nbret=n55)
            call getvr8('AFFE', 'PHI_INF_X', iocc=iocc, scal=valv(56), nbret=n56)
            call getvr8('AFFE', 'PHI_SUP_X', iocc=iocc, scal=valv(57), nbret=n57)
            call getvr8('AFFE', 'PHI_INF_Y', iocc=iocc, scal=valv(58), nbret=n58)
            call getvr8('AFFE', 'PHI_SUP_Y', iocc=iocc, scal=valv(59), nbret=n59)
            call getvr8('AFFE', 'PHI_INF_Z', iocc=iocc, scal=valv(60), nbret=n60)
            call getvr8('AFFE', 'PHI_SUP_Z', iocc=iocc, scal=valv(61), nbret=n61)
            call getvr8('AFFE', 'RHO_ACIER', iocc=iocc, scal=valv(41), nbret=n41)
            call getvr8('AFFE', 'ALPHA_REINF', iocc=iocc, scal=valv(42), nbret=n42)
            call getvr8('AFFE', 'ALPHA_SHEAR', iocc=iocc, scal=valv(43), nbret=n43)
            call getvr8('AFFE', 'ALPHA_STIRRUPS', iocc=iocc, scal=valv(44), nbret=n44)
            call getvr8('AFFE', 'RHO_CRIT', iocc=iocc, scal=valv(45), nbret=n45)
            call getvr8('AFFE', 'DNSTRA_CRIT', iocc=iocc, scal=valv(46), nbret=n46)
            call getvr8('AFFE', 'L_CRIT', iocc=iocc, scal=valv(47), nbret=n47)
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(31), nbret=n31)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(32), nbret=n32)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(33), nbret=n33)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(35), nbret=n35)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n36)
            if (typdiag .eq. 'B1') valv(36) = 1.d0
            if (typdiag .eq. 'B2') valv(36) = 2.d0
        else if (typco .eq. 'EC2') then
!           RECUPERATION DES MOTS CLES POUR CODIFICATION = 'EC2'
            call getvtx('AFFE', 'TYPE_STRUCTURE', iocc=iocc, scal=typstru, nbret=n8)
            if (typstru .eq. '2D') valv(8) = 0.d0
            if (typstru .eq. '1D') valv(8) = 1.d0
            call getvtx('AFFE', 'FERR_SYME', iocc=iocc, scal=ferrsyme, nbret=n9)
            if (ferrsyme .eq. 'NON') valv(9) = 0.d0
            if (ferrsyme .eq. 'OUI') valv(9) = 1.d0
            call getvr8('AFFE', 'SEUIL_SYME', iocc=iocc, scal=valv(10), nbret=n10)
            call getvtx('AFFE', 'FERR_COMP', iocc=iocc, scal=ferrcomp, nbret=n11)
            if (ferrcomp .eq. 'NON') valv(11) = 0.d0
            if (ferrcomp .eq. 'OUI') valv(11) = 1.d0
            call getvtx('AFFE', 'EPURE_CISA', iocc=iocc, scal=epucisa, nbret=n12)
            if (epucisa .eq. 'NON') valv(12) = 0.d0
            if (epucisa .eq. 'OUI') valv(12) = 1.d0
            call getvtx('AFFE', 'FERR_MIN', iocc=iocc, scal=ferrmin, nbret=n13)
            if (ferrmin .eq. 'NON') valv(13) = 0.d0
            if (ferrmin .eq. 'OUI') valv(13) = 1.d0
            if (ferrmin .eq. 'CODE') valv(13) = 2.d0
            call getvr8('AFFE', 'RHO_LONGI_MIN', iocc=iocc, scal=valv(14), nbret=n14)
            call getvr8('AFFE', 'RHO_TRNSV_MIN', iocc=iocc, scal=valv(15), nbret=n15)
            call getvtx('AFFE', 'UTIL_COMPR', iocc=iocc, scal=compress, nbret=n16)
            if (compress .eq. 'NON') valv(16) = 0.d0
            if (compress .eq. 'OUI') valv(16) = 1.d0
            call getvr8('AFFE', 'ALPHA_E', iocc=iocc, scal=valv(17), nbret=n17)
            call getvr8('AFFE', 'C_INF', iocc=iocc, scal=valv(18), nbret=n18)
            call getvr8('AFFE', 'C_SUP', iocc=iocc, scal=valv(19), nbret=n19)
            call getvr8('AFFE', 'C_INF_Y', iocc=iocc, scal=valv(20), nbret=n20)
            call getvr8('AFFE', 'C_SUP_Y', iocc=iocc, scal=valv(21), nbret=n21)
            call getvr8('AFFE', 'C_INF_Z', iocc=iocc, scal=valv(22), nbret=n22)
            call getvr8('AFFE', 'C_SUP_Z', iocc=iocc, scal=valv(23), nbret=n23)
            call getvr8('AFFE', 'SIGS_ELS', iocc=iocc, scal=valv(24), nbret=n24)
            call getvr8('AFFE', 'SIGC_INF_ELS', iocc=iocc, scal=valv(25), nbret=n25)
            call getvr8('AFFE', 'SIGC_SUP_ELS', iocc=iocc, scal=valv(26), nbret=n26)
            call getvr8('AFFE', 'SIGC_INF_Y_ELS', iocc=iocc, scal=valv(27), nbret=n27)
            call getvr8('AFFE', 'SIGC_SUP_Y_ELS', iocc=iocc, scal=valv(28), nbret=n28)
            call getvr8('AFFE', 'SIGC_INF_Z_ELS', iocc=iocc, scal=valv(29), nbret=n29)
            call getvr8('AFFE', 'SIGC_SUP_Z_ELS', iocc=iocc, scal=valv(30), nbret=n30)
            call getvr8('AFFE', 'FYK', iocc=iocc, scal=valv(34), nbret=n34)
            call getvr8('AFFE', 'FCK', iocc=iocc, scal=valv(37), nbret=n37)
            call getvtx('AFFE', 'CLASSE_ACIER', iocc=iocc, scal=clacier, nbret=n38)
            if (clacier .eq. 'A') valv(38) = 0.d0
            if (clacier .eq. 'B') valv(38) = 1.d0
            if (clacier .eq. 'C') valv(38) = 2.d0
            call getvr8('AFFE', 'WMAX_INF', iocc=iocc, scal=valv(48), nbret=n48)
            call getvr8('AFFE', 'WMAX_SUP', iocc=iocc, scal=valv(49), nbret=n49)
            call getvr8('AFFE', 'WMAX_INF_Y', iocc=iocc, scal=valv(50), nbret=n50)
            call getvr8('AFFE', 'WMAX_SUP_Y', iocc=iocc, scal=valv(51), nbret=n51)
            call getvr8('AFFE', 'WMAX_INF_Z', iocc=iocc, scal=valv(52), nbret=n52)
            call getvr8('AFFE', 'WMAX_SUP_Z', iocc=iocc, scal=valv(53), nbret=n53)
            call getvr8('AFFE', 'SIGC_ELS_QP', iocc=iocc, scal=valv(54), nbret=n54)
            call getvr8('AFFE', 'KT', iocc=iocc, scal=valv(55), nbret=n55)
            call getvr8('AFFE', 'PHI_INF_X', iocc=iocc, scal=valv(56), nbret=n56)
            call getvr8('AFFE', 'PHI_SUP_X', iocc=iocc, scal=valv(57), nbret=n57)
            call getvr8('AFFE', 'PHI_INF_Y', iocc=iocc, scal=valv(58), nbret=n58)
            call getvr8('AFFE', 'PHI_SUP_Y', iocc=iocc, scal=valv(59), nbret=n59)
            call getvr8('AFFE', 'PHI_INF_Z', iocc=iocc, scal=valv(60), nbret=n60)
            call getvr8('AFFE', 'PHI_SUP_Z', iocc=iocc, scal=valv(61), nbret=n61)
            call getvr8('AFFE', 'RHO_ACIER', iocc=iocc, scal=valv(41), nbret=n41)
            call getvr8('AFFE', 'ALPHA_REINF', iocc=iocc, scal=valv(42), nbret=n42)
            call getvr8('AFFE', 'ALPHA_SHEAR', iocc=iocc, scal=valv(43), nbret=n43)
            call getvr8('AFFE', 'ALPHA_STIRRUPS', iocc=iocc, scal=valv(44), nbret=n44)
            call getvr8('AFFE', 'RHO_CRIT', iocc=iocc, scal=valv(45), nbret=n45)
            call getvr8('AFFE', 'DNSTRA_CRIT', iocc=iocc, scal=valv(46), nbret=n46)
            call getvr8('AFFE', 'L_CRIT', iocc=iocc, scal=valv(47), nbret=n47)
            call getvr8('AFFE', 'ALPHA_CC', iocc=iocc, scal=valv(31), nbret=n31)
            call getvr8('AFFE', 'GAMMA_S', iocc=iocc, scal=valv(32), nbret=n32)
            call getvr8('AFFE', 'GAMMA_C', iocc=iocc, scal=valv(33), nbret=n33)
            call getvr8('AFFE', 'EYS', iocc=iocc, scal=valv(35), nbret=n35)
            call getvtx('AFFE', 'TYPE_DIAGRAMME', iocc=iocc, scal=typdiag, nbret=n36)
            if (typdiag .eq. 'B1') valv(36) = 1.d0
            if (typdiag .eq. 'B2') valv(36) = 2.d0
        end if

!
!       VERIFICATION DE LA COHERENCE DES PARAMETRES

        if (typstru .eq. '2D') then
            if ((valv(4) .le. 0) .or. (valv(4) .gt. 10)) then
                call utmess('F', 'CALCULEL7_31')
            end if
            if ((valv(3) .eq. 2.d0) .and. (valv(1) .ne. 0.d0)) then
                call utmess('F', 'CALCULEL7_37')
            end if
            if ((valv(3) .eq. 2.d0) .and. (valv(1) .eq. 0.d0)) then
                if ((valv(6) .le. 0) .or. (valv(6) .gt. 0.2)) then
                    call utmess('F', 'CALCULEL7_33')
                end if
            end if
        end if

        if (typstru .eq. '1D') then
!           VERIFICATION DES ENROBAGES 1D
            if (n20 .eq. 0 .or. n21 .eq. 0 .or. n22 .eq. 0 &
                & .or. n23 .eq. 0) then
                call utmess('F', 'CALCULEL7_18')
            end if
        elseif (typstru .eq. '2D') then
!           VERIFICATION DES ENROBAGES 2D
            if (n18 .eq. 0 .or. n19 .eq. 0) then
                call utmess('F', 'CALCULEL7_19')
            end if
        end if

        if (ferrsyme .eq. 'OUI') then
!           VERIFICATION DES SYMETRIES
            if (n10 .eq. 0) then
                call utmess('F', 'CALCULEL7_30')
            end if
            if (typstru .eq. '1D') then
                if (valv(20) .ne. valv(21) &
                    & .or. valv(22) .ne. valv(23)) then
                    call utmess('F', 'CALCULEL7_20')
                end if
                if (typcb .eq. 'ELS') then
                    if (valv(27) .ne. valv(28) &
                        & .or. valv(29) .ne. valv(30)) then
                        call utmess('F', 'CALCULEL7_21')
                    end if
                elseif (typcb .eq. 'ELS_QP') then
                    if (valv(50) .ne. valv(51) &
                        & .or. valv(52) .ne. valv(53)) then
                        call utmess('F', 'CALCULEL7_22')
                    end if
                    if (valv(58) .ne. valv(59) &
                        & .or. valv(60) .ne. valv(61)) then
                        call utmess('F', 'CALCULEL7_23')
                    end if
                end if
            elseif (typstru .eq. '2D') then
                if (valv(18) .ne. valv(19)) then
                    call utmess('F', 'CALCULEL7_24')
                end if
                if (typcb .eq. 'ELS') then
                    if (valv(25) .ne. valv(26)) then
                        call utmess('F', 'CALCULEL7_25')
                    end if
                elseif (typcb .eq. 'ELS_QP') then
                    if (valv(48) .ne. valv(49)) then
                        call utmess('F', 'CALCULEL7_26')
                    end if
                    if (valv(56) .ne. valv(57) &
                        & .or. valv(58) .ne. valv(59)) then
                        call utmess('F', 'CALCULEL7_27')
                    end if
                end if
            end if
        end if

        if (valv(13) .eq. (1.d0)) then
            if ((n14 .eq. 0) .or. (n15 .eq. 0)) then
                call utmess('F', 'CALCULEL7_11')
            end if
        end if

        if (valv(41) .lt. 0.d0) then
!           MASSE VOLUMIQUE NEGATIVE
            call utmess('I', 'CALCULEL_89')
        end if
!
        if (typcb .eq. 'ELU') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELU
            if (n32 .eq. 0 .or. n33 .eq. 0 .or. n34 .eq. 0 &
                & .or. n35 .eq. 0 .or. n36 .eq. 0 .or. n37 .eq. 0) then
                call utmess('F', 'CALCULEL_74')
            end if

        else if (typcb .eq. 'ELS') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELS
            if (n17 .eq. 0 .or. n24 .eq. 0) then
                call utmess('F', 'CALCULEL_82')
            end if
            if (typstru .eq. '2D') then
                if (n24 .eq. 0 .or. n26 .eq. 0) then
                    call utmess('F', 'CALCULEL_82')
                end if
            elseif (typstru .eq. '1D') then
                if (n27 .eq. 0 .or. n28 .eq. 0 .or. n29 .eq. 0 .or. n30 .eq. 0) then
                    call utmess('F', 'CALCULEL_82')
                end if
            end if
            if (typco .eq. 'EC2' .and. n37 .eq. 0) then
                call utmess('F', 'CALCULEL_82')
            end if
            if (typco .eq. 'BAEL91') then
!               MESSAGE D'INFORMATION : PAS DE CALCUL DES ACIERS
!               TRANSVERSAUX POUR LA CODIFICATION BAEL91
                call utmess('I', 'CALCULEL_80')
            end if

        else if (typcb .eq. 'ELS_QP') then
!           MOTS-CLE OBLIGATOIRES POUR UN CALCUL A L'ELS QP
            if (n17 .eq. 0 .or. n34 .eq. 0 .or. n35 .eq. 0 .or. n37 .eq. 0 &
                & .or. n54 .eq. 0 .or. n55 .eq. 0) then
                call utmess('F', 'CALCULEL7_8')
            end if
            if (typstru .eq. '2D') then
                if (n48 .eq. 0 .or. n49 .eq. 0 .or. n56 .eq. 0 &
                    & .or. n57 .eq. 0 .or. n58 .eq. 0 .or. n59 .eq. 0) then
                    call utmess('F', 'CALCULEL7_8')
                end if
            elseif (typstru .eq. '1D') then
                if (n50 .eq. 0 .or. n51 .eq. 0 .or. n52 .eq. 0 &
                    & .or. n53 .eq. 0 .or. n58 .eq. 0 .or. n59 .eq. 0 &
                    & .or. n60 .eq. 0 .or. n61 .eq. 0) then
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
