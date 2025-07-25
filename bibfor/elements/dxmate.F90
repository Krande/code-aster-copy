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

subroutine dxmate(fami, df, dm, dmf, dc, &
                  dci, dmc, dfc, nno, pgl, &
                  multic, coupmf, t2iu, t2ui, t1ve)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/codent.h"
#include "asterfort/coqrep.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moyte2.h"
#include "asterfort/moytem.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcadlv.h"
#include "asterfort/rcvalt.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utbtab.h"
#include "asterfort/utdtab.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
    aster_logical :: coupmf
    integer(kind=8) :: nno, multic
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: pgl(3, 3), t2iu(4), t2ui(4), t1ve(9)
    character(len=4) :: fami
!     CALCUL DES MATRICES DE RIGIDITE DE FLEXION, MEMBRANE , COUPLAGE
!     MEMBRANE-FLEXION ET CISAILLEMENT POUR UN MATERIAU ISOTROPE OU
!     MULTICOUCHE
!     OUT MULTIC :
!        1 POUR UN MATERIAU MULTICOUCHE SANS COUPLAGE MEMBRANE-FLEXION
!        0 DANS LES AUTRES CAS
!     OUT COUPMF :
!        .TRUE. POUR UN MATERIAU AVEC COUPLAGE MEMBRANE-FLEXION
!  ------------------------------------------------------------------
! aslint: disable=W1501
    integer(kind=8) :: jcoqu, jmate, nbv, i, j, k, nbpar, elasco, c_elasco(45)
    integer(kind=8) :: iazi, iazk24, npg, jcou, ncou, iret, npgh, iret1
    integer(kind=8) :: ndim, nnos, ipoids, ivf, idfde, jgano
    integer(kind=8) :: icodre(33), icodre1, jadr, n1
    real(kind=8) :: kcis, cdf, cdm, cdc, gcis, valres(33), v_elasco(45)
    real(kind=8) :: young, nu, epais, temp, excent
    real(kind=8) :: xab1(3, 3), xab2(2, 2), xab3(3, 2)
    real(kind=8) :: s, c
    real(kind=8) :: alpha, beta, det
    real(kind=8) :: zero, deux
    real(kind=8) :: em, ef, num, nuf
    character(len=3) :: nume
    character(len=8) :: nompar
    character(len=16) :: nomres(33)
    character(len=32) :: phenom
    character(len=16) :: nomte
!
!---------------------------------------------------------------------
    zero = 0.0d0
    deux = 2.0d0
    elasco = 0
    multic = 0
    coupmf = .false.
    call r8inir(9, zero, dmf, 1)
    call r8inir(6, zero, dmc, 1)
    call r8inir(6, zero, dfc, 1)
    dc = 0.d0
    dci = 0.d0
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    call tecael(iazi, iazk24, noms=0)
    nomte = zk24(iazk24-1+3+nno+1) (1:16)
!
    call jevech('PCACOQU', 'L', jcoqu)
    epais = zr(jcoqu)
    alpha = zr(jcoqu+1)*r8dgrd()
    beta = zr(jcoqu+2)*r8dgrd()
    excent = zr(jcoqu+4)
    call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jcou)
    if (iret .eq. 0) then
        ncou = zi(jcou)
        npgh = 3
    else
        npgh = 1
        ncou = 1
    end if

    call tecach('NNO', 'PMATERC', 'L', iret, iad=jmate)
    if (iret .ne. 0) then
        multic = 0
        goto 999
    end if
    call rccoma(zi(jmate), 'ELAS', 1, phenom, icodre(1))

!   -- calcul de t1ve :
!   --------------------
    if (phenom .ne. 'ELAS') then
        call coqrep(pgl, alpha, beta, t2iu, t2ui, &
                    c, s)
!       CALCUL DE LA MATRICE T1VE DE PASSAGE D'UNE MATRICE
!       (3,3) DU REPERE DE LA VARIETE AU REPERE ELEMENT
!
        t1ve(1) = c*c
        t1ve(4) = s*s
        t1ve(7) = c*s
        t1ve(2) = t1ve(4)
        t1ve(5) = t1ve(1)
        t1ve(8) = -t1ve(7)
        t1ve(3) = -t1ve(7)-t1ve(7)
        t1ve(6) = t1ve(7)+t1ve(7)
        t1ve(9) = t1ve(1)-t1ve(4)
    end if

!   -- calcul de nbv et nomres :
!   ----------------------------
    if (phenom .eq. 'ELAS_COQMU') then
        nbv = 26
        do i = 1, nbv
            call codent(i, 'G', nume)
            nomres(i) = 'HOM_'//nume
        end do

    else if (phenom .eq. 'ELAS') then
        nbv = 2
        nomres(1) = 'E'
        nomres(2) = 'NU'

    else if (phenom .eq. 'ELAS_GLRC') then
        nbv = 6
        nomres(1) = 'E_M'
        nomres(2) = 'NU_M'
        nomres(3) = 'E_F'
        nomres(4) = 'NU_F'
        nomres(5) = 'BT1'
        nomres(6) = 'BT2'

    else if (phenom .eq. 'ELAS_COQUE') then
!        call utmess('A', 'ELEMENTS_93', sk=phenom)
!       -- on remplit nomres plus tard ...

    else if (phenom .eq. 'ELAS_DHRC') then
!      -- pour ELAS_DHRC, on n'utilise pas nomres

    else if (phenom .eq. 'ELAS_ORTH') then
        call utmess('F', 'ELEMENTS_91', sk=phenom)

    else if (phenom .eq. 'ELAS_ISTR') then
        call utmess('F', 'ELEMENTS_92', sk=phenom)

    else
        call utmess('F', 'ELEMENTS_42', sk=phenom)
    end if

    if (nomte .eq. 'MEDKQG4' .or. nomte .eq. 'MEDKTG3') then
        call moyte2(fami, npg, '+', temp, iret1)
    else
        call moytem(fami, npg, npgh*ncou, '+', temp, &
                    iret1)
    end if

    nbpar = 1
    nompar = 'TEMP'

    if (phenom .eq. 'ELAS') then
!
!        ------ MATERIAU ISOTROPE AVEC DECOUPLAGE MEMBRANE FLEXION----
!
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    nbv, nomres, valres, icodre, 1)
!
        young = valres(1)
        nu = valres(2)
!
        multic = 0
        kcis = 5.d0/6.d0
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = young*epais*epais*epais/12.d0/(1.d0-nu*nu)
        df(:, :) = 0.d0
        dmf(:, :) = 0.d0
        df(1, 1) = cdf
        df(1, 2) = cdf*nu
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
        df(3, 3) = cdf*(1.d0-nu)/2.d0
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = epais*young/(1.d0-nu*nu)
        dm(:, :) = 0.d0
        dm(1, 1) = cdm
        dm(1, 2) = cdm*nu
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
        dm(3, 3) = cdm*(1.d0-nu)/2.d0
!      --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
        gcis = young/2.d0/(1.d0+nu)
        cdc = gcis*kcis*epais
        dc(1, 1) = cdc
        dc(2, 2) = dc(1, 1)
        dc(1, 2) = 0.d0
        dc(2, 1) = 0.d0
!      --- CALCUL DE SON INVERSE ------------------------------------
        dci(1, 1) = 1.d0/dc(1, 1)
        dci(2, 2) = dci(1, 1)
        dci(1, 2) = 0.d0
        dci(2, 1) = 0.d0
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do
!
    else if (phenom .eq. 'ELAS_GLRC') then
!
!        ------ MATERIAU ISOTROPE --------------------------------------
!
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    2, nomres, valres, icodre, 1)
!
        em = valres(1)
        num = valres(2)
!
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    2, nomres(3), valres, icodre, 0)
!
        if (icodre(1) .eq. 0) then
            ef = valres(1)
        else
            ef = em
        end if
!
        if (icodre(2) .eq. 0) then
            nuf = valres(2)
        else
            nuf = num
        end if
!
        multic = 0
        kcis = 5.d0/6.d0
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = ef*epais*epais*epais/12.d0/(1.d0-nuf*nuf)
        df(:, :) = 0.d0
        dmf(:, :) = 0.d0
        df(1, 1) = cdf
        df(1, 2) = cdf*nuf
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
        df(3, 3) = cdf*(1.d0-nuf)/2.d0
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = epais*em/(1.d0-num*num)
        dm(:, :) = 0.d0
        dm(1, 1) = cdm
        dm(1, 2) = cdm*num
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
        dm(3, 3) = cdm*(1.d0-num)/2.d0
!      --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    2, nomres(5), valres, icodre, 0)
!
        if (icodre(1) .eq. 0) then
            dc(1, 1) = valres(1)
            dc(2, 2) = valres(2)
        else
            gcis = em/2.d0/(1.d0+num)
            cdc = gcis*kcis*epais
            dc(1, 1) = cdc
            dc(2, 2) = dc(1, 1)
        end if
!
        dc(1, 2) = 0.d0
        dc(2, 1) = 0.d0
!      --- CALCUL DE SON INVERSE ------------------------------------
        dci(1, 1) = 1.d0/dc(1, 1)
        dci(2, 2) = dci(1, 1)
        dci(1, 2) = 0.d0
        dci(2, 1) = 0.d0
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do
!
    else if (phenom .eq. 'ELAS_COQUE') then
!       call rcvalb(fami, 1, 1, '+', zi(jmate),&
!                   ' ', phenom, nbpar, nompar, [temp],&
!                   nbv, nomres, valres, icodre, 1)

!       -- on recupere TOUS les parametres de ELAS_COQUE :
        call rcvalt(fami, 1, 1, '+', zi(jmate), ' ', &
                    'ELAS_COQUE', nbpar, nompar, [temp], &
                    45, v_elasco, c_elasco, 1)

!       -- selon le type d'elasticite :
        if (c_elasco(1) .eq. 0) then
            elasco = 1
            nbv = 10
            nomres(1) = 'MEMB_L  '
            nomres(2) = 'MEMB_LT '
            nomres(3) = 'MEMB_T  '
            nomres(4) = 'MEMB_G_LT'
            nomres(5) = 'FLEX_L  '
            nomres(6) = 'FLEX_LT '
            nomres(7) = 'FLEX_T  '
            nomres(8) = 'FLEX_G_LT'
            nomres(9) = 'CISA_L  '
            nomres(10) = 'CISA_T  '
            do k = 1, nbv
                ASSERT(c_elasco(k) .eq. 0)
                valres(k) = v_elasco(k)
            end do
        else
            ASSERT(c_elasco(11) .eq. 0)
            elasco = 2
            nbv = 33
            coupmf = .true.
            nomres(1) = 'M_LLLL  '
            nomres(2) = 'M_LLTT  '
            nomres(3) = 'M_LLLT  '
            nomres(4) = 'M_TTTT  '
            nomres(5) = 'M_TTLT  '
            nomres(6) = 'M_LTLT  '
            nomres(7) = 'F_LLLL  '
            nomres(8) = 'F_LLTT  '
            nomres(9) = 'F_LLLT  '
            nomres(10) = 'F_TTTT  '
            nomres(11) = 'F_TTLT  '
            nomres(12) = 'F_LTLT  '
            nomres(13) = 'MF_LLLL '
            nomres(14) = 'MF_LLTT '
            nomres(15) = 'MF_LLLT '
            nomres(16) = 'MF_TTTT '
            nomres(17) = 'MF_TTLT '
            nomres(18) = 'MF_LTLT '
            nomres(19) = 'MC_LLLZ '
            nomres(20) = 'MC_LLTZ '
            nomres(21) = 'MC_TTLZ '
            nomres(22) = 'MC_TTTZ '
            nomres(23) = 'MC_LTLZ '
            nomres(24) = 'MC_LTTZ '
            nomres(25) = 'FC_LLLZ '
            nomres(26) = 'FC_LLTZ '
            nomres(27) = 'FC_TTLZ '
            nomres(28) = 'FC_TTTZ '
            nomres(29) = 'FC_LTLZ '
            nomres(30) = 'FC_LTTZ '
            nomres(31) = 'C_LZLZ  '
            nomres(32) = 'C_LZTZ  '
            nomres(33) = 'C_TZTZ  '
            do k = 1, nbv
                ASSERT(c_elasco(10+k) .eq. 0)
                valres(k) = v_elasco(10+k)
            end do
        end if

        if (elasco .eq. 1) then
            multic = 0
!
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = valres(1)
            dm(1, 2) = valres(2)
            dm(1, 3) = 0.d0
            dm(2, 1) = dm(1, 2)
            dm(2, 2) = valres(3)
            dm(2, 3) = 0.d0
            dm(3, 1) = 0.d0
            dm(3, 2) = 0.d0
            dm(3, 3) = valres(4)
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = valres(5)
            df(1, 2) = valres(6)
            df(1, 3) = 0.d0
            df(2, 1) = df(1, 2)
            df(2, 2) = valres(7)
            df(2, 3) = 0.d0
            df(3, 1) = 0.d0
            df(3, 2) = 0.d0
            df(3, 3) = valres(8)
!        --- COUPLAGE  MEMBRANE FLEXION --------------------------------
            dmf(1, 1) = 0.d0
            dmf(1, 2) = 0.d0
            dmf(1, 3) = 0.d0
            dmf(2, 1) = 0.d0
            dmf(2, 2) = 0.d0
            dmf(2, 3) = 0.d0
            dmf(3, 1) = 0.d0
            dmf(3, 2) = 0.d0
            dmf(3, 3) = 0.d0
!        --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
            dc(1, 1) = valres(9)
            dc(1, 2) = 0.d0
            dc(2, 1) = 0.d0
            dc(2, 2) = valres(10)
!        --- CALCUL DE SON INVERSE -------------------------------------
            dci(1, 1) = 1/valres(9)
            dci(1, 2) = 0.d0
            dci(2, 1) = 0.d0
            dci(2, 2) = 1/valres(10)
!
        else if (elasco .eq. 2) then
            multic = 0
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = valres(1)
            dm(1, 2) = valres(2)
            dm(1, 3) = valres(3)
            dm(2, 1) = dm(1, 2)
            dm(2, 2) = valres(4)
            dm(2, 3) = valres(5)
            dm(3, 1) = dm(1, 3)
            dm(3, 2) = dm(2, 3)
            dm(3, 3) = valres(6)
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = valres(7)
            df(1, 2) = valres(8)
            df(1, 3) = valres(9)
            df(2, 1) = df(1, 2)
            df(2, 2) = valres(10)
            df(2, 3) = valres(11)
            df(3, 1) = df(1, 3)
            df(3, 2) = df(2, 3)
            df(3, 3) = valres(12)
!        --- COUPLAGE  MEMBRANE FLEXION --------------------------------
            dmf(1, 1) = valres(13)
            dmf(1, 2) = valres(14)
            dmf(1, 3) = valres(15)
            dmf(2, 1) = dmf(1, 2)
            dmf(2, 2) = valres(16)
            dmf(2, 3) = valres(17)
            dmf(3, 1) = dmf(1, 3)
            dmf(3, 2) = dmf(2, 3)
            dmf(3, 3) = valres(18)
!        --- COUPLAGE  MEMBRANE CISAILLEMENT ---------------------------
            dmc(1, 1) = valres(19)
            dmc(1, 2) = valres(20)
            dmc(2, 1) = valres(21)
            dmc(2, 2) = valres(22)
            dmc(3, 1) = valres(23)
            dmc(3, 2) = valres(24)
!        --- COUPLAGE  FLEXION CISAILLEMENT ---------------------------
            dfc(1, 1) = valres(25)
            dfc(1, 2) = valres(26)
            dfc(2, 1) = valres(27)
            dfc(2, 2) = valres(28)
            dfc(3, 1) = valres(29)
            dfc(3, 2) = valres(30)
!        --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
            dc(1, 1) = valres(31)
            dc(1, 2) = valres(32)
            dc(2, 1) = dc(1, 2)
            dc(2, 2) = valres(33)
!        --- CALCUL DE SON INVERSE -------------------------------------
            det = dc(1, 1)*dc(2, 2)-dc(1, 2)*dc(2, 1)
            if (det .gt. r8prem()) then
                dci(1, 1) = dc(2, 2)/det
                dci(1, 2) = -dc(1, 2)/det
                dci(2, 1) = -dc(2, 1)/det
                dci(2, 2) = dc(1, 1)/det
            else
                call utmess('F', 'ELEMENTS_43')
            end if
        end if
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION ET DE LA --------
!        --- MATRICE DE COUPLAGE FLEXION-CISAILLEMENT          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
        do i = 1, 3
            do j = 1, 2
                dfc(i, j) = dfc(i, j)+excent*dmc(i, j)
            end do
        end do
!
        if (nomte .ne. 'MEDKQG4' .and. nomte .ne. 'MEDKTG3') then
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
            call utbtab('ZERO', 3, 3, dm, t1ve, &
                        xab1, dm)
            call utbtab('ZERO', 3, 3, df, t1ve, &
                        xab1, df)
            call utbtab('ZERO', 3, 3, dmf, t1ve, &
                        xab1, dmf)
            call utbtab('ZERO', 2, 2, dc, t2ui, &
                        xab2, dc)
            call utbtab('ZERO', 2, 2, dci, t2ui, &
                        xab2, dci)
            if (elasco .eq. 2) then
                call utdtab('ZERO', 3, 2, 2, 3, &
                            dmc, t2ui, t1ve, xab3, dmc)
                call utdtab('ZERO', 3, 2, 2, 3, &
                            dfc, t2ui, t1ve, xab3, dfc)
            end if
        end if
!
    else if (phenom .eq. 'ELAS_COQMU') then
!        ------ MATERIAU MULTICOUCHE -----------------------------------
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    18, nomres, valres, icodre, 1)
        dm(1, 1) = valres(1)
        dm(1, 2) = valres(2)
        dm(1, 3) = valres(3)
        dm(2, 2) = valres(4)
        dm(2, 3) = valres(5)
        dm(3, 3) = valres(6)
        dm(2, 1) = dm(1, 2)
        dm(3, 1) = dm(1, 3)
        dm(3, 2) = dm(2, 3)
        dmf(1, 1) = valres(7)
        dmf(1, 2) = valres(8)
        dmf(1, 3) = valres(9)
        dmf(2, 2) = valres(10)
        dmf(2, 3) = valres(11)
        dmf(3, 3) = valres(12)
        dmf(2, 1) = dmf(1, 2)
        dmf(3, 1) = dmf(1, 3)
        dmf(3, 2) = dmf(2, 3)
        df(1, 1) = valres(13)
        df(1, 2) = valres(14)
        df(1, 3) = valres(15)
        df(2, 2) = valres(16)
        df(2, 3) = valres(17)
        df(3, 3) = valres(18)
        df(2, 1) = df(1, 2)
        df(3, 1) = df(1, 3)
        df(3, 2) = df(2, 3)
        call rcvalb(fami, 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, [temp], &
                    6, nomres(21), valres(21), icodre(21), 1)
        dci(1, 1) = valres(21)
        dci(2, 2) = valres(22)
        dci(1, 2) = valres(23)
        dci(2, 1) = dci(1, 2)
        dc(1, 1) = valres(24)
        dc(2, 2) = valres(25)
        dc(1, 2) = valres(26)
        dc(2, 1) = dc(1, 2)
!
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
        call utbtab('ZERO', 3, 3, dm, t1ve, &
                    xab1, dm)
        call utbtab('ZERO', 3, 3, df, t1ve, &
                    xab1, df)
        call utbtab('ZERO', 3, 3, dmf, t1ve, &
                    xab1, dmf)
        call utbtab('ZERO', 2, 2, dc, t2ui, &
                    xab2, dc)
        call utbtab('ZERO', 2, 2, dci, t2ui, &
                    xab2, dci)
!
        multic = 1
!
    else if (phenom .eq. 'ELAS_DHRC') then
        multic = 0
        coupmf = .true.
!        ------ MATERIAU ELAS_DHRC -----------------------------------
        call rcadlv(fami, 1, 1, '+', zi(jmate), ' ', &
                    'ELAS_DHRC', 'A0', nbpar, nompar, [temp], jadr, n1, icodre1, 1)
        ASSERT(icodre1 .eq. 0 .and. n1 .eq. 21)

        dm(1, 1) = zr(jadr-1+1)
        dm(1, 2) = zr(jadr-1+2)
        dm(1, 3) = zr(jadr-1+3)
        dm(2, 2) = zr(jadr-1+7)
        dm(2, 3) = zr(jadr-1+8)
        dm(3, 3) = zr(jadr-1+12)
        dm(2, 1) = dm(1, 2)
        dm(3, 1) = dm(1, 3)
        dm(3, 2) = dm(2, 3)
        dmf(1, 1) = zr(jadr-1+4)
        dmf(1, 2) = zr(jadr-1+5)
        dmf(1, 3) = zr(jadr-1+6)
        dmf(2, 1) = zr(jadr-1+9)
        dmf(2, 2) = zr(jadr-1+10)
        dmf(2, 3) = zr(jadr-1+11)
        dmf(3, 1) = zr(jadr-1+13)
        dmf(3, 2) = zr(jadr-1+14)
        dmf(3, 3) = zr(jadr-1+15)
        df(1, 1) = zr(jadr-1+16)
        df(1, 2) = zr(jadr-1+17)
        df(1, 3) = zr(jadr-1+18)
        df(2, 2) = zr(jadr-1+19)
        df(2, 3) = zr(jadr-1+20)
        df(3, 3) = zr(jadr-1+21)
        df(2, 1) = df(1, 2)
        df(3, 1) = df(1, 3)
        df(3, 2) = df(2, 3)
!
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
        call utbtab('ZERO', 3, 3, dm, t1ve, &
                    xab1, dm)
        call utbtab('ZERO', 3, 3, df, t1ve, &
                    xab1, df)
        call utbtab('ZERO', 3, 3, dmf, t1ve, &
                    xab1, dmf)
    end if
!
    do i = 1, 3
        do j = 1, 3
            if (abs(dmf(i, j)) .gt. 1.d-10) coupmf = .true.
        end do
    end do
!
999 continue
!
end subroutine
