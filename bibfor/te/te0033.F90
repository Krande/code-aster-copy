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

subroutine te0033(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/coqrep.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqedg.h"
#include "asterfort/dkqsie.h"
#include "asterfort/dktedg.h"
#include "asterfort/dktsie.h"
#include "asterfort/dsqedg.h"
#include "asterfort/dsqsie.h"
#include "asterfort/dstedg.h"
#include "asterfort/dstsie.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxsiro.h"
#include "asterfort/dxsit2.h"
#include "asterfort/dxsit3.h"
#include "asterfort/dxsith.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/q4gedg.h"
#include "asterfort/q4gsie.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/t3gedg.h"
#include "asterfort/t3gsie.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
    character(len=16) :: option, nomte
!     CALCUL DE CONTRAINTES, DEFORMATIONS, EFFORTS ET DEFORMATIONS
!     GENERALISES POUR LES ELEMENTS DKT, DKTG, DST, DKQ, DSQ ET Q4G
!     POUR UN MATERIAU ISOTROPE OU MULTICOUCHE
!         OPTIONS TRAITEES  ==>  SIEF_ELGA
!                                EPSI_ELGA
!                                DEGE_ELGA
!                                DEGE_ELNO
!     IN   K16   OPTION : NOM DE L'OPTION A CALCULER
!     IN   K16   NOMTE  : NOM DU TYPE_ELEMENT
!     ------------------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano
    integer(kind=8) :: iret, jcara, vali(2)
    integer(kind=8) :: jdepg, jeffg, jgeom, jmate, jsigm
    integer(kind=8) :: np, multic
    integer(kind=8) :: jnbspi, nbcou, icou
    integer(kind=8) :: icodre(1), kpg, spt
!
    real(kind=8) :: epi(1), epais, eptot, alpha, beta
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), r8bid, valr(2)
    real(kind=8) :: depl(24)
    real(kind=8) :: effgt(32), effpg(32)
    real(kind=8) :: t2iu(4), t2ui(4), c, s
!
    aster_logical :: dkg
!
    character(len=2) :: val
    character(len=3) :: num
    character(len=4) :: fami
    character(len=8) :: famil, poum
    character(len=16) :: nomres
    character(len=32) :: phenom
!     ------------------------------------------------------------------
!
    r8bid = 0.d0
!
    if (option(6:9) .eq. 'ELNO') then
        fami = 'NOEU'
    else
        fami = 'RIGI'
    end if
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
    if (option .ne. 'SIEF_ELGA' .and. option .ne. 'EPSI_ELGA' .and. option .ne. 'DEGE_ELNO' &
        .and. option .ne. 'DEGE_ELGA') then
! OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    end if
!
    dkg = .false.
    if ((nomte .eq. 'MEDKTG3') .or. (nomte .eq. 'MEDKQG4')) then
        dkg = .true.
    end if
!
    call r8inir(32, 0.d0, effgt, 1)
!
    call jevech('PGEOMER', 'L', jgeom)
!
! --- VERIFICATION DE LA COHERENCE DES INFORMATIONS
! --- PROVENANT DE DEFI_COQU_MULT ET DE AFFE_CARA_ELEM
!     ----------------------------------
    jnbspi = 0
    if (option .eq. 'SIEF_ELGA' .or. option .eq. 'EPSI_ELGA') then
        call jevech('PMATERC', 'L', jmate)
        call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jnbspi)
        if (iret .eq. 0) then
            nbcou = zi(jnbspi)
            icou = 0
            eptot = 0.d0
            epi(1) = 0.d0
            call jevech('PCACOQU', 'L', jcara)
            epais = zr(jcara)
5           continue
            icou = icou+1
            call codent(icou, 'G', num)
            call codent(1, 'G', val)
            nomres = 'C'//num//'_V'//val
            famil = 'FPG1'
            kpg = 1
            spt = 1
            poum = '+'
            call rcvalb(famil, kpg, spt, poum, zi(jmate), &
                        ' ', 'ELAS_COQMU', 0, ' ', [r8bid], &
                        1, nomres, epi, icodre(1), 0)
            if (icodre(1) .eq. 0) then
                eptot = eptot+epi(1)
                goto 5
            end if
            if (eptot .ne. 0.d0) then
                if ((icou-1) .ne. nbcou) then
                    vali(1) = icou-1
                    vali(2) = nbcou
                    call utmess('F', 'ELEMENTS3_51', ni=2, vali=vali)
                end if
                if (abs(epais-eptot)/epais .gt. 1.d-2) then
                    valr(1) = eptot
                    valr(2) = epais
                    call utmess('F', 'ELEMENTS3_52', nr=2, valr=valr)
                end if
            end if
        end if
    end if
!
    if (option(8:9) .eq. 'GA') then
        np = npg
    else if (option(8:9) .eq. 'NO') then
        np = nno
    end if
!
    if (nno .eq. 3) then
        call dxtpgl(zr(jgeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jgeom), pgl)
    end if
!
    call utpvgl(nno, 3, pgl, zr(jgeom), xyzl)
!
    call jevech('PCACOQU', 'L', jcara)
    alpha = zr(jcara+1)*r8dgrd()
    beta = zr(jcara+2)*r8dgrd()
    call coqrep(pgl, alpha, beta, t2iu, t2ui, &
                c, s)
!
    call jevech('PDEPLAR', 'L', jdepg)
    call utpvgl(nno, 6, pgl, zr(jdepg), depl)
!
!     ---------- CONTRAINTES ET DEFORMATIONS --------------------------
    if (option(1:9) .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jsigm)
!
        if (dkg) then
            nbcou = 1
        else
            call jevech('PNBSP_I', 'L', jnbspi)
            nbcou = zi(jnbspi)
            if (nbcou .le. 0) then
                call utmess('F', 'ELEMENTS_46')
            end if
        end if
!
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else
! TYPE D ELEMENT INVALIDE
            ASSERT(.false.)
        end if
!
        call rccoma(zi(jmate), 'ELAS', 1, phenom, icodre(1))
        if (phenom .eq. 'ELAS' .or. phenom .eq. 'ELAS_ORTH' .or. phenom .eq. 'ELAS_ISTR') then
            call dxsith(nomte, zi(jmate), zr(jsigm))
        else if (phenom .eq. 'ELAS_COQMU') then
            call dxsit2(nomte, pgl, zr(jsigm))
        elseif (phenom .eq. 'ELAS_COQUE') then
            call dxsit3(nomte, zi(jmate), pgl, zr(jsigm))
        else
            call utmess('F', 'PLATE1_1', nk=2, valk=[option, phenom])
        end if
!     ----------------------------
    else if (option(1:9) .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jsigm)
!
        if (dkg) then
            nbcou = 1
        else
            call jevech('PNBSP_I', 'L', jnbspi)
            nbcou = zi(jnbspi)
            if (nbcou .le. 0) then
                call utmess('F', 'ELEMENTS_46')
            end if
        end if
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(option, fami, xyzl, pgl, depl, &
                        nbcou, zr(jsigm))
        end if
        call dxsiro(np*nbcou*3, t2iu, zr(jsigm), zr(jsigm))
!     ----------------------------
    else if (option(1:9) .eq. 'DEGE_ELNO') then
        call jevech('PDEFOGR', 'E', jeffg)
!
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(xyzl, option, pgl, depl, effgt, &
                        multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(xyzl, option, pgl, depl, effgt)
        end if
!
! ---    PASSAGE DES DEFORMATIONS GENERALISEES DU REPERE INTRINSEQUE
! ---    A L'ELEMENT AU REPERE LOCAL DE LA COQUE
        call dxefro(np, t2iu, effgt, zr(jeffg))
!     ----------------------------
    else if (option(1:9) .eq. 'DEGE_ELGA') then
        call jevech('PDEFOPG', 'E', jeffg)
!
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(xyzl, option, pgl, depl, effpg, &
                        multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(xyzl, option, pgl, depl, effpg)
        end if
! ---    PASSAGE DES DEFORMATIONS GENERALISEES DU REPERE INTRINSEQUE
! ---    A L'ELEMENT AU REPERE LOCAL DE LA COQUE
        call dxefro(np, t2iu, effpg, zr(jeffg))
    end if
!
    if (option .eq. 'SIEF_ELGA') then
! ---    PASSAGE DES CONTRAINTES DANS LE REPERE UTILISATEUR :
        call cosiro(nomte, 'PCONTRR', 'E', 'IU', 'G', &
                    jsigm, 'S')
    end if
!
end subroutine
