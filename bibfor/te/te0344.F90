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

subroutine te0344(option, nomte)
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DU VECTEUR ELEMENTAIRE EFFORT GENERALISE,
!     POUR LES ELEMENTS DE POUTRE DE TIMOSHENKO AVEC GAUCHISSEMENT.
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        SIPM_ELNO
!        SIPO_ELNO
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        MECA_POU_D_TG : POUTRE DROITE DE TIMOSHENKO AVEC GAUCHISSEMENT
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmavec.h"
#include "asterfort/porigi.h"
#include "asterfort/posigr.h"
#include "asterfort/posipr.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ptforp.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/verift.h"
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lmater, jmat, nbmat, imat, icomp, nbpar, i, npg, nno, nc
    integer(kind=8) :: ncc, jeffo, iret, itype
    integer(kind=8) :: lorien, jdepl, lforcr, lforcf
    real(kind=8) :: valpar, zero, e, g, a, xl, epsith
    real(kind=8) :: nu, fe(12), fi(12), flr(14), klv(105)
    real(kind=8) :: ulr(14), ugr(14), pgl(14, 14), klc(14, 14)
    character(len=32) :: messk(2)
    aster_logical :: okopt
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbres
    parameter(nbres=2)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
    character(len=8) :: nompar
    data nomres/'E', 'NU'/
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 2
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'TVAR'/
! --------------------------------------------------------------------------------------------------
!
    okopt = (option .eq. 'SIPM_ELNO') .or. (option .eq. 'SIPO_ELNO')
    ASSERT(okopt)
!
!   recuperation des caracteristiques materiaux
    call jevech('PMATERC', 'L', lmater)
!
!   option valide avec un seul phenomene : elas
    jmat = zi(lmater)
    nbmat = zi(jmat)
!   un seul materiau
    if (nbmat .ne. 1) then
        messk(1) = option
        call utmess('F', 'ELEMENTS4_59', sk=messk(1))
    end if
!   le 1er materiau
    imat = jmat+zi(jmat+nbmat+1)
!   seul elas est autorise
    do icomp = 1, zi(imat+1)
        if (zk32(zi(imat)+icomp-1) .ne. 'ELAS') then
            messk(1) = option
            messk(2) = zk32(zi(imat)+icomp-1)
            call utmess('F', 'ELEMENTS4_64', nk=2, valk=messk)
        end if
    end do
!
    nbpar = 0
    nompar = '  '
    valpar = 0.d0
    zero = 0.d0
    valres(:) = zero
!
    npg = 3
    call moytem('RIGI', npg, 1, '+', valpar, iret)
    nbpar = 1
    nompar = 'TEMP'
    call rcvalb('RIGI', 1, 1, '+', zi(lmater), ' ', 'ELAS', nbpar, nompar, [valpar], &
                2, nomres, valres, codres, 1)
!
    e = valres(1)
    nu = valres(2)
    g = e/(2.d0*(1.d0+nu))
!   recuperation des carac des sections utiles, longueur et pgl
    nno = 2
    nc = 7
    ncc = 6
    call jevech('PCAORIE', 'L', lorien)
    xl = lonele()
    call matrot(zr(lorien), pgl)
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
    a = vale_cara(1)
    itype = nint(vale_cara(2))
!   calcul de la matrice de rigidite locale
    call porigi(nomte, e, nu, xl, klv)
!   matrice rigidite ligne > matrice rigidite carre
    call vecma(klv, 105, klc, 14)
!
    call jevech('PDEPLAR', 'L', jdepl)
    do i = 1, 14
        ugr(i) = zr(jdepl+i-1)
    end do
!   vecteur deplacement local  ULR = PGL * UGR
    call utpvgl(nno, nc, pgl, ugr, ulr)
!   vecteur effort local  FLR = KLC * ULR
    call pmavec('ZERO', 14, klc, ulr, flr)
!   tenir compte des efforts dus a la dilatation
    call verift('RIGI', npg, 1, '+', zi(lmater), epsth_=epsith)
    ugr(:) = 0.d0
    ugr(1) = -epsith*xl
    ugr(8) = -ugr(1)
!   calcul des forces induites
    do i = 1, 7
        flr(i) = flr(i)-klc(i, 1)*ugr(1)
        flr(i+7) = flr(i+7)-klc(i+7, 1+7)*ugr(1+7)
    end do
!   prise en compte des efforts repartis
    call tecach('ONO', 'PFR1D1D', 'L', iret, iad=lforcr)
    if (lforcr .ne. 0) then
        call ptforp(itype, 'CHAR_MECA_FR1D1D', nomte, a, a, xl, 1, nno, ncc, pgl, fe, fi)
        do i = 1, 6
            flr(i) = flr(i)-fe(i)
            flr(i+7) = flr(i+7)-fe(i+6)
        end do
    end if
!   prise en compte des efforts repartis (sous forme de fonction)
    call tecach('ONO', 'PFF1D1D', 'L', iret, iad=lforcf)
    if (lforcf .ne. 0) then
        call ptforp(itype, 'CHAR_MECA_FF1D1D', nomte, a, a, xl, 1, nno, ncc, pgl, fe, fi)
        do i = 1, 6
            flr(i) = flr(i)-fe(i)
            flr(i+7) = flr(i+7)-fe(i+6)
        end do
    end if
!   archivage
    if (option .eq. 'SIPM_ELNO') then
        call jevech('PSIMXRR', 'E', jeffo)
        do i = 1, 6
            fe(i) = flr(i)
            fe(i+6) = flr(i+7)
        end do
        call posigr(nomte, fe, zr(jeffo))
    else if (option .eq. 'SIPO_ELNO') then
        call jevech('PCONTPO', 'E', jeffo)
        do i = 1, 6
            fe(i) = flr(i)
            fe(i+6) = flr(i+7)
        end do
        call posipr(nomte, fe, zr(jeffo))
    end if
!
end subroutine
