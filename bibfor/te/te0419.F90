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
subroutine te0419(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btldth.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrth.h"
#include "asterfort/trnflg.h"
#include "asterfort/vectan.h"
#include "asterfort/vexpan.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          POUR LES ELEMENTS MEC3QU9H, MEC3TR7H
!                          OPTIONS : 'CHAR_MECA_TEMP_R'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
    integer(kind=8) :: nb1, nb2, nddle, npge, npgsr, npgsn
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), vecpt(9, 3, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: forthi(42), forcth(42), vecl(51)
    real(kind=8) :: young, nu, alpha
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ib, indic, indith, inte, intsn, intsr
    integer(kind=8) :: j, jcara, jgeom, jvecg, kwgt, lzi, lzr
!
    real(kind=8) :: epais, temper
!-----------------------------------------------------------------------
    parameter(npge=2)
    real(kind=8) :: epsval(npge), ksi3s2, ksi3
    data epsval/-0.577350269189626d0, 0.577350269189626d0/
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PVECTUR', 'E', jvecg)
!
!     RECUPERATION DES OBJETS
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
    call jevech('PCACOQU', 'L', jcara)
    epais = zr(jcara)
!
    nddle = 5*nb1+2
!
!
!     RECUPERATION DE LA TEMPERATURE DE REFERENCE
!
!
    do i = 1, nddle
        forcth(i) = 0.d0
    end do
!
    call vectan(nb1, nb2, zr(jgeom), zr(lzr), vecta, &
                vectn, vectpt)
!
    kwgt = 0
    do inte = 1, npge
        ksi3s2 = epsval(inte)/2.d0
!
!     CALCUL DE BTDMR, BTDSR : M=MEMBRANE , S=CISAILLEMENT , R=REDUIT
!
        do intsr = 1, npgsr
            call mahsms(0, nb1, zr(jgeom), ksi3s2, intsr, &
                        zr(lzr), epais, vectn, vectg, vectt, &
                        hsfm, hss)
!
            call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                        hsj1m, hsj1s)
!
            call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                        epais, vectpt, hsj1m, hsj1s, btdm, &
                        btds)
        end do
!
        do intsn = 1, npgsn
!
!     CALCUL DE BTDFN : F=FLEXION , N=NORMAL
!     ET DEFINITION DE WGT=PRODUIT DES POIDS ASSOCIES AUX PTS DE GAUSS
!                          (NORMAL) ET DU DETERMINANT DU JACOBIEN
!
            call mahsf(1, nb1, zr(jgeom), ksi3s2, intsn, &
                       zr(lzr), epais, vectn, vectg, vectt, &
                       hsf)
!
            call hsj1f(intsn, zr(lzr), epais, vectg, vectt, &
                       hsf, kwgt, hsj1fx, wgt)
!
            call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                       zr(lzr), epais, vectpt, hsj1fx, btdf)
!
!     CALCUL DE BTDMN, BTDSN
!     ET
!     FORMATION DE BTILD
!
            call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                        btdm, btdf, btds, btild)
!
            call matrth('MASS', npgsn, young, nu, alpha, &
                        indith)
!
!     CALCUL DU CHAMP DE TEMPERATURE ET(OU) DES EFFORTS THERMIQUES
!     INDIC=1 : TEMPERATURE ET EFFORTS THERMIQUES
!     INDIC=0 : TEMPERATURE
!
            indic = 1
            ksi3 = epsval(inte)
            call btldth('MASS', ksi3, nb1, intsn, btild, &
                        wgt, indic, young, nu, alpha, &
                        temper, forthi)
!
            do i = 1, nddle
                forcth(i) = forcth(i)+forthi(i)
            end do
!
        end do
    end do
!
    call vexpan(nb1, forcth, vecl)
    do i = 1, 3
        vecl(6*nb1+i) = 0.d0
    end do
!
    do ib = 1, nb2
        do i = 1, 2
            do j = 1, 3
                vecpt(ib, i, j) = vectpt(ib, i, j)
            end do
        end do
        vecpt(ib, 3, 1) = vectn(ib, 1)
        vecpt(ib, 3, 2) = vectn(ib, 2)
        vecpt(ib, 3, 3) = vectn(ib, 3)
    end do
!
    call trnflg(nb2, vecpt, vecl, zr(jvecg))
!
end subroutine
