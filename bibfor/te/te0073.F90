! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine te0073(option, nomte)
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_TEXT_F/R'
!                                   'CHAR_THER_RAYO_F/R'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
    integer :: nbres
    parameter (nbres=3)
    character(len=8) :: nompar(nbres), elrefe, alias8
    real(kind=8) :: valpar(nbres), poids, r, z, nx, ny, tpg, coen, coenp1, texn
    real(kind=8) :: texnp1, coorse(18), vectt(9), theta, sigm1, sigmn, eps1
    real(kind=8) :: epsn, tpf1, tpfn, tz0
    integer :: nno, nnos, jgano, ndim, kp, npg, ipoids, ivf, idfde, igeom
    integer :: itemps, ivectt, i, l, li, itex, icoefh, iray, itemp, nnop2
    integer :: c(6, 9), ise, nse, j, ier, icode, ibid
    aster_logical :: laxi, ltext
!
!
    call elref1(elrefe)
!
    if (lteatt('LUMPE','OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'SE3') elrefe='SE2'
    endif
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos,&
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    tz0 = r8t0()
    laxi = .false.
    if (lteatt('AXIS','OUI')) laxi = .true.
!
    if (option(11:14) .eq. 'TEXT') then
        ltext = .true.
    else if (option(11:14).eq.'RAYO') then
        ltext = .false.
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    endif
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
!
    if (ltext) then
! CHAR_.._TEXT
        call jevech('PTEMPER', 'L', itemp)
        call jevech('PCOEFHF', 'L', icoefh)
        call jevech('PT_EXTF', 'L', itex)
    else
! CHAR_..._RAYO
        call jevech('PRAYONF', 'L', iray)
        call jevech('PTEMPER', 'L', itemp)
! FIN DU IF LTEXT
    endif
!
! TRONC COMMUN
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PVECTTR', 'E', ivectt)
!
!====
! 1.3 PREALABLES LIES AUX CALCULS
!====
    theta = zr(itemps+2)
    call connec(nomte, nse, nnop2, c)
    do i = 1, nnop2
        vectt(i) = 0.d0
    end do
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'INST'
!
!====
! 2. CALCULS TERMES DE MASSE
!====
!
! BOUCLE SUR LES SOUS-ELEMENTS
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2* (i-1)+j) = zr(igeom-1+2* (c(ise,i)-1)+j)
            end do
        end do
!
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde,&
                        coorse, nx, ny, poids)
            r = 0.d0
            z = 0.d0
            tpg = 0.d0
            if (itemp .ne. 0) then
                do i = 1, nno
! CALCUL DE T-
                    l = (kp-1)*nno + i
                    tpg = tpg + zr(itemp-1+c(ise,i))*zr(ivf+l-1)
                end do
            endif
!
            do i = 1, nno
                l = (kp-1)*nno + i
                r = r + coorse(2* (i-1)+1)*zr(ivf+l-1)
                z = z + coorse(2* (i-1)+2)*zr(ivf+l-1)
            end do
            if (laxi) poids = poids*r
            valpar(1) = r
            valpar(2) = z
            valpar(3) = zr(itemps)
!
!====
! 2.1 OPTION CHAR_THER_TEXT_F/R
!====
            if (ltext) then
!
                call fointe('FM', zk8(icoefh), 3, nompar, valpar,&
                            coenp1, icode)
                ASSERT(icode.eq.0)
                if (theta .ne. 1.0d0) then
                    valpar(3) = zr(itemps) - zr(itemps+1)
                    call fointe('FM', zk8(icoefh), 3, nompar, valpar,&
                                coen, icode)
                    ASSERT(icode.eq.0)
                else
                    coen = 0.d0
                endif
!
                valpar(3) = zr(itemps)
                call fointe('FM', zk8(itex), 3, nompar, valpar,&
                            texnp1, icode)
                ASSERT(icode.eq.0)
                if (theta .ne. 1.0d0) then
                    valpar(3) = zr(itemps) - zr(itemps+1)
                    call fointe('FM', zk8(itex), 3, nompar, valpar,&
                                texn, icode)
                    ASSERT(icode.eq.0)
                else
                    texn = 0.d0
                endif
                do i = 1, nno
                    li = ivf + (kp-1)*nno + i - 1
                    vectt(c(ise,i)) = vectt(&
                                      c(ise,i)) + poids*zr(li)* (theta*coenp1*texnp1+ (1.0d0-thet&
                                      &a)*coen* (texn- tpg)&
                                      )
                end do
!====
! 2.2 OPTION CHAR_THER_RAYO_F/R
!====
            else
!
                call fointe('FM', zk8(iray), 3, nompar, valpar,&
                            sigm1, ier)
                ASSERT(ier.eq.0)
                if (theta .ne. 1.0d0) then
                    valpar(3) = zr(itemps) - zr(itemps+1)
                    call fointe('FM', zk8(iray), 3, nompar, valpar,&
                                sigmn, ier)
                    ASSERT(ier.eq.0)
                else
                    sigmn = 0.d0
                endif
!
                valpar(3) = zr(itemps)
                call fointe('FM', zk8(iray+1), 3, nompar, valpar,&
                            eps1, ier)
                ASSERT(ier.eq.0)
                if (theta .ne. 1.0d0) then
                    valpar(3) = zr(itemps) - zr(itemps+1)
                    call fointe('FM', zk8(iray+1), 3, nompar, valpar,&
                                epsn, ier)
                    ASSERT(ier.eq.0)
                else
                    epsn = 0.d0
                endif
!
                valpar(3) = zr(itemps)
                call fointe('FM', zk8(iray+2), 3, nompar, valpar,&
                            tpf1, ier)
                ASSERT(ier.eq.0)
                if (theta .ne. 1.0d0) then
                    valpar(3) = zr(itemps) - zr(itemps+1)
                    call fointe('FM', zk8(iray+2), 3, nompar, valpar,&
                                tpfn, ier)
                    ASSERT(ier.eq.0)
                else
                    tpfn = 0.d0
                endif
                do i = 1, nno
                    li = ivf + (kp-1)*nno + i - 1
                    vectt(c(ise,i)) = vectt(&
                                      c(ise,i)) + poids*zr(li)* (theta*sigm1*eps1* (tpf1+tz0)**4+&
                                      & (1.0d0-theta)* sigmn* epsn* ((tpfn+tz0)**4- (tpg+tz0)**4)&
                                      )
                end do
!
! FIN DU IF LTEXT
            endif
!
! FIN DE BOUCLE SUR LES PTS DE GAUSS
        end do
! FIN DE BOUCLE SUR LES SOUS-ELEMENTS
    end do
!
    do i = 1, nnop2
        zr(ivectt-1+i) = vectt(i)
    end do
end subroutine
