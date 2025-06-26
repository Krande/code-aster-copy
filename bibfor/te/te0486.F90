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
subroutine te0486(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/antisy.h"
#include "asterfort/assert.h"
#include "asterfort/b1tdb2.h"
#include "asterfort/btsig.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matdn.h"
#include "asterfort/provec.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
    character(len=16) :: option, nomte
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT EN PRESSION SUIVEUSE
!          POUR LES PLAQUES ET COQUES
!          MATRICE DE RIGIDITE TANGENTE POUR LES COQUES 3D
!          (TOUJOURS EN PRESSION SUIVEUSE)
!
!          OPTIONS : 'CHAR_MECA_PRSU_R '
!                    'RIGI_MECA_PRSU_R '
!                    'CHAR_MECA_PRSU_F '
!                    'RIGI_MECA_PRSU_F '
!                    'CHAR_MECA_SRCO3D '
!                    'RIGI_MECA_SRCO3D' (COQUES 3D SEULEMENT)
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
! ----------------------------------------------------------------------
    character(len=24) :: valk
! ----------------------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano
    integer(kind=8) :: igeom, ium, iup, ipres, jpres
    integer(kind=8) :: ino, lzi, lzr, iadzi, iazk24, iret
    integer(kind=8) :: i, j, in, kn, ii, komptn, nb1, nb2
    integer(kind=8) :: ivectu, imatun, intsn, npgsn
    integer(kind=8) :: irco3d, ifco3d, itemps, ierz, mxnoeu
    parameter(mxnoeu=9)
    real(kind=8) :: pres, presno(9), madn(3, 51), nks1(3, 51), nks2(3, 51)
    real(kind=8) :: a1(3), a2(3), anta1(3, 3), anta2(3, 3), surf(3)
    real(kind=8) :: rigns(2601), pr
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), valpar(7)
    real(kind=8) :: geom_reac(3*mxnoeu)
    aster_logical :: locapr
    character(len=8) :: nomail, nompar(7)
!     ------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLMR', 'L', ium)
    call jevech('PDEPLPR', 'L', iup)
!
!
!     --- AUTRES ELEMENTS QUE LA COQUE 3 D ---
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
    call tecach('NNO', 'PPRESSR', 'L', iret, iad=ipres)
!
!
    if (option .eq. 'CHAR_MECA_PRSU_R' .or. option .eq. 'RIGI_MECA_PRSU_R') then
        call jevech('PPRESSR', 'L', jpres)
    else if (option .eq. 'CHAR_MECA_PRSU_F' .or. option .eq. 'RIGI_MECA_PRSU_F') then
        call jevech('PPRESSF', 'L', jpres)
    end if
!
    if (nomte .ne. 'MEC3QU9H' .and. nomte .ne. 'MEC3TR7H') then
        if (ipres .eq. 0) then
            call jevech('PFRCO3D', 'L', ipres)
        end if
!
        do ino = 0, nno-1
            if (zr(ipres+ino) .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk = nomail
                call utmess('F', 'ELEMENTS4_93', sk=valk)
            end if
        end do
!
    else
!
!     --- ELEMENTS COQUES 3D ---
!
        call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
!        --- NOMBRE DE NOEUDS ( NB1 : SERENDIP , NB2 : LAGRANGE )
        nb1 = zi(lzi-1+1)
        nb2 = zi(lzi-1+2)
        ASSERT(nb2-1 .lt. mxnoeu)
!        --- NBRE POINTS INTEGRATIONS (NPGSN : NORMALE )
        npgsn = zi(lzi-1+4)
!
!        ---  ( FONCTIONS DE FORMES, DERIVEES ET POIDS )
        call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
!        ---- RECUPERATION DES POINTEURS ( E : ECRITURE ) SELON OPTION
        if (option(1:16) .eq. 'CHAR_MECA_SRCO3D' .or. option(1:16) .eq. &
            'CHAR_MECA_SFCO3D' .or. option(1:16) .eq. 'CHAR_MECA_PRSU_R' .or. &
            option(1:16) .eq. 'CHAR_MECA_PRSU_F') then
!           --- VECTEUR DES FORCES INTERNES
            call jevech('PVECTUR', 'E', ivectu)
        end if
!
        if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. &
            'RIGI_MECA_SFCO3D' .or. option(1:16) .eq. 'RIGI_MECA_PRSU_R' .or. &
            option(1:16) .eq. 'RIGI_MECA_PRSU_F') then
!         --- MATRICE TANGENTE DE RIGIDITE ET INITIALISATION
!
!           CALL JEVECH ( 'PMATUUR' , 'E' , IMATUU )
!            POUR UNE MATRICE NON SYMETRIQUE (VOIR AUSSI MECGME)
            call jevech('PMATUNS', 'E', imatun)
            call r8inir(51*51, 0.d0, rigns, 1)
        end if
!
!        ---- RECUPERATION DE L ADRESSE DES VARIABLES NODALES TOTALES
        call jevech('PDEPLMR', 'L', ium)
        call jevech('PDEPLPR', 'L', iup)
!
!        ---  CALCUL DE LA GEOMETRIE REACTUALISEE
!
        do in = 1, nb2-1
            do ii = 1, 3
                geom_reac(3*(in-1)+ii) = zr( &
                                                 igeom-1+3*(in-1)+ii)+zr(ium-1 &
                                                 &+6*(in-1)+ii)+zr(iup-1+6*(i&
                                                 &n-1)+ii &
                                                 )
            end do
        end do
!
!        ---- RECUPERATION DU VECTEUR DE PRESSION NODALE A INTERPOLER
        if (option(10:16) .eq. '_SRCO3D') then
            call jevech('PFRCO3D', 'L', irco3d)
            locapr = abs(zr(irco3d-1+7)-3.d0) .lt. 1.d-3
            if (locapr) then
                do kn = 1, nb2
                    presno(kn) = zr(irco3d-1+(kn-1)*7+3)
                end do
            else
                do kn = 1, nb2
                    do i = 1, 6
                        if (abs(zr(irco3d-1+(kn-1)*7+i)) .gt. r8prem()) then
                            call utmess('F', 'CHARGES_11')
                        end if
                    end do
                end do
                presno(1:nb2) = 0.d0
            end if
!
        else if (option(10:16) .eq. '_SFCO3D') then
            call jevech('PFFCO3D', 'L', ifco3d)
            call jevech('PINSTR', 'L', itemps)
            valpar(4) = zr(itemps)
            nompar(4) = 'INST'
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            nompar(5) = 'XF'
            nompar(6) = 'YF'
            nompar(7) = 'ZF'
!
            do i = 1, nb2
                presno(i) = 0.d0
            end do
!
            locapr = zk8(ifco3d+6) .eq. 'LOCAL_PR'
!
            do in = 0, nb2-1
                valpar(1) = zr(igeom+3*in)
                valpar(2) = zr(igeom+3*in+1)
                valpar(3) = zr(igeom+3*in+2)
                valpar(5) = geom_reac(3*in+1)
                valpar(6) = geom_reac(3*in+2)
                valpar(7) = geom_reac(3*in+3)
!
                if (locapr) then

                    call fointe('FM', zk8(ifco3d+2), 4, nompar, valpar, &
                                pr, ierz)
                    presno(in+1) = pr
                    if (ierz .ne. 0) then
                        call utmess('F', 'ELEMENTS4_1')
                    end if
                else
                    do i = 1, 6
                        call fointe('FM', zk8(ifco3d-1+i), 4, nompar, valpar, &
                                    pr, ierz)
                        if (abs(pr) .gt. r8prem()) call utmess('F', 'CHARGES_11')
                    end do
                    presno(in+1) = 0.d0
                end if
!
            end do
!
!
        else if (option(10:16) .eq. '_PRSU_R') then
            do kn = 1, nb2
                presno(kn) = zr(jpres-1+(kn-1)*1+1)
            end do
!
        else if (option(10:16) .eq. '_PRSU_F') then
            call jevech('PINSTR', 'L', itemps)
            valpar(4) = zr(itemps)
            nompar(4) = 'INST'
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            do j = 0, nb1-1
                valpar(1) = zr(igeom+3*j)
                valpar(2) = zr(igeom+3*j+1)
                valpar(3) = zr(igeom+3*j+2)
                call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                            pr, ierz)
                if (ierz .ne. 0) then
                    call utmess('F', 'ELEMENTS4_1')
                end if
                if (pr .ne. 0.d0) then
                    call tecael(iadzi, iazk24)
                    nomail = zk24(iazk24-1+3) (1:8)
                    valk = nomail
                    call utmess('F', 'ELEMENTS4_92', sk=valk)
                end if
            end do
        end if
!
!        ---- VECTEURS TANGENTS A1 ET A2 AUX NOEUDS NON NORMALISES
!
        call vectan(nb1, nb2, geom_reac, zr(lzr), vecta, &
                    vectn, vectpt)
!
!        ---- BOUCLE SUR LES POINTS D INTEGRATION NORMALE
        do intsn = 1, npgsn
!
!          ---- VECTEURS DE BASE AUX POINTS DE GAUSS A KSI3 = 0.D0
            call r8inir(3, 0.d0, a1, 1)
            call r8inir(3, 0.d0, a2, 1)
!
!          ---- INTERPOLATIONS PRESSION
!               VECTEURS TANGENTS NON NORMES A1 A2
            pres = 0.d0
!
            do kn = 1, nb2
!
                pres = pres+zr(lzr-1+459+9*(intsn-1)+kn)*presno(kn)
!
                do ii = 1, 3
!
                    a1(ii) = a1(ii)+zr(lzr-1+459+9*( &
                                       intsn-1)+kn)*vecta(kn, 1, ii)
!
                    a2(ii) = a2(ii)+zr(lzr-1+459+9*( &
                                       intsn-1)+kn)*vecta(kn, 2, ii)
!
                end do
!
            end do
!
!          ---- A1 VECTORIEL A2
            call provec(a1, a2, surf)
!
!          --- MATRICE D INTERPOLATION POUR LES DEPLACEMENTS
            call matdn(nb1, zr(lzr), intsn, madn, nks1, &
                       nks2)
!
            if (option(1:16) .eq. 'CHAR_MECA_SRCO3D' .or. option(1:16) .eq. &
                'CHAR_MECA_SFCO3D' .or. option(1:16) .eq. 'CHAR_MECA_PRSU_R') then
!
!             --- FORCE EXTERNE NODALE AU SIGNE DU TE0423 ET NMPR3D
                call btsig(6*nb1+3, 3, -pres*zr(lzr-1+127+intsn-1), madn, surf, &
                           zr(ivectu))
            end if
!
            if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. &
                'RIGI_MECA_SFCO3D' .or. option(1:16) .eq. 'RIGI_MECA_PRSU_R') then
!
!             --- MATRICE ANTISYM DE A1 ET DE A2
                call antisy(a1, 1.d0, anta1)
                call antisy(a2, 1.d0, anta2)
!
!            --- PREMIER TERME
                call b1tdb2(madn, nks2, anta1, pres*zr(lzr-1+127+intsn-1), 3, &
                            6*nb1+3, rigns)
!
!           --- DEUXIEME TERME
                call b1tdb2(madn, nks1, anta2, -pres*zr(lzr-1+127+intsn-1), 3, &
                            6*nb1+3, rigns)
            end if
!
        end do
!
!
        if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. &
            'RIGI_MECA_SFCO3D' .or. option(1:16) .eq. 'RIGI_MECA_PRSU_R') then
!
!        --- PARTIE SYMETRIQUE DE LA MATRICE TANGENTE
!
            komptn = 0
!            KOMPTU = 0
            do j = 1, 6*nb1+3
!
!            POUR UNE MATRICE NON SYMETRIQUE (VOIR AUSSI MECGME)
                do i = 1, 6*nb1+3
                    zr(imatun+komptn) = -rigns((6*nb1+3)*(i-1)+ &
                                               j)
                    komptn = komptn+1
                end do
!
!              DO 100  I = 1 , J
!                 KOMPTU = KOMPTU + 1
!                 ZR ( IMATUU - 1 + KOMPTU ) = -0.5D0 *
!     &              (    RIGNS ( ( 6 * NB1 + 3 ) * ( J - 1 ) + I )
!     &              +  RIGNS ( ( 6 * NB1 + 3 ) * ( I - 1 ) + J )  )
!
!  100         CONTINUE
!
            end do
!
        end if
!
    end if
!
end subroutine
