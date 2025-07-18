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

subroutine te0150(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!     CALCULE LE CHARGEMENT INDUIT PAR UNE ELEVATION UNIFORME DE
!     TEMPERATURE DANS LES POUTRES D'EULER ET DE TIMOSHENKO
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!       'FC1D1D_MECA'       : FORCES LINEIQUES (COMP)
!       'FR1D1D_MECA'       : FORCES LINEIQUES (REEL)
!       'FF1D1D_MECA'       : FORCES LINEIQUES (FONCTION)
!       'SR1D1D_MECA'       : FORCES LINEIQUES SUIVEUSES (REEL)
!       'SF1D1D_MECA'       : FORCES LINEIQUES SUIVEUSES (FONCTION)
!       'CHAR_MECA_PESA_R'  : CHARGES DE PESANTEUR
!       'CHAR_MECA_ROTA_R'  : CHARGES DE ROTATION
!       'CHAR_MECA_TEMP_R'  : DEFORMATIONS THERMIQUES
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!       'MECA_POU_D_E'  : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!       'MECA_POU_D_T'  : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!       'MECA_POU_D_EM' : POUTRE DROITE MULTIFIBRE D EULER (SECT. CONST)
!       'MECA_POU_D_TG' : POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!       'MECA_POU_D_TGM': POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!                         MULTI-FIBRES SECTION CONSTANTE
!       'MECA_POU_D_SQUE' : SQUELETTE D'ASSEMBLAGE COMBUSTIBLE
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmfrig.h"
#include "asterfort/porigi.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ptfocp.h"
#include "asterfort/ptforp.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "asterfort/verifm.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbpar, lmater, iret
    integer(kind=8) :: lorien, lvect
    integer(kind=8) :: itype, nc, ind, i, j, ncf
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids
    integer(kind=8) :: ivf, idfdx, jgano
    real(kind=8) :: valpar(3)
    real(kind=8) :: e, nu, g
    real(kind=8) :: a1, a2, xl
    real(kind=8) :: pgl(3, 3), de(18), ffe(18)
    real(kind=8) :: bsm(18, 18), matk(171)
    real(kind=8) :: epsith
    real(kind=8) :: fr(18), fi(18), fgr(18), fgi(18)
    real(kind=8) :: fer(18), fei(18)
!
    character(len=4) :: fami
    character(len=8) :: nompar(3), materi
    character(len=16) :: ch16
    aster_logical :: lrho
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbres = 2
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
    data nomres/'E', 'NU'/
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter          :: nb_cara = 3
    real(kind=8)                :: vale_cara(nb_cara)
    character(len=8), parameter :: noms_cara(nb_cara) = (/'A1  ', 'A2  ', 'TVAR'/)
! --------------------------------------------------------------------------------------------------
!
!   POUR LA PESANTEUR ET LA ROTATION, ON N'A BESOIN DE RHO QUI EST CONSTANT DANS LA MAILLE
    lrho = (option .eq. 'CHAR_MECA_PESA_R') .or. (option .eq. 'CHAR_MECA_ROTA_R')
!
    fami = 'RIGI'
    nno = 2
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
    a1 = vale_cara(1)
    a2 = vale_cara(2)
    itype = nint(vale_cara(3))
!
    if (nomte .eq. 'MECA_POU_D_E') then
!       poutre droite d'euler à 6 DDL
        nc = 6
        ncf = 6
    else if (nomte .eq. 'MECA_POU_D_T') then
!       poutre droite de timoskenko à 6 ddl
        nc = 6
        ncf = 6
    else if (nomte .eq. 'MECA_POU_D_EM') then
!       poutre multifibre droite d'euler à 6 DDL
        nc = 6
        ncf = 6
        itype = 20
        if (lrho) itype = 0
    else if (nomte .eq. 'MECA_POU_D_TG') then
!       poutre droite de timoskenko à 7 ddl (gauchissement)
        nc = 7
        ncf = 6
        itype = 30
    else if (nomte .eq. 'MECA_POU_D_TGM') then
!       poutre droite de timoskenko à 7 ddl (gauchissement, multifibres)
        nc = 7
        ncf = 6
        itype = 30
    else if (nomte .eq. 'MECA_POU_D_SQUE') then
!       element squelette 9 DDL
        nc = 9
        ncf = 9
        itype = 20
        if (lrho) itype = 0
    else
        ASSERT(.FALSE.)
    end if
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
!   RECUPERATION DES CARACTERISTIQUES MATERIAUX (MOYENNE)
    nbpar = 1
    nompar(1) = 'TEMP'
    if ((option(13:16) .ne. '1D1D') .and. (.not. lrho)) then
        call jevech('PMATERC', 'L', lmater)
        call moytem(fami, npg, 1, '+', valpar(1), iret)
    end if
!
!
    materi = ' '
    e = 0.d0
    if ((nomte .ne. 'MECA_POU_D_EM') .and. (nomte .ne. 'MECA_POU_D_TGM') .and. &
        (nomte .ne. 'MECA_POU_D_SQUE')) then
!       poutres classiques
        if ((option(13:16) .ne. '1D1D') .and. (.not. lrho)) then
            valres(:) = 0.d0
            call rcvalb(fami, 1, 1, '+', zi(lmater), materi, 'ELAS', nbpar, nompar, valpar, &
                        nbres, nomres, valres, codres, 1)
!
            e = valres(1)
            nu = valres(2)
            g = e/(2.d0*(1.d0+nu))
        end if
    end if
!
!   RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
    xl = lonele()
!
!   passage du repere local au repere global
    if (option .eq. 'CHAR_MECA_FC1D1D') then
        call jevech('PVECTUC', 'E', lvect)
        if (nomte .eq. 'MECA_POU_D_TG' .or. nomte .eq. 'MECA_POU_D_TGM') then
            call ptfocp(itype, option, xl, nno, 6, pgl, fr, fi)
            call utpvlg(nno, 6, pgl, fr, fgr)
            call utpvlg(nno, 6, pgl, fi, fgi)
            do i = 1, 6
                zc(lvect+i-1) = dcmplx(fgr(i), fgi(i))
                zc(lvect+i-1+7) = dcmplx(fgr(i+6), fgi(i+6))
            end do
            zc(lvect+7-1) = dcmplx(0.d0, 0.d0)
            zc(lvect+14-1) = dcmplx(0.d0, 0.d0)
        else
            call ptfocp(itype, option, xl, nno, nc, pgl, fr, fi)
            call utpvlg(nno, nc, pgl, fr, fgr)
            call utpvlg(nno, nc, pgl, fi, fgi)
            do i = 1, 12
                zc(lvect+i-1) = dcmplx(fgr(i), fgi(i))
            end do
        end if
    else if ((option .eq. 'CHAR_MECA_FR1D1D') .or. &
             (option .eq. 'CHAR_MECA_FF1D1D') .or. &
             (option .eq. 'CHAR_MECA_SR1D1D') .or. &
             (option .eq. 'CHAR_MECA_SF1D1D') .or. &
             (option .eq. 'CHAR_MECA_ROTA_R') .or. &
             (option .eq. 'CHAR_MECA_PESA_R')) then
        if ((nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
            call ptforp(0, option, nomte, a1, a2, xl, 1, nno, ncf, pgl, fer, fei)
        else
            call ptforp(itype, option, nomte, a1, a2, xl, 1, nno, ncf, pgl, fer, fei)
        end if
        ffe(:) = 0.0
        do i = 1, ncf
            ffe(i) = fer(i)
            ffe(i+nc) = fer(i+ncf)
        end do
    else
        matk(:) = 0.d0
        if ((nomte .eq. 'MECA_POU_D_EM') .or. (nomte .eq. 'MECA_POU_D_TGM') .or. &
            (nomte .eq. 'MECA_POU_D_SQUE')) then
!           poutre droite multifibre a section constante
            call pmfrig(nomte, zi(lmater), matk)
        else
            call porigi(nomte, e, nu, xl, matk)
        end if
!       remplissage de la matrice carree
        ind = 0
        do i = 1, nc*2
            de(i) = 0.d0
            do j = 1, i-1
                ind = ind+1
                bsm(i, j) = matk(ind)
                bsm(j, i) = matk(ind)
            end do
            ind = ind+1
            bsm(i, i) = matk(ind)
        end do
!
        if (option .eq. 'CHAR_MECA_TEMP_R') then
!           calcul du deplacement local induit par l'elevation de temp.
            call verifm(fami, npg, 1, '+', zi(lmater), epsith, iret)
        else
            ch16 = option
            call utmess('F', 'ELEMENTS2_47', sk=ch16)
        end if
        de(1) = -epsith*xl
        de(1+nc) = -de(1)
!       calcul des forces induites
        ffe(:) = 0.0
        do i = 1, nc
            do j = 1, nc
                ffe(i) = ffe(i)+bsm(i, j)*de(j)
                ffe(i+nc) = ffe(i+nc)+bsm(i+nc, j+nc)*de(j+nc)
            end do
        end do
    end if
!
    if (option .ne. 'CHAR_MECA_FC1D1D') then
        call jevech('PVECTUR', 'E', lvect)
!       matrice de passage du repere global au repere local : PGL
        call utpvlg(nno, nc, pgl, ffe, zr(lvect))
    end if
!
end subroutine
