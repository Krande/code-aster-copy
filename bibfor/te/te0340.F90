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
subroutine te0340(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cgfint.h"
#include "asterfort/cginit.h"
#include "asterfort/cgtang.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: CABLE_GAINE
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: comporKit(COMPOR_SIZE)
    character(len=8) :: typmod(2), lielrf(10)
    integer(kind=8) :: nno1, nno2, npg, imatuu, lgpg, lgpg1, lgpg2
    integer(kind=8) :: iw, ivf1, idf1, igeom, imate
    integer(kind=8) :: npgn, idf1n
    integer(kind=8) :: ivf2, ino, i, nddl1
    integer(kind=8) :: ivarim, ivarip, iinstm, iinstp
    integer(kind=8) :: iddlm, iddld, icarcr
    integer(kind=8) :: ivectu, icontp
    integer(kind=8) :: ivarix
    integer(kind=8) :: jtab(7), jcret, codret
    integer(kind=8) :: ndim, iret, ntrou
    integer(kind=8) :: iu(3, 3), iuc(3), im(3), isect, icontm
    real(kind=8) :: tang(3, 3), a, geom(3, 3)
    character(len=16) :: defo_comp, rela_comp, rela_cpla
    aster_logical :: lVect, lMatr, lVari, lSigm
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    imatuu = 1
    ivectu = 1
    icontp = 1
    ivarip = 1
!
! - FONCTIONS DE FORME
!
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', jvf=ivf1, jdfde=idf1)
    call elrefe_info(elrefe=lielrf(1), fami='NOEU', nno=nno1, npg=npgn, jdfde=idf1n)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', nno=nno2, npg=npg, jpoids=iw, &
                     jvf=ivf2)
    ndim = 3
    nddl1 = 5
!
! - DECALAGE D'INDICE POUR LES ELEMENTS D'INTERFACE
    call cginit(nomte, iu, iuc, im)
!
! - TYPE DE MODELISATION
!
    typmod(1) = '1D'
    typmod(2) = ' '
!
!
! - Get input fields
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', iddlm)
    call jevech('PDEPLPR', 'L', iddld)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PCAGNBA', 'L', isect)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, lMatr, lVect, lVari, &
                         lSigm, codret)
!
! - Properties of behaviour
!
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    rela_cpla = compor(PLANESTRESS)
!
!     MISE A JOUR EVENTUELLE DE LA GEOMETRIE
!
    if (defo_comp .eq. 'PETIT') then
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr(igeom-1+(ino-1)*ndim+i)
            end do
        end do
    else if (defo_comp .eq. 'PETIT_REAC') then
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr( &
                               igeom-1+(ino-1)*ndim+i)+zr(iddlm-1+(ino-1)*nddl1+i)+zr(iddld-1+(&
                               &ino-1)*nddl1+i &
                               )
            end do
        end do
    else
        call utmess('F', 'CABLE0_6', sk=defo_comp)
    end if
!
!     DEFINITION DES TANGENTES
!
    call cgtang(3, nno1, npgn, geom, zr(idf1n), &
                tang)
!
!     SECTION DE LA BARRE
    a = zr(isect)
!
! - ON VERIFIE QUE PVARIMR ET PVARIPR ONT LE MEME NOMBRE DE V.I. :
!
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg1 = max(jtab(6), 1)*jtab(7)
!
    if (lVari) then
        call tecach('OOO', 'PVARIPR', 'E', iret, nval=7, &
                    itab=jtab)
        lgpg2 = max(jtab(6), 1)*jtab(7)
        ASSERT(lgpg1 .eq. lgpg2)
    end if
    lgpg = lgpg1
!
! - Get output fields
!
    if (lMatr) then
        call jevech('PMATUNS', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
!
! - FORCES INTERIEURES ET MATRICE TANGENTE
!
    do i = 1, COMPOR_SIZE
        comporKit(i) = compor(i)
    end do
    call cgfint(ndim, nno1, nno2, npg, zr(iw), &
                zr(ivf1), zr(ivf2), zr(idf1), geom, tang, &
                typmod, option, zi(imate), comporKit, lgpg, &
                zr(icarcr), zr(iinstm), zr(iinstp), zr(iddlm), zr(iddld), &
                iu, iuc, im, a, zr(icontm), &
                zr(ivarim), zr(icontp), zr(ivarip), zr(imatuu), zr(ivectu), &
                codret)
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
end subroutine
