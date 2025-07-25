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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine te0590(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nifilg.h"
#include "asterfort/nifipd.h"
#include "asterfort/nifism.h"
#include "asterfort/niinit.h"
#include "asterfort/nmtstm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/tecach.h"
#include "asterfort/tgverm.h"
#include "asterfort/tgveri_use.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_INCO_UPG, AXIS_INCO_UPG, D_PLAN_INCO_UPG
! Options: FULL_MECA, FULL_MECA_ELAS
!          RAPH_MECA
!          RIGI_MECA_ELAS, RIGI_MECA_TANG
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno1, nno2, nno3, npg, nb_elrefe
    integer(kind=8) :: icoret, codret, iret
    integer(kind=8) :: iw, ivf1, ivf2, ivf3, idf1, idf2
    integer(kind=8) :: jtab(7), lgpg
    integer(kind=8) :: vu(3, 27), vg(27), vp(27), vpi(3, 27)
    integer(kind=8) :: igeom, imate, icontm, ivarim, iuse
    integer(kind=8) :: iinstm, iinstp, iddlm, iddld, icarcr, ivarix
    integer(kind=8) :: ivectu, icontp, ivarip, imatuu
    integer(kind=8) :: idbg, nddl, ia, ja
    real(kind=8) :: angl_naut(3)
    character(len=8) :: list_elrefe(10), typmod(2)
    aster_logical :: matsym
    character(len=16) :: defo_comp, rela_comp, type_comp
    aster_logical :: lVect, lMatr, lVari, lSigm, lMatrPred
!     POUR TGVERI
    real(kind=8) :: sdepl(135), svect(135), scont(6*27)
    real(kind=8) :: epsilo, epsilp, epsilg
    real(kind=8) :: tab_out(27*3*27*3)
    real(kind=8), pointer :: varia(:) => null(), smatr(:) => null()
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: na, os, nb, ib, kk
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    idbg = 0
    codret = 0
    matsym = ASTER_TRUE
    typmod = ' '
!
! - List of ELREFE
!
    call elref2(nomte, 10, list_elrefe, nb_elrefe)
    ASSERT(nb_elrefe .ge. 3)
!
! - Get shape functions
!
! - PRES
    call elrefe_info(elrefe=list_elrefe(3), fami='RIGI', nno=nno3, jvf=ivf3)
! - GONF
    call elrefe_info(elrefe=list_elrefe(2), fami='RIGI', nno=nno2, jvf=ivf2, jdfde=idf2)
! - DX, DY, DZ
    call elrefe_info(elrefe=list_elrefe(1), fami='RIGI', ndim=ndim, nno=nno1, npg=npg, &
                     jpoids=iw, jvf=ivf1, jdfde=idf1)
!
    nddl = nno1*ndim+nno2+nno3
!
! - Modelling
!
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Get index of dof
!
    call niinit(typmod, ndim, nno1, nno2, nno3, &
                0, vu, vg, vp, vpi)
!
! - Get input fields
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', iddlm)
    call jevech('PDEPLPR', 'L', iddld)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
!
! - Get orientation
!
    call getElemOrientation(ndim, nno1, igeom, angl_naut)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, lMatr, lVect, lVari, &
                         lSigm, codret)
    lMatrPred = option .eq. 'RIGI_MECA_TANG'
!
! - Properties of behaviour
!
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    type_comp = compor(INCRELAS)
!
! - Get output fields
!
    if (lMatr) then
        call nmtstm(zr(icarcr), imatuu, matsym)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', icoret)
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
    call tgveri_use(option, zr(icarcr), compor, iuse)
    if (iuse == 1) then
        allocate (varia(2*135*135))
        allocate (smatr(135*135))
    end if
!
100 continue
! - PETITES DEFORMATIONS
    if (defo_comp(1:6) .eq. 'PETIT ') then
        call nifipd(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    vu, vg, vp, zr(igeom), typmod, &
                    option, zi(imate), compor, lgpg, zr(icarcr), &
                    zr(iinstm), zr(iinstp), zr(iddlm), zr(iddld), angl_naut, &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), lMatr, &
                    lVect, zr(ivectu), zr(imatuu), codret)
    else if (defo_comp .eq. 'GDEF_LOG') then
!
! - PARAMETRES EN SORTIE
        call nifilg(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    vu, vg, vp, zr(igeom), typmod, &
                    option, zi(imate), compor, lgpg, zr(icarcr), &
                    zr(iinstm), zr(iinstp), zr(iddlm), zr(iddld), angl_naut, &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), lMatr, &
                    lVect, lSigm, lVari, zr(ivectu), zr(imatuu), &
                    matsym, codret)
    else if (defo_comp .eq. 'SIMO_MIEHE') then
!
! - PARAMETRES EN SORTIE
        typmod(2) = 'INCO'
        call nifism(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    idf2, vu, vg, vp, zr(igeom), &
                    typmod, option, zi(imate), compor, lgpg, &
                    zr(icarcr), zr(iinstm), zr(iinstp), zr(iddlm), zr(iddld), &
                    angl_naut, zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), &
                    lMatr, lVect, lMatrPred, zr(ivectu), zr(imatuu), &
                    codret)
    else
        call utmess('F', 'ELEMENTS3_16', sk=defo_comp)
    end if
!
    if (codret .ne. 0) goto 200
!       Calcul eventuel de la matrice TGTE par PERTURBATION
    call tgverm(option, zr(icarcr), compor, nno1, nno2, &
                nno3, zr(igeom), ndim, nddl, zr(iddld), &
                sdepl, vu, vg, vp, zr(ivectu), &
                svect, ndim*2*npg, zr(icontp), scont, npg*lgpg, &
                zr(ivarip), zr(ivarix), zr(imatuu), smatr, matsym, &
                epsilo, epsilp, epsilg, varia, iret)
    if (iret .ne. 0) goto 100
!
200 continue
!
!
! - Save return code
!
    if (lSigm) then
        zi(icoret) = codret
    end if
!
    if (idbg .eq. 1) then
        if (lMatr) then
            write (6, *) 'MATRICE TANGENTE'
            if (matsym) then
                do ia = 1, nddl
                    write (6, '(108(1X,E11.4))') (zr(imatuu+(ia*(ia-1)/2)+ja-1), ja=1, ia)
                end do
            else
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
!
                write (6, *) 'KUU'
                ja = 1
                do na = 1, nno1
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nno1
                            do ib = 1, ndim
                                kk = os+vu(ib, nb)
                                tab_out(ja) = zr(imatuu+kk-1)
                                ja = ja+1
                            end do
                        end do
                    end do
                end do
                do ia = 1, nno1*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno1*ndim+ja), ja=1, nno1*ndim)
                end do
!
! - TERME K:GG      KGG(NNO2,NNO2)
!
                write (6, *) 'KGG'
                ja = 1
                do na = 1, nno2
                    os = (vg(na)-1)*nddl
                    do ia = 1, nno2
                        kk = os+vg(ia)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nno2
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno2+ja), ja=1, nno2)
                end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO3)
!
                write (6, *) 'KUP'
                ja = 1
                do na = 1, nno1
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nno3
                            kk = os+vp(nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nno1*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno3+ja), ja=1, nno3)
                end do
!
! - TERME K:PU      KPU(NDIM,NNO3,NNO1)
!
                write (6, *) 'KPU'
                ja = 1
                do ia = 1, nno3
                    os = (vp(ia)-1)*nddl
                    do nb = 1, nno1
                        do ib = 1, ndim
                            kk = os+vu(ib, nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nno3
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno1*ndim+ja), ja=1, nno1*ndim)
                end do
!
! - TERME K:UG      KUG(NDIM,NNO1,NNO2)
!
                write (6, *) 'KUG'
                ja = 1
                do na = 1, nno1
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nno2
                            kk = os+vg(nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nno1*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno2+ja), ja=1, nno2)
                end do
!
! - TERME K:GU      KGU(NDIM,NNO2,NNO1)
!
                write (6, *) 'KGU'
                ja = 1
                do ia = 1, nno2
                    os = (vg(ia)-1)*nddl
                    do nb = 1, nno1
                        do ib = 1, ndim
                            kk = os+vu(ib, nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nno2
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno1*ndim+ja), ja=1, nno1*ndim)
                end do
!
! - TERME K:PG      KPG(NNO3,NNO2)
!
                write (6, *) 'KPG'
                ja = 1
                do ia = 1, nno3
                    os = (vp(ia)-1)*nddl
                    do ib = 1, nno2
                        kk = os+vg(ib)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nno3
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno2+ja), ja=1, nno2)
                end do
!
! - TERME K:GP      KPG(NNO2,NNO3)
!
                write (6, *) 'KGP'
                ja = 1
                do ia = 1, nno2
                    os = (vg(ia)-1)*nddl
                    do ib = 1, nno3
                        kk = os+vp(ib)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nno2
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nno3+ja), ja=1, nno3)
                end do
            end if
        end if
        if (lVect) then
            write (6, *) 'FORCE INTERNE'
            write (6, '(108(1X,E11.4))') (zr(ivectu+ja-1), ja=1, nddl)
        end if
    end if
!
! - Free large arrays
    if (iuse == 1) then
        deallocate (smatr)
        deallocate (varia)
    end if
!
end subroutine
