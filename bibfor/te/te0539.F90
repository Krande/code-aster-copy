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
! person_in_charge: samuel.geniaut at edf.fr
!
subroutine te0539(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ltequa.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmtstm.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/xnmel.h"
#include "asterfort/xnmpl.h"
#include "asterfort/xteddl.h"
#include "asterfort/xteini.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dcopy.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: XFEM
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          RIGI_MECA (linear)
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: typmod(2), enr, lag
    character(len=16) :: compor(4), elref
    integer(kind=8) :: jgano, nno, npg, i, imatuu, lgpg, ndim, iret, nfiss
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: iinstm, iinstp, ideplm, ideplp, icarcr
    integer(kind=8) :: ivectu, icontp, ivarip, li, jcret, codret
    integer(kind=8) :: ivarix
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jbaslo, jlsn, jlst, jstno, jpmilt, jheavn
    integer(kind=8) :: jtab(7), nnos, idim, jfisno
    integer(kind=8) :: nfh, ddlc, nddl, nnom, nfe, ibid, ddls, ddlm
    aster_logical :: matsym, l_nonlin, l_line
    real(kind=8) :: bary(3), crit(1), sig(1), vi(1)
    character(len=16), pointer :: ccompor(:) => null()
    character(len=16) :: defo_comp, rela_comp, type_comp
    aster_logical :: lVect, lMatr, lVari, lSigm, lMatrPred
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    codret = 0
!
! - Get element parameters
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    call elref1(elref)
    ASSERT(nno .le. 27)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, ibid, ddlc, &
                nnom, ddls, nddl, ddlm, nfiss, &
                ibid)
!
! - Type of finite element
!
    if (ndim .eq. 3) then
        typmod(1) = '3D'
        typmod(2) = ' '
    else
        if (lteatt('AXIS', 'OUI')) then
            typmod(1) = 'AXIS'
        else if (lteatt('C_PLAN', 'OUI')) then
            typmod(1) = 'C_PLAN'
        else if (lteatt('D_PLAN', 'OUI')) then
            typmod(1) = 'D_PLAN'
        else
            ASSERT(lteatt('C_PLAN', 'OUI'))
        end if
        typmod(2) = ' '
    end if
!
! - General options
!
    l_nonlin = (option(1:9) .eq. 'FULL_MECA') .or. (option .eq. 'RAPH_MECA') .or. &
               (option(1:10) .eq. 'RIGI_MECA_')
    l_line = option .eq. 'RIGI_MECA'
    ASSERT(l_nonlin .or. l_line)
    lVect = ASTER_FALSE
    lMatr = ASTER_FALSE
    lVari = ASTER_FALSE
    lSigm = ASTER_FALSE
    if (l_line) then
        lMatr = ASTER_TRUE
    end if
    lMatrPred = option(1:9) .eq. 'RIGI_MECA'
!
! - Get input fields
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
    if (l_nonlin) then
        ASSERT(.not. l_line)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PCOMPOR', 'L', vk16=ccompor)
        call jevech('PCARCRI', 'L', icarcr)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                    itab=jtab)
        lgpg = max(jtab(6), 1)*jtab(7)
        npg = jtab(2)
    end if
!
! - Get input fields for XFEM
!
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PBASLOR', 'L', jbaslo)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call jevech('PSTANO', 'L', jstno)
    call teattr('S', 'XFEM', enr, ibid)
    if (enr(1:2) .eq. 'XH') then
        call jevech('PHEA_NO', 'L', jheavn)
    end if
    if ((ibid .eq. 0) .and. ltequa(elref, enr)) then
        call jevech('PPMILTO', 'L', jpmilt)
    end if
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
    end if
!
! - Rigidity matrix (linear)
!
    if (l_line) then
        call jevech('PMATUUR', 'E', imatuu)
        lgpg = 0
        compor = ' '
        call xnmel(nno, nfh, nfe, ddlc, ddlm, &
                   igeom, typmod, option, zi(imate), compor, &
                   lgpg, crit, jpintt, zi(jcnset), zi(jheavt), &
                   zi(jlonch), zr(jbaslo), zr(iinstm), zr(iinstp), ibid, &
                   zr(jlsn), zr(jlst), sig, vi, zr(imatuu), &
                   ibid, codret, jpmilt, nfiss, jheavn, &
                   jstno, l_line, l_nonlin, lMatr, lVect, &
                   lSigm)
        matsym = .true.
        goto 999
    end if
!
! - Compute barycentric center
!
    bary(:) = 0.d0
    do i = 1, nno
        do idim = 1, ndim
            bary(idim) = bary(idim)+zr(igeom+idim+ndim*(i-1)-1)/nno
        end do
    end do
!
! - Select objects to construct from option name (non-linear case)
!
    call behaviourOption(option, ccompor, lMatr, lVect, lVari, &
                         lSigm, codret)
!
! - Properties of behaviour
!
    rela_comp = ccompor(RELA_NAME)
    defo_comp = ccompor(DEFO)
    type_comp = ccompor(INCRELAS)
    if (defo_comp .ne. 'PETIT') call utmess('F', 'XFEM_50', sk=defo_comp)
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
! - HYPER-ELASTICITE
!
    if (type_comp .eq. 'COMP_ELAS') then
        if (lMatrPred) then
            call xnmel(nno, nfh, nfe, ddlc, ddlm, &
                       igeom, typmod, option, zi(imate), ccompor, &
                       lgpg, zr(icarcr), jpintt, zi(jcnset), zi(jheavt), &
                       zi(jlonch), zr(jbaslo), zr(iinstm), zr(iinstp), ideplm, &
                       zr(jlsn), zr(jlst), zr(icontm), zr(ivarim), zr(imatuu), &
                       ivectu, codret, jpmilt, nfiss, jheavn, &
                       jstno, l_line, l_nonlin, lMatr, lVect, &
                       lSigm)
        else
            do li = 1, nddl
                zr(ideplp+li-1) = zr(ideplm+li-1)+zr(ideplp+li-1)
            end do
            call xnmel(nno, nfh, nfe, ddlc, ddlm, &
                       igeom, typmod, option, zi(imate), ccompor, &
                       lgpg, zr(icarcr), jpintt, zi(jcnset), zi(jheavt), &
                       zi(jlonch), zr(jbaslo), zr(iinstm), zr(iinstp), ideplp, &
                       zr(jlsn), zr(jlst), zr(icontp), zr(ivarip), zr(imatuu), &
                       ivectu, codret, jpmilt, nfiss, jheavn, &
                       jstno, l_line, l_nonlin, lMatr, lVect, &
                       lSigm)
        end if
    else
!
! - HYPO-ELASTICITE
!
        if (defo_comp .eq. 'PETIT') then
            call xnmpl(nno, nfh, nfe, ddlc, ddlm, &
                       igeom, zr(iinstm), zr(iinstp), ideplp, zr(icontm), &
                       zr(ivarip), typmod, option, zi(imate), ccompor, &
                       lgpg, zr(icarcr), jpintt, zi(jcnset), zi(jheavt), &
                       zi(jlonch), zr(jbaslo), ideplm, zr(jlsn), zr(jlst), &
                       zr(icontp), zr(ivarim), zr(imatuu), ivectu, codret, &
                       jpmilt, nfiss, jheavn, jstno, lMatr, &
                       lVect, lSigm)
        else
            ASSERT(.false.)
        end if
    end if
!
999 continue
!
!     SUPPRESSION DES DDLS SUPERFLUS
    call teattr('C', 'XLAG', lag, ibid)
    if (ibid .eq. 0 .and. lag .eq. 'ARETE') then
        nno = nnos
    end if
!
!   OPTIONS RELATIVES A UNE MATRICE UNIQUEMENT
    if (option .eq. 'RIGI_MECA') then
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), .false._1, matsym, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    mat=zr(imatuu))
!   OPTIONS RELATIVES A UN VECTEUR UNIQUEMENT
    else if (option .eq. 'RAPH_MECA') then
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), .false._1, matsym, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    vect=zr(ivectu))
!   OPTIONS RELATIVES A UNE MATRICE ET UN VECTEUR
    else if (option .eq. 'FULL_MECA' .or. option .eq. 'RIGI_MECA_TANG') then
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), .false._1, matsym, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    mat=zr(imatuu), vect=zr(ivectu))
    else
        ASSERT(.false.)
    end if
!
! - Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = codret
    end if
!
end subroutine
