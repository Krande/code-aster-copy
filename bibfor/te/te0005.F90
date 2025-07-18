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
subroutine te0005(option, nomte)
!
    use dil_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dilcar.h"
#include "asterfort/dilele.h"
#include "asterfort/dilini.h"
#include "asterfort/fnodil.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/terefe.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN_DIL and 3D_DIL
! Formulations: DIL and DIL_INCO
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          FORC_NODA, REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lSigm, lMatr, lVect
    integer(kind=8) :: ivf, ivf2, idfde, idfde2, jgano, ndim, lgpg, ipoids, npi
    integer(kind=8) :: ipoid2, dimdef, ichg, ichn
    integer(kind=8) :: icontm, ideplm, ideplp, igeom, imate, jcret, nddls, nddlm
    integer(kind=8) :: imatuu, ivectu, icontp, nnos, nnom, dimuel, dimcon
    integer(kind=8) :: ivarim, ivarip, icarcr, iinstm, iinstp
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: typmod(2)
    real(kind=8) :: sigref, lagref, epsref
    real(kind=8), allocatable:: sref(:)
    type(dil_modelisation) :: ds_dil
!
! --------------------------------------------------------------------------------------------------
!
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lVect = L_VECT(option)
!
! - Get adresses for fields
!

    call dilcar(option, compor, icontm, ivarim, ideplm, ideplp, &
                igeom, imate, imatuu, ivectu, icontp, &
                ivarip, ichg, ichn, jcret, icarcr, iinstm, iinstp)

! ======================================================================
! --- INITIALISATION DES VARIABLES DE L'ELEMENT ------------------------
! ======================================================================
    call dilini(option, ivf, ivf2, idfde, &
                idfde2, jgano, ndim, ipoids, ipoid2, &
                npi, dimdef, nddls, nddlm, lgpg, &
                dimcon, typmod, dimuel, nnom, nnos, ds_dil)

! ======================================================================
! --- CALCUL DES OPTIONS -----------------------------------------------
! ======================================================================
    if (option(1:9) .eq. 'RIGI_MECA') then
        call dilele(option, typmod, ds_dil, ndim, nnos, &
                    nnom, npi, dimuel, dimdef, ipoids, zr(ivf), &
                    zr(ivf2), idfde, idfde2, zr(igeom), compor, &
                    zi(imate), lgpg, zr(icarcr), zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), zr(icontm), zr(ivarim), &
                    zr(icontm), zr(ivarim), &
                    zr(ivectu), zr(imatuu), lMatr, lVect, lSigm, zi(jcret))
    else if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) &
             .eq. 'FULL_MECA') then
        call dilele(option, typmod, ds_dil, ndim, nnos, &
                    nnom, npi, dimuel, dimdef, ipoids, zr(ivf), &
                    zr(ivf2), idfde, idfde2, zr(igeom), compor, &
                    zi(imate), lgpg, zr(icarcr), zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), zr(icontm), zr(ivarim), &
                    zr(icontp), zr(ivarip), &
                    zr(ivectu), zr(imatuu), lMatr, lVect, lSigm, zi(jcret))
    else if (option .eq. 'FORC_NODA') then
        call fnodil(option, typmod, ds_dil, ndim, nnos, &
                    nnom, npi, dimuel, dimdef, ipoids, zr(ivf), &
                    zr(ivf2), idfde, idfde2, zr(igeom), compor, &
                    zr(icontm), zr(ivectu))
    else if (option .eq. 'REFE_FORC_NODA') then

        allocate (sref(dimdef))

        call terefe('SIGM_REFE', 'MECA_DIL', sigref)
        call terefe('LAGR_REFE', 'MECA_DIL', lagref)
        call terefe('EPSI_REFE', 'MECA_DIL', epsref)

        if (ndim .eq. 2) then
            sref(1:dimdef) = [sigref, sigref, sigref, sigref, lagref, &
                              lagref, sigref, sigref, epsref]
        else if (ndim .eq. 3) then
            sref(1:dimdef) = [sigref, sigref, sigref, sigref, sigref, &
                              sigref, lagref, lagref, sigref, &
                              sigref, sigref, epsref]
        end if

        call fnodil(option, typmod, ds_dil, ndim, nnos, &
                    nnom, npi, dimuel, dimdef, ipoids, zr(ivf), &
                    zr(ivf2), idfde, idfde2, zr(igeom), compor, &
                    transpose(spread(sref, 1, npi)), &
                    zr(ivectu))

        deallocate (sref)

    else
        ASSERT(ASTER_FALSE)
    end if

!
end subroutine
