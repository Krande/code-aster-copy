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
subroutine te0313(option, nomte)
!
    use THM_type
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterfort/aseihm.h"
#include "asterfort/caeihm.h"
#include "asterfort/eiangl.h"
#include "asterfort/fneihm.h"
#include "asterfort/jevech.h"
#include "asterfort/poeihm.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_JOINT_HYME
!           PLAN_JOINT_HYME
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          FORC_NODA, VARI_ELNO, SIEF_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jgano, imatuu, ndim, imate, iinstm, jcret, nb_strain_meca
    integer(kind=8) :: iret, ichg, ichn, itabin(7), itabou(7)
    integer(kind=8) :: ivf2
    integer(kind=8) :: idf2, npi, npg
    integer(kind=8) :: codret
    integer(kind=8) :: ipoids, ivf1, idf1, igeom
    integer(kind=8) :: iinstp, ideplm, ideplp, icamas
    integer(kind=8) :: icontm, ivarip, ivarim, ivectu, icontp, jvSief
    integer(kind=8) :: mecani(8), press1(9), press2(9), tempe(5), dimuel
    integer(kind=8) :: dimdef, dimcon, nbvari, nb_vari_meca
    integer(kind=8) :: nno1, nno2
    integer(kind=8) :: iadzi, iazk24
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: iu(3, 18), ip(2, 9), ipf(2, 2, 9), iq(2, 2, 9)
    real(kind=8) :: r(22)
    real(kind=8) :: ang(24)
    character(len=3) :: modint
    character(len=8) :: nomail
    type(THM_DS) :: ds_thm
    integer(kind=8) :: li
    aster_logical :: axi
    aster_logical :: fnoevo
    real(kind=8) :: dt
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: compor_copy(COMPOR_SIZE)
    integer(kind=8) :: iCompor
!
! --------------------------------------------------------------------------------------------------
!

!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm)
!
! - Preparation
!
    call caeihm(ds_thm, nomte, axi, mecani, press1, &
                press2, tempe, dimdef, dimcon, ndim, &
                nno1, nno2, npi, npg, dimuel, &
                ipoids, ivf1, idf1, ivf2, idf2, &
                jgano, iu, ip, ipf, iq, &
                modint)
!
    call tecael(iadzi, iazk24)
    nomail = zk24(iazk24-1+3) (1:8)
!
! RECUPERATION DES ANGLES NAUTIQUES DEFINIS PAR AFFE_CARA_ELEM
    if ((option .eq. 'FORC_NODA') .or. (option(1:9) .eq. 'RIGI_MECA') .or. &
        (option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
!
        call jevech('PCAMASS', 'L', icamas)
        if (zr(icamas) .lt. 0.d0) then
            call utmess('F', 'ELEMENTS5_48')
        end if
!
! DEFINITION DES ANGLES NAUTIQUES AUX NOEUDS SOMMETS
        call eiangl(ndim, nno2, zr(icamas+1), ang)
    end if
! =====================================================================
! --- DEBUT DES DIFFERENTES OPTIONS -----------------------------------
! =====================================================================
! --- 2. OPTIONS : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA -------------
! =====================================================================
    if ((option(1:9) .eq. 'RIGI_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'FULL_MECA')) then
! ----- Input fields
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PMATERC', 'L', imate)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PCONTMR', 'L', icontm)

! ----- Make copy of COMPOR map
        do iCompor = 1, COMPOR_SIZE
            compor_copy(iCompor) = compor(iCompor)
        end do

! ----- Force DEFO_LDC="MECANIQUE" for THM
        compor_copy(DEFO_LDC) = "MECANIQUE"

        read (compor_copy(NVAR), '(I16)') nbvari

! ----- Select objects to construct from option name
        call behaviourOption(option, compor_copy, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)
! ----- Output fields
        imatuu = ismaem()
        ivectu = ismaem()
        icontp = ismaem()
        ivarip = ismaem()
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVari) then
            call jevech('PVARIPR', 'E', ivarip)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
            call jevech('PCODRET', 'E', jcret)
            zi(jcret) = 0
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
! ----- Integration
        codret = 0
        if (option(1:9) .eq. 'RIGI_MECA') then
            call aseihm(ds_thm, option, &
                        lSigm, lVari, lMatr, lVect, &
                        axi, ndim, nno1, nno2, &
                        npi, npg, dimuel, dimdef, dimcon, &
                        nbvari, zi(imate), iu, ip, ipf, &
                        iq, mecani, press1, press2, tempe, &
                        zr(ivf1), zr(ivf2), zr(idf2), zr(iinstm), zr(iinstp), &
                        zr(ideplm), zr(ideplm), zr(icontm), zr(icontp), zr(ivarim), &
                        zr(ivarim), nomail, zr(ipoids), zr(igeom), ang, &
                        compor_copy, zr(ivectu), zr(imatuu), &
                        codret)
        else
            do li = 1, dimuel
                zr(ideplp+li-1) = zr(ideplm+li-1)+zr(ideplp+li-1)
            end do
            call aseihm(ds_thm, option, &
                        lSigm, lVari, lMatr, lVect, &
                        axi, ndim, nno1, nno2, &
                        npi, npg, dimuel, dimdef, dimcon, &
                        nbvari, zi(imate), iu, ip, ipf, &
                        iq, mecani, press1, press2, tempe, &
                        zr(ivf1), zr(ivf2), zr(idf2), zr(iinstm), zr(iinstp), &
                        zr(ideplm), zr(ideplp), zr(icontm), zr(icontp), zr(ivarim), &
                        zr(ivarip), nomail, zr(ipoids), zr(igeom), ang, &
                        compor_copy, zr(ivectu), zr(imatuu), &
                        codret)
            if (lSigm) then
                zi(jcret) = codret
            end if
        end if
    end if
!
! ======================================================================
! --- 3. OPTION : FORC_NODA --------------------------------------------
! ======================================================================
    if (option .eq. 'FORC_NODA') then
! ======================================================================
! --- PARAMETRES EN ENTREE ---------------------------------------------
! ======================================================================
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PMATERC', 'L', imate)

        fnoevo = .false.
        dt = 0.d0
!
! ======================================================================
! --- PARAMETRES EN SORTIE ---------------------------------------------
! ======================================================================
        call jevech('PVECTUR', 'E', ivectu)
!
        call fneihm(ds_thm, fnoevo, dt, nno1, nno2, &
                    npi, npg, zr(ipoids), iu, ip, &
                    ipf, iq, zr(ivf1), zr(ivf2), zr(idf2), &
                    zr(igeom), ang, zr(jvSief), r, zr(ivectu), &
                    mecani, press1, press2, dimdef, &
                    dimcon, dimuel, ndim, axi)
!
    end if
!
! ======================================================================
! --- 4. OPTION : SIEF_ELNO ---------------------------------------
! ======================================================================
    if (option .eq. 'SIEF_ELNO') then
        call jevech('PCONTRR', 'L', ichg)
        call jevech('PSIEFNOR', 'E', ichn)
        nb_strain_meca = mecani(6)
        call poeihm(nomte, option, modint, jgano, nno1, &
                    nno2, dimcon, nb_strain_meca, zr(ichg), zr(ichn))
    end if
!
! ======================================================================
! --- 5. OPTION : VARI_ELNO ---------------------------------------
! ======================================================================
    if (option .eq. 'VARI_ELNO') then
        call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, itab=itabin)
        call tecach('OOO', 'PVARINR', 'E', iret, nval=7, itab=itabou)
        ichg = itabin(1)
        ichn = itabou(1)
!
        call jevech('PCOMPOR', 'L', vk16=compor)
        read (compor(NVAR), '(I16)') nbvari
        read (compor(MECA_NVAR), '(I16)') nb_vari_meca
        call poeihm(nomte, option, modint, jgano, nno1, &
                    nno2, nbvari, nb_vari_meca, zr(ichg), zr(ichn))
    end if
!
! ======================================================================
end subroutine
