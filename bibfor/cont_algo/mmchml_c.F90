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
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
subroutine mmchml_c(ds_contact, ligrcf, chmlcf, sddyna, time_incr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfr.h"
#include "asterfort/ndynlo.h"
#include "Contact_type.h"
#include "jeveux.h"
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: ligrcf
    character(len=19), intent(in) :: chmlcf
    character(len=19), intent(in) :: sddyna
    real(kind=8), intent(in) :: time_incr
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Create and fill input field
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  ligrcf           : name of LIGREL for contact element
! In  chmlcf           : name of CHAM_LEM for input field
! In  sddyna           : datastructure for dynamic
! In  time_incr        : time increment
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ncmp = 60
    integer(kind=8), parameter :: nceld1 = 4
    integer(kind=8), parameter :: nceld2 = 4
    integer(kind=8), parameter :: nceld3 = 4
    integer(kind=8) :: ztabf
    integer(kind=8) :: i_cont_poin, i_zone, nt_cont_poin
    integer(kind=8) :: vale_indx, decal
    aster_logical :: l_dyna
    integer(kind=8) :: dyna_form
    real(kind=8) :: coef_fric, glis_maxi
    integer(kind=8) :: i_algo_cont, i_algo_fric, i_reso_fric, i_reso_geom
    integer(kind=8) :: nt_liel, nb_grel, nb_liel, i_grel, i_liel
    character(len=24) :: chmlcf_celv
    integer(kind=8) :: jv_chmlcf_celv
    character(len=24) :: chmlcf_celd
    integer(kind=8), pointer :: v_chmlcf_celd(:) => null()
    integer(kind=8), pointer :: v_ligrcf_liel(:) => null()
    character(len=24) :: sdcont_tabfin, sdcont_jsupco
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    real(kind=8), pointer :: v_sdcont_jsupco(:) => null()
    character(len=24) :: sdcont_cychis
    real(kind=8), pointer :: v_sdcont_cychis(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    glis_maxi = 0.

! - Active functionnalities
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')

! - Access to contact objects
    sdcont_jsupco = ds_contact%sdcont_solv(1:14)//'.JSUPCO'
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    sdcont_cychis = ds_contact%sdcont_solv(1:14)//'.CYCHIS'
    call jeveuo(sdcont_jsupco, 'L', vr=v_sdcont_jsupco)
    call jeveuo(sdcont_tabfin, 'L', vr=v_sdcont_tabfin)
    call jeveuo(sdcont_cychis, 'L', vr=v_sdcont_cychis)
    ztabf = cfmmvd('ZTABF')
!
! - Get_contact parameters
!
    i_reso_fric = cfdisi(ds_contact%sdcont_defi, 'ALGO_RESO_FROT')
    i_reso_geom = cfdisi(ds_contact%sdcont_defi, 'ALGO_RESO_GEOM')
    nt_cont_poin = nint(v_sdcont_tabfin(1))
!
! - Get dynamic parameters
!
    dyna_form = 0
    if (l_dyna) then
        dyna_form = 1
    end if
!
! - Access to input field
!
    chmlcf_celd = chmlcf//'.CELD'
    chmlcf_celv = chmlcf//'.CELV'
    call jeveuo(chmlcf_celd, 'L', vi=v_chmlcf_celd)
    call jeveuo(chmlcf_celv, 'E', jv_chmlcf_celv)
    nb_grel = v_chmlcf_celd(2)
!
! - Fill input field
!
    nt_liel = 0
    do i_grel = 1, nb_grel
        decal = v_chmlcf_celd(nceld1+i_grel)
        nb_liel = v_chmlcf_celd(decal+1)
        ASSERT(v_chmlcf_celd(decal+3) .eq. ncmp)
        call jeveuo(jexnum(ligrcf//'.LIEL', i_grel), 'L', vi=v_ligrcf_liel)
        do i_liel = 1, nb_liel
            i_cont_poin = -v_ligrcf_liel(i_liel)
            i_zone = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+14))
            coef_fric = mminfr(ds_contact%sdcont_defi, 'COEF_COULOMB', i_zone)
            i_algo_cont = mminfi(ds_contact%sdcont_defi, 'ALGO_CONT', i_zone)
            i_algo_fric = mminfi(ds_contact%sdcont_defi, 'ALGO_FROT', i_zone)
! --------- Adress in CHAM_ELEM
            vale_indx = jv_chmlcf_celv-1+v_chmlcf_celd(decal+nceld2+nceld3*(i_liel-1)+4)
! --------- Set values in CHAM_ELEM
            zr(vale_indx-1+1) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+4)
            zr(vale_indx-1+2) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+5)
            zr(vale_indx-1+3) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+6)
            zr(vale_indx-1+4) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+7)
            zr(vale_indx-1+5) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+8)
            zr(vale_indx-1+6) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+9)
            zr(vale_indx-1+7) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+10)
            zr(vale_indx-1+8) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+11)
            zr(vale_indx-1+9) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+12)
            zr(vale_indx-1+10) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+13)
            zr(vale_indx-1+11) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+15)
            zr(vale_indx-1+12) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+23)
            zr(vale_indx-1+13) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+17)
            zr(vale_indx-1+14) = v_sdcont_jsupco(i_cont_poin)
            zr(vale_indx-1+15) = i_algo_cont
            zr(vale_indx-1+16) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2)
!            A la premiere iteration on ne passe pas par mmalgo
!            On prend directement la valeur de coef*_cont venant de nmprma
            if (nint(ds_contact%update_init_coefficient) .eq. 1 .and. &
                ds_contact%iteration_newton .le. 1) then
                zr(vale_indx-1+16) = max(ds_contact%estimated_coefficient, &
                                         zr(vale_indx-1+16))
                if (i_algo_cont .ne. 3) zr(vale_indx-1+16) = zr(vale_indx-1+16)*1.d-6
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = zr(vale_indx-1+16)
            end if
            zr(vale_indx-1+17) = i_reso_fric
            zr(vale_indx-1+25) = i_reso_geom
            zr(vale_indx-1+18) = i_algo_fric
            zr(vale_indx-1+19) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6)
            if ((i_algo_cont .ne. 3) .and. (i_algo_fric .eq. 3) .and. &
                (ds_contact%iteration_newton .eq. 1)) then
                glis_maxi = mminfr(ds_contact%sdcont_defi, 'GLIS_MAXI', i_zone)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6) = &
                    1.d-3*ds_contact%estimated_coefficient*ds_contact%arete_min/glis_maxi
            end if
            zr(vale_indx-1+20) = coef_fric
            zr(vale_indx-1+21) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+20)
            zr(vale_indx-1+22) = dyna_form
            zr(vale_indx-1+23) = time_incr
            zr(vale_indx-1+24) = 0.d0
!           Previous iteration state
            ! previous pressure
            zr(vale_indx-1+26) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+3)
            ! previous contact status
            zr(vale_indx-1+27) = nint(v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+1))
            ! alpha_cont_matr
            zr(vale_indx-1+28) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+59)
            ! alpha_frot_matr
            zr(vale_indx-1+42) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+55)
            ! alpha_frot_vect
            zr(vale_indx-1+43) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+54)
            ! previous gap
            zr(vale_indx-1+29) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+4)
            ! treatment of cycling or not cycling
            !    contact cycling
            zr(vale_indx-1+30) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+57)
            !    glis_av-glis_ar cycling
            zr(vale_indx-1+44) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+50)
            ! alpha_cont_vect
            zr(vale_indx-1+31) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+56)
            ! Previous tangentials
            zr(vale_indx-1+32) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+13)
            zr(vale_indx-1+33) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+14)
            zr(vale_indx-1+34) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+15)
            zr(vale_indx-1+35) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+16)
            zr(vale_indx-1+36) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+17)
            zr(vale_indx-1+37) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+18)
            ! Previous contact points coordinates and projections
            zr(vale_indx-1+38) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+19)
            zr(vale_indx-1+39) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+20)
            zr(vale_indx-1+40) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+21)
            zr(vale_indx-1+41) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+22)
            !mode robuste contact
            zr(vale_indx-1+45) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+51)
            !mode robuste frottement
            zr(vale_indx-1+46) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+52)
            !mode adaptatif : frottement penalise
            if ((i_algo_cont .eq. 3) .and. (i_algo_fric .eq. 1)) then
                zr(vale_indx-1+47) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+5)
            end if
            ! and ALGO_RESO_POINT_FIXE="POINT_FIXE"
            !if (ds_contact%iteration_newton .le. 0) then
            !    zr(vale_indx-1+48) = 0
            !else
            zr(vale_indx-1+48) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+73)
            !endif
            !wpg old
            zr(vale_indx-1+49) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+75)
        end do
        nt_liel = nt_liel+nb_liel
    end do
    ASSERT(nt_liel .eq. nt_cont_poin)
!
    call jedema()
end subroutine
