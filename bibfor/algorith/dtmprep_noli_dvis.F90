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

subroutine dtmprep_noli_dvis(sd_dtm_, sd_nl_, icomp)
    implicit none
! dtmprep_noli_dvis : prepare the calculations for a localized nonlinearity
!                     of type : generalized Zener damper. This routine adds a
!                     single occurence to sd_nl and increments NB_NOLI in sd_dtm
!
!             icomp : an integer giving the index of occurence of the
!                     nonlinearity to be treated under the factor kw
!                     COMPORTEMENT of the command DYNA_VIBRA.
!
!
!                     MODELE D'AMORTISSEUR DE ZENZER GENERALISE
!
!                                         e2
!                         e1     |-----======-----|
!                     ---=====---|                |-----
!                                |--=====----=]---|
!                                     e3    n3,a3
!
!
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/angvx.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/gloloc.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdchdl.h"
#include "asterfort/nlget.h"
#include "asterfort/nlinivec.h"
#include "asterfort/nlsav.h"
#include "asterfort/nltype.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), intent(in) :: icomp
!
!   -0.2- Local variables
    aster_logical     :: lnoeu2
    integer(kind=8)           :: i, n1, ibid, nbdisv, nbnoli
    integer(kind=8)           :: nbmcl, ier, nbno1, nbno2, ino1
    integer(kind=8)           :: ino2, ind1, ind2, nbmode, info
    integer(kind=8)           :: vali, j, neq, mxlevel, nexcit
!
    real(kind=8)      :: dv_c, dv_a, dv_k1, dv_k2, dv_k3
    real(kind=8)      :: r8bid, alpha, beta, sina, cosa
    real(kind=8)      :: sinb, cosb, sing, cosg, valr(10)
    real(kind=8)      :: ddpilo(3), dpiglo(6), dpiloc(6), one, axe(3)
    real(kind=8)      :: res_inte
!
    character(len=8)  :: sd_dtm, sd_nl, mesh, mesh1, mesh2
    character(len=8)  :: nume, nume1, nume2, no1_name, no2_name
    character(len=8)  :: monmot, intk
    character(len=16) :: typnum, typem
    character(len=16) :: limocl(2), tymocl(2), valk(2), obst_typ, motfac
    character(len=19) :: nomres
    character(len=24) :: nl_title
!
    integer(kind=8), pointer       :: ddlcho(:) => null()
    real(kind=8), pointer       :: coor_no1(:) => null()
    real(kind=8), pointer       :: coor_no2(:) => null()
    real(kind=8), pointer       :: vale(:) => null()
    real(kind=8), pointer       :: sincos_angle_a(:) => null()
    real(kind=8), pointer       :: sincos_angle_b(:) => null()
    real(kind=8), pointer       :: sincos_angle_g(:) => null()
    real(kind=8), pointer       :: defmod1(:) => null()
    real(kind=8), pointer       :: defmod2(:) => null()
    real(kind=8), pointer       :: ps2del1(:) => null()
    real(kind=8), pointer       :: ps2del2(:) => null()
    real(kind=8), pointer       :: origob(:) => null()
    real(kind=8), pointer       :: bmodal_v(:) => null()
    real(kind=8), pointer       :: ps1del_v(:) => null()
!
    character(len=8), pointer  :: noeud(:) => null()
!
#define ps1del(m,n) ps1del_v((n-1)*neq+m)
#define bmodal(m,n) bmodal_v((n-1)*neq+m)
!
!   --- 0. Diverse initializations
    call jemarq()
!
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    lnoeu2 = .false.
    one = 1.d0
    !
    motfac = 'COMPORTEMENT'
    call nlget(sd_nl, _MAX_LEVEL, iscal=mxlevel)
    nbnoli = mxlevel+1
    i = mxlevel+1
!
    call infmaj()
    call infniv(ibid, info)
!
!   --- 1 - Basic information about the mesh and numbering
!
    call dtmget(sd_dtm, _NUM_DDL, kscal=nume)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    call gettco(nume, typnum)

!
!   --- 1.1 - Case with a simple modal projection (direct calculation)
    if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
        call dismoi('NOM_MAILLA', nume, 'NUME_DDL', repk=mesh)
        mesh1 = mesh
        nume1 = nume
        mesh2 = mesh
        nume2 = nume
        call nlsav(sd_nl, _NUMDDL_1, 1, iocc=i, kscal=nume1(1:8))
        call nlsav(sd_nl, _MESH_1, 1, iocc=i, kscal=mesh1)

!   --- 1.2 - Case with double (or triple) projections (sub-structuring case)
!             Not supported for buckling non linearities
    else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
        call utmess('F', 'ALGORITH5_36')
    else
        ASSERT(.false.)
    end if
!
!
!   --- 2 - Localisation (support nodes) of the buckling non linearity
!
!   --- 2.1 - Definition using nodes or nodal groups (NOEUD/GROUP_NO)
!             Unlike the chocs case, here only a single nonlinearity
!             can be defined per occurence
    typem = 'NO_NOEUD'
    nbmcl = 2
    limocl(1) = 'GROUP_NO_1'
    limocl(2) = 'NOEUD_1'
    tymocl(1) = 'GROUP_NO'
    tymocl(2) = 'NOEUD'
    call reliem(' ', mesh1, typem, motfac, icomp, &
                nbmcl, limocl, tymocl, sd_nl//'.INDI_NO1.TEMP', nbno1)
!
    if (nbno1 .gt. 0) then
        ASSERT(nbno1 .eq. 1)
        call jeveuo(sd_nl//'.INDI_NO1.TEMP', 'L', vk8=noeud)
        no1_name = noeud(1)
        call nlsav(sd_nl, _NO1_NAME, 1, iocc=i, kscal=no1_name)
        call jedetr(sd_nl//'.INDI_NO1.TEMP')

        typem = 'NO_NOEUD'
        nbmcl = 2
        limocl(1) = 'GROUP_NO_2'
        limocl(2) = 'NOEUD_2'
        call reliem(' ', mesh2, typem, motfac, icomp, &
                    nbmcl, limocl, tymocl, sd_nl//'.INDI_NO2.TEMP', nbno2)

        if (nbno2 .gt. 0) then
            ASSERT(nbno2 .eq. 1)
            call jeveuo(sd_nl//'.INDI_NO2.TEMP', 'L', vk8=noeud)
            no2_name = noeud(1)
            call nlsav(sd_nl, _NO2_NAME, 1, iocc=i, kscal=no2_name)
            call jedetr(sd_nl//'.INDI_NO2.TEMP')
            lnoeu2 = .true.
        else
            ASSERT(.false.)
        end if

        call nlsav(sd_nl, _NUMDDL_2, 1, iocc=i, kscal=nume2(1:8))
        call nlsav(sd_nl, _MESH_2, 1, iocc=i, kscal=mesh2)
    end if
!
!   --- 3 - Filling up the sd_nl with further information regarding the
!           nonlinearity(ies)
!
    AS_ALLOCATE(vi=ddlcho, size=6)
!
!   --- Loop over the detected nonlineary in the current COMPORTEMENT
!       occurence
!
    call nlsav(sd_nl, _NL_TYPE, 1, iocc=i, iscal=NL_DIS_VISC)
!
!   --- 3.1 - DOF numbering localisation index for the concerned nodes
    call mdchdl(lnoeu2, i, ddlcho, ier)
!
!   --- 3.2 - Coordinates of the nodes
    call jeveuo(mesh1//'.COORDO    .VALE', 'L', vr=vale)
    ino1 = char8_to_int(no1_name)
    ind1 = 1+3*(ino1-1)
    ind2 = ind1+3
    call nlsav(sd_nl, _COOR_NO1, 3, iocc=i, rvect=vale(ind1:ind2))
    if (lnoeu2) then
        if (mesh2 .ne. mesh1) then
            call jeveuo(mesh2//'.COORDO    .VALE', 'L', vr=vale)
        end if
        ino2 = char8_to_int(no2_name)
        ind1 = 1+3*(ino2-1)
        ind2 = ind1+3
        call nlsav(sd_nl, _COOR_NO2, 3, iocc=i, rvect=vale(ind1:ind2))
    end if
!
!   --- 3.3 - Other information are read from the user input
    call codent(i, 'D0', intk)
    nl_title = nltype(NL_DIS_VISC)//intk
    call nlsav(sd_nl, _NL_TITLE, 1, iocc=i, kscal=nl_title)

    call getvr8(motfac, 'K1', iocc=icomp, scal=r8bid, nbret=n1)
    if (n1 .eq. 1) then
        call nlsav(sd_nl, _DISVISC_K1, 1, iocc=i, rscal=1.0d0/r8bid)
    else
        call getvr8(motfac, 'UNSUR_K1', iocc=icomp, scal=r8bid, nbret=n1)
        ASSERT(n1 .eq. 1)
        call nlsav(sd_nl, _DISVISC_K1, 1, iocc=i, rscal=r8bid)
    end if

    call getvr8(motfac, 'K2', iocc=icomp, scal=r8bid, nbret=n1)
    if (n1 .eq. 1) then
        call nlsav(sd_nl, _DISVISC_K2, 1, iocc=i, rscal=r8bid)
    else
        call getvr8(motfac, 'UNSUR_K2', iocc=icomp, scal=r8bid, nbret=n1)
        ASSERT(n1 .eq. 1)
        call nlsav(sd_nl, _DISVISC_K2, 1, iocc=i, rscal=1.0d0/r8bid)
    end if

    call getvr8(motfac, 'K3', iocc=icomp, scal=r8bid, nbret=n1)
    if (n1 .eq. 1) then
        call nlsav(sd_nl, _DISVISC_K3, 1, iocc=i, rscal=1.0d0/r8bid)
    else
        call getvr8(motfac, 'UNSUR_K3', iocc=icomp, scal=r8bid, nbret=n1)
        ASSERT(n1 .eq. 1)
        call nlsav(sd_nl, _DISVISC_K3, 1, iocc=i, rscal=r8bid)
    end if

    call getvr8(motfac, 'C', iocc=icomp, scal=dv_c, nbret=n1)
    call nlsav(sd_nl, _DISVISC_C, 1, iocc=i, rscal=dv_c)

    call getvr8(motfac, 'PUIS_ALPHA', iocc=icomp, scal=dv_a, nbret=n1)
    call nlsav(sd_nl, _DISVISC_A, 1, iocc=i, rscal=dv_a)

!   --- Calculation of : 1/K1 + K2/(K1*K3) + 1/K3
    call nlget(sd_nl, _DISVISC_K1, iocc=i, rscal=dv_k1)
    call nlget(sd_nl, _DISVISC_K2, iocc=i, rscal=dv_k2)
    call nlget(sd_nl, _DISVISC_K3, iocc=i, rscal=dv_k3)
    r8bid = (dv_k1+dv_k2*dv_k1*dv_k3+dv_k3)
    if (r8bid .le. r8miem()) then
        call utmess('F', 'DISCRETS_41')
    end if
    call getvis(motfac, 'ITER_INTE_MAXI', iocc=icomp, scal=vali, nbret=n1)
    call nlsav(sd_nl, _MAX_INTE, 1, iocc=i, iscal=vali)

    call getvr8(motfac, 'RESI_INTE', iocc=icomp, scal=res_inte, nbret=n1)
    call nlsav(sd_nl, _RES_INTE, 1, iocc=i, rscal=res_inte)
!
!
!   --- 3.4 - Calculation of geometrical properties :
!             play, orientation, local coordinates, distances
!             Vector x1x2 must not be zero
    call nlget(sd_nl, _COOR_NO1, iocc=i, vr=coor_no1)
    call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
    axe(1) = (coor_no2(1)-coor_no1(1))
    axe(2) = (coor_no2(2)-coor_no1(2))
    axe(3) = (coor_no2(3)-coor_no1(3))

    r8bid = axe(1)**2+axe(2)**2+axe(3)**2
    if (r8bid .le. r8miem()) then
        call utmess('F', 'DISCRETS_43')
    end if

    call angvx(axe, alpha, beta)
    call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=i, rvect=[sin(alpha), cos(alpha)])
    call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=i, rvect=[sin(beta), cos(beta)])
    call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=i, rvect=[0.0d0, 1.0d0])

    obst_typ = 'BI_PLANY'
    call nlsav(sd_nl, _OBST_TYP, 1, iocc=i, kscal=obst_typ)
    call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=i, &
               rvect=[0.5d0*(coor_no1(1)+coor_no2(1)), &
                      0.5d0*(coor_no1(2)+coor_no2(2)), &
                      0.5d0*(coor_no1(3)+coor_no2(3))])

    call nlget(sd_nl, _COOR_ORIGIN_OBSTACLE, iocc=i, vr=origob)
    call nlget(sd_nl, _SINCOS_ANGLE_A, iocc=i, vr=sincos_angle_a)
    call nlget(sd_nl, _SINCOS_ANGLE_B, iocc=i, vr=sincos_angle_b)
    call nlget(sd_nl, _SINCOS_ANGLE_G, iocc=i, vr=sincos_angle_g)

    sina = sincos_angle_a(1)
    cosa = sincos_angle_a(2)
    sinb = sincos_angle_b(1)
    cosb = sincos_angle_b(2)
    sing = sincos_angle_g(1)
    cosg = sincos_angle_g(2)
!
!   --- Global to local (obstacle) coordinates
    dpiglo(1) = coor_no1(1)
    dpiglo(2) = coor_no1(2)
    dpiglo(3) = coor_no1(3)
    call gloloc(dpiglo, origob, sina, cosa, sinb, &
                cosb, sing, cosg, dpiloc)
!   --- Differential distance given a single node
    ddpilo(1) = dpiloc(1)
    ddpilo(2) = dpiloc(2)
    ddpilo(3) = dpiloc(3)
!
    if (obst_typ(1:2) .eq. 'BI') then
!       --- Initial position of node 2 in the local (obstacle) reference
        call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
        dpiglo(4) = coor_no2(1)
        dpiglo(5) = coor_no2(2)
        dpiglo(6) = coor_no2(3)
        call gloloc(dpiglo(4), origob, sina, cosa, sinb, &
                    cosb, sing, cosg, dpiloc(4))
!   --- Differential coordinates (distance) for the binodal system
        ddpilo(1) = dpiloc(1)-dpiloc(4)
        ddpilo(2) = dpiloc(2)-dpiloc(5)
        ddpilo(3) = dpiloc(3)-dpiloc(6)
    end if
    call nlsav(sd_nl, _SIGN_DYZ, 2, iocc=i, rvect=[-sign(one, ddpilo(2)), -sign(one, ddpilo(3))])
!
!
!   --- 3.5 - Printing out user information
    if (info .eq. 2) then
        vali = i
        call nlget(sd_nl, _NO1_NAME, iocc=i, kscal=valk(1))
        call utmess('I', 'ALGORITH16_2', sk=valk(1), si=vali)

        call nlget(sd_nl, _COOR_NO1, iocc=i, rvect=valr)
        call utmess('I', 'ALGORITH16_4', nr=3, valr=valr)

        call nlget(sd_nl, _NO2_NAME, iocc=i, kscal=valk(1))
        call utmess('I', 'ALGORITH16_5', sk=valk(1))
        if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
            call nlget(sd_nl, _SS2_NAME, iocc=i, kscal=valk(1))
            call utmess('I', 'ALGORITH16_3', sk=valk(1))
        end if
        call nlget(sd_nl, _COOR_NO2, iocc=i, rvect=valr)
        call utmess('I', 'ALGORITH16_4', nr=3, valr=valr)

        valr(1) = 0.d0
        valr(2) = origob(1)
        valr(3) = origob(2)
        valr(4) = origob(3)
        valr(5) = sincos_angle_a(1)
        valr(6) = sincos_angle_a(2)
        valr(7) = sincos_angle_b(1)
        valr(8) = sincos_angle_b(2)
        valr(9) = sincos_angle_g(1)
        valr(10) = sincos_angle_g(2)
        call utmess('I', 'ALGORITH16_8', nr=10, valr=valr)

        valr(1) = 0.d0
        call utmess('I', 'ALGORITH16_9', sr=valr(1))

        call utmess('I', 'VIDE_1')
    end if
!
!   -- 3.6 - Modal displacements of the node(s)
!            Note : if a single node is used, we fill with zeros the
!                   deformations for node_2, this simplifies the
!                   case treatments for calculating the forces
    call dtmget(sd_dtm, _BASE_VEC, vr=bmodal_v)
    call dtmget(sd_dtm, _NB_PHYEQ, iscal=neq)
    call nlinivec(sd_nl, _MODAL_DEPL_NO1, 3*nbmode, iocc=i, vr=defmod1)
    call nlinivec(sd_nl, _MODAL_DEPL_NO2, 3*nbmode, iocc=i, vr=defmod2)

    do j = 1, nbmode
        defmod1(3*(j-1)+1) = bmodal(ddlcho(1), j)
        defmod1(3*(j-1)+2) = bmodal(ddlcho(2), j)
        defmod1(3*(j-1)+3) = bmodal(ddlcho(3), j)

        if (obst_typ(1:2) .eq. 'BI') then
            defmod2(3*(j-1)+1) = bmodal(ddlcho(4), j)
            defmod2(3*(j-1)+2) = bmodal(ddlcho(5), j)
            defmod2(3*(j-1)+3) = bmodal(ddlcho(6), j)
        else
            defmod2(3*(j-1)+1) = 0.d0
            defmod2(3*(j-1)+2) = 0.d0
            defmod2(3*(j-1)+3) = 0.d0
        end if
    end do
!
!   --- 3.7 - Multi supported systems, filling up of psixdelta for the
!             concerned nodes
    call dtmget(sd_dtm, _MULTI_AP, kscal=monmot)
    if (monmot(1:3) .eq. 'OUI') then
        call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
        call dtmget(sd_dtm, _NB_EXC_T, iscal=nexcit)
        call jeveuo(nomres//'.IPSD', 'E', vr=ps1del_v)

        call nlinivec(sd_nl, _PSI_DELT_NO1, 3*nexcit, iocc=i, vr=ps2del1)
        do j = 1, nexcit
            ps2del1(3*(j-1)+1) = ps1del(ddlcho(1), j)
            ps2del1(3*(j-1)+2) = ps1del(ddlcho(2), j)
            ps2del1(3*(j-1)+3) = ps1del(ddlcho(3), j)
        end do
        if (obst_typ(1:2) .eq. 'BI') then
            call nlinivec(sd_nl, _PSI_DELT_NO2, 3*nexcit, iocc=i, vr=ps2del2)
            do j = 1, nexcit
                ps2del2(3*(j-1)+1) = ps1del(ddlcho(4), j)
                ps2del2(3*(j-1)+2) = ps1del(ddlcho(5), j)
                ps2del2(3*(j-1)+3) = ps1del(ddlcho(6), j)
            end do
        end if
    end if
!
!   --- 4 - Updating indices for sd_nl and sd_dtm
    call nlsav(sd_nl, _MAX_LEVEL, 1, iscal=nbnoli)
    call dtmsav(sd_dtm, _NB_NONLI, 1, iscal=nbnoli)
!
    call nlget(sd_nl, _NB_DIS_VISC, iscal=nbdisv)
    nbdisv = nbdisv+1
    call nlsav(sd_nl, _NB_DIS_VISC, 1, iscal=nbdisv)
!
    AS_DEALLOCATE(vi=ddlcho)
!
    call jedema()
end subroutine
