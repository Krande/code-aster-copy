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

subroutine dtmprep_noli_flam(sd_dtm_, sd_nl_, icomp)
    implicit none
! dtmprep_noli_flam : prepare the calculations for a localized nonlinearity
!                     of type : buckling. This routine adds a single
!                     occurence to sd_nl and increments NB_NOLI in sd_dtm
!
!             icomp : an integer giving the index of occurence of the
!                     nonlinearity to be treated under the factor kw
!                     COMPORTEMENT of the command DYNA_VIBRA.
!
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/gloloc.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdchan.h"
#include "asterfort/mdchdl.h"
#include "asterfort/mdchre.h"
#include "asterfort/nlget.h"
#include "asterfort/nlinivec.h"
#include "nlinc.h"
#include "asterfort/nlsav.h"
#include "asterfort/nltype.h"
#include "asterfort/reliem.h"
#include "asterfort/tbliva.h"
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
    integer(kind=8)           :: i, n1, ibid, nbbuck, nbnoli
    integer(kind=8)           :: iret, nbmcl, ier, nbno1, nbno2
    integer(kind=8)           :: ino1, ino2, ind1, ind2, nbmode
    integer(kind=8)           :: info, vali, j, neq, mxlevel
    integer(kind=8)           :: nbchoc, nexcit, nl_type, tomove, k
    integer(kind=8)           :: cntr, amorin
!
    real(kind=8)      :: r8bid, gap, xjeu, sina, cosa
    real(kind=8)      :: sinb, cosb, sing, cosg, valr(10)
    real(kind=8)      :: dist_no1, dist_no2, ddpilo(3), dpiglo(6)
    real(kind=8)      :: dpiloc(6), one
!
    complex(kind=8)   :: cbid
!
    character(len=8)  :: sd_dtm, sd_nl, mesh, mesh1, mesh2
    character(len=8)  :: nume, nume1, nume2, no1_name, no2_name
    character(len=8)  :: monmot, k8typ, kbid, repere, node
    character(len=8)  :: intk
    character(len=8)  :: fx, rigi_nor, amor_nor
    character(len=16) :: typnum, typem, refo, limocl(2), tymocl(2)
    character(len=16) :: valk(2), obst_typ, motfac
    character(len=19) :: nomres
    character(len=24) :: nl_title, mdgene, jvname, jvname0
    character(len=3)  :: typamor
!
    integer(kind=8), pointer  :: ddlcho(:) => null()
    real(kind=8), pointer  :: coor_no1(:) => null()
    real(kind=8), pointer  :: coor_no2(:) => null()
    real(kind=8), pointer  :: vale(:) => null()
    real(kind=8), pointer  :: sincos_angle_a(:) => null()
    real(kind=8), pointer  :: sincos_angle_b(:) => null()
    real(kind=8), pointer  :: sincos_angle_g(:) => null()
    real(kind=8), pointer  :: defmod1(:) => null()
    real(kind=8), pointer  :: defmod2(:) => null()
    real(kind=8), pointer  :: ps2del1(:) => null()
    real(kind=8), pointer  :: ps2del2(:) => null()
    real(kind=8), pointer  :: origob(:) => null()
    real(kind=8), pointer  :: bmodal_v(:) => null()
    real(kind=8), pointer  :: ps1del_v(:) => null()

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
    dist_no1 = 0.d0
    dist_no2 = 0.d0
!
    lnoeu2 = .false.
    one = 1.d0
    !
    motfac = 'COMPORTEMENT'
    call nlget(sd_nl, _MAX_LEVEL, iscal=mxlevel)
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
!             Unlike the previous case, here only a single nonlinearity
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
        end if

        call nlsav(sd_nl, _NUMDDL_2, 1, iocc=i, kscal=nume2(1:8))
        call nlsav(sd_nl, _MESH_2, 1, iocc=i, kscal=mesh2)
    end if

!   --- 2.2 - Check whether a/several previous stop/choc nonlinearity(ies)
!             have been defined for the same support node, if yes, remove
!             these stop/choc nonlinearity(ies)
    call nlget(sd_nl, _NB_CHOC, iscal=nbchoc)
    nbnoli = mxlevel+1
    cntr = 0
    if (nbchoc .gt. 0) then
        do j = 1, mxlevel
            call nlget(sd_nl, _NL_TYPE, iocc=j-cntr, iscal=nl_type)
            if (nl_type .eq. NL_CHOC) then
                call nlget(sd_nl, _NO1_NAME, iocc=j-cntr, kscal=node)
                if (node .eq. no1_name) then
                    call utmess('A', 'ALGORITH5_38')
                    nbchoc = nbchoc-1
                    nbnoli = nbnoli-1
                    call nlget(sd_nl, 1, iocc=j-cntr, savejv=jvname)
                    call detrsd(' ', jvname(1:15))
                    if (i .eq. (mxlevel+1)) then
                        tomove = i
                        i = j-cntr
                    else
!                       --- In this case, several chocs are detected and
!                           are removed, this causes empty entries in the
!                           sd_nl that need to be filled up
                        tomove = mxlevel-cntr
                    end if

                    if (tomove .gt. (j-cntr)) then
                        do k = 1, _NL_NBPAR
                            if (parind(k) .gt. 0) then
                                call nlget(sd_nl, k, iocc=tomove, savejv=jvname)
                                call jeexin(jvname, iret)
                                if (iret .ne. 0) then
                                    call nlget(sd_nl, k, iocc=j-cntr, savejv=jvname0)
                                    call jedupo(jvname, 'V', jvname0, .false._1)
                                end if
                            end if
                        end do
                        call detrsd(' ', jvname(1:15))
                        if (i .ne. (j-cntr)) cntr = cntr+1
                    end if
                    call dtmsav(sd_dtm, _NB_NONLI, 1, iscal=nbnoli)
                    call nlsav(sd_nl, _NB_CHOC, 1, iscal=nbchoc)
                end if
            end if
        end do
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
    call nlsav(sd_nl, _NL_TYPE, 1, iocc=i, iscal=NL_BUCKLING)
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
    call getvtx(motfac, 'INTITULE', iocc=icomp, scal=nl_title, nbret=n1)
    if (n1 .ne. 0) then
        call nlsav(sd_nl, _NL_TITLE, 1, iocc=i, kscal=nl_title)
    else
        call codent(i, 'D0', intk)
        nl_title = nltype(NL_BUCKLING)//intk
        call nlsav(sd_nl, _NL_TITLE, 1, iocc=i, kscal=nl_title)
    end if

    call nlsav(sd_nl, _GAP, 1, iocc=i, rscal=0.d0)
    call getvr8(motfac, 'JEU', iocc=icomp, scal=gap, nbret=n1)
    if (n1 .gt. 0) call nlsav(sd_nl, _GAP, 1, iocc=i, rscal=gap)

    call nlsav(sd_nl, _DIST_NO1, 1, iocc=i, rscal=0.d0)
    call getvr8(motfac, 'DIST_1', iocc=icomp, scal=dist_no1, nbret=n1)
    if (n1 .gt. 0) call nlsav(sd_nl, _DIST_NO1, 1, iocc=i, rscal=dist_no1)

    call nlsav(sd_nl, _DIST_NO2, 1, iocc=i, rscal=0.d0)
    call getvr8(motfac, 'DIST_2', iocc=icomp, scal=dist_no2, nbret=n1)
    if (n1 .gt. 0) call nlsav(sd_nl, _DIST_NO2, 1, iocc=i, rscal=dist_no2)

    call getvtx(motfac, 'CRIT_AMOR', iocc=icomp, scal=typamor, nbret=n1)
    if (typamor(1:3) .eq. 'INC') then
        amorin = 1
        call nlsav(sd_nl, _BUCKLING_AMOR_IN, 1, iocc=i, iscal=amorin)
    else if (typamor(1:3) .eq. 'EXC') then
        amorin = 0
        call nlsav(sd_nl, _BUCKLING_AMOR_IN, 1, iocc=i, iscal=amorin)
    end if

    call getvid(motfac, 'FX', iocc=icomp, scal=fx)
    call getvid(motfac, 'RIGI_NOR', iocc=icomp, scal=rigi_nor)
    call getvid(motfac, 'AMOR_NOR', iocc=icomp, scal=amor_nor)

    call nlsav(sd_nl, _BUCKLING_F_CRIT, 1, iocc=i, kscal=fx)
    call nlsav(sd_nl, _BUCKLING_F_RIGI, 1, iocc=i, kscal=rigi_nor)
    call nlsav(sd_nl, _BUCKLING_F_AMOR, 1, iocc=i, kscal=amor_nor)

    call getvid(motfac, 'OBSTACLE', iocc=icomp, scal=obst_typ, nbret=n1)

!   --- 3.4 - Obstacle type
    call tbliva(obst_typ, 1, 'LIEU', [ibid], [r8bid], &
                [cbid], 'DEFIOBST', kbid, [r8bid], 'TYPE', &
                k8typ, ibid, r8bid, cbid, refo, &
                ier)
    ASSERT(ier .eq. 0)
    if (refo(1:9) .eq. 'BI_PLAN_Y') then
        obst_typ = 'BI_PLANY'
    else if (refo(1:9) .eq. 'BI_PLAN_Z') then
        obst_typ = 'BI_PLANZ'
    else if (refo(1:11) .eq. 'BI_CERC_INT') then
        obst_typ = 'BI_CERCI'
    else if (refo(1:7) .ne. 'DISCRET') then
        obst_typ = refo(1:8)
    end if
    call nlsav(sd_nl, _OBST_TYP, 1, iocc=i, kscal=obst_typ)

    if ((obst_typ(1:8) .eq. 'BI_CERCI') .and. (dist_no2 .lt. dist_no1)) then
        call utmess('F', 'ALGORITH5_35')
    end if

!
!   --- 3.4 - Calculation of geometrical properties :
!             play, orientation, local coordinates, distances
    call nlget(sd_nl, _COOR_NO1, iocc=i, vr=coor_no1)
    xjeu = 0.d0
    if (obst_typ(1:2) .eq. 'BI') then
        call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
        xjeu = (coor_no2(1)-coor_no1(1))**2+ &
               (coor_no2(2)-coor_no1(2))**2+ &
               (coor_no2(3)-coor_no1(3))**2
    end if
!
    call mdchre('FLAMBAGE  ', icomp, i, mdgene, typnum, &
                repere, lnoeu2)
!
    call mdchan('FLAMBAGE  ', icomp, i, mdgene, typnum, &
                repere, xjeu)

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

        if (obst_typ(1:2) .eq. 'BI') then
            call nlget(sd_nl, _NO2_NAME, iocc=i, kscal=valk(1))
            call utmess('I', 'ALGORITH16_5', sk=valk(1))
            if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
                call nlget(sd_nl, _SS2_NAME, iocc=i, kscal=valk(1))
                call utmess('I', 'ALGORITH16_3', sk=valk(1))
            end if
            call nlget(sd_nl, _COOR_NO2, iocc=i, rvect=valr)
            call utmess('I', 'ALGORITH16_4', nr=3, valr=valr)
        end if
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
        if (obst_typ(1:2) .eq. 'BI') then
            xjeu = sqrt(xjeu)-(dist_no1+dist_no2)
            valr(1) = xjeu
            call utmess('I', 'ALGORITH16_9', sr=valr(1))
        end if
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
    call nlget(sd_nl, _NB_FLAMB, iscal=nbbuck)
    nbbuck = nbbuck+1
    call nlsav(sd_nl, _NB_FLAMB, 1, iscal=nbbuck)
!
    AS_DEALLOCATE(vi=ddlcho)

!
    call jedema()
end subroutine
