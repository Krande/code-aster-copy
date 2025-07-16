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

subroutine dtmprep_noli_choc(sd_dtm_, sd_nl_, icomp)
    implicit none
! dtmprep_noli_choc : prepare the calculations for a localized nonlinearity
!                     of type : stop/choc. This routine adds one or more
!                     occurences to sd_nl and increments NB_NOLI in sd_dtm
!
!             icomp : an integer giving the index of occurence of the
!                     nonlinearity to be treated under the factor kw
!                     COMPORTEMENT of the command DYNA_VIBRA.
!
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
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
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdchan.h"
#include "asterfort/mdchdl.h"
#include "asterfort/mdchre.h"
#include "asterfort/mgutdm.h"
#include "asterfort/nlget.h"
#include "asterfort/nlinivec.h"
#include "asterfort/nlsav.h"
#include "asterfort/nltype.h"
#include "asterfort/orient.h"
#include "asterfort/reliem.h"
#include "asterfort/resmod.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/locglo.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), intent(in) :: icomp
!
!   -0.2- Local variables
    aster_logical     :: lnoeu2, friction
    integer(kind=8)           :: i, n1, ibid, ind, nbbuck
    integer(kind=8)           :: jcoor, iret, nbmcl, nbma, im
    integer(kind=8)           :: jmama, nbnma, ier, nbno1, nbno2
    integer(kind=8)           :: ino1, ino2, ind1, ind2, irett
    integer(kind=8)           :: namtan, nbmode, ind_mmax, info, vali(10)
    integer(kind=8)           :: j, neq, start, finish, mxlevel, nbchoc
    integer(kind=8)           :: nl_type, nexcit, unidir
!
    real(kind=8)      :: ctang, ktang, fric_static, fric_dynamic, r8bid
    real(kind=8)      :: gap, xjeu, k, amor, mmax, tdepl2, gdepl(3), ldepl(3)
    real(kind=8)      :: sina, cosa, sinb, cosb, sing
    real(kind=8)      :: cosg, valr(10), damp_normal, stif_normal, dist_no1
    real(kind=8)      :: dist_no2, ddpilo(3), dpiglo(6), dpiloc(6), one
    real(kind=8)      :: coor(3), VectIN(3), VectOUT(3)
    real(kind=8)      :: ldepl_norm

!
    character(len=3)  :: unidirk
    character(len=8)  :: sd_dtm, sd_nl, mesh, mesh1, mesh2
    character(len=8)  :: nume, sst1, sst2, nomma, no1_name
    character(len=8)  :: no2_name, monmot, k8typ, kbid, repere, nomfon
    character(len=8)  :: node, intk
    character(len=14) :: nume1, nume2
    character(len=16) :: typnum, typem, refo, limocl(2), tymocl(2)
    character(len=16) :: valk(2), obst_typ, motfac
    character(len=19) :: nomres
    character(len=24) :: typfro, nl_title, numero, mdgene, jvname, mdssno
!
    complex(kind=8)   :: cbid

    integer(kind=8), pointer       :: ddlcho(:) => null()
    integer(kind=8), pointer       :: elems(:) => null()
    real(kind=8), pointer       :: coor_no1(:) => null()
    real(kind=8), pointer       :: coor_no2(:) => null()
    real(kind=8), pointer       :: bc_norm(:) => null()
    real(kind=8), pointer       :: vale(:) => null()
    real(kind=8), pointer       :: omega2(:) => null()
    real(kind=8), pointer       :: amogen(:) => null()
    real(kind=8), pointer       :: masgen(:) => null()
    real(kind=8), pointer       :: sincos_angle_a(:) => null()
    real(kind=8), pointer       :: sincos_angle_b(:) => null()
    real(kind=8), pointer       :: sincos_angle_g(:) => null()
    real(kind=8), pointer       :: defmod1(:) => null()
    real(kind=8), pointer       :: defmod2(:) => null()
    real(kind=8), pointer       :: dplcho(:) => null()
    real(kind=8), pointer       :: ps2del1(:) => null()
    real(kind=8), pointer       :: ps2del2(:) => null()
    real(kind=8), pointer       :: origob(:) => null()
    real(kind=8), pointer       :: bmodal_v(:) => null()
    real(kind=8), pointer       :: ps1del_v(:) => null()
    real(kind=8), pointer       :: mass_normed(:) => null()
    character(len=8), pointer  :: noeud(:) => null()
    character(len=24), pointer  :: refn(:) => null()
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
    start = i
    finish = i
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
    else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
        call jeveuo(nume//'      .NUME.REFN', 'L', vk24=refn)
        mdgene = refn(1)

        call getvtx(motfac, 'SOUS_STRUC_1', iocc=icomp, scal=sst1, nbret=n1)
        if (n1 .eq. 0) call utmess('F', 'ALGORITH5_31', sk='SOUS_STRUC_1')

        mdssno = mdgene(1:14)//'.MODG.SSNO'
        call jenonu(jexnom(mdssno, sst1), iret)
        if (iret .eq. 0) call utmess('F', 'ALGORITH5_32')

        call mgutdm(mdgene, sst1, ibid, 'NOM_NUME_DDL', ibid, nume1)
        call mgutdm(mdgene, sst1, ibid, 'NOM_MAILLAGE', ibid, mesh1)
!
        call nlsav(sd_nl, _SS1_NAME, 1, iocc=i, kscal=sst1)
        call nlsav(sd_nl, _NUMDDL_1, 1, iocc=i, kscal=nume1(1:8))
        call nlsav(sd_nl, _MESH_1, 1, iocc=i, kscal=mesh1)

        mesh2 = mesh1
        nume2 = nume1
        sst2 = sst1
        call getvtx(motfac, 'SOUS_STRUC_2', iocc=icomp, scal=sst2, nbret=n1)
        if (n1 .ne. 0) then
            call jenonu(jexnom(mdssno, sst2), iret)
            if (iret .eq. 0) call utmess('F', 'ALGORITH5_32', sk='SOUS_STRUC_2')

            call mgutdm(mdgene, sst2, ibid, 'NOM_NUME_DDL', ibid, nume2)
            call mgutdm(mdgene, sst2, ibid, 'NOM_MAILLAGE', ibid, mesh2)
            !
        end if
        call nlsav(sd_nl, _SS2_NAME, 1, iocc=i, kscal=sst2)
        call nlsav(sd_nl, _NUMDDL_2, 1, iocc=i, kscal=nume2(1:8))
        call nlsav(sd_nl, _MESH_2, 1, iocc=i, kscal=mesh2)
    else
        ASSERT(.false.)
    end if
!
!
!   --- 2 - Localisation (support nodes) of the stop/choc non linearity
!
!   --- 2.1 - Definition using elements or element groups (MAILLE/GROUP_MA)
!             Caution : GROUP_MA may contain several segments and thus
!                       several non linearities
    typem = 'NU_MAILLE'
    nbmcl = 2
    limocl(1) = 'GROUP_MA'
    limocl(2) = 'MAILLE'
    call reliem(' ', mesh1, typem, motfac, icomp, &
                nbmcl, limocl, limocl, sd_nl//'.INDI_SUPP.TEMP', nbma)
!
    if (nbma .gt. 0) then
        ASSERT(mesh2 .eq. mesh1)
        call jeveuo(sd_nl//'.INDI_SUPP.TEMP', 'L', vi=elems)
        ind = i
        do im = 1, nbma
            call jeveuo(jexnum(mesh1//'.CONNEX', elems(im)), 'L', jmama)
            call jelira(jexnum(mesh1//'.CONNEX', elems(im)), 'LONMAX', nbnma)
            if (nbnma .ne. 2) then
                nomma = int_to_char8(elems(im))
                valk(1) = nomma
                valk(2) = 'SEG2'
                call utmess('F', 'ALGORITH13_39', nk=2, valk=valk)
            end if
            no1_name = int_to_char8(zi(jmama))
            no2_name = int_to_char8(zi(jmama+1))
            call nlsav(sd_nl, _NO1_NAME, 1, iocc=ind, kscal=no1_name)
            call nlsav(sd_nl, _NO2_NAME, 1, iocc=ind, kscal=no2_name)

            call nlsav(sd_nl, _NUMDDL_1, 1, iocc=ind, kscal=nume1(1:8))
            call nlsav(sd_nl, _MESH_1, 1, iocc=ind, kscal=mesh1)

            call nlsav(sd_nl, _NUMDDL_2, 1, iocc=ind, kscal=nume1(1:8))
            call nlsav(sd_nl, _MESH_2, 1, iocc=ind, kscal=mesh1)
!
            ind = ind+1
        end do
        call jedetr(sd_nl//'.INDI_SUPP.TEMP')
        finish = ind-1
        lnoeu2 = .true.
    end if
!
!   --- 2.2 - Definition using nodes or nodal groups (NOEUD/GROUP_NO)
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

!   --- 2.3 - Check whether this stop/choc nonlinearity(ies) conflicts
!             with a buckling one for the same support node, if yes, stop
!             the treatment immediately
    call nlget(sd_nl, _NB_FLAMB, iscal=nbbuck)
    if (nbbuck .gt. 0) then
        do j = 1, mxlevel
            call nlget(sd_nl, _NL_TYPE, iocc=j, iscal=nl_type)
            if (nl_type .eq. NL_BUCKLING) then
                call nlget(sd_nl, _NO1_NAME, iocc=j, kscal=node)
                if (node .eq. no1_name) then
                    call utmess('A', 'ALGORITH5_38')
                    call nlget(sd_nl, 1, iocc=i, savejv=jvname)
                    call detrsd(' ', jvname(1:15))
                    goto 999
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
!   --- Loop over the detected nonlinearitiesi in the current COMPORTEMENT
!       occurence
    do i = start, finish
!
        call nlsav(sd_nl, _NL_TYPE, 1, iocc=i, iscal=NL_CHOC)
!
!       --- 3.1 - DOF numbering localisation index for the concerned nodes
        call mdchdl(lnoeu2, i, ddlcho, ier)
!
!       --- 3.2 - Coordinates of the nodes
        if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
            call jeveuo(mesh1//'.COORDO    .VALE', 'L', vr=vale)
            call nlget(sd_nl, _NO1_NAME, iocc=i, kscal=no1_name)
            ino1 = char8_to_int(no1_name)
            ind1 = 1+3*(ino1-1)
            ind2 = ind1+2
            call nlsav(sd_nl, _COOR_NO1, 3, iocc=i, rvect=vale(ind1:ind2))
            if (lnoeu2) then
                if (mesh2 .ne. mesh1) then
                    call jeveuo(mesh2//'.COORDO    .VALE', 'L', vr=vale)
                end if
                call nlget(sd_nl, _NO2_NAME, iocc=i, kscal=no2_name)
                ino2 = char8_to_int(no2_name)
                ind1 = 1+3*(ino2-1)
                ind2 = ind1+2
                call nlsav(sd_nl, _COOR_NO2, 3, iocc=i, rvect=vale(ind1:ind2))
            end if
        else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
            call jeveuo(mesh1//'.COORDO    .VALE', 'L', jcoor)
            call nlget(sd_nl, _NO1_NAME, iocc=i, kscal=no1_name)
            ino1 = char8_to_int(no1_name)
            call orient(mdgene, sst1, jcoor, ino1, coor, 1)
            call nlsav(sd_nl, _COOR_NO1, 3, iocc=i, rvect=coor)
            if (lnoeu2) then
                if (mesh2 .ne. mesh1) then
                    call jeveuo(mesh2//'.COORDO    .VALE', 'L', jcoor)
                end if
                call nlget(sd_nl, _NO2_NAME, iocc=i, kscal=no2_name)
                ino2 = char8_to_int(no2_name)
                call orient(mdgene, sst2, jcoor, ino2, coor, 1)
                call nlsav(sd_nl, _COOR_NO2, 3, iocc=i, rvect=coor)
            end if
        end if
!
!       --- 3.3 - Other information are read from the user input
        call getvtx(motfac, 'INTITULE', iocc=icomp, scal=nl_title, nbret=n1)
        if (n1 .ne. 0) then
            call nlsav(sd_nl, _NL_TITLE, 1, iocc=i, kscal=nl_title)
        else
            call codent(i, 'D0', intk)
            nl_title = nltype(NL_CHOC)//intk
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

!       Normal Stiffness or Function
        call nlsav(sd_nl, _STIF_NORMAL, 1, iocc=i, rscal=0.d0)
        call nlsav(sd_nl, _NL_FUNC_NAME, 1, iocc=i, kscal='.')
!       Dans le catalogue c'est l'un ou l'autre
        call getvr8(motfac, 'RIGI_NOR', iocc=icomp, scal=stif_normal, nbret=n1)
        if (n1 .gt. 0) call nlsav(sd_nl, _STIF_NORMAL, 1, iocc=i, rscal=stif_normal)
        call getvid(motfac, 'FX', iocc=icomp, scal=nomfon, nbret=n1)
        if (n1 .gt. 0) call nlsav(sd_nl, _NL_FUNC_NAME, 1, iocc=i, kscal=nomfon)
!
        call nlsav(sd_nl, _DAMP_NORMAL, 1, iocc=i, rscal=0.d0)
        call getvr8(motfac, 'AMOR_NOR', iocc=icomp, scal=damp_normal, nbret=n1)
        if (n1 .gt. 0) call nlsav(sd_nl, _DAMP_NORMAL, 1, iocc=i, rscal=damp_normal)

        call nlsav(sd_nl, _RIGI_TANGENTIAL, 1, iocc=i, rscal=0.d0)
        call getvr8(motfac, 'RIGI_TAN', iocc=icomp, scal=ktang, nbret=n1)
        if (n1 .gt. 0) call nlsav(sd_nl, _RIGI_TANGENTIAL, 1, iocc=i, rscal=ktang)

        friction = .true.
        unidir = 0
        call getvtx(motfac, 'FROTTEMENT', iocc=icomp, scal=typfro, nbret=n1)
        if (typfro(1:10) .eq. 'COULOMB   ') then
            call getvr8(motfac, 'COULOMB', iocc=icomp, scal=fric_dynamic, nbret=n1)
            call getvtx(motfac, 'UNIDIRECTIONNEL', iocc=icomp, scal=unidirk, nbret=n1)
            if (unidirk .eq. 'OUI') unidir = 1
            call nlsav(sd_nl, _FRIC_DYNAMIC, 1, iocc=i, rscal=fric_dynamic)
            call nlsav(sd_nl, _FRIC_STATIC, 1, iocc=i, rscal=fric_dynamic)
            call nlsav(sd_nl, _FRIC_UNIDIR, 1, iocc=i, iscal=unidir)
        else if (typfro(1:16) .eq. 'COULOMB_STAT_DYN') then
            call getvr8(motfac, 'COULOMB_DYNA', iocc=icomp, scal=fric_dynamic, nbret=n1)
            call getvr8(motfac, 'COULOMB_STAT', iocc=icomp, scal=fric_static, nbret=n1)
            call getvtx(motfac, 'UNIDIRECTIONNEL', iocc=icomp, scal=unidirk, nbret=n1)
            if (unidirk .eq. 'OUI') unidir = 1
            call nlsav(sd_nl, _FRIC_DYNAMIC, 1, iocc=i, rscal=fric_dynamic)
            call nlsav(sd_nl, _FRIC_STATIC, 1, iocc=i, rscal=fric_static)
            call nlsav(sd_nl, _FRIC_UNIDIR, 1, iocc=i, iscal=unidir)
        else
            call nlsav(sd_nl, _FRIC_DYNAMIC, 1, iocc=i, rscal=0.d0)
            call nlsav(sd_nl, _FRIC_STATIC, 1, iocc=i, rscal=0.d0)
            friction = .false.
        end if
!
        ctang = 0.d0
        call getvr8(motfac, 'AMOR_TAN', iocc=icomp, nbret=namtan)
        if (namtan .gt. 0) call getvr8(motfac, 'AMOR_TAN', iocc=icomp, scal=ctang)

        call getvid(motfac, 'OBSTACLE', iocc=icomp, scal=obst_typ, nbret=n1)
!       --- 3.4 - Obstacle type
        call tbliva(obst_typ, 1, 'LIEU', [ibid], [r8bid], &
                    [cbid], 'DEFIOBST', kbid, [r8bid], 'TYPE', &
                    k8typ, ibid, r8bid, cbid, refo, &
                    irett)
        ASSERT(irett .eq. 0)
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
!       --- 3.5 - Calculation of geometrical properties :
!                 play, orientation, local coordinates, distances
        call nlget(sd_nl, _COOR_NO1, iocc=i, vr=coor_no1)
        xjeu = 0.d0
        if (obst_typ(1:2) .eq. 'BI') then
            call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
            xjeu = (coor_no2(1)-coor_no1(1))**2+ &
                   (coor_no2(2)-coor_no1(2))**2+ &
                   (coor_no2(3)-coor_no1(3))**2
        end if
!
        call mdchre('DIS_CHOC  ', icomp, i, mdgene, typnum, &
                    repere, lnoeu2)
!
        call mdchan('DIS_CHOC  ', icomp, i, mdgene, typnum, &
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
!       --- Global to local (obstacle) coordinates
        dpiglo(1) = coor_no1(1)
        dpiglo(2) = coor_no1(2)
        dpiglo(3) = coor_no1(3)
        call gloloc(dpiglo, origob, sina, cosa, sinb, &
                    cosb, sing, cosg, dpiloc)
!       --- Differential distance given a single node
        ddpilo(1) = dpiloc(1)
        ddpilo(2) = dpiloc(2)
        ddpilo(3) = dpiloc(3)
!
        if (obst_typ(1:2) .eq. 'BI') then
!           --- Initial position of node 2 in the local (obstacle) reference
            call nlget(sd_nl, _COOR_NO2, iocc=i, vr=coor_no2)
            dpiglo(4) = coor_no2(1)
            dpiglo(5) = coor_no2(2)
            dpiglo(6) = coor_no2(3)
            call gloloc(dpiglo(4), origob, sina, cosa, sinb, &
                        cosb, sing, cosg, dpiloc(4))
!       --- Differential coordinates (distance) for the binodal system
            ddpilo(1) = dpiloc(1)-dpiloc(4)
            ddpilo(2) = dpiloc(2)-dpiloc(5)
            ddpilo(3) = dpiloc(3)-dpiloc(6)
        end if
        call nlsav(sd_nl, _SIGN_DYZ, 2, iocc=i, &
                   rvect=[-sign(one, ddpilo(2)), -sign(one, ddpilo(3))])
!
!
!       --- 3.6 - Printing out user information
        if (info .eq. 2) then
            call nlget(sd_nl, _NO1_NAME, iocc=i, kscal=valk(1))
            call utmess('I', 'ALGORITH16_2', sk=valk(1), si=i)
            if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
                call nlget(sd_nl, _SS1_NAME, iocc=i, kscal=valk(1))
                call utmess('I', 'ALGORITH16_3', sk=valk(1))
            end if
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
            valr(1) = ctang
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

            VectIN(1) = 1.0d0
            VectIN(2) = 0.0d0
            VectIN(3) = 0.0d0
            call locglo(VectIN, sina, cosa, sinb, &
                        cosb, sing, cosg, VectOUT)
            if (unidir .eq. 1) then
                call utmess('I', 'ALGORITH16_97', nr=3, valr=VectOUT)
            end if

            if (obst_typ(1:2) .eq. 'BI') then
                xjeu = sqrt(xjeu)-(dist_no1+dist_no2)
                valr(1) = xjeu
                call utmess('I', 'ALGORITH16_9', sr=valr(1))
            end if
            call utmess('I', 'VIDE_1')
        end if
!
!       -- 3.7 - Modal displacements of the node(s)
!                Note : if a single node is used, we fill with zeros the
!                       deformations for node_2, this simplifies the
!                       case treatments for calculating the forces
        call dtmget(sd_dtm, _BASE_VEC, vr=bmodal_v)
        call dtmget(sd_dtm, _NB_PHYEQ, iscal=neq)
        call nlinivec(sd_nl, _MODAL_DEPL_NO1, 3*nbmode, iocc=i, vr=defmod1)
        if (obst_typ(1:2) .eq. 'BI') then
            call nlinivec(sd_nl, _MODAL_DEPL_NO2, 3*nbmode, iocc=i, vr=defmod2)
        end if

        if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
            do j = 1, nbmode
                defmod1(3*(j-1)+1) = bmodal(ddlcho(1), j)
                defmod1(3*(j-1)+2) = bmodal(ddlcho(2), j)
                defmod1(3*(j-1)+3) = bmodal(ddlcho(3), j)
            end do
            if (obst_typ(1:2) .eq. 'BI') then
                do j = 1, nbmode
                    defmod2(3*(j-1)+1) = bmodal(ddlcho(4), j)
                    defmod2(3*(j-1)+2) = bmodal(ddlcho(5), j)
                    defmod2(3*(j-1)+3) = bmodal(ddlcho(6), j)
                end do
            end if
!
        else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
            numero = nume
            AS_ALLOCATE(vr=dplcho, size=nbmode*6)
            AS_ALLOCATE(vk8=noeud, size=3)
            call nlget(sd_nl, _NO1_NAME, iocc=i, kscal=noeud(1))
            call nlget(sd_nl, _SS1_NAME, iocc=i, kscal=noeud(2))
            call nlget(sd_nl, _NUMDDL_1, iocc=i, kscal=noeud(3))
            call resmod(bmodal_v, nbmode, neq, numero, mdgene, &
                        noeud, dplcho)
            do j = 1, nbmode
                defmod1(3*(j-1)+1) = dplcho(j)
                defmod1(3*(j-1)+2) = dplcho(j+nbmode)
                defmod1(3*(j-1)+3) = dplcho(j+2*nbmode)
            end do

            if (obst_typ(1:2) .eq. 'BI') then
                call nlget(sd_nl, _NO2_NAME, iocc=i, kscal=noeud(1))
                call nlget(sd_nl, _SS2_NAME, iocc=i, kscal=noeud(2))
                call nlget(sd_nl, _NUMDDL_2, iocc=i, kscal=noeud(3))
                call resmod(bmodal_v, nbmode, neq, numero, mdgene, &
                            noeud, dplcho)
                do j = 1, nbmode
                    defmod2(3*(j-1)+1) = dplcho(j)
                    defmod2(3*(j-1)+2) = dplcho(j+nbmode)
                    defmod2(3*(j-1)+3) = dplcho(j+2*nbmode)
                end do
            end if
            AS_DEALLOCATE(vr=dplcho)
            AS_DEALLOCATE(vk8=noeud)
        end if
!
!       --- 3.8 - Multi supported systems, filling up of psixdelta for the
!                 concerned nodes
        call dtmget(sd_dtm, _MULTI_AP, kscal=monmot)
        if (monmot(1:3) .eq. 'OUI') then

            call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
            call dtmget(sd_dtm, _NB_EXC_T, iscal=nexcit)
            call jeveuo(nomres//'.IPSD', 'E', vr=ps1del_v)

            if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
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
            else if (typnum(1:13) .eq. 'NUME_DDL_GENE') then
                call utmess('E', 'ALGORITH5_37')
            end if
        end if

!       --- 3.9 - Special treatment for the tangential damping
!           If no value is given, an optimized value is calculated and
!           based on the maximum modal mass and the tangential stiffness
        if (namtan .eq. 0 .and. ktang .gt. r8prem()) then
            call dtmget(sd_dtm, _MASS_DIA, vr=masgen)
            call dtmget(sd_dtm, _OMEGA_SQ, vr=omega2)
!
!           --- Calculate a normalizing factor based on the effective
!               *tangential* local displacement for each mode related to
!               the present non linearity
            AS_ALLOCATE(vr=mass_normed, size=nbmode)
!
!           --- 1) Determine the local tangential (in plane) displacement
!                  for each mode
            do j = 1, nbmode
!               --- Detect (and skip) non dynamic modes => fictional mass
                if ((omega2(j)*omega2(j)) .lt. 1.d-10) then
                    ! write(*,*) 'Rigid body movement => not considered'
                    mass_normed(j) = 0.d0
                else
                    gdepl(1) = defmod1(3*(j-1)+1)
                    gdepl(2) = defmod1(3*(j-1)+2)
                    gdepl(3) = defmod1(3*(j-1)+3)
                    if (obst_typ(1:2) .eq. 'BI') then
                        gdepl(1) = defmod2(3*(j-1)+1)-gdepl(1)
                        gdepl(2) = defmod2(3*(j-1)+2)-gdepl(2)
                        gdepl(3) = defmod2(3*(j-1)+3)-gdepl(3)
                    end if
                    call gloloc(gdepl, [0.d0, 0.d0, 0.d0], sina, cosa, sinb, &
                                cosb, sing, cosg, ldepl)
                    if (obst_typ .eq. 'BI_PLANY' .or. obst_typ .eq. 'PLAN_Y') then
                        tdepl2 = ldepl(1)*ldepl(1)+ldepl(3)*ldepl(3)
                    else if (obst_typ .eq. 'BI_PLANZ' .or. obst_typ .eq. 'PLAN_Z') then
                        tdepl2 = ldepl(1)*ldepl(1)+ldepl(2)*ldepl(2)
                    else if (obst_typ .eq. 'BI_CERCLE') then
                        bc_norm = coor_no1-coor_no2
                        bc_norm = bc_norm/sqrt(dot_product(bc_norm, bc_norm))
                        ldepl_norm = dot_product(ldepl, bc_norm)
                        tdepl2 = dot_product(ldepl, ldepl)-ldepl_norm**2
                    else
                        tdepl2 = ldepl(1)*ldepl(1)
                    end if

!                   --- Detect (and skip) non dynamic modes => fictional mass
                    mass_normed(j) = 0.d0
!                    if ((tdepl2*tdepl2) .gt. 1.d-10) then
                    mass_normed(j) = masgen(j)*sqrt(tdepl2)
!                    endif
                end if
            end do

!           --- Determine the mode with the maximum mass, relevant to the
!               tangential motion at the contact interface
!               mmax     : value of the max modal mass
!               ind_mmax : index of the mode with the greatest mass
            mmax = mass_normed(1)
            ind_mmax = 1
            if (nbmode .gt. 1) then
                do im = 2, nbmode
                    if (mass_normed(im) .gt. mmax) then
                        mmax = mass_normed(im)
                        ind_mmax = im
                    end if
                end do
            end if
            mmax = masgen(ind_mmax)
            AS_DEALLOCATE(vr=mass_normed)

!           --- Access index of the corresponding modal damping (ksi)
!               (depending whether the damping matrix is diagonal or no)
            call dtmget(sd_dtm, _AMOR_FUL, lonvec=iret)
            if (iret .gt. 0) then
!               --- amogen = C (full matrix)
!               --- Acrit (critical damping) = 2*sqrt(k*m) = 2*omega*m
                call dtmget(sd_dtm, _AMOR_FUL, vr=amogen)
                amor = amogen(ind_mmax+nbmode*(ind_mmax-1))
                amor = amor/(2.d0*sqrt(omega2(ind_mmax))*masgen(ind_mmax))
            else
!               --- amogen = 2*omega*ksi (diagonal matrix)
                call dtmget(sd_dtm, _AMOR_DIA, vr=amogen)
                amor = 0.5*amogen(ind_mmax)/sqrt(omega2(ind_mmax))
            end if
!
!
!           --- The stiffness (k) is calculated based on mmax which is normed.
!               Hence (k) itself is also normed to local choc displacement = 1
            k = omega2(ind_mmax)*mmax
            ctang = 2.d0*sqrt(mmax*(k+ktang))-2.d0*amor*sqrt(k*mmax)

!           --- Insure that the optimized tangential damping is not negative
            ctang = max(0.d0, ctang)
!
            if (ctang .gt. r8prem()) then
                vali(1) = i
                vali(2) = ind_mmax
                valr(1) = ctang
                valr(2) = mmax
                valr(3) = k
                valr(4) = amor
                valr(5) = ktang
                call utmess('I', 'ALGORITH16_10', ni=2, vali=vali, nr=5, &
                            valr=valr)
            end if
        end if
        call nlsav(sd_nl, _DAMP_TANGENTIAL, 1, iocc=i, rscal=ctang)
    end do
!
!   --- 4 - Updating indices for sd_nl and sd_dtm
    call nlget(sd_nl, _MAX_LEVEL, iscal=mxlevel)
    mxlevel = mxlevel+finish-start+1
    call nlsav(sd_nl, _MAX_LEVEL, 1, iscal=mxlevel)
    call dtmsav(sd_dtm, _NB_NONLI, 1, iscal=mxlevel)
!
    call nlget(sd_nl, _NB_CHOC, iscal=nbchoc)
    nbchoc = nbchoc+finish-start+1
    call nlsav(sd_nl, _NB_CHOC, 1, iscal=nbchoc)
!
    AS_DEALLOCATE(vi=ddlcho)
!
999 continue
    call jedema()
end subroutine
