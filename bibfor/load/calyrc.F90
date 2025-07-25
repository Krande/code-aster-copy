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
subroutine calyrc(load, mesh, model, geomDime)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/calirg.h"
#include "asterfort/canort.h"
#include "asterfort/char_read_tran.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprel.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nbnlma.h"
#include "asterfort/orilma.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load LIAISON_CYCL
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  model            : model
! In  geomDime         : space dimension
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: k, kk, nuno1, nuno2, ino1, ino2, nocc, iocc
    integer(kind=8) :: nnomx, idmax
    integer(kind=8) :: nbma1, nbma2
    integer(kind=8) :: nbno2, idcal1, idcal2, nul
    integer(kind=8) :: nno11, nno12, i, indire, lno
    integer(kind=8) :: nbtyp, nddl2, jlistk, jdim, ndim1
    integer(kind=8) :: idim, ij, norien, ntrait
    integer(kind=8) :: icoef1, icoef2, icoef3, iagno3, nbno3, nbma3
    aster_logical :: lrota, dnor, lreori
    real(kind=8) :: beta, coef1, mrota(3, 3), zero, normal(3)
    real(kind=8) :: r8b
    real(kind=8) :: coef11, coef12, coef3
    complex(kind=8) :: betac
    character(len=4) :: valeType
    character(len=4) :: typcoe, typlia
    character(len=8) :: kbeta, nono1, nono2, cmp, ddl2, listyp(11)
    character(len=16) :: motfac, cores1, cores2, tymocl(4), motcle(4)
    character(len=19) :: lisrel
    character(len=24) :: geom3
    character(len=24) :: valk(2)
    aster_logical :: l_tran
    real(kind=8) :: tran(3)
    aster_logical :: l_cent
    real(kind=8) :: cent(3)
    aster_logical :: l_angl_naut
    real(kind=8) :: angl_naut(3)
    character(len=24) :: list_node
    real(kind=8), pointer :: coef(:) => null()
    complex(kind=8), pointer :: coemuc(:) => null()
    integer(kind=8), pointer :: dim(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: nomddl(:) => null()
    character(len=8), pointer :: nomnoe(:) => null()
    integer(kind=8), pointer :: conu1(:) => null()
    integer(kind=8), pointer :: conu2(:) => null()
    integer(kind=8), pointer :: limanu3(:) => null()
    integer(kind=8), pointer :: ln(:) => null()
    integer(kind=8), pointer :: limanu2(:) => null()
    real(kind=8), pointer :: cocf1(:) => null()
    real(kind=8), pointer :: cocf2(:) => null()
    real(kind=8), pointer :: normale(:) => null()
    integer(kind=8), pointer :: conb1(:) => null()
    integer(kind=8), pointer :: conb2(:) => null()
    integer(kind=8), pointer :: limanu1(:) => null()
    aster_logical :: l_error
    character(len=8) :: elem_error
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    nul = 0
    motfac = 'LIAISON_CYCL'
    call getfac(motfac, nocc)
    if (nocc .eq. 0) goto 310
!
    typlia = 'DEPL'
!
!
    valeType = 'REEL'
    typcoe = 'REEL'
!
    lisrel = '&&CALYRC.RLLISTE'
    zero = 0.0d0
    beta = 0.0d0
    betac = (0.0d0, 0.0d0)
    kbeta = ' '
    ndim1 = 3
    lreori = .false.
!
    if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
        call utmess('F', 'CHARGES2_6')
    end if
!
    if (geomDime .eq. 2) then
        nbtyp = 3
        listyp(1) = 'SEG2'
        listyp(2) = 'SEG3'
        listyp(3) = 'SEG4'
    else if (geomDime .eq. 3) then
        nbtyp = 11
        listyp(1) = 'TRIA3'
        listyp(2) = 'TRIA6'
        listyp(3) = 'TRIA9'
        listyp(4) = 'QUAD4'
        listyp(5) = 'QUAD8'
        listyp(6) = 'QUAD9'
        listyp(7) = 'QUAD12'
        listyp(8) = 'SEG2'
        listyp(9) = 'SEG3'
        listyp(10) = 'SEG4'
        listyp(11) = 'TRIA7'
    end if
!
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nnomx)
    idmax = nnomx+3
    AS_ALLOCATE(vk8=nomnoe, size=idmax)
    AS_ALLOCATE(vk8=nomddl, size=idmax)
    AS_ALLOCATE(vr=coef, size=idmax)
    AS_ALLOCATE(vc=coemuc, size=1)
    AS_ALLOCATE(vr=direct, size=idmax*3)
    AS_ALLOCATE(vi=dim, size=idmax)
!
    cores1 = '&&CALYRC.CORES1'
    cores2 = '&&CALYRC.CORES2'
!
    do iocc = 1, nocc
!
        dnor = .false.
        if (typlia .eq. 'DEPL') then
            call getvtx(motfac, 'DDL_ESCL', iocc=iocc, scal=ddl2, nbret=nddl2)
            if (nddl2 .gt. 0) dnor = .true.
        end if
!
!        1.1 RECUPERATION DE LA LISTE DES MAILLE_MAIT :
!        ----------------------------------------------
!        -- 1er groupe maitre --
        motcle(1) = 'MAILLE_MAIT1'
        tymocl(1) = 'MAILLE'
        motcle(2) = 'GROUP_MA_MAIT1'
        tymocl(2) = 'GROUP_MA'
        call reliem(model, mesh, 'NU_MAILLE', motfac, iocc, &
                    2, motcle, tymocl, '&&CALYRC.LIMANU1', nbma1)
        call jeveuo('&&CALYRC.LIMANU1', 'L', vi=limanu1)
!        -- 2eme groupe maitre --
        motcle(1) = 'MAILLE_MAIT2'
        tymocl(1) = 'MAILLE'
        motcle(2) = 'GROUP_MA_MAIT2'
        tymocl(2) = 'GROUP_MA'
        call reliem(model, mesh, 'NU_MAILLE', motfac, iocc, &
                    2, motcle, tymocl, '&&CALYRC.LIMANU2', nbma2)
        if (nbma2 .gt. 0) then
            call jeveuo('&&CALYRC.LIMANU2', 'L', vi=limanu2)
        end if
!
!        1.2 RECUPERATION DES NOEUD_ESCL
!        -------------------------------
        if (.not. dnor) then
!
!        -- RECUPERATION DE LA LISTE DES NOEUD_ESCL :
!        --------------------------------------------
            motcle(1) = 'NOEUD_ESCL'
            tymocl(1) = 'NOEUD'
            motcle(2) = 'GROUP_NO_ESCL'
            tymocl(2) = 'GROUP_NO'
            motcle(3) = 'MAILLE_ESCL'
            tymocl(3) = 'MAILLE'
            motcle(4) = 'GROUP_MA_ESCL'
            tymocl(4) = 'GROUP_MA'
            call reliem(' ', mesh, 'NU_NOEUD', motfac, iocc, &
                        4, motcle, tymocl, '&&CALYRC.LINONU', nbno3)
            call jeveuo('&&CALYRC.LINONU', 'L', iagno3)
!
        else
!
!        -- RECUPERATION DE LA LISTE DES MAILLE_ESCL :
!        ---------------------------------------------
            motcle(1) = 'MAILLE_ESCL'
            tymocl(1) = 'MAILLE'
            motcle(2) = 'GROUP_MA_ESCL'
            tymocl(2) = 'GROUP_MA'
            call reliem(model, mesh, 'NU_MAILLE', motfac, iocc, &
                        2, motcle, tymocl, '&&CALYRC.LIMANU3', nbma3)
            if (nbma3 .eq. 0) then
                valk(1) = motcle(1)
                valk(2) = motcle(2)
                call utmess('F', 'MODELISA8_49', nk=2, valk=valk)
            end if
            call jeveuo('&&CALYRC.LIMANU3', 'L', vi=limanu3)
!
            norien = 0
            call orilma(mesh, geomDime, limanu3, nbma3, norien, &
                        ntrait, lreori, 0, [0])
            if (norien .ne. 0) then
                call utmess('F', 'MODELISA3_19')
            end if
!
! ---        CREATION DU TABLEAU DES NUMEROS DES NOEUDS '&&NBNLMA.LN'
! ---        ET DES NOMBRES D'OCCURENCES DE CES NOEUDS '&&NBNLMA.NBN'
! ---        DES MAILLES DE PEAU MAILLE_ESCL :
!            -------------------------------
            call nbnlma(mesh, nbma3, limanu3, nbtyp, listyp, &
                        nbno3, l_error, elem_error)
            if (l_error) then
                call utmess('F', 'CHARGES6_5', sk=elem_error)
            end if
!
! ---        CALCUL DES NORMALES EN CHAQUE NOEUD :
!            -----------------------------------
            call wkvect('&&CALYRC.LISTK', 'V V K8', 1, jlistk)
            call jeveuo('&&NBNLMA.LN', 'L', vi=ln)
!
! ---        CREATION DU TABLEAU D'INDIRECTION ENTRE LES INDICES
! ---        DU TABLEAU DES NORMALES ET LES NUMEROS DES NOEUDS :
!            -------------------------------------------------
            call wkvect('&&CALYRC.INDIRE', 'V V I', nnomx, indire)
            call jelira('&&NBNLMA.LN', 'LONUTI', lno)
!
            do i = 1, lno
                zi(indire+ln(i)-1) = i
            end do
!
            call canort(mesh, nbma3, limanu3, geomDime, nbno3, &
                        ln, 1)
            call jeveuo('&&CANORT.NORMALE', 'L', vr=normale)
            call jedupo('&&NBNLMA.LN3', 'V', '&&CALYRC.LINONU', .false._1)
            call jeveuo('&&CALYRC.LINONU', 'L', iagno3)
        end if
!
!
!       1.3 TRANSFORMATION DE LA GEOMETRIE DE GRNO2 :
!       ------------------------------------------
        call char_read_tran(motfac, iocc, geomDime, l_tran, tran, &
                            l_cent, cent, l_angl_naut, angl_naut)
!
        geom3 = '&&CALYRC.GEOM_TRANSF'
        list_node = '&&CALYRC.LINONU'
        call calirg(mesh, nbno3, list_node, tran, cent, &
                    l_angl_naut, angl_naut, geom3, lrota, mrota)
!
!
!       2. CALCUL DE CORRES :
!       -------------------
        if (geomDime .eq. 2) then
!        -- 1er groupe esclave / 1er groupe maitre --
            call pj2dco('PARTIE', model, model, nbma1, limanu1, &
                        nbno3, zi(iagno3), ' ', geom3, cores1, &
                        .false._1, r8b, 0.d0)
            if (nbma2 .gt. 0) then
!        -- 1er groupe esclave  / 2eme groupe maitre --
                call pj2dco('PARTIE', model, model, nbma2, limanu2, &
                            nbno3, zi(iagno3), ' ', geom3, cores2, &
                            .false._1, r8b, 0.d0)
            end if
        else if (geomDime .eq. 3) then
!        -- 1er groupe esclave / 1er groupe maitre --
            call pj3dco('PARTIE', model, model, nbma1, limanu1, &
                        nbno3, zi(iagno3), ' ', geom3, cores1, &
                        .false._1, r8b, 0.d0)
            if (nbma2 .gt. 0) then
!        -- 1er groupe esclave  / 2eme groupe maitre --
                call pj3dco('PARTIE', model, model, nbma2, limanu2, &
                            nbno3, zi(iagno3), ' ', geom3, cores2, &
                            .false._1, r8b, 0.d0)
            end if
        end if
!
!        -- 1er groupe maitre --
        call jeveuo(cores1//'.PJEF_NB', 'L', vi=conb1)
        call jeveuo(cores1//'.PJEF_NU', 'L', vi=conu1)
        call jeveuo(cores1//'.PJEF_CF', 'L', vr=cocf1)
        call jelira(cores1//'.PJEF_NB', 'LONMAX', nbno2)
!
        if (nbma2 .gt. 0) then
!        -- 2eme groupe maitre --
            call jeveuo(cores2//'.PJEF_NB', 'L', vi=conb2)
            call jeveuo(cores2//'.PJEF_NU', 'L', vi=conu2)
            call jeveuo(cores2//'.PJEF_CF', 'L', vr=cocf2)
!         CALL JELIRA(CORES2//'.PJEF_NB','LONMAX',NBNO2,KB)
        end if
!
!
!       3. ECRITURE DES RELATIONS LINEAIRES :
!       =====================================
!
!       3.0 RECUPERATION D'UN FACTEUR :
!       =================
!        -- 1er groupe maitre --
        call getvr8(motfac, 'COEF_MAIT1', iocc=iocc, scal=coef11, nbret=icoef1)
        if (icoef1 .le. 0) then
            coef11 = 1.d0
        end if
        if (nbma2 .gt. 0) then
!        -- 2eme groupe maitre --
            call getvr8(motfac, 'COEF_MAIT2', iocc=iocc, scal=coef12, nbret=icoef2)
            if (icoef2 .le. 0) then
                coef12 = 1.d0
            end if
        end if
!        -- 1er groupe esclave --
        call getvr8(motfac, 'COEF_ESCL', iocc=iocc, scal=coef3, nbret=icoef3)
        if (icoef3 .le. 0) then
            coef3 = 1.d0
        end if
!
!       3.1 CAS "DEPL" :
!       =================
        if (typlia .eq. 'DEPL') then
!
!       -- 3.1.1 S'IL N'Y A PAS DE ROTATION :
!       -------------------------------------
            if (.not. lrota) then
                idcal1 = 0
                idcal2 = 0
                do ino2 = 1, nbno3
!           NNO11: NB DE NOEUD_MAIT LIES A INO2 SELON CORES1
                    nno11 = conb1(ino2)
!           NNO12: NB DE NOEUD_MAIT LIES A INO2 SELON CORES2
                    if (nbma2 .gt. 0) then
                        nno12 = conb2(ino2)
                    else
                        nno12 = 0
                    end if
                    if ((nno11 .eq. 0) .and. (nno12 .eq. 0)) goto 90
!
                    nuno2 = ino2
                    nono2 = int_to_char8(nuno2)
!
                    nomnoe(1) = nono2
                    coef(1) = -1.d0*coef3
!
                    do ino1 = 1, nno11
                        nuno1 = conu1(1+idcal1-1+ino1)
                        coef1 = cocf1(1+idcal1-1+ino1)
                        nono1 = int_to_char8(nuno1)
                        nomnoe(ino1+1) = nono1
                        coef(ino1+1) = coef1*coef11
                    end do
                    do ino1 = 1, nno12
                        nuno1 = conu2(1+idcal2-1+ino1)
                        coef1 = cocf2(1+idcal2-1+ino1)
                        nono1 = int_to_char8(nuno1)
                        nomnoe(1+nno11+ino1) = nono1
                        coef(1+nno11+ino1) = coef1*coef12
                    end do
!
!           -- AFFECTATION DES RELATIONS CONCERNANT LE NOEUD INO2 :
!           -----------------------------------------------------
                    if (dnor) then
                        do ino1 = 1, nno11+nno12+1
                            dim(ino1) = geomDime
                            nomddl(ino1) = 'DEPL'
                            do idim = 1, geomDime
                                direct(1+(ino1-1)*ndim1+idim-1) = &
                                    normale(1+(zi(indire+ino2-1)-1)*geomDime+ &
                                            idim-1)
                            end do
                        end do
                        call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                                    direct, nno11+nno12+1, beta, betac, kbeta, &
                                    typcoe, valeType, 1.d-6, lisrel)
                    else
                        do k = 1, geomDime
                            if (k .eq. 1) cmp = 'DX'
                            if (k .eq. 2) cmp = 'DY'
                            if (k .eq. 3) cmp = 'DZ'
                            do ino1 = 1, nno11+nno12+1
                                nomddl(ino1) = cmp
                            end do
                            call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                                        direct, nno11+nno12+1, beta, betac, kbeta, &
                                        typcoe, valeType, 1.d-6, lisrel)
                            call imprel(motfac, nno11+nno12+1, coef, nomddl, nomnoe, &
                                        beta)
                        end do
                    end if
                    idcal1 = idcal1+nno11
                    idcal2 = idcal2+nno12
90                  continue
                end do
!
!       -- 3.1.2  S'IL Y A UNE ROTATION :
!       ---------------------------------
            else
                idcal1 = 0
                idcal2 = 0
!
                do ino2 = 1, nbno2
!
! ---       NNO1: NB DE NOEUD_MAIT LIES A INO2 :
!           ------------------------------------
                    nno11 = conb1(ino2)
                    if (nbma2 .gt. 0) then
                        nno12 = conb2(ino2)
                    else
                        nno12 = 0
                    end if
                    if ((nno11 .eq. 0) .and. (nno12 .eq. 0)) goto 250
                    do k = 1, idmax
                        nomnoe(k) = ' '
                        nomddl(k) = ' '
                        coef(k) = zero
                        dim(k) = 0
                        do kk = 1, 3
                            direct(1+3*(k-1)+kk-1) = zero
                        end do
                    end do
!
                    normal(1) = zero
                    normal(2) = zero
                    normal(3) = zero
!
                    nuno2 = ino2
                    nono2 = int_to_char8(nuno2)
!
                    if (dnor) then
                        ij = 1
                    else
                        ij = geomDime
                    end if
!
                    do ino1 = 1, nno11
                        nuno1 = conu1(1+idcal1-1+ino1)
                        coef1 = cocf1(1+idcal1-1+ino1)
                        nono1 = int_to_char8(nuno1)
                        nomnoe(1+ij+ino1-1) = nono1
                        coef(1+ij+ino1-1) = coef1*coef11
                    end do
                    do ino1 = 1, nno12
                        nuno1 = conu2(1+idcal2-1+ino1)
                        coef1 = cocf2(1+idcal2-1+ino1)
                        nono1 = int_to_char8(nuno1)
                        nomnoe(1+nno11+ij+ino1-1) = nono1
                        coef(1+nno11+ij+ino1-1) = coef1*coef12
                    end do
!
!
!           -- AFFECTATION DES RELATIONS CONCERNANT LE NOEUD INO2 :
!           -----------------------------------------------------
                    if (dnor) then
                        do idim = 1, geomDime
                            do jdim = 1, geomDime
                                normal(idim) = normal(idim)+mrota(jdim, idim)*normale(1+(zi(in&
                                               &dire+ino2-1)-1)*geomDime+jdim-1)
                            end do
                        end do
                        coef(1) = 1.0d0*coef3
                        nomnoe(1) = nono2
                        nomddl(1) = 'DEPL'
                        dim(1) = geomDime
                        do idim = 1, geomDime
                            direct(idim) = normale(1+(zi(indire+ino2-1)-1)*geomDime+idim-1)
                        end do
                        do ino1 = 2, nno11+1
                            dim(ino1) = geomDime
                            nomddl(ino1) = 'DEPL'
                            do idim = 1, geomDime
                                direct(1+(ino1-1)*ndim1+idim-1) = -normal(idim)
                            end do
                        end do
                        do ino1 = 2, nno12+1
                            dim(1+nno11+ino1-1) = geomDime
                            nomddl(nno11+ino1) = 'DEPL'
                            do idim = 1, geomDime
                                direct(1+nno11+(ino1-1)*ndim1+idim- &
                                       1) = -normal(idim)
                            end do
                        end do
                        call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                                    direct, nno11+nno12+1, beta, betac, kbeta, &
                                    typcoe, valeType, 1.d-6, lisrel)
                    else
                        do k = 1, geomDime
                            if (k .eq. 1) cmp = 'DX'
                            if (k .eq. 2) cmp = 'DY'
                            if (k .eq. 3) cmp = 'DZ'
                            do ino1 = 1, nno11
                                nomddl(1+geomDime+ino1-1) = cmp
                            end do
                            do ino1 = 1, nno12
                                nomddl(1+nno11+geomDime+ino1-1) = cmp
                            end do
                            do kk = 1, geomDime
                                if (kk .eq. 1) cmp = 'DX'
                                if (kk .eq. 2) cmp = 'DY'
                                if (kk .eq. 3) cmp = 'DZ'
                                nomnoe(kk) = nono2
                                nomddl(kk) = cmp
                                coef(kk) = -mrota(kk, k)*coef3
                            end do
                            call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                                        direct, nno11+nno12+geomDime, beta, betac, kbeta, &
                                        typcoe, valeType, 1.d-6, lisrel)
                            call imprel(motfac, nno11+nno12+geomDime, coef, nomddl, nomnoe, &
                                        beta)
                        end do
                    end if
                    idcal1 = idcal1+nno11
                    idcal2 = idcal2+nno12
250                 continue
                end do
            end if
!
!
!       3.2 CAS "TEMP" :
!       =================
        else if (typlia .eq. 'TEMP') then
            idcal1 = 0
            idcal2 = 0
            do ino2 = 1, nbno2
!           NNO1: NB DE NOEUD_MAIT LIES A INO2
                nno11 = conb1(ino2)
                nno12 = conb2(ino2)
                if ((nno11 .eq. 0) .and. (nno12 .eq. 0)) goto 290
!
                nuno2 = ino2
                nono2 = int_to_char8(nuno2)
!
                nomnoe(1) = nono2
                coef(1) = -1.d0*coef3
!
                do ino1 = 1, nno11
                    nuno1 = conu1(1+idcal1-1+ino1)
                    coef1 = cocf1(1+idcal1-1+ino1)
                    nono1 = int_to_char8(nuno1)
                    nomnoe(ino1+1) = nono1
                    coef(ino1+1) = coef1*coef11
                end do
                do ino1 = 1, nno12
                    nuno1 = conu2(1+idcal2-1+ino1)
                    coef1 = cocf2(1+idcal2-1+ino1)
                    nono1 = int_to_char8(nuno1)
                    nomnoe(1+nno11+ino1) = nono1
                    coef(1+nno11+ino1) = coef1*coef12
                end do
!
!           -- AFFECTATION DE LA RELATION CONCERNANT LE NOEUD INO2 :
!           -----------------------------------------------------
                cmp = 'TEMP'
                do ino1 = 1, nno11+nno12+1
                    nomddl(ino1) = cmp
                end do
                call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                            direct, nno11+nno12+1, beta, betac, kbeta, &
                            typcoe, valeType, 1.d-6, lisrel)
                call imprel(motfac, nno11+nno12+1, coef, nomddl, nomnoe, &
                            beta)
                idcal1 = idcal1+nno11
                idcal2 = idcal2+nno12
290             continue
            end do
        else
            ASSERT(.false.)
        end if
!
        call detrsd('CORRESP_2_MAILLA', cores1)
        call detrsd('CORRESP_2_MAILLA', cores2)
        call jedetr(geom3)
        call jedetr('&&CALYRC.LIMANU1')
        call jedetr('&&CALYRC.LIMANU2')
        call jedetr('&&CALYRC.LIMANU3')
        call jedetr('&&CALYRC.LINONU')
        call jedetr('&&CALYRC.LISTK')
        call jedetr('&&CALYRC.INDIRE')
        call jedetr('&&NBNLMA.LN')
        call jedetr('&&NBNLMA.NBN')
        call jedetr('&&CANORT.NORMALE')
!
    end do
!
    AS_DEALLOCATE(vc=coemuc)
    AS_DEALLOCATE(vk8=nomnoe)
    AS_DEALLOCATE(vk8=nomddl)
    AS_DEALLOCATE(vr=coef)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dim)
!
! --- AFFECTATION DE LA LISTE DE RELATIONS A LA load :
!     ------------------------------------------------
    call aflrch(lisrel, load, 'LIN')
!
310 continue
    call jedema()
end subroutine
