! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine calipj(chargz)
!
! --------------------------------------------------------------------------------------------------
!
!     Création des relations entre les ddls des noeuds d'un maillage esclave et les ddls des
!     noeuds des mailles d'un maillage maître.
!
!     les relations sont fabriquées à partir de la sd_corresp_2_mailla sortant de proj_champ
!
! --------------------------------------------------------------------------------------------------
!
!   chargz : nom utilisateur du résultat de charge
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    implicit none
!
    character(len=*) :: chargz
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/indik8.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprel.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: zeror, coenul
    parameter (zeror=0.0d0, coenul=1.0d-6)
!
    integer :: k, nuno1, ino1, ino2, iret, nocc, iocc
    integer :: ibid, jnoma, nnomx, idmax, igeom
    integer :: nbno2t, nbno1, jconb, jconu, jcocf
    integer :: jcom1, idecal, ima1, nddls, jcok1, jprnm, jtyma0
    integer :: inomg, nbcmpg, icmpg, nbec, ityp0
    integer :: jcoord, kddl, iexc
    aster_logical :: excent, lbid
    real(kind=8)        :: xyzom(3), coeffi
    complex(kind=8)     :: zeroc
    character(len=1)    :: k1bid
    character(len=4)    :: fonree, typcoe
    character(len=8)    :: noma0, noma1, noma2, model, k8bid, relats, kbeta, nomg
    character(len=8)    :: nono1, nono2, charge, cmp, nomtyp, cmpddl
    character(len=16)   :: motfac, nomcmd, corres, k16bid
    character(len=19)   :: ligrmo, lisrel
!
    character(len=30)   :: valmk(4)
!
    aster_logical               :: okddlrot(3)
    character(len=3), parameter :: ddlrot(3)=['DRX','DRY','DRZ']
!
! --------------------------------------------------------------------------------------------------
!
    integer, pointer ::             idimen(:) => null()
    real(kind=8), pointer ::        idirec(:) => null()
    real(kind=8), pointer ::        idcoef(:) => null()
    character(len=8), pointer ::    idnomn(:) => null()
    character(len=8), pointer ::    idnomd(:) => null()
    character(len=8), pointer ::    jddls(:)  => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    motfac='LIAISON_PROJ'
    call getfac(motfac, nocc)
    if (nocc .eq. 0) goto 999
!
    call getres(k1bid, k1bid, nomcmd)
    if (nomcmd .ne. 'AFFE_CHAR_MECA') then
        ASSERT(ASTER_FALSE)
    endif
!
    fonree = 'REEL'
    typcoe = 'REEL'
    charge = chargz
!
    lisrel = '&&CALIPJ.RLLISTE'
    zeroc = (0.0d0,0.0d0)
    kbeta = ' '
!
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', ibid, model)
    ligrmo=model//'.MODELE'
!   maillage associé au modèle
    call jeveuo(ligrmo//'.LGRF', 'L', jnoma)
    noma0=zk8(jnoma)
    call jeveuo(ligrmo//'.PRNM', 'L', jprnm)
!   coordonnées des noeuds
    call jeveuo(noma0//'.COORDO    .VALE', 'L', jcoord)
!   les types de mailles
    call jeveuo(noma0//'.TYPMAIL', 'L', jtyma0)
!
    call dismoi('DIM_GEOM', model, 'MODELE', igeom, k8bid)
    if ( (igeom.ne.2) .and. (igeom.ne.3) ) then
        call utmess('F', 'MODELISA2_75')
    endif
!
    call dismoi('NB_NO_MAILLA', noma0, 'MAILLAGE', nnomx, k1bid)
    ! idmax : nombre max de termes d'une relation linéaire
    !   = au maximum 27 noeuds dans un élément
    !   U1 = UI2 + DRI2 ^ O2M
    idmax = 1 + 3*27
    AS_ALLOCATE(size=idmax,   vi=idimen)
    AS_ALLOCATE(size=idmax*3, vr=idirec)
    AS_ALLOCATE(size=idmax,   vr=idcoef)
    AS_ALLOCATE(size=idmax,   vk8=idnomn)
    AS_ALLOCATE(size=idmax,   vk8=idnomd)
!
!   Récupération dans les catalogues des composantes de la grandeur
    nomg='DEPL_R'
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', inomg)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmpg)
    call dismoi('NB_EC', nomg, 'GRANDEUR', nbec, k8bid)
!
    do iocc = 1, nocc
        ! il faut remettre à zéro ces 2 objets entre 2 occurrences
        idimen(:) = 0
        idirec(:) = 0.0d0
        ! Récupération des ddls
        call getvtx(motfac, 'DDL', iocc=iocc, nbval=0, nbret=nddls)
        nddls = -nddls
        AS_ALLOCATE(size=nddls,   vk8=jddls)
        call getvtx(motfac, 'DDL', iocc=iocc, nbval=nddls, vect=jddls)
        ! Récupération du concept "corresp_2_mailla"
        call getvid(motfac, 'MATR_PROJECTION', iocc=iocc, scal=relats)
        corres = relats
        ! Récupération du type de relation
        !   TYPE = IDENTITE         par défaut
        !   TYPE = EXCENTREMENT
        call getvtx(motfac, 'TYPE', iocc=iocc, scal=k16bid, nbret=iret)
        excent = ASTER_FALSE
        if (iret .ne. 0) then
            excent = k16bid .eq. 'EXCENTREMENT'
        endif
        !
        call jeveuo(corres//'.PJEF_NB', 'L', jconb)
        call jeveuo(corres//'.PJEF_M1', 'L', jcom1)
        call jeveuo(corres//'.PJEF_NU', 'L', jconu)
        call jeveuo(corres//'.PJEF_CF', 'L', jcocf)
        ! Vérification que les maillages sont tous les mêmes
        call jeveuo(corres//'.PJXX_K1', 'L', jcok1)
        noma1 = zk24(jcok1)(1:8)
        noma2 = zk24(jcok1+1)(1:8)
        if ((noma0.ne.noma1) .or. (noma0.ne.noma2)) then
            valmk(1) = charge
            valmk(2) = noma0
            valmk(3) = noma1
            valmk(4) = noma2
            call utmess('F', 'MODELISA2_74', nk=4, valk=valmk)
        endif
        ! Nombre de noeuds du maillage M2
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbno2t)
        idecal = 0
        ! On boucle sur tous les noeuds du maillage M2
        cino2: do ino2 = 1, nbno2t
            ! IMA1: maille à connecter à ino2
            ima1 = zi(jcom1-1+ino2)
            ! Si ima1=0, c'est que ino2 n'est dans aucune maille
            if (ima1 .eq. 0) cycle cino2
            if (excent) then
                ! Le type de la maille de ma0
                ityp0 = zi(jtyma0-1+ima1)
                call jenuno(jexnum('&CATA.TM.NOMTM', ityp0), nomtyp)
                lbid = (nomtyp(1:5).eq.'QUAD4') .or. (nomtyp(1:5) .eq.'TRIA3')
                if (.not. lbid) then
                    valmk(1) = 'QUAD4 TRIA3'
                    call utmess('F', 'MODELISA2_73', nk=1, valk=valmk)
                endif
            endif
            ! Nbno1: nb noeuds de ima1 dans la relation
            nbno1 = zi(jconb-1+ino2)
            call jenuno(jexnum(noma0//'.NOMNOE', ino2), nono2)
            ! Est-ce que les composantes données par l'utilisateur existent sur le noeud ?
            ibid = jprnm-1+ (ino2-1)*nbec+1
            do kddl = 1, nddls
                cmpddl= jddls(kddl)
                icmpg = indik8(zk8(inomg),cmpddl,1,nbcmpg)
                if (.not.exisdg(zi(ibid),icmpg)) then
                    valmk(1) = nono2
                    valmk(2) = cmpddl
                    call utmess('F', 'MODELISA2_71', nk=2, valk=valmk)
                endif
            enddo
            ! Équations à écrire sur les DDL de nono2
            !   Sans excentrement
            !       U1 = UI2                ==> -UI2 + U1*al(i)
            !   Avec excentrement
            !       U1 = UI2 + DRI2^O2M     ==> -UI2 + U1*al(i) - DRI2^O2M*al(i)
            idnomn(1) = nono2
            idcoef(1) = -1.d0
            do ino1 = 1, nbno1
                nuno1 = zi(jconu+idecal-1+ino1)
                ! Si la relation est une tautologie, on ne l'écrit pas
                if (nuno1 .eq. ino2) then
                    call utmess('A', 'MODELISA2_72', sk=nono1)
                    goto 130
                endif
                !
                call jenuno(jexnum(noma0//'.NOMNOE', nuno1), nono1)
                ! Est-ce que les composantes existent sur le noeud du maillage en face ?
                ! Si le coefficient est <= coenul , la relation sera éliminée. Il n'y a donc pas
                ! besoin de vérifier si les ddl existent.
                coeffi = zr(jcocf+idecal-1+ino1)
                if ( abs(coeffi) .gt. coenul ) then
                    ibid = jprnm-1+ (nuno1-1)*nbec+1
                    do kddl = 1, nddls
                        cmpddl= jddls(kddl)
                        icmpg = indik8(zk8(inomg),cmpddl,1,nbcmpg)
                        if (.not.exisdg(zi(ibid),icmpg)) then
                            valmk(1) = nono1
                            valmk(2) = cmpddl
                            call utmess('F', 'MODELISA2_71', nk=2, valk=valmk)
                        endif
                    enddo
                endif
                idnomn(1+ino1) = nono1
                idcoef(1+ino1) = coeffi
            enddo
            !
            ! Relations concernant les ddls
            do k = 1, nddls
                cmp = jddls(k)
                idnomd(1) = cmp
                do  ino1 = 1, nbno1
                    idnomd(1+ino1) = cmp
                enddo
                iexc = 0
                if (excent) then
                    do ino1 = 1, nbno1
                        nuno1 = zi(jconu+idecal-1+ino1)
                        call jenuno(jexnum(noma0//'.NOMNOE', nuno1), nono1)
                        ! Coefficients intervenant dans la relation
                        coeffi   = zr(jcocf+idecal-1+ino1)
                        xyzom(1) = zr(jcoord+3*(ino2-1))     - zr(jcoord+3*(nuno1-1))
                        xyzom(2) = zr(jcoord+3*(ino2-1) + 1) - zr(jcoord+3*(nuno1-1) + 1)
                        xyzom(3) = zr(jcoord+3*(ino2-1) + 2) - zr(jcoord+3*(nuno1-1) + 2)
                        ! Est-ce que les ddl de rotation existent sur le noeud d'en face ?
                        okddlrot(:) = ASTER_TRUE
                        ibid = jprnm-1+ (nuno1-1)*nbec+1
                        do kddl = 1, 3
                            cmpddl= ddlrot(kddl)
                            icmpg = indik8(zk8(inomg),cmpddl,1,nbcmpg)
                            if (.not.exisdg(zi(ibid),icmpg)) then
                                okddlrot(kddl) = ASTER_FALSE
                            endif
                        enddo
                        if (cmp .eq. 'DX') then
                            call TestDDL(2, 3,  coeffi, iexc)
                            call TestDDL(3, 2, -coeffi, iexc)
                        else if (cmp.eq.'DY') then
                            call TestDDL(1, 3, -coeffi, iexc)
                            call TestDDL(3, 1,  coeffi, iexc)
                        else if (cmp.eq.'DZ') then
                            call TestDDL(1, 2,  coeffi, iexc)
                            call TestDDL(2, 1, -coeffi, iexc)
                        endif
                    enddo
                endif
                call afrela(idcoef, [zeroc],      idnomd, idnomn, idimen, &
                            idirec, nbno1+1+iexc, zeror,  zeroc,  kbeta,  &
                            typcoe, fonree,       coenul, lisrel)
                call imprel(motfac, nbno1+1+iexc, idcoef, idnomd, idnomn, zeror, epsi=coenul)
            enddo
130         continue
            idecal=idecal+nbno1
        enddo cino2
!
        AS_DEALLOCATE(vk8=jddls)
    enddo
!
    AS_DEALLOCATE(vi=idimen)
    AS_DEALLOCATE(vr=idirec)
    AS_DEALLOCATE(vr=idcoef)
    AS_DEALLOCATE(vk8=idnomn)
    AS_DEALLOCATE(vk8=idnomd)
!   affectation de la liste de relations à la charge
    call aflrch(lisrel, charge, 'NLIN')
!
999 continue
    call jedema()
!
! ==================================================================================================
!
contains
! Variables globales : évite de passer par un common !!! NE PAS MODIFIER LES NOMS !!!
!   idnomn, idnomd, idcoef, nbno1, ddlrot, xyzom, nono1, coenul, okddlrot
    subroutine TestDDL(kd, ig, xc, incr)
        integer       :: ig, kd, incr
        real(kind=8)  :: xc
!
#include "asterfort/utmess.h"
!
        if ( abs(xc*xyzom(ig)) .gt. coenul) then
            if ( okddlrot(kd) ) then
                incr= incr+1
                idnomn(1+nbno1+incr) = nono1
                idnomd(1+nbno1+incr) = ddlrot(kd)
                idcoef(1+nbno1+incr) = xc*xyzom(ig)
            else
                valmk(1) = nono1
                valmk(2) = ddlrot(kd)
                call utmess('F', 'MODELISA2_71', nk=2, valk=valmk)
            endif
        endif
!
    end subroutine TestDDL


end subroutine calipj
