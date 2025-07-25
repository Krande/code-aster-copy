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

subroutine tran77(nomres, typres, nomin, basemo)

    use DynaGene_module

    implicit none
!
!
! IN  : NOMRES : NOM UTILISATEUR POUR LA COMMANDE REST_SOUS_STRUC
! IN  : TYPRES : TYPE DE RESULTAT : 'DYNA_TRANS'
! IN  : NOMIN  : NOM UTILISATEUR DU CONCEPT TRAN_GENE AMONT
! IN  : BASEMO : NOM UTILISATEUR DU CONCEPT MODE_MECA AMONT
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/cnocre.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/extrac.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdgeph.h"
#include "asterfort/rbph01.h"
#include "asterfort/rbph02.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rstran.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/vtcrec.h"
#include "asterfort/vtdefs.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=24) :: valk(2)
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, itresu(8), i_cham(8)
    integer(kind=8) :: foci, focf, fomi, fomf, fomo
    real(kind=8) :: r8b, epsi
    complex(kind=8) :: cbid
    character(len=8) :: k8b, blanc, basemo, crit, interp, basem2, mailla, nomres
    character(len=8) :: nomin, mode, nomma, matgen, nomgd
    character(len=14) :: numddl
    character(len=16) :: typres, type(8), typcha, typbas(8), concep
    character(len=19) :: kinst, knume, trange, typref(8), prof
    character(len=24) :: matric, chamno, nomcha, objve1
    character(len=24) :: objve2, objve3, objve4
    aster_logical :: tousno, multap, leffor, prems
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iarchi, ibid, ich, iadrif
    integer(kind=8) :: idec, idefm, inocmp
    integer(kind=8) :: inoecp, inuddl, inumno, iret, iretou, isk
    integer(kind=8) :: jc, jinst, jnume, linst
    integer(kind=8) :: lvale, n1, n2, n3, n4, nbcham, nbinsg
    integer(kind=8) :: nbinst, nbmode, nbnoeu, ncmp, neq, nfonct, shift, i_bloc
    real(kind=8), pointer :: base(:) => null()
    real(kind=8), pointer :: vectgene(:) => null()
    character(len=24), pointer :: refn(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: resu(:) => null()
    type(DynaGene) :: dyna_gene
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
!      DATA CHAMN2   /'&&TRAN77.CHAMN2'/
!      DATA NOMCMP   /'DX      ','DY      ','DZ      ',
!     &               'DRX     ','DRY     ','DRZ     '/
!     ------------------------------------------------------------------
    call jemarq()
!
    blanc = '        '
    mode = basemo
    trange = nomin
    call gettco(nomin, concep)
    nomcha = ' '
    numddl = ' '
!
!     --- RECUPERATION DES ENTITES DU MAILLAGE SUR LESQUELLES ---
!     ---                PORTE LA RESTITUTION                 ---
    tousno = .true.
    prems = .true.
    call getvtx(' ', 'GROUP_NO', nbval=0, nbret=n1)
    call getvtx(' ', 'NOEUD', nbval=0, nbret=n2)
    call getvtx(' ', 'GROUP_MA', nbval=0, nbret=n3)
    call getvtx(' ', 'MAILLE', nbval=0, nbret=n4)
    if (n1+n2+n3+n4 .ne. 0) tousno = .false.
!
!     --- RECUPERATION DE LA BASE MODALE ---
!
    call jeveuo(trange//'.DESC', 'L', vi=desc)
    nbmode = desc(2)
!
!
    if (mode .eq. blanc) then
        call dismoi('REF_RIGI_PREM', trange, 'RESU_DYNA', repk=matgen, arret='C')
        call dismoi('REF_INTD_PREM', trange, 'RESU_DYNA', repk=basemo)
        if (matgen(1:8) .ne. blanc) then
            call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric, arret='C')
            if (matric .ne. blanc) then
                call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=numddl)
                call dismoi('NOM_MAILLA', matric, 'MATR_ASSE', repk=mailla)
                if (tousno) call dismoi('NB_EQUA', matric, 'MATR_ASSE', repi=neq)
            else
                call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=numddl)
                call dismoi('NOM_GD', numddl, 'NUME_DDL', repk=nomgd)
                call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
                if (tousno) call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
            end if
        else
!  POUR LES CALCULS SANS MATRICE GENERALISEE (PROJ_MESU_MODAL)
            call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=matric)
            if (matric(1:8) .eq. blanc) then
                call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric)
                call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=numddl)
            else
                numddl = matric(1:8)
            end if
            call jeveuo(numddl//'.NUME.REFN', 'L', vk24=refn)
            matric = refn(1)
            mailla = matric(1:8)
            call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric)
            if (tousno) call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
        end if
!
        basem2 = basemo
!
!
    else
!         --- BASE MODALE CALCULEE PAR SOUS-STRUCTURATION
!
        call rsexch('F', basemo, 'DEPL', 1, nomcha, &
                    iret)
        nomcha = nomcha(1:19)//'.REFE'
        call dismoi('NOM_GD', nomcha, 'CHAM_NO', repk=nomgd)
        call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=mailla)
!
!------ON VERIFIE QUE L'UTILISATEUR A RENSEIGNE LE MEME SUPPORT DE
!------RESTITUTION DANS LE FICHIER DE COMMANDE
        call getvid(' ', 'SQUELETTE', scal=nomma, nbret=isk)
        if (isk .ne. 0) then
            if (nomma .ne. mailla) then
                valk(1) = nomma
                valk(2) = mailla
                call utmess('F', 'SOUSTRUC2_9', nk=2, valk=valk)
            end if
        end if
!
        if (tousno) then
            call dismoi('NB_EQUA', nomcha, 'CHAM_NO', repi=neq)
        end if
        basem2 = ' '
        call jeveuo(nomcha, 'L', iadrif)
        matric = zk24(iadrif+1)
        numddl = matric(1:14)
    end if
!
    multap = .false.
!
!     ---   RECUPERATION DES VECTEURS DEPLACEMENT, VITESSE ET   ---
!     --- ACCELERATION GENERALISES SUIVANT LES CHAMPS SOUHAITES ---
    nfonct = 0
    call rbph01(trange, nbcham, type, itresu, nfonct, &
                basem2, typref, typbas, tousno, multap, i_cham)
!
!     --- RECUPERATION DES NUMEROS DES NOEUDS ET DES DDLS ASSOCIES ---
!     ---         DANS LE CAS D'UNE RESTITUTION PARTIELLE          ---
!
    if (.not. tousno) then
        objve1 = '&&TRAN77.NUME_NOEUD  '
        objve2 = '&&TRAN77.NOM_CMP     '
        objve3 = '&&TRAN77.NB_NEQ      '
        objve4 = '&&TRAN77.NUME_DDL    '
        call rbph02(mailla, numddl, nomcha, nomgd, neq, &
                    nbnoeu, objve1, ncmp, objve2, objve3, &
                    objve4)
        call jeveuo(objve1, 'L', inumno)
        call jeveuo(objve2, 'L', inocmp)
        call jeveuo(objve3, 'L', inoecp)
        call jeveuo(objve4, 'L', inuddl)
    end if
!
!     --- RECUPERATION DES INSTANTS ---
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERPOL', scal=interp, nbret=n1)
!
    knume = '&&TRAN77.NUM_RANG'
    kinst = '&&TRAN77.INSTANT'
    call rstran(interp, trange, ' ', 1, kinst, &
                knume, nbinst, iretou)
    if (iretou .ne. 0) then
        call utmess('F', 'UTILITAI4_24')
    end if
    call jeexin(kinst, iret)
    if (iret .gt. 0) then
        call jeveuo(kinst, 'L', jinst)
        call jeveuo(knume, 'L', jnume)
    end if
!
!     --- CREATION DE LA SD RESULTAT ---
    call rscrsd('G', nomres, typres, nbinst)
!
!     --- RESTITUTION SUR LA BASE REELLE ---
!
! VERIFICATION QU'IL Y UN DE CES MOTS CLEFS :
!  'LIST_INST', 'LIST_FREQ', 'INST' ou 'FREQ'
! A MOINS QUE L'ON NE SOIT DANS UN CAS DE DOUBLE RESTITUTION
! APRES UNE DOUBLE PROJECTION (PRESENCE DU MOT CLEF 'MODE_MECA')
    foci = 0
    focf = 0
    fomi = 0
    fomf = 0
    fomo = 0
    call getvid(' ', 'LIST_INST', scal=k8b, nbret=foci)
    call getvid(' ', 'LIST_FREQ', scal=k8b, nbret=focf)
    call getvr8(' ', 'INST', scal=r8b, nbret=fomi)
    call getvr8(' ', 'FREQ', scal=r8b, nbret=fomf)
    call getvid(' ', 'MODE_MECA', scal=k8b, nbret=fomo)
    if ((interp(1:3) .ne. 'NON') .and. &
        ( &
        foci .eq. 0 .and. focf .eq. 0 .and. fomi .eq. 0 .and. fomf .eq. 0 .and. fomo .eq. 0 &
        )) then
        call utmess('F', 'ALGORITH10_95')
    end if
!

    call dyna_gene%init(trange(1:8))

    if (interp(1:3) .ne. 'NON') then
        AS_ALLOCATE(vr=vectgene, size=nbmode)
    end if

    do ich = 1, nbcham
        leffor = .true.
        if (type(ich) .eq. 'DEPL' .or. type(ich) .eq. 'VITE' .or. type(ich) .eq. 'ACCE' .or. &
            type(ich) .eq. 'ACCE_ABSOLU') leffor = .false.
!
!            --- RECUPERATION DES DEFORMEES MODALES ---
!
        typcha = typbas(ich)
        call rsexch('F', basemo, typcha, 1, nomcha, &
                    iret)
        nomcha = nomcha(1:19)//'.VALE'
        call jeexin(nomcha, ibid)
        if (ibid .gt. 0) then
            nomcha(20:24) = '.VALE'
        else
            nomcha(20:24) = '.CELV'
        end if
        if (leffor) call jelira(nomcha, 'LONMAX', neq)
        AS_ALLOCATE(vr=base, size=nbmode*neq)
        if (tousno) then
            call copmod(basemo, champ=typcha, numer=numddl, bmodr=base, nequa=neq)
        else
            do j = 1, nbmode
                call rsexch('F', basemo, typcha, j, nomcha, &
                            iret)
                call jeexin(nomcha(1:19)//'.VALE', ibid)
                if (ibid .gt. 0) then
                    nomcha(20:24) = '.VALE'
                else
                    nomcha(20:24) = '.CELV'
                end if
                call jeveuo(nomcha, 'L', idefm)
                idec = 0
                do i = 1, nbnoeu
                    do jc = 1, ncmp
                        if (zi(inoecp-1+(i-1)*ncmp+jc) .eq. 1) then
                            idec = idec+1
                            base(1+(j-1)*neq+idec-1) = zr(idefm+zi(inuddl+idec-1)-1)
                        end if
                    end do
                end do
            end do
        end if
        iarchi = 0
        if (interp(1:3) .eq. 'NON') then
            call dyna_gene%has_field(dyna_gene%ordr, iret)
            if (iret .ne. 0 .and. zi(jnume) .eq. 1) iarchi = -1
        end if
        do i = 0, nbinst-1
            iarchi = iarchi+1
            call rsexch(' ', nomres, type(ich), iarchi, chamno, &
                        iret)
            if (iret .eq. 0) then
                call utmess('A', 'ALGORITH2_64', sk=chamno)
            else if (iret .eq. 100) then
                if (tousno) then
                    if (mode .eq. blanc) then
                        if (leffor) then
                            call vtdefs(chamno, typref(ich), 'G', 'R')
                        else
                            call vtcreb(chamno, 'G', 'R', &
                                        nume_ddlz=numddl, &
                                        nb_equa_outz=neq)
                        end if
                    else
                        call vtcrec(chamno, nomcha, 'G', 'R', neq)
                    end if
                else
                    if (prems) then
                        prems = .false.
                        call cnocre(mailla, nomgd, nbnoeu, zi(inumno), ncmp, &
                                    zk8(inocmp), zi(inoecp), 'G', ' ', chamno)
                        call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=prof)
                    else
                        call cnocre(mailla, nomgd, nbnoeu, zi(inumno), ncmp, &
                                    zk8(inocmp), zi(inoecp), 'G', prof, chamno)
                    end if
                end if
            else
                ASSERT(.false.)
            end if
            chamno(20:24) = '.VALE'
            call jeexin(chamno, ibid)
            if (ibid .gt. 0) then
                chamno(20:24) = '.VALE'
            else
                chamno(20:24) = '.CELV'
            end if
            call jeveuo(chamno, 'E', lvale)
!
            if (leffor .or. .not. tousno) call jelira(chamno, 'LONMAX', neq)
            if (interp(1:3) .ne. 'NON') then
                call dyna_gene%get_values_by_disc(i_cham(ich), zr(jinst+i), length=nbinsg, vr=resu)
                call dyna_gene%get_current_bloc(i_cham(ich), i_bloc)
                call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)
                call extrac(interp, epsi, crit, nbinsg, disc, &
                            zr(jinst+i), resu, nbmode, vectgene, ibid)
            else
                call dyna_gene%get_values_by_index(i_cham(ich), zi(jnume+i), shift, vr=resu)
                vectgene => resu((zi(jnume+i)-1-shift)*nbmode+1:(zi(jnume+i)-shift)*nbmode)
            end if
            call mdgeph(neq, nbmode, base, vectgene, zr(lvale))
!
            call rsnoch(nomres, type(ich), iarchi)
            call rsadpa(nomres, 'E', 1, 'INST', iarchi, &
                        0, sjv=linst, styp=k8b)
            zr(linst) = zr(jinst+i)
        end do
        AS_DEALLOCATE(vr=base)
    end do
!

    if (interp(1:3) .ne. 'NON') then
        AS_DEALLOCATE(vr=vectgene)
    end if

    call dyna_gene%free()

!
!
    if (mode .eq. blanc) then
        call refdcp(basemo, nomres)
    end if
!
    call jedetr('&&TRAN77.NUME_NOEUD  ')
    call jedetr('&&TRAN77.NOM_CMP     ')
    call jedetr('&&TRAN77.NB_NEQ      ')
    call jedetr('&&TRAN77.NUME_DDL    ')
    call jedetr('&&TRAN77.NUM_RANG')
    call jedetr('&&TRAN77.INSTANT')
!
    call titre()
!
    call jedema()
end subroutine
