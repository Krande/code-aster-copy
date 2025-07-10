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
! aslint: disable=W1501
!
subroutine tran75(nomres, typres, nomin, basemo)
!
    use DynaGene_module
!
    implicit none
!
!     ------------------------------------------------------------------
!     OPERATEUR DE RETOUR A LA BASE PHYSIQUE A PARTIR DE DONNEES
!     GENERALISEES DANS LE CAS D'UN CALCUL TRANSITOIRE
!     ------------------------------------------------------------------
! IN  : NOMRES : NOM UTILISATEUR POUR LA COMMANDE REST_BASE_PHYS
! IN  : TYPRES : TYPE DE RESULTAT : 'DYNA_TRANS'
! IN  : NOMIN  : NOM UTILISATEUR DU CONCEPT TRAN_GENE AMONT
! IN  : NOMCMD : NOM DE LA COMMANDE : 'REST_BASE_PHYS'
! IN  : BASEMO : NOM UTILISATEUR DU CONCEPT MODE_MECA AMONT
!                (SI CALCUL MODAL PAR SOUS-STRUCTURATION)
!                BLANC SINON
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "blas/daxpy.h"
#include "asterfort/gettco.h"
#include "asterc/r8prem.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/cnocre.h"
#include "asterfort/copmod.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/extrac.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/idensd.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdgep3.h"
#include "asterfort/mdgeph.h"
#include "asterfort/normev.h"
#include "asterfort/pteddl.h"
#include "asterfort/rbph01.h"
#include "asterfort/rbph02.h"
#include "asterfort/refdaj.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rstran.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vecini.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcrea.h"
#include "asterfort/vtcreb.h"
#include "asterfort/vtcrec.h"
#include "asterfort/vtdefs.h"
#include "asterfort/wkvect.h"
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, itresu(8), i_cham(8)
    integer(kind=8) :: foci, focf, fomi, fomf, fomo
    real(kind=8) :: r8b, epsi, alpha, xnorm, coef(3), direction(9), zero
    character(len=1) :: typ1
    character(len=8) :: k8b, blanc, basemo, crit, interp, basem2, mailla, nomres
    character(len=8) :: nomin, nomcmp(6), mode, monmot(2), matgen, nomgd
    character(len=14) :: numddl
    character(len=16) :: typres, type(8), fieldType, typbas(8), concep
    character(len=19) :: fonct(3), kinst, knume, numeq, numeq1, trange
    character(len=19) :: typref(8), prof
    character(len=24) :: matric, chamno, crefe(2), nomcha, chamn2, objve1
    character(len=24) :: objve2, objve3, objve4, chmod, tmpcha, valk
    aster_logical :: tousno, multap, leffor, prems, l_corr_stat, l_multi_app
    integer(kind=8) :: iexi
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iarchi, ibid, ich, id
    integer(kind=8) :: idec, idefm, ie
    integer(kind=8) :: ier, inocmp, inoecp, inuddl, inumno, ipsdel, ir
    integer(kind=8) :: irou, jc, jinst
    integer(kind=8) :: jnume, jpsdel, jvec, linst
    integer(kind=8) :: lpsdel, lval2, lvale, n1, n2, n3, borne1, borne2
    integer(kind=8) :: n4, nbcham, nbd, nbexci, nbinsg, nbinst
    integer(kind=8) :: nbmode, nbnoeu, ncmp, neq, nfonct, neq0, ifonct, vali(2), neq1
    integer(kind=8) :: shift, i_bloc
    integer(kind=8) :: indice, taille, inst, indice1, taille1, ifm, niv, jddl
    complex(kind=8) :: cbid
    real(kind=8), pointer :: base(:) => null()
    real(kind=8), pointer :: vectgene(:) => null()
    character(len=8), pointer :: fvit(:) => null()
    character(len=24), pointer :: refn(:) => null()
    character(len=8), pointer :: facc(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: resu(:) => null()
    character(len=8), pointer :: fdep(:) => null()
    type(DynaGene) :: dyna_gene
    character(len=4) :: kacce, kprof, kprofold
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    data blanc/'        '/
    data chamn2/'&&TRAN75.CHAMN2'/
    data nomcmp/'DX      ', 'DY      ', 'DZ      ',&
     &               'DRX     ', 'DRY     ', 'DRZ     '/
!     ------------------------------------------------------------------
    call jemarq()
! Init. pour routines mdgep...
    kacce = 'CAS3'
    kprof = 'NON'
    kprofold = kprof
    zero = 0.d0
    call infniv(ifm, niv)
!
    cbid = dcmplx(0.d0, 0.d0)
    mode = basemo
    trange = nomin
    call gettco(nomin, concep)
!
    nomcha = ' '
    numddl = ' '
    numeq = ' '
!
!
!     --- RECHERCHE SI UNE ACCELERATION D'ENTRAINEMENT EXISTE ---
    nfonct = 0
    call getvid(' ', 'ACCE_MONO_APPUI', nbret=nfonct)
    if (nfonct .ne. 0) then
        nfonct = abs(nfonct)
        call getvid(' ', 'ACCE_MONO_APPUI', nbval=nfonct, vect=fonct)
        call getvr8(' ', 'DIRECTION', nbval=3*nfonct, vect=direction, nbret=nbd)
        if (nbd .lt. 3*nfonct) then
            vali(1) = nfonct
            vali(2) = 3*nfonct
            call utmess('F', 'ALGORITH9_85', ni=2, vali=vali)
        end if
        do ifonct = 1, nfonct
            call normev(direction(3*(ifonct-1)+1:3*ifonct), xnorm)
            if (xnorm .lt. r8prem()) then
                call utmess('F', 'ALGORITH9_81')
            end if
        end do
    end if
!
!     --- RECUPERATION DES ENTITES DU MAILLAGE SUR LESQUELLES ---
!     ---                PORTE LA RESTITUTION                 ---
    tousno = .true.
    call getvtx(' ', 'GROUP_NO', nbval=0, nbret=n1)
    call getvtx(' ', 'NOEUD', nbval=0, nbret=n2)
    call getvtx(' ', 'GROUP_MA', nbval=0, nbret=n3)
    call getvtx(' ', 'MAILLE', nbval=0, nbret=n4)
    if (n1+n2+n3+n4 .ne. 0) tousno = .false.
!
!     --- RECUPERATION DE LA BASE MODALE ---
!
!
    call jeveuo(trange//'.DESC', 'L', vi=desc)
    nbmode = desc(2)
!   correction statique dans le résultat d'entrée
    l_corr_stat = ASTER_FALSE
    call dismoi('CORR_STAT', trange, 'RESU_DYNA', repk=k8b)
    if (k8b .eq. 'OUI') then
        l_corr_stat = ASTER_TRUE
    end if
!   multi-appui dans le résultat d'entrée
    l_multi_app = ASTER_FALSE
    call dismoi('MULT_APPUI', trange, 'RESU_DYNA', repk=k8b)
    if (k8b .eq. 'OUI') then
        l_multi_app = ASTER_TRUE
    end if
!
    if (nfonct .ne. 0 .and. l_multi_app) then
        call utmess('F', 'UTILITAI4_2')
    end if
!
!
    if (mode .eq. blanc) then
!
        call dismoi('BASE_MODALE', trange, 'RESU_DYNA', repk=basemo, arret='C', ier=ir)
        call dismoi('REF_RIGI_PREM', trange, 'RESU_DYNA', repk=matgen, arret='C', ier=ir)
!
        if (matgen(1:8) .ne. blanc) then
            call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=numddl, arret='C', ier=ir)
            call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric, arret='C', ier=ir)
            if (numddl .eq. blanc) then
                if (matric .ne. blanc) then
                    call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=numddl)
                end if
            end if
            numeq = numddl//'.NUME'
            call dismoi('NOM_GD', numddl, 'NUME_DDL', repk=nomgd)
            call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
            if (tousno) call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
        else
!          -- POUR LES CALCULS SANS MATRICE GENERALISEE
!             (PROJ_MESU_MODAL)
            call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=matric, arret='C', &
                        ier=ir)
            if (matric(1:8) .eq. blanc) then
                call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric)
                call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=numddl)
            else
                numddl = matric(1:8)
            end if
            numeq = numddl//'.NUME'
            call jeveuo(numddl//'.NUME.REFN', 'L', vk24=refn)
            matric = refn(1)
            mailla = matric(1:8)
            call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric, arret='C')
            if (tousno) call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
        end if
!
        basem2 = basemo
!
!
    else
!         --- BASE MODALE CALCULEE PAR SOUS-STRUCTURATION
!
        call rsexch('F', basemo, 'DEPL', 1, chmod, &
                    ir)
        call dismoi('NOM_GD', chmod, 'CHAM_NO', repk=nomgd)
        call dismoi('NUME_EQUA', chmod, 'CHAM_NO', repk=numeq)
        call dismoi('NOM_MAILLA', chmod, 'CHAM_NO', repk=mailla)
        crefe(1) = " "
        crefe(2) = numeq
        if (tousno) then
            call dismoi('NB_EQUA', numeq, 'NUME_EQUA', repi=neq)
        end if
        basem2 = ' '
    end if
!
! --- MULTI-APPUIS
!
    multap = .false.
    monmot(1) = blanc
    monmot(2) = blanc
    call getvtx(' ', 'MULT_APPUI', scal=monmot(1), nbret=n1)
    call getvtx(' ', 'CORR_STAT', scal=monmot(2), nbret=n2)
!
    if (monmot(1) .eq. 'OUI') then
        multap = .true.
        valk = 'MULT_APPUI'
        if (.not. l_multi_app) call utmess('F', 'UTILITAI4_1', sk=valk)
    else if (monmot(2) .eq. 'OUI') then
        multap = .true.
        valk = 'CORR_STAT'
        if (.not. l_corr_stat) call utmess('F', 'UTILITAI4_1', sk=valk)
    end if
    if (l_multi_app .and. monmot(1) .ne. 'OUI') then
        call utmess('A', 'UTILITAI4_3')
    end if
!
!     ---   RECUPERATION DES VECTEURS DEPLACEMENT, VITESSE ET   ---
!     --- ACCELERATION GENERALISES SUIVANT LES CHAMPS SOUHAITES ---
    call rbph01(trange, nbcham, type, itresu, nfonct, &
                basem2, typref, typbas, tousno, multap, &
                i_cham)
!
!     --- RECUPERATION DES NUMEROS DES NOEUDS ET DES DDLS ASSOCIES ---
!     ---         DANS LE CAS D'UNE RESTITUTION PARTIELLE          ---
!
    if (.not. tousno) then
        objve1 = '&&TRAN75.NUME_NOEUD  '
        objve2 = '&&TRAN75.NOM_CMP     '
        objve3 = '&&TRAN75.NB_NEQ      '
        objve4 = '&&TRAN75.NUME_DDL    '
        call rbph02(mailla, numddl, chmod, nomgd, neq, &
                    nbnoeu, objve1, ncmp, objve2, objve3, &
                    objve4)
        call jeveuo(objve1, 'L', inumno)
        call jeveuo(objve2, 'L', inocmp)
        call jeveuo(objve3, 'L', inoecp)
        call jeveuo(objve4, 'L', inuddl)
    end if
!
! --- MULTI-APPUIS : RECUP DE L EXCITATION ET DE PSI*DELTA
!
    if (multap) then
        call jeveuo(trange//'.FDEP', 'L', vk8=fdep)
        call jeveuo(trange//'.FVIT', 'L', vk8=fvit)
        call jeveuo(trange//'.FACC', 'L', vk8=facc)
        call jeveuo(trange//'.IPSD', 'L', ipsdel)
        call jelira(trange//'.FDEP', 'LONMAX', nbexci)
        nbexci = nbexci/2
        if (tousno) then
            call vtcreb(chamn2, 'V', 'R', nume_ddlz=numddl, nb_equa_outz=neq)
            chamn2(20:24) = '.VALE'
            call jeveuo(chamn2, 'E', lval2)
            lpsdel = ipsdel
        else
            call wkvect('&&TRAN75.PSI_DELTA', 'V V R', neq, jpsdel)
            idec = 0
            do i = 1, nbnoeu
                do j = 1, ncmp
                    if (zi(inoecp-1+(i-1)*ncmp+j) .eq. 1) then
                        idec = idec+1
                        zr(jpsdel+idec-1) = zr(ipsdel+zi(inuddl+idec-1)-1)
                    end if
                end do
            end do
            call wkvect('&&TRAN75.VAL2', 'V V R', neq, lval2)
            lpsdel = jpsdel
        end if
    end if
!
!     --- RECUPERATION DES INSTANTS ---
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERPOL', scal=interp, nbret=n1)
!
    knume = '&&TRAN75.NUM_RANG'
    kinst = '&&TRAN75.INSTANT'
    call rstran(interp, trange, ' ', 1, kinst, &
                knume, nbinst, irou)
    if (irou .ne. 0) call utmess('F', 'UTILITAI4_24')
    call jeexin(kinst, ir)
    if (ir .gt. 0) then
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
!
    if (interp(1:3) .ne. 'NON') then
        AS_ALLOCATE(vr=vectgene, size=nbmode)
    end if
!
    neq0 = neq
!
!******** Boucle sur les champs
    do ich = 1, nbcham
        leffor = .true.
        if (type(ich) .eq. 'DEPL' .or. type(ich) .eq. 'VITE' .or. type(ich) .eq. 'ACCE' .or. &
            type(ich) .eq. 'ACCE_ABSOLU') leffor = .false.
!
!            --- RECUPERATION DES DEFORMEES MODALES ---
!
        fieldType = typbas(ich)
        call rsexch('F', basemo, fieldType, 1, nomcha, ir)
        nomcha = nomcha(1:19)//'.VALE'
        call jeexin(nomcha, ibid)
        if (ibid .gt. 0) then
            nomcha(20:24) = '.VALE'
        else
            nomcha(20:24) = '.CELV'
        end if
!
        if (leffor) then
            call jelira(nomcha, 'LONMAX', neq)
        else
            neq = neq0
        end if
!
        AS_ALLOCATE(vr=base, size=nbmode*neq)
        if (tousno) then
            call copmod(basemo, champ=fieldType, numer=numeq(1:14), bmodr=base, nequa=neq, &
                        nbmodes=nbmode)
        else
            crefe(1) = " "
            crefe(2) = numeq
            tmpcha = '&&TRAN75.CHAMP'
            call dismoi('NB_EQUA', numeq, 'NUME_EQUA', repi=neq1)
            do j = 1, nbmode
                call rsexch('F', basemo, fieldType, j, nomcha, ir)
                call jeexin(nomcha(1:19)//'.VALE', iexi)
                ASSERT(iexi .gt. 0)
                call jelira(nomcha(1:19)//'.VALE', 'TYPE', cval=typ1)
                ASSERT(typ1 .eq. 'R')
!
!              SI NOMCHA N'A PAS LA BONNE NUMEROTATION, ON ARRETE TOUT :
                ASSERT(numeq .ne. ' ')
                call dismoi('NUME_EQUA', nomcha, 'CHAM_NO', repk=numeq1)
                if (.not. idensd('NUME_EQUA', numeq, numeq1)) then
                    call vtcrea(tmpcha, crefe, 'V', 'R', neq1)
                    call vtcopy(nomcha, tmpcha, ir)
                    if (ir .ne. 0) then
                        call utmess("A", "FIELD0_3", sk=fieldType)
                    end if
                    nomcha = tmpcha
                end if
                nomcha(20:24) = '.VALE'
                call jelira(nomcha(1:19)//'.VALE', 'LONMAX', neq1)
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
                call jeexin(tmpcha(1:19)//'.VALE', iexi)
                if (iexi .gt. 0) then
                    call detrsd('CHAM_NO', tmpcha)
                end if
            end do
        end if
        iarchi = 0
        if (interp(1:3) .eq. 'NON') then
            call dyna_gene%has_field(dyna_gene%ordr, ir)
            if (ir .ne. 0 .and. zi(jnume) .eq. 1) iarchi = -1
        end if
        prems = .true.
!
!******** Boucle sur les instants
        indice = 1
        indice1 = 1
        taille = neq
        taille1 = neq
        kprof = kprofold
!       --- GESTTION DES INCOMPATIBILITES SUIVANT LES ACCELERATIONS
!       --- CAS DE FIGURES RESTANT A ADAPTER SI BESOIN EN PERFORMANCE
        if (kacce(1:4) .eq. 'CAS3') then
            if ((interp(1:3) .ne. 'NON') .or. (multap) .or. &
                ((type(ich) .eq. 'ACCE_ABSOLU') .and. (nfonct .gt. 0))) kacce = 'CAS2'
        end if
        do i = 0, nbinst-1
! Pour MDGEP*
            if (i .eq. (nbinst-1)) then
                if (kprof(1:3) .eq. 'OUI') kprof = 'OUIA'
                inst = -1*nbinst
            else
                inst = i+1
            end if
!
            iarchi = iarchi+1
            call rsexch(' ', nomres, type(ich), iarchi, chamno, &
                        ir)
            if (ir .eq. 0) then
                call utmess('A', 'ALGORITH2_64', sk=chamno)
            else if (ir .eq. 100) then
                if (tousno) then
                    if (mode .eq. blanc) then
                        if (leffor) then
                            call vtdefs(chamno, typref(ich), 'G', 'R')
                        else
                            call vtcreb(chamno, 'G', 'R', nume_ddlz=numddl, nb_equa_outz=neq)
                        end if
                    else
                        call vtcrec(chamno, chmod, 'G', 'R', neq)
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
!
            chamno(20:24) = '.VALE'
            call jeexin(chamno, ibid)
            if (ibid .le. 0) chamno(20:24) = '.CELV'
!
            call jeveuo(chamno, 'E', lvale)
!
            if (leffor .or. .not. tousno) call jelira(chamno, 'LONMAX', neq)
!
            if (interp(1:3) .ne. 'NON') then
                call dyna_gene%get_values_by_disc(i_cham(ich), zr(jinst+i), length=nbinsg, &
                                                  vr=resu)
                call dyna_gene%get_current_bloc(i_cham(ich), i_bloc)
                call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)
                call extrac(interp, epsi, crit, nbinsg, disc, &
                            zr(jinst+i), resu, nbmode, vectgene, ibid)
            else
                call dyna_gene%get_values_by_index(i_cham(ich), zi(jnume+i), shift, vr=resu)
                vectgene => resu((zi(jnume+i)-1-shift)*nbmode+1:(zi(jnume+i)-shift)*nbmode)
            end if
            call mdgeph(neq, nbmode, base, vectgene, zr(lvale), &
                        kacce=kacce, kprof=kprof, inst=inst, instt=nbinst, indice=indice, &
                        taille=taille, kcham=chamno)
            if (multap) then
                ASSERT(kacce .ne. 'CAS3')
                if (type(ich) .eq. 'DEPL') call mdgep3(neq, nbexci, zr(lpsdel), zr(jinst+i), &
                                                       fdep, zr(lval2), kacce=kacce, kprof=kprof, &
                                                       inst=inst, indice=indice1, taille=taille1)
                if (type(ich) .eq. 'VITE') call mdgep3(neq, nbexci, zr(lpsdel), zr(jinst+i), &
                                                       fvit, zr(lval2), kacce=kacce, kprof=kprof, &
                                                       inst=inst, indice=indice1, taille=taille1)
                if (type(ich) .eq. 'ACCE') call mdgep3(neq, nbexci, zr(lpsdel), zr(jinst+i), &
                                                       facc, zr(lval2), kacce=kacce, kprof=kprof, &
                                                       inst=inst, indice=indice1, taille=taille1)
                if (type(ich) .eq. 'ACCE_ABSOLU') call mdgep3(neq, nbexci, zr(lpsdel), &
                                                              zr(jinst+i), facc, zr(lval2), &
                                                              kacce=kacce, kprof=kprof, &
                                                              inst=inst, indice=indice1, &
                                                              taille=taille1)
!               --- VERIFICATION AU CAS OU
                ASSERT((indice .eq. indice1) .and. (taille .eq. taille1))
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, 1.d0, zr(lval2), b_incx, zr(lvale), &
                           b_incy)
            end if
!            --- PRISE EN COMPTE D'UNE ACCELERATION D'ENTRAINEMENT
!            --- PRISE EN COMPTE D'UNE ACCELERATION D'ENTRAINEMENT
            if (type(ich) .eq. 'ACCE_ABSOLU' .and. nfonct .gt. 0) then
                ASSERT(kacce .ne. 'CAS3')
                ir = 0
                coef = 0
                do ifonct = 1, nfonct
                    call fointe('F', fonct(ifonct), 1, 'INST', zr(jinst+i), &
                                alpha, ier)
                    coef = coef+alpha*direction(3*(ifonct-1)+1:3*ifonct)
                end do
!               --- ACCELERATION ABSOLUE = RELATIVE + ENTRAINEMENT
                if (i .eq. 0) then
                    call wkvect('&&TRAN75.VECTEUR', 'V V R', neq, jvec)
                    call wkvect('&&TRAN75.VECTEUI', 'V V I', 3*neq, jddl)
                    if (tousno) then
                        call pteddl('NUME_DDL', numddl, 3, nomcmp, neq, &
                                    tabl_equa=zi(jddl))
                    else
                        call pteddl('CHAM_NO', chamno, 3, nomcmp, neq, &
                                    tabl_equa=zi(jddl))
                    end if
                else
                    call vecini(neq, zero, zr(jvec))
                end if
! Preparation du vecteur auxiliaire zr(jvec)
                borne1 = indice-1
                borne2 = indice+taille-2
                do id = 1, 3
                    do ie = borne1, borne2
                        zr(jvec+ie) = zr(jvec+ie)+zi(jddl+neq*(id-1)+ie)*coef(id)
                    end do
                end do
! COM MPI sur le vecteur auxiliaire zr(jvec) que si au moins 1 ddl par processus
                if (taille .ne. neq) call asmpi_comm_vect('MPI_SUM', 'R', nbval=neq, vr=zr(jvec))
! Somme du vecteur principal, zr(lvale) et du vecteur auxiliaire, zr(jvec)
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, 1.d0, zr(jvec), b_incx, zr(lvale), &
                           b_incy)
!
                if (i .eq. nbinst-1) then
                    call jedetr('&&TRAN75.VECTEUR')
                    call jedetr('&&TRAN75.VECTEUI')
                end if
!            --- FIN PRISE EN COMPTE D'UNE ACCELERATION D'ENTRAINEMENT
            end if
            call rsnoch(nomres, type(ich), iarchi)
            call rsadpa(nomres, 'E', 1, 'INST', iarchi, &
                        0, sjv=linst, styp=k8b)
            zr(linst) = zr(jinst+i)
!
!******** Fin boucle sur les instants
        end do
        AS_DEALLOCATE(vr=base)
!
!******** Fin boucle sur les champs
    end do
!
    if (interp(1:3) .ne. 'NON') then
        AS_DEALLOCATE(vr=vectgene)
    end if
!
    call dyna_gene%free()
!
!
!
    if (mode .eq. blanc) then
        call refdcp(basemo, nomres)
    else
        call refdaj(' ', nomres, -1, ' ', 'INIT', &
                    ' ', ir)
    end if
!
! --- MENAGE
    call detrsd('CHAM_NO', '&&TRAN75.CHAMN2')
    call jedetr('&&TRAN75.NUME_NOEUD  ')
    call jedetr('&&TRAN75.NOM_CMP     ')
    call jedetr('&&TRAN75.NB_NEQ      ')
    call jedetr('&&TRAN75.NUME_DDL    ')
    call jedetr('&&TRAN75.PSI_DELTA')
    call jedetr('&&TRAN75.VAL2')
    call jedetr('&&TRAN75.NUM_RANG')
    call jedetr('&&TRAN75.INSTANT')
!
    call titre()
!
    call jedema()
end subroutine
