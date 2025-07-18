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

subroutine caliag(fonrez, chargz, phenom)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/indik8.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/caexno.h"
#include "asterfort/calemn.h"
#include "asterfort/calinn.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxcadr.h"
#include "asterfort/lxcaps.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=*), intent(in) :: fonrez
    character(len=*), intent(in) :: chargz
    character(len=4), intent(in) :: phenom
!
! ----------------------------------------------------------------------
!
!     CREER LES CARTES CHAR.CHME.CMULT ET CHAR.CHME.CIMPO
!          ET REMPLIR LIGRCH, POUR LE MOT-CLE LIAISON_GROUP
!
! IN  : FONREE : 'REEL' OU 'FONC'
! IN  : CHARGE : NOM UTILISATEUR DU RESULTAT DE CHARGE
!-----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, iret, iocc, ifm, niv, nmocl, ier
    integer(kind=8) :: vali(2)
!-----------------------------------------------------------------------
    integer(kind=8) :: icmpz, idco1, idco2, idconi, idconr, iddl1
    integer(kind=8) :: iddl2, idg1, idg2, idim, idmax
    integer(kind=8) :: iec, iexcm1, iexcm2, imult1
    integer(kind=8) :: imult2, ino1, ino2, inom, jprnm
    integer(kind=8) :: lonli1, lonli2, nb, nbcmp, nbec, nbno, nbterm
    integer(kind=8) :: nddl1, nddl2, nddla, nliag, nmult1, nmult2
!-----------------------------------------------------------------------
    parameter(nmocl=300)
    real(kind=8) :: beta
    complex(kind=8) :: betac
    aster_logical :: dnor, lcolle
    character(len=4) :: fonree, typcoe
    character(len=7) :: typcha
    character(len=8) :: nomno1, nomno2, charge, nomg, noma, mod, nomdep
    character(len=8) :: kbeta, cmp, nomcmp(nmocl), valk(2)
    character(len=16) :: motfac, mcgrex, mcex
    character(len=19) :: prefix, ligrmo, lisrel
    character(len=24) :: coni, conr, nomdd1, nomdd2, coef1, coef2, lisin1
    character(len=24) :: lisin2
    real(kind=8), pointer :: coef(:) => null()
    complex(kind=8), pointer :: coemuc(:) => null()
    integer(kind=8), pointer :: dim(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    integer(kind=8), pointer :: nbnor(:) => null()
    character(len=8), pointer :: nomddl(:) => null()
    character(len=8), pointer :: nomnoe(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    motfac = 'LIAISON_GROUP'
    call getfac(motfac, nliag)
    if (nliag .eq. 0) goto 999
!
    fonree = fonrez
    charge = chargz
!
    typcoe = 'REEL'
    if (fonree .eq. 'COMP') typcoe = 'COMP'
!
    lisrel = '&&CALIAG.RLLISTE'
    lisin1 = '&&CALIAG.LISNO1'
    lisin2 = '&&CALIAG.LISNO2'
    betac = (1.0d0, 0.0d0)
    nomdep = 'DEPL'
!
! --- MODELE ASSOCIE AU LIGREL DE CHARGE ---
!
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', repk=mod)
    ligrmo = mod(1:8)//'.MODELE'
    call jeveuo(ligrmo//'.LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)

    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
    mcgrex = 'SANS_GROUP_NO'
    mcex = 'SANS_NOEUD'
    prefix = charge//'.LIAG.COUPL'
    coni = prefix//'.CONI'
    call jecrec(coni, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
    conr = prefix//'.CONR'
    call jecrec(conr, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
    nomdd1 = prefix(1:13)//'.NOMDDL1'
    call jecrec(nomdd1, 'V V K8', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
    nomdd2 = prefix(1:13)//'.NOMDDL2'
    call jecrec(nomdd2, 'V V K8', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
    coef1 = prefix(1:13)//'.CMULT1'
    call jecrec(coef1, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
    coef2 = prefix(1:13)//'.CMULT2'
    call jecrec(coef2, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nliag)
!
    do iocc = 1, nliag
!
! --- LECTURE DES MOTS CLES GROUP_MA_1 OU 2 OU MAILLE_1 OU 2 OU ---
! --- GROUP_NO_1 OU 2 OU NOEUD_1 OU 2)                          ---
!
        call calemn(motfac, noma, iocc, lisin1, lonli1, &
                    lisin2, lonli2)
!
! --- CONSTRUCTION DES VIS A VIS DES LISTES DE NOEUDS         ---
! --- LISIN1 ET LISIN2 DANS L'OJB CONI(IOCC)                  ---
!
        call calinn(prefix, noma, motfac, iocc, lisin1, &
                    lonli1, lisin2, lonli2, mod)
!
! --- LECTURE DES MOTS CLES SANS_GROUP_NO ET SANS_NOEUD ---
! --- MISE A JOUR DE CONI ET CONR SI IL EXISTE          ---
!
        call caexno(coni, noma, motfac, mcgrex, mcex, &
                    iocc)
        call jeveuo(jexnum(coni, iocc), 'L', idconi)
        nbno = zi(idconi)
        call jeexin(jexnum(conr, iocc), iret)
        if (iret .eq. 0) then
            dnor = .false.
        else
            dnor = .true.
        end if
!
! --- LECTURE DES DDLS IMPOSES SUR LA LISTE 1 ---
!
        call getvtx(motfac, 'DDL_1', iocc=iocc, nbval=0, nbret=nddl1)
        nddl1 = -nddl1
!
! --- LECTURE DES COEF. MULT. ASSOCIES AUX DDLS IMPOSES ---
! --- SUR LA LISTE 1                                    ---
!
        call getvr8(motfac, 'COEF_MULT_1', iocc=iocc, nbval=0, nbret=nmult1)
        nmult1 = -nmult1
        if (nddl1 .ne. nmult1) then
            vali(1) = nddl1
            vali(2) = nmult1
            call utmess('F', 'MODELISA8_43', ni=2, vali=vali)
        end if
!
        call jecroc(jexnum(nomdd1, iocc))
        call jeecra(jexnum(nomdd1, iocc), 'LONMAX', nddl1)
        call jeveuo(jexnum(nomdd1, iocc), 'E', iddl1)
        call jeecra(jexnum(nomdd1, iocc), 'LONUTI', nddl1)
        call getvtx(motfac, 'DDL_1', iocc=iocc, nbval=nddl1, vect=zk8(iddl1))
        do k = 1, nddl1
            call lxcaps(zk8(iddl1-1+k))
            call lxcadr(zk8(iddl1-1+k))
        end do
!
        call jecroc(jexnum(coef1, iocc))
        call jeecra(jexnum(coef1, iocc), 'LONMAX', nddl1)
        call jeveuo(jexnum(coef1, iocc), 'E', imult1)
        call getvr8(motfac, 'COEF_MULT_1', iocc=iocc, nbval=nddl1, vect=zr(imult1))
!
! --- CAS DE DNOR : ON VA GENERE UNE LIAISON SUR DX,DY DZ POUR ---
! --- CHAQUE COUPLE DE LA LIST(UN GREL PAR COUPLE)             ---
!
        if ((nddl1 .eq. 1) .and. (zk8(iddl1) .eq. 'DNOR')) then
            if (.not. dnor) then
                call utmess('F', 'MODELISA2_94')
            end if
        end if
!
! --- LECTURE DES DDLS IMPOSES SUR LA LISTE 2 ---
!
        call getvtx(motfac, 'DDL_2', iocc=iocc, nbval=0, nbret=nddl2)
        nddl2 = -nddl2
!
! --- LECTURE DES COEF. MULT. ASSOCIES AUX DDLS IMPOSES ---
! --- SUR LA LISTE 2                                    ---
!
        call getvr8(motfac, 'COEF_MULT_2', iocc=iocc, nbval=0, nbret=nmult2)
        nmult2 = -nmult2
        if (nddl2 .ne. nmult2) then
            vali(1) = nddl2
            vali(2) = nmult2
            call utmess('F', 'MODELISA8_44', ni=2, vali=vali)
        end if
!
        call jecroc(jexnum(nomdd2, iocc))
        call jeecra(jexnum(nomdd2, iocc), 'LONMAX', nddl2)
        call jeveuo(jexnum(nomdd2, iocc), 'E', iddl2)
        call jeecra(jexnum(nomdd2, iocc), 'LONUTI', nddl2)
        call getvtx(motfac, 'DDL_2', iocc=iocc, nbval=nddl2, vect=zk8(iddl2))
        do k = 1, nddl2
            call lxcaps(zk8(iddl2-1+k))
            call lxcadr(zk8(iddl2-1+k))
        end do
!
        call jecroc(jexnum(coef2, iocc))
        call jeecra(jexnum(coef2, iocc), 'LONMAX', nddl2)
        call jeveuo(jexnum(coef2, iocc), 'E', imult2)
        call getvr8(motfac, 'COEF_MULT_2', iocc=iocc, nbval=nddl2, vect=zr(imult2))
        if ((nddl2 .eq. 1) .and. (zk8(iddl2) .eq. 'DNOR')) then
            if (.not. dnor) then
                call utmess('F', 'MODELISA2_94')
            end if
        end if
    end do
!
! --- TYPE DE LA CHARGE ---
!
    call dismoi('TYPE_CHARGE', charge(1:8), 'CHARGE', repk=typcha)
!
    if (typcha(1:2) .eq. 'TH') then
        nomg = 'TEMP_R'
    else
        nomg = 'DEPL_R'
    end if
!
! --- NOMBRE D'ENTIERS CODES ASSOCIE A LA GRANDEUR ---
!
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 10)
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', inom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmp)
    nddla = nbcmp-1
    if (nddla .gt. nmocl) then
        vali(1) = nmocl
        vali(2) = nddla
        call utmess('F', 'MODELISA8_29', ni=2, vali=vali)
    end if
    do i = 1, nddla
        nomcmp(i) = zk8(inom-1+i)
    end do
!
    call jeveuo(ligrmo//'.PRNM', 'L', jprnm)
!
! --- CREATION ET AFFECTATION DES RELATIONS A LA LISTE DE ---
! --- RELATIONS                                           ---
!
    icmpz = indik8(nomcmp, 'DZ', 1, nddla)
!
    do iocc = 1, nliag
        if (fonree .eq. 'REEL') then
            call getvr8(motfac, 'COEF_IMPO', iocc=iocc, scal=beta, nbret=nb)
        else
            call getvid(motfac, 'COEF_IMPO', iocc=iocc, scal=kbeta, nbret=nb)
        end if
        call jeveuo(jexnum(coni, iocc), 'L', idconi)
! --- NOMBRE DE NOEUDS DE CHACUNES DES LISTES EN VIS A VIS
        nbno = zi(idconi)
!
        call jeveuo(jexnum(nomdd1, iocc), 'L', iddl1)
        call jeveuo(jexnum(nomdd2, iocc), 'L', iddl2)
        call jeveuo(jexnum(coef1, iocc), 'L', idco1)
        call jeveuo(jexnum(coef2, iocc), 'L', idco2)
!
! --- NOMBRE DE DDL IMPOSES POUR LES NOEUDS DE LA 1ERE LISTE ---
!
        call jelira(jexnum(nomdd1, iocc), 'LONUTI', nddl1)
!
! --- NOMBRE DE DDL IMPOSES POUR LES NOEUDS DE LA 2EME LISTE ---
!
        call jelira(jexnum(nomdd2, iocc), 'LONUTI', nddl2)
!
        idmax = 3*(nddl1+nddl2)
!
! ---  ALLOCATION D'UN TABLEAU BIDON POUR AFRELA ---
!
        AS_ALLOCATE(vc=coemuc, size=idmax)
!
! ---  ALLOCATION DU TABLEAU DES NOMS DES NOEUDS DE LA RELATION ---
!
        AS_ALLOCATE(vk8=nomnoe, size=idmax)
!
! ---  ALLOCATION DU TABLEAU DES NOMS DES DDLS DE LA RELATION ---
!
        AS_ALLOCATE(vk8=nomddl, size=idmax)
!
! ---  ALLOCATION DU TABLEAU DES COEFFICIENTS DE LA RELATION ---
!
        AS_ALLOCATE(vr=coef, size=idmax)
!
! ---  ALLOCATION DU TABLEAU DES DIRECTIONS DES COMPOSANTES ---
! ---  DE LA RELATION                                       ---
!
        AS_ALLOCATE(vr=direct, size=3*idmax)
!
! ---  ALLOCATION DU TABLEAU DE LA DIMENSION DU PROBLEME  ---
! ---  RELATIVE A CHAQUE COMPOSANTE DE LA RELATION        ---
!
        AS_ALLOCATE(vi=dim, size=idmax)
!
! ---  ALLOCATION DU TABLEAU DES DIMENSIONS DES VECTEURS NORMAUX ---
! ---  EN CHAQUE NOEUD POUR TOUTES LES RELATIONS                 ---
!
        AS_ALLOCATE(vi=nbnor, size=2*nbno)
!
! ---  AFFECTATION DE CE VECTEUR ---
!
        do j = 1, nbno
            ino1 = zi(idconi+2*(j-1)+1)
            ino2 = zi(idconi+2*(j-1)+2)
!
            iexcm1 = 0
            iexcm2 = 0
            do iec = 1, nbec
                if (zi(jprnm-1+(ino1-1)*nbec+iec) .ne. 0) then
                    iexcm1 = 1
                    goto 60
                end if
            end do
60          continue
!
            do iec = 1, nbec
                if (zi(jprnm-1+(ino2-1)*nbec+iec) .ne. 0) then
                    iexcm2 = 1
                    goto 80
                end if
            end do
80          continue
            idg1 = jprnm-1+(ino1-1)*nbec+1
            idg2 = jprnm-1+(ino2-1)*nbec+1
!
            if (iexcm1 .eq. 0) then
                nomno1 = int_to_char8(ino1, lcolle, noma, "NOEUD")
                call utmess('F', 'CHARGES2_33', sk=nomno1)
            end if
            nbnor(2*(j-1)+1) = 3
            if ((icmpz .eq. 0) .or. (.not. exisdg(zi(idg1), icmpz))) then
                nbnor(2*(j-1)+1) = 2
            end if
!
            if (iexcm2 .eq. 0) then
                nomno2 = int_to_char8(ino2, lcolle, noma, "NOEUD")
                call utmess('F', 'CHARGES2_33', sk=nomno2)
            end if
            nbnor(2*(j-1)+2) = 3
            if ((icmpz .eq. 0) .or. (.not. exisdg(zi(idg2), icmpz))) then
                nbnor(2*(j-1)+2) = 2
            end if
        end do
!
! ---  AFFECTATION DES RELATIONS ---
!
        do j = 1, nbno
            k = 0
!
! --- PREMIER NOEUD DE LA RELATION ---
!
            ino1 = zi(idconi+2*(j-1)+1)
            nomno1 = int_to_char8(ino1, lcolle, noma, "NOEUD")
            cmp = zk8(iddl1)
            if (cmp .eq. 'DNOR') then
                call jeveuo(jexnum(conr, iocc), 'L', idconr)
                idim = nbnor(2*(j-1)+1)
                k = k+1
                coef(k) = zr(idco1)
                nomnoe(k) = nomno1
                nomddl(k) = nomdep
                dim(k) = idim
                do i = 1, idim
                    direct(1+3*(k-1)+i-1) = zr(idconr-1+(2*idim+1)*(j-1)+i)
                end do
            else
                do i = 1, nddl1
                    k = k+1
                    nomnoe(k) = nomno1
                    nomddl(k) = zk8(iddl1+i-1)
                    coef(k) = zr(idco1+i-1)
                end do
            end if
!
! --- DEUXIEME NOEUD DE LA RELATION ---
!
            ino2 = zi(idconi+2*(j-1)+2)
            nomno2 = int_to_char8(ino2, lcolle, noma, "NOEUD")
            cmp = zk8(iddl2)
            if (cmp .eq. 'DNOR') then
                call jeveuo(jexnum(conr, iocc), 'L', idconr)
                idim = nbnor(2*(j-1)+2)
                k = k+1
                coef(k) = zr(idco2)
                nomnoe(k) = nomno2
                nomddl(k) = nomdep
                dim(k) = idim
                do i = 1, idim
                    direct(1+3*(k-1)+i-1) = zr(idconr-1+(2*idim+1)*(j-1)+idim+i)
                end do
            else
                do i = 1, nddl2
                    k = k+1
                    nomnoe(k) = nomno2
                    nomddl(k) = zk8(iddl2+i-1)
                    coef(k) = zr(idco2+i-1)
                end do
            end if
!
! --- NOMBRE DE TERMES DE LA RELATION ---
!
            nbterm = k
!
! --- AFFECTATION DE LA RELATION ---
!
            call afrela(coef, coemuc, nomddl, nomnoe, dim, &
                        direct, nbterm, beta, betac, kbeta, &
                        typcoe, fonree, 0.d0, lisrel)
!
! --- FIN DE LA BOUCLE SUR LES RELATIONS                       ---
! --- (I.E. LES COUPLES DE NOEUDS EN VIS A VIS POUR LE MOT-CLE ---
! --- LIAISON-GROUP COURANT)                                   ---
!
        end do
!
! --- IMPRESSION DES COUPLES DE NOEUDS EN VIS-A-VIS ---
!
        call infniv(ifm, niv)
        if (niv .eq. 2) then
            call utmess('I', 'CHARGES2_35', si=iocc)
            do j = 1, nbno
                ino1 = zi(idconi+2*(j-1)+1)
                nomno1 = int_to_char8(ino1, lcolle, noma, "NOEUD")
                ino2 = zi(idconi+2*(j-1)+2)
                nomno2 = int_to_char8(ino2, lcolle, noma, "NOEUD")
                valk(1) = nomno1
                valk(2) = nomno2
                call utmess('I', 'CHARGES2_36', nk=2, valk=valk)
            end do
        end if
!
! --- DESTRUCTION DES TABLEAUX DE TRAVAIL  ---
!
        AS_DEALLOCATE(vc=coemuc)
        AS_DEALLOCATE(vk8=nomnoe)
        AS_DEALLOCATE(vk8=nomddl)
        AS_DEALLOCATE(vr=coef)
        AS_DEALLOCATE(vr=direct)
        AS_DEALLOCATE(vi=dim)
        AS_DEALLOCATE(vi=nbnor)
        call jedetr(lisin1)
        call jedetr(lisin2)
!
! --- FIN DE LA BOUCLE SUR LES OCCURENCES DU MOT-CLE  ---
! --- LIAISON-GROUP                                   ---
!
    end do
!
! --- AFFECTATION DE LA LISTE DE RELATIONS A LA CHARGE  ---
!
    if (phenom .eq. 'MECA') then
    end if
    call aflrch(lisrel, charge, 'LIN')
!
999 continue
    call jedema()
end subroutine
