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
subroutine aceapo(noma, nomo, lmax, npoutr, nbocc, &
                  mclf, nbepo, ntyele, ivr, zjdlm)
    implicit none
!
    integer(kind=8) :: lmax, npoutr, nbocc, nbepo, zjdlm(*)
    integer(kind=8) :: ntyele(*), ivr(*)
    character(len=8) :: noma, nomo
    character(len=*) :: mclf
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!
!     AFFECTATION DES CARACTÉRISTIQUES POUR L'ÉLEMENT POUTRE
!
! --------------------------------------------------------------------------------------------------
!
! IN  : NOMA   : NOM DU MAILLAGE
!       NOMO   : NOM DU MODELE
!       LMAX   : NOMBRE MAX DE MAILLE OU GROUPE DE MAILLE
!       NPOUTR : NOMBRE DE POUTRE DU MODELE
!       NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE POUTRE
!       NBEPO  : NOMBRE D'ELEMENT DE TYPE POUTRE
!       NTYELE : TABLEAU DES TYPES D'ELEMENTS
!       IVR    : TABLEAU DES INDICES DE VERIFICATION
!       JDLM   : ADRESSE DES MAILLES
!
! --------------------------------------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/acedat.h"
#include "asterfort/affdef.h"
#include "asterfort/affgen.h"
#include "asterfort/affpou.h"
#include "asterfort/alcart.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calc_cara_homo.h"
#include "asterfort/codent.h"
#include "asterfort/coecis.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/modiMetriVeri.h"
#include "asterfort/nocart.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6)  :: kioc
    character(len=8)  :: k8b, nomu, fcx, nomsec
    character(len=8)  :: caram(4)
    character(len=16) :: k16b, sec, concep, cmd, varsec
    character(len=16) :: nunoel
    character(len=19) :: cartpo, cartge, cartpf, tabcar, napcis, foncis
    character(len=24) :: tmpnpo, tmpvpo, tmpgen, tmpnge, tmpvge, typca, nommai
    character(len=24) :: tmpnpf, tmpvpf, tmpgef, modmai, mlggma, mlgnma
    character(len=24) :: vmessk(2)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii, idw, ier, iisec, iivar, ioc, isec, ifm, jj, kk
    integer(kind=8) :: itabl, ivar, ivect, ixma, j
    integer(kind=8) :: jdcge, jdcpo, jdcpof, jdge, jdgef
    integer(kind=8) :: jdgm, jdme, jdvge, jdvpo, jdvpof
    integer(kind=8) :: nbcar, nbcolo, nblign, nbmagr, nbmail, nbo, nbval
    integer(kind=8) :: ncar, ncarac, ndim, nfcx, ng, nm, nnosec
    integer(kind=8) :: npoaff, nsec, nsecpo, ntab, ntypse, nummai, nutyel
    integer(kind=8) :: nval, nvsec, nutyptu(3)
    real(kind=8) :: epy1, hy1
!
    integer(kind=8), pointer :: ncp(:) => null()
    integer(kind=8), pointer :: tab_para(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
!
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: valem(:) => null()
!
    character(len=8), pointer :: cara(:) => null()
    character(len=8), pointer :: carpou(:) => null()
    character(len=8), pointer :: exppou(:) => null()
    character(len=8), pointer :: tabpou(:) => null()
    character(len=16), pointer :: typ_sect(:) => null()
    character(len=24), pointer :: poutre(:) => null()
    character(len=24), pointer :: tblp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    if (npoutr .le. 0) then
        if (isParallelMesh(noma)) goto 999
        ASSERT(.false.)
    end if
    call getres(nomu, concep, cmd)
!
    AS_ALLOCATE(vi=tab_para, size=10)
    call acedat('POUTRE', 0, tab_para, k16b, k8b, k8b, k8b)
    nsecpo = tab_para(1)
    ntypse = tab_para(1+1)
    nbo = tab_para(1+2)
    nbcar = tab_para(1+3)
    nbval = tab_para(1+4)
!   Ceinture et bretelle
    ASSERT(nbcar .eq. nbval)
!
    AS_ALLOCATE(vi=ncp, size=ntypse)
    do ii = 1, ntypse
        ncp(ii) = tab_para(1+4+ii)
    end do
    ndim = ncp(1)*(nsecpo+1)
    AS_ALLOCATE(vk16=typ_sect, size=ntypse)
    AS_ALLOCATE(vk8=exppou, size=nbo)
    AS_ALLOCATE(vk8=tabpou, size=nbo)
    AS_ALLOCATE(vk8=carpou, size=ndim*ntypse)
    call acedat('POUTRE', 1, tab_para, typ_sect, exppou, tabpou, carpou)
    AS_ALLOCATE(vk8=cara, size=nbcar)
    AS_ALLOCATE(vr=vale, size=nbval)
!
    modmai = nomo//'.MAILLE'
    mlgnma = noma//'.TYPMAIL'
    mlggma = noma//'.GROUPEMA'
    ier = 0
    call jelira(mlgnma, 'LONMAX', nbmail)
    call jeexin(modmai, ixma)
    if (ixma .ne. 0) call jeveuo(modmai, 'L', jdme)
!
!   Construction des cartes
    tmpgen = nomu//'.POUTRE'
    cartpo = nomu//'.CARGENPO'
    cartge = nomu//'.CARGEOPO'
    tmpnpo = cartpo//'.NCMP'
    tmpvpo = cartpo//'.VALV'
    tmpnge = cartge//'.NCMP'
    tmpvge = cartge//'.VALV'
!
    tmpgef = nomu//'.VENT'
    cartpf = nomu//'.CVENTCXF'
    tmpnpf = cartpf//'.NCMP'
    tmpvpf = cartpf//'.VALV'
!
!   Création d un objet tampon (surdimensionné a nbo*npoutr)  :
    call jecrec(tmpgen, 'V V R', 'NO', 'CONTIG', 'CONSTANT', npoutr)
    call jeecra(tmpgen, 'LONMAX', nbo)
    call jecrec(tmpgef, 'V V K8', 'NO', 'CONTIG', 'CONSTANT', npoutr)
    call jeecra(tmpgef, 'LONMAX', 1)
    AS_ALLOCATE(vk24=poutre, size=lmax)
!
!   Récupération des numéros des types d'éléments tuyaux
    kk = 0
    do j = 1, nbepo
        call jenuno(jexnum('&CATA.TE.NOMTE', ntyele(j)), nunoel)
        if ((nunoel .eq. 'MET3SEG3') .or. (nunoel .eq. 'MET6SEG3') .or. &
            (nunoel .eq. 'MET3SEG4')) then
            kk = kk+1
            nutyptu(kk) = ntyele(j)
        end if
    end do
    ASSERT(kk .eq. 3)
!
!   Lecture et stockage des données  dans l objet tampon
    do ioc = 1, nbocc

        call modiMetriVeri(noma, ioc, modmai, nutyptu)

        call getvtx('POUTRE', 'SECTION', iocc=ioc, scal=sec, nbret=nsec)
        if (sec .eq. 'COUDE') cycle
        !
        call codent(ioc, 'G', kioc)
        call getvem(noma, 'GROUP_MA', 'POUTRE', 'GROUP_MA', ioc, lmax, poutre, ng)
        call getvem(noma, 'MAILLE', 'POUTRE', 'MAILLE', ioc, lmax, poutre, nm)
        call getvtx('POUTRE', 'VARI_SECT', iocc=ioc, scal=varsec, nbret=nvsec)
        !
        call getvid('POUTRE', 'TABLE_CARA', iocc=ioc, scal=tabcar, nbret=ntab)
        if (ntab .eq. 1) then
            call getvtx('POUTRE', 'NOM_SEC', iocc=ioc, scal=nomsec, nbret=nnosec)
            ASSERT(nnosec .eq. 1)
            ! on recherche nomsec dans la 1ere colonne
            call jeveuo(tabcar//'.TBNP', 'L', vi=tbnp)
            !  nombre de colonnes, lignes
            nbcolo = tbnp(1)
            nblign = tbnp(2)
            !
            call jeveuo(tabcar//'.TBLP', 'L', vk24=tblp)
            typca = tblp(2)
            if (typca(1:2) .ne. 'K8' .and. typca(1:3) .ne. 'K24') then
                call utmess('F', 'MODELISA8_17', sk=tabcar)
            end if
            call jeveuo(tblp(3), 'L', itabl)
            iisec = 0
            if (typca .eq. 'K8') then
                do ii = 1, nblign
                    if (zk8(itabl-1+ii) .eq. nomsec) then
                        iisec = ii
                        goto 97
                    end if
                end do
            else
                do ii = 1, nblign
                    if (zk24(itabl-1+ii) (1:8) .eq. nomsec) then
                        iisec = ii
                        goto 97
                    end if
                end do
            end if
            vmessk(1) = tabcar
            vmessk(2) = nomsec
            call utmess('F', 'MODELISA8_18', nk=2, valk=vmessk)
97          continue
            !
            jj = 0
            cii1: do ii = 1, nbcolo-1
                if (tblp(1+4*ii+1) .ne. 'R') cycle cii1
                do kk = 1, nbo
                    if (tblp(1+4*ii) .eq. exppou(kk)) then
                        jj = jj+1
                        ! Ceinture
                        ASSERT(jj .le. nbcar)
                        cara(jj) = tblp(1+4*ii) (1:8)
                        call jeveuo(tblp(1+4*ii+2), 'L', ivect)
                        vale(jj) = zr(ivect-1+iisec)
                        cycle cii1
                    end if
                end do
            end do cii1
            ncarac = jj
        else
            call getvtx('POUTRE', 'CARA', iocc=ioc, nbval=nbcar, vect=cara, nbret=ncar)
            call getvr8('POUTRE', 'VALE', iocc=ioc, nbval=nbval, vect=vale, nbret=nval)
            ASSERT(ncar .gt. 0)
            ncarac = ncar
        end if
        !
        fcx = '.'
        call getvid('POUTRE', 'FCX', iocc=ioc, scal=fcx, nbret=nfcx)
        !
        ivar = 2
        ! Type de section et de variation de section pour cette occurence
        !   test de zk8(jcara) seul > vérification d'homogénéité déjà faite
        if (varsec(1:4) .eq. 'AFFI') then
            ivar = 1
            hy1 = 0.d0
            epy1 = 0.d0
            do kk = 1, nbcar
                if (cara(kk) (1:3) .eq. 'HY ') then
                    hy1 = vale(kk)
                    cara(kk) = 'HY1'
                end if
                if (cara(kk) (1:4) .eq. 'EPY ') then
                    epy1 = vale(kk)
                    cara(kk) = 'EPY1'
                end if
            end do
            ncar = ncar+1
            cara(ncar) = 'HY2'
            vale(ncar) = hy1
            if (epy1 .ne. 0.d0) then
                ncar = ncar+1
                cara(ncar) = 'EPY2'
                vale(ncar) = epy1
            end if
            ncarac = ncar
        end if
        !
        if (ntab .eq. 0) then
            do ii = 1, ntypse
                if (sec .eq. typ_sect(ii)) then
                    isec = ii-1
                    do j = 1, ncp(ii)
                        if (cara(1) .eq. carpou(1+j+ndim*(ii-1)-1)) then
                            ivar = 0
                            goto 24
                        end if
                    end do
                end if
            end do
        else
            ! si on a donné TABLE_CARA la section est constante
            ivar = 0
            isec = 0
        end if
24      continue
        iivar = ivar
        ! "GROUP_MA" = toutes les mailles possibles de la liste des groupes de mailles
        if (ng .gt. 0) then
            do ii = 1, ng
                call jeveuo(jexnom(mlggma, poutre(ii)), 'L', jdgm)
                call jelira(jexnom(mlggma, poutre(ii)), 'LONUTI', nbmagr)
                ! traitement des affectations CERCLE HOMOTHETIQUE GROUP_MA
                if (isec .eq. 2 .and. ivar .eq. 2) then
                    AS_ALLOCATE(vr=valem, size=ncarac*nbmagr)
                    call calc_cara_homo(noma, poutre(ii), zi(jdgm), nbmagr, ncarac, &
                                        cara, vale, caram, valem)
                end if
                do j = 1, nbmagr
                    nummai = zi(jdgm+j-1)
                    nommai = int_to_char8(nummai)
                    nutyel = zi(jdme+nummai-1)
                    do kk = 1, nbepo
                        if (nutyel .eq. ntyele(kk)) then
                            if (isec .eq. 2 .and. ivar .eq. 2) then
                                call affpou(tmpgen, tmpgef, fcx, nommai, isec, &
                                            iivar, caram, ncarac, &
                                            valem(ncarac*(j-1)+1:ncarac*j), tabpou, &
                                            exppou, nbo, kioc, ier)
                            else
                                call affpou(tmpgen, tmpgef, fcx, nommai, isec, &
                                            iivar, cara, ncarac, &
                                            vale, tabpou, &
                                            exppou, nbo, kioc, ier)
                            end if
                            iivar = ivar
                            goto 42
                        end if
                    end do
                    vmessk(1) = mclf
                    vmessk(2) = nommai
                    call utmess('F', 'MODELISA_8', nk=2, valk=vmessk)
42                  continue
                end do
                if (isec .eq. 2 .and. ivar .eq. 2) then
                    AS_DEALLOCATE(vr=valem)
                end if
            end do
        end if
        !
        ! "MAILLE" = Toutes les mailles possibles de la liste de mailles
        if (nm .gt. 0) then
            do ii = 1, nm
                nommai = poutre(ii)
                nummai = char8_to_int(nommai)
                nutyel = zi(jdme+nummai-1)
                do j = 1, nbepo
                    if (nutyel .eq. ntyele(j)) then
                        call affpou(tmpgen, tmpgef, fcx, nommai, isec, &
                                    iivar, cara, ncarac, vale, tabpou, &
                                    exppou, nbo, kioc, ier)
                        iivar = ivar
                        goto 50
                    end if
                end do
                vmessk(1) = mclf
                vmessk(2) = nommai
                call utmess('F', 'MODELISA_8', nk=2, valk=vmessk)
50              continue
            end do
        end if
    end do
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA_14')
    end if
    call jelira(tmpgen, 'NUTIOC', npoaff)
!
!   Vérification des obligations et affectation des défauts
    do ii = 1, npoaff
        call jenuno(jexnum(tmpgen, ii), nommai)
        nummai = char8_to_int(nommai)
        nutyel = zi(jdme+nummai-1)
        call affdef(tmpgen, nommai, nutyel, ntyele, tabpou, ier)
    end do
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA_15')
    end if
!
!   Calcul des donnees generales a partir des donnees geometriques (cerc+rect) et
!   affectations dans le tampon
!
!   Utilisation de la nappe pour les sections rectangulaires et de la fonction pour les
!   sections circulaires afin d'interpoler les coefficients de cisaillement
    call coecis(napcis, foncis)
    !
    do ii = 1, npoaff
        call jenuno(jexnum(tmpgen, ii), nommai)
        nummai = char8_to_int(nommai)
        nutyel = zi(jdme+nummai-1)
        call affgen(tmpgen, nommai, nutyel, ntyele, napcis, foncis)
    end do
!
!   Impression des valeurs affectées dans le tampon si demande
    if (ivr(3) .eq. 2) then
        ifm = ivr(4)
        ! Impression des données générales
        write (ifm, 200)
        do ii = 1, npoaff
            call jenuno(jexnum(tmpgen, ii), nommai)
            call jeveuo(jexnum(tmpgen, ii), 'L', jdge)
            ivar = nint(zr(jdge+22))
            isec = nint(zr(jdge+35))
            write (ifm, 201) nommai, (zr(jdge+j-1), j=1, 22), (zr(jdge+j-1), j=37, 44), ivar, isec
            call jenuno(jexnum(tmpgef, ii), nommai)
            call jeveuo(jexnum(tmpgef, ii), 'L', jdgef)
            write (ifm, *) 'CX : ', zk8(jdgef)
        end do
        ! Impression des données géométriques
        idw = 0
        do ii = 1, npoaff
            call jenuno(jexnum(tmpgen, ii), nommai)
            call jeveuo(jexnum(tmpgen, ii), 'L', jdge)
            isec = nint(zr(jdge+35))
            if (isec .gt. 0) then
                if (idw .eq. 0) then
                    write (ifm, 202)
                    idw = 1
                end if
                write (ifm, 203) nommai, (zr(jdge+j-1), j=24, 35), isec
            end if
        end do
    end if
!
200 format(/, 3x, '<SECTION> ', &
            'VALEURS DE TYPE GENERALE AFFECTEES AUX POUTRES', //, 3x, &
            'MAILLE   ', &
            'A1  ', 8x, 'IY1  ', 7x, 'IZ1  ', 7x, 'AY1  ', 7x, 'AZ1  ', /, &
            12x, 'EY1 ', 8x, 'EZ1  ', 7x, 'JX1  ', 7x, 'RY1  ', 7x, 'RZ1  ', /, &
            12x, 'RT1 ', 8x, 'A2   ', 7x, 'IY2  ', 7x, 'IZ2  ', 7x, 'AY2  ', /, &
            12x, 'AZ2 ', 8x, 'EY2  ', 7x, 'EZ2  ', 7x, 'JX2  ', 7x, 'RY2  ', /, &
            12x, 'RZ2 ', 8x, 'RT2  ', 7x, 'AI1  ', 7x, 'AI2  ', 7x, 'JG1  ', /, &
            12x, 'JG2 ', 8x, 'IYR21', 7x, 'IYR22', 7x, 'IZR21', 7x, 'IZR22', /, &
            12x, 'TVAR', 8x, 'TSEC ')
201 format(/, 1p, 3x, a8, 1x, 5(1pd12.5, 1x), 5(/, 12x, 5(1pd12.5, 1x)), &
            /, 12x, i6, 6x, i6)
202 format(/, 3x, '<SECTION> ', &
            'VALEURS DE TYPE GEOMETRIQUE AFFECTEES AUX POUTRES', //, 3x, &
            'MAILLE   HY1         HZ1         EPY1        EPZ1' &
            , /, 12x, 'HY2         HZ2         EPY2        EPZ2' &
            , /, 12x, 'R1          EP1         R2          EP2', 9x, 'TSEC')
203 format(/, 1p, 3x, a8, 1x, 4(1pd12.5, 1x), 2(/, 12x, 4(1pd12.5, 1x)), i6)
!
!   Allocation des cartes
    call alcart('G', cartpo, noma, 'CAGNPO_R')
    call alcart('G', cartge, noma, 'CAGEPO_R')
    call jeveuo(tmpnpo, 'E', jdcpo)
    call jeveuo(tmpvpo, 'E', jdvpo)
    call jeveuo(tmpnge, 'E', jdcge)
    call jeveuo(tmpvge, 'E', jdvge)
    call jeveuo(tmpnpf, 'E', jdcpof)
    call jeveuo(tmpvpf, 'E', jdvpof)
!   Affectations des données générales
    do ii = 1, 23
        zk8(jdcpo+ii-1) = tabpou(ii)
    end do
    do ii = 24, 31
        zk8(jdcpo+ii-1) = tabpou(1+ii+13-1)
    end do
!   POUR LA CARTE DE VENT ==> FCXP
    zk8(jdcpof) = 'FCXP'
!
    do ii = 1, npoaff
        call jenuno(jexnum(tmpgen, ii), nommai)
        nummai = char8_to_int(nommai)
        zjdlm(nummai) = -1
        call jeveuo(jexnum(tmpgen, ii), 'L', jdge)
        do j = 1, 23
            zr(jdvpo+j-1) = zr(jdge+j-1)
        end do
        do j = 24, 31
            zr(jdvpo+j-1) = zr(jdge+j+13-1)
        end do
        call jeveuo(jexnum(tmpgef, ii), 'L', jdgef)
        zk8(jdvpof) = zk8(jdgef)
        call nocart(cartpo, 3, 31, mode='NOM', nma=1, limano=[nommai])
        call nocart(cartpf, 3, 1, mode='NOM', nma=1, limano=[nommai])
    end do
!
!   Affectations données géométriques (on affecte toutes les cmps)
    do ii = 1, 13
        zk8(jdcge+ii-1) = tabpou(1+ii+23-1)
    end do
    do j = 1, npoaff
        call jenuno(jexnum(tmpgen, j), nommai)
        call jeveuo(jexnum(tmpgen, j), 'L', jdge)
        isec = nint(zr(jdge+35))
        if (isec .eq. 0) then
            ! GENERALE
            do ii = 1, 13
                zr(jdvge+ii-1) = 0.d0
            end do
            call nocart(cartge, 3, 13, mode='NOM', nma=1, limano=[nommai])
        else
            ! RECTANGLE OU CERCLE
            do ii = 1, 13
                zr(jdvge+ii-1) = zr(jdge+ii+22)
            end do
            call nocart(cartge, 3, 13, mode='NOM', nma=1, limano=[nommai])
        end if
    end do
!
!   Compactage des cartes (mais on ne cherche pas de rémanence)
    call tecart(cartpo)
    call tecart(cartge)
!
!   Nettoyage
    AS_DEALLOCATE(vi=tab_para)
    AS_DEALLOCATE(vi=ncp)
    AS_DEALLOCATE(vk16=typ_sect)
    AS_DEALLOCATE(vk8=exppou)
    AS_DEALLOCATE(vk8=tabpou)
    AS_DEALLOCATE(vk8=carpou)
    AS_DEALLOCATE(vk8=cara)
    AS_DEALLOCATE(vr=vale)
    AS_DEALLOCATE(vk24=poutre)
    call jedetr(tmpgen)
    call jedetr(tmpgef)
    call jedetr(tmpnpo)
    call jedetr(tmpvpo)
    call jedetr(tmpnge)
    call jedetr(tmpvge)
!
999 continue
    call jedema()
end subroutine
