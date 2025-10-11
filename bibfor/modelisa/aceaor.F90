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

subroutine aceaor(noma, nomo, lmax, nbepo, ntyele, nomele, ivr, nbocc)
!
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES ORIENTATIONS
!
! --------------------------------------------------------------------------------------------------
!
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : NOMO   : NOM DU MODELE
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use cara_elem_parameter_module
    implicit none
    integer(kind=8) :: lmax, nbepo, ntyele(*), ivr(*), nbocc(*)
    character(len=8) :: noma, nomo
    character(len=16) :: nomele(*)
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8rddg.h"
#include "asterfort/aceatu.h"
#include "asterfort/affori.h"
#include "asterfort/alcart.h"
#include "asterfort/angvx.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/tbcarapou.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ii, ifm, ioc, ixma
    integer(kind=8) :: jj, jad, jin, jdcmpo, jdco, jdgm, nbid
    integer(kind=8) :: jdls, jdme, jdno, jdori, jdtm, jinit
    integer(kind=8) :: jdvlvo, nbmagr, nbmail
    integer(kind=8) :: nbval, ncar, ng
    integer(kind=8) :: no1, no2, nocaor, ntpoi, ntseg, ntseg3, ntseg4
    integer(kind=8) :: nummai, nutyel, nutyma, nbalarme
    integer(kind=8) :: nval
    parameter(nbval=6)
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: iret
    real(kind=8)        :: valcara(2)
    integer(kind=8)             :: okcara(2)
    character(len=8)    :: nomsec, nomcara(2)
    character(len=19)   :: tabcar
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: val(nbval), x1(3), x2(3), x3(3), longseg, longseuil
    real(kind=8) :: rddg, alpha, AlphaGeom(2), beta, gamma
    character(len=4) :: exituy
    character(len=8) :: nomu
    character(len=16) :: oricara
    character(len=16) :: concep, cmd, nunomel
    character(len=19) :: cartor, ligrel
    character(len=24) :: tmpnor, tmpvor, tmpori, tmpini
    character(len=24) :: mlgtma, mlggno, mlggma, mlgcoo, mlgcnx
    character(len=24) :: modmai, nommai
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    rddg = r8rddg()
    call getres(nomu, concep, cmd)
    tmpori = nomu//'.ORIENTATION'
    tmpini = nomu//'.ORIENTINIT'
! --------------------------------------------------------------------------------------------------
!   CONSTRUCTION DES CARTES
    cartor = nomu//'.CARORIEN'
    tmpnor = cartor//'.NCMP'
    tmpvor = cartor//'.VALV'
!
!   RECONSTRUCTION DES NOMS JEVEUX DU CONCEPT MODELE
    call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=ligrel)
    modmai = ligrel//'.TYFE'
!
!   RECONSTRUCTION DES NOMS JEVEUX DU CONCEPT MAILLAGE ASSOCIE
    mlgtma = noma//'.TYPMAIL'
    mlgcnx = noma//'.CONNEX'
    mlggno = noma//'.GROUPENO'
    mlggma = noma//'.GROUPEMA'
    mlgcoo = noma//'.COORDO    .VALE'
!
    call jelira(mlgtma, 'LONMAX', nbmail)
    call jeexin(modmai, ixma)
    if (ixma .ne. 0) call jeveuo(modmai, 'L', jdme)
! --------------------------------------------------------------------------------------------------
!   Récupération des adresses jeveux utiles
    call jeveuo(mlgtma, 'L', jdtm)
    call jeveuo(mlgcoo, 'L', jdco)
!
!   Récupération des numéros des types mailles poi1/seg2
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ntpoi)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), ntseg)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG3'), ntseg3)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG4'), ntseg4)
!
!
    call wkvect('&&TMPORIEN', 'V V K24', lmax, jdls)
    call wkvect(tmpori, 'V V R', nbmail*3, jdori)
    call wkvect(tmpini, 'V V I', nbmail*3, jinit)
!
!   Initialisation des angles nautiques sur toutes les mailles non nulles. Repère local par défaut
    do ii = 1, nbmail*3
        zr(jdori+ii-1) = 0.d0
        zi(jinit+ii-1) = 0
    end do
!
    do nummai = 1, nbmail
        nutyma = zi(jdtm+nummai-1)
        jad = jdori+(nummai-1)*3
        jin = jinit+(nummai-1)*3
        if (nutyma .eq. ntseg) then
            call jeveuo(jexnum(mlgcnx, nummai), 'L', jdno)
            no1 = zi(jdno)
            no2 = zi(jdno+1)
            do ii = 1, 3
                x1(ii) = zr(jdco+(no1-1)*3+ii-1)
                x2(ii) = zr(jdco+(no2-1)*3+ii-1)
                x3(ii) = x2(ii)-x1(ii)
            end do
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            longseg = sqrt(ddot(b_n, x3, b_incx, x3, b_incy))
            if (longseg .gt. 0.0d0) then
                call angvx(x3, alpha, beta)
                zr(jad) = zr(jad)+alpha
                zr(jad+1) = zr(jad+1)+beta
                zi(jin) = 1
                zi(jin+1) = 1
            end if
        end if
    end do
!
! --------------------------------------------------------------------------------------------------
!   Affectation des valeurs lues dans l'objet tampon
    if (nbocc(ACE_ORIENTATION) .ne. 0) then
        nbalarme = 0
        do ioc = 1, nbocc(ACE_ORIENTATION)
            ! Pour les MAILLES
            call getvem(noma, 'GROUP_MA', 'ORIENTATION', 'GROUP_MA', ioc, lmax, zk24(jdls), ng)
            ! Seuil correspondant à la longueur nulle pour une maille :
            !   si seglong .LT. longseuil ==> maille de taille nulle
            call getvr8('ORIENTATION', 'PRECISION', iocc=ioc, scal=longseuil, nbret=nbid)
            if (nbid .ne. 1) longseuil = -1.0d0
            call getvtx('ORIENTATION', 'CARA', iocc=ioc, scal=oricara, nbret=ncar)
            call getvr8('ORIENTATION', 'VALE', iocc=ioc, nbval=nbval, vect=val, nbret=nval)
            ! Dans le catalogue, si oricara == GENE_TUYAU c'est GROUP_NO ou NOEUD ==> Tuyaux
            if (ng .gt. 0) then
                AlphaGeom(:) = 0.0d0
                ! Si oricara = VECT_MAIL_Y VECT_MAIL_Z, on va lire TABLE_CARA, NOM_SEC
                if ((oricara .eq. "VECT_MAIL_Y") .or. (oricara .eq. "VECT_MAIL_Z")) then
                    ! Obligatoire dans le catalogue : TABLE_CARA et NOM_SEC
                    call getvid('ORIENTATION', 'TABLE_CARA', iocc=ioc, scal=tabcar, nbret=iret)
                    call getvtx('ORIENTATION', 'NOM_SEC', iocc=ioc, scal=nomsec, nbret=iret)
                    !
                    nomcara(1) = 'ALPHA'
                    call tbcarapou(tabcar, nomsec, 1, nomcara, valcara, okcara)
                    AlphaGeom(1) = valcara(1)
                    !
                end if
                !
                ! GROUP_MA = toutes les mailles possibles de la liste des groupes de mailles
                do ii = 1, ng
                    call jeveuo(jexnom(mlggma, zk24(jdls+ii-1)), 'L', jdgm)
                    call jelira(jexnom(mlggma, zk24(jdls+ii-1)), 'LONUTI', nbmagr)
                    do jj = 1, nbmagr
                        nummai = zi(jdgm+jj-1)
                        nommai = int_to_char8(nummai)
                        call jeveuo(jexnum(mlgcnx, nummai), 'L', jdno)
                        nutyma = zi(jdtm+nummai-1)
                        jad = jdori+(nummai-1)*3
                        jin = jinit+(nummai-1)*3
                        if ((nutyma .ne. ntseg3) .and. (nutyma .ne. ntseg4)) then
                            call affori('MAILLE', nommai, oricara, val, jad, jin, &
                                        jdno, jdco, nutyma, ntseg, &
                                        lseuil=longseuil, nbseuil=nbalarme, alphayz=AlphaGeom)
                        end if
                    end do
                end do
            end if
        end do
        if (nbalarme .gt. 0) call utmess('A', 'MODELISA_95', si=nbalarme, sr=longseuil)
    end if
!
! --------------------------------------------------------------------------------------------------
!   Impression des valeurs des orientations si demandé
    nocaor = 0
    ifm = ivr(4)
    if (ivr(3) .eq. 2) write (ifm, 100)
    cnum1: do nummai = 1, nbmail
        nutyel = zi(jdme+nummai-1)
        do jj = 1, ACE_NB_TYPE_ELEM
            if (nutyel .eq. ntyele(jj)) then
                nocaor = nocaor+1
                if (ivr(3) .eq. 2) then
                    nommai = int_to_char8(nummai)
                    jad = jdori+(nummai-1)*3
                    alpha = rddg*zr(jad)
                    beta = rddg*zr(jad+1)
                    gamma = rddg*zr(jad+2)
                    write (ifm, 110) nommai, nomele(jj), alpha, beta, gamma
                end if
                cycle cnum1
            end if
        end do
    end do cnum1
!
100 format(/, 3x, '<ANGL> ORIENTATIONS SUR LES MAILLES DE TYPE POUTRE BARRE OU DISCRET', //, 3x, &
            'NOM      TYPE             ALPHA         BETA          GAMMA')
110 format(3x, a8, 1x, a16, 1x, 3(1pd13.6, 2x))
!
! --------------------------------------------------------------------------------------------------
!   Affectation des valeurs du tampon dans la carte orientation :
    if (nocaor .gt. 0) then
        call alcart('G', cartor, noma, 'CAORIE_R')
        call jeveuo(tmpnor, 'E', jdcmpo)
        call jeveuo(tmpvor, 'E', jdvlvo)
        zk8(jdcmpo) = 'ALPHA'
        zk8(jdcmpo+1) = 'BETA'
        zk8(jdcmpo+2) = 'GAMMA'
!
!       Affectation des mailles du maillage (poutre, barre ou discret)
        cnum2: do nummai = 1, nbmail
            nutyel = zi(jdme+nummai-1)
            do jj = 1, ACE_NB_TYPE_ELEM
                if (nutyel .eq. ntyele(jj)) then
!                   Récupération des numéros des noms des éléments
                    call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), nunomel)
                    ! Pas de carte d'orientation sur les :
                    !   TUYAUX                  : MET3SEG3 MET3SEG4 MET6SEG3
                    !   "meca_plate_skin"       : MEBODKT  MEBODST  MEBOQ4G
                    !   "meca_coque_3d_skin"    : MEBOCQ3
                    if ((nunomel .ne. 'MET3SEG3') .and. (nunomel .ne. 'MET3SEG4') .and. &
                        (nunomel .ne. 'MET6SEG3') .and. &
                        (nunomel .ne. 'MEBODKT') .and. (nunomel .ne. 'MEBODST') .and. &
                        (nunomel .ne. 'MEBOQ4G') .and. &
                        (nunomel .ne. 'MEBOCQ3')) then
                        jad = jdori+(nummai-1)*3
                        zr(jdvlvo) = zr(jad)
                        zr(jdvlvo+1) = zr(jad+1)
                        zr(jdvlvo+2) = zr(jad+2)
                        call nocart(cartor, 3, 3, mode='NUM', nma=1, limanu=[nummai])
                        cycle cnum2
                    end if
                end if
            end do
        end do cnum2
    end if
!
! --------------------------------------------------------------------------------------------------
!   Affectation des elements tuyaux
    call dismoi('EXI_TUYAU', nomo, 'MODELE', repk=exituy)
    if (exituy .eq. 'OUI') then
        call aceatu(noma, nomo, nbepo, ntyele, ivr, nbocc)
    end if
! --------------------------------------------------------------------------------------------------
!
    call jedetr('&&TMPORIEN')
    call jedetr(tmpnor)
    call jedetr(tmpvor)
    call jedetr(tmpori)
    call jedetr(tmpini)
!
    call jedema()
end subroutine
