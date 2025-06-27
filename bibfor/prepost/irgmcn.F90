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
subroutine irgmcn(chamsy, partie, ifi, nomcon, ordr, &
                  nbordr, coord, connx, point, nobj, &
                  nbel, nbcmpi, nomcmp, lresu, para, &
                  versio, tycha)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cnocns.h"
#include "asterfort/codent.h"
#include "asterfort/irgmor.h"
#include "asterfort/irgmpv.h"
#include "asterfort/irgnal.h"
#include "asterfort/irgnte.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ifi, nbordr, nbcmpi, versio
    integer(kind=8) :: ordr(*), connx(*), point(*)
    real(kind=8) :: coord(*), para(*)
    aster_logical :: lresu
    character(len=*) :: nomcon, chamsy, nomcmp(*), partie
!     NBRE, NOM D'OBJET POUR CHAQUE TYPE D'ELEMENT
    integer(kind=8) :: neletr
    parameter(neletr=8)
    integer(kind=8) :: tord(neletr)
    integer(kind=8) :: nbel(*)
    character(len=24) :: nobj(*)
    character(len=8) :: tycha
!
!        IMPRESSION D'UN CHAM_NO AU FORMAT GMSH
!
!        CHAMSY : NOM SYMBOLIQUE DU CHAM_NO A ECRIRE
!        IFI    : NUMERO D'UNITE LOGIQUE DU FICHIER DE SORTIE GMSH
!        NOMCON : NOM DU CONCEPT A IMPRIMER
!        PARTIE : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
!        ORDR   : LISTE DES NUMEROS D'ORDRE A IMPRIMER
!        NBORDR : NOMBRE DE NUMEROS D'ORDRE DANS LE TABLEAU ORDR
!        COORD  : VECTEUR COORDONNEES DES NOEUDS DU MAILLAGE
!        CONNX  : VECTEUR CONNECTIVITES DES NOEUDS DU MAILLAGE
!        POINT  : VECTEUR DU NOMBRE DE NOEUDS DES MAILLES DU MAILLAGE
!        NOBJ(i): NOM JEVEUX DEFINISSANT LES ELEMENTS DU MAILLAGE
!        NBEL(i): NOMBRE D'ELEMENTS DU MAILLAGE DE TYPE i
!        NBCMPI : NOMBRE DE COMPOSANTES DEMANDEES A IMPRIMER
!        NOMCMP : NOMS DES COMPOSANTES DEMANDEES A IMPRIMER
!        LRESU  : LOGIQUE INDIQUANT SI NOMCON EST UNE SD RESULTAT
!        PARA   : VALEURS DES VARIABLES D'ACCES (INST, FREQ)
!        TYCHA  : TYPE DE CHAMP A IMPRIMER (VERSION >= 1.2)
!                 = SCALAIRE/VECT_2D/VECT_3D/TENS_2D/TENS_3D
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: i, ine
    integer(kind=8) :: ior, k, ncmp, iret, nbord2, ncmpu
    integer(kind=8) :: jcnsk, jtype
    aster_logical :: scal, vect, tens
    character(len=8) :: k8b, nocmp, tbcmp(3)
    character(len=19) :: noch19, champs
    integer(kind=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    integer(kind=8), pointer :: cnsl(:) => null()
    integer(kind=8), pointer :: cnsv(:) => null()
    character(len=8), pointer :: vnocmp(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- ORDRE D'IMPRESSION DES VALEURS
    call irgmor(tord, versio)
!
    nbord2 = max(1, nbordr)
!
    AS_ALLOCATE(vi=cnsd, size=nbord2)
    AS_ALLOCATE(vi=cnsc, size=nbord2)
    AS_ALLOCATE(vi=cnsv, size=nbord2)
    AS_ALLOCATE(vi=cnsl, size=nbord2)
    call wkvect('&&IRGMCN.TYPE', 'V V K8', nbord2, jtype)
!
!
    do ior = 1, nbord2
        if (lresu) then
            call rsexch(' ', nomcon, chamsy, ordr(ior), noch19, &
                        iret)
            if (iret .ne. 0) goto 100
        else
            noch19 = nomcon
        end if
        call codent(ior, 'D0', k8b)
        champs = '&&IRGMCN.CH'//k8b
        call cnocns(noch19, 'V', champs)
        call jeveuo(champs//'.CNSK', 'L', jcnsk)
        call jeveuo(champs//'.CNSD', 'L', cnsd(ior))
        call jeveuo(champs//'.CNSC', 'L', cnsc(ior))
        call jeveuo(champs//'.CNSV', 'L', cnsv(ior))
        call jeveuo(champs//'.CNSL', 'L', cnsl(ior))
        call jelira(champs//'.CNSV', 'TYPE', cval=zk8(jtype+ior-1))
!
!
100     continue
    end do
!
! --- RECUPERATION DES COMPOSANTES POUR L'IMPRESSION
!     D'UN CHAMP SCALAIRE PAR COMPOSANTE
!
    ncmp = zi(cnsd(1)-1+2)
    ncmpu = 0
    AS_ALLOCATE(vk8=vnocmp, size=ncmp)
    if (nbcmpi .eq. 0) then
        do k = 1, ncmp
            nocmp = zk8(cnsc(1)-1+k)
            ncmpu = ncmpu+1
            vnocmp(ncmpu) = nocmp
        end do
    else
        do k = 1, nbcmpi
            ncmpu = ncmpu+1
            vnocmp(ncmpu) = nomcmp(k)
        end do
    end if
!
! -- VERSION GMSH = 1.0 :
!    LA DETERMINATION DU TYPE DE CHAMP A IMPRIMER
!    EST FONCTION DES COMPOSANTES FOURNIES OU TROUVEES
!     1/ ON RECHERCHE LES COMPOSANTES DX, DY, DZ
!        ==> IMPRESSION D'UN CHAMP VECTORIEL
!     2/ POUR LES AUTRES COMPOSANTES
!        ==> IMPRESSION D'UN CHAMP SCALAIRE PAR COMPOSANTE
! -- VERSION GMSH = 1.2 :
!    ON UTILISE TYCHA POUR DETERMINER LE TYPE DE CHAMP A IMPRIMER
!    SI NOMCMP ABSENT : => TYCHA='SCALAIRE'
!    SI NOMCMP PRESENT: => TYCHA='SCALAIRE'/'VECT_xD'/'TENS_xD'
!
    scal = .false.
    vect = .false.
    tens = .false.
!
    if (versio .eq. 1) then
        ncmp = zi(cnsd(1)-1+2)
        if (nbcmpi .eq. 0) then
            do k = 1, ncmp
                nocmp = zk8(cnsc(1)-1+k)
                if (nocmp .eq. 'DX' .or. nocmp .eq. 'DY' .or. nocmp .eq. 'DZ') then
                    vect = .true.
                else
                    scal = .true.
                end if
            end do
        else
            scal = .true.
        end if
    else if (versio .ge. 2) then
        if (tycha(1:4) .eq. 'SCAL') then
            scal = .true.
        else if (tycha(1:4) .eq. 'VECT') then
            vect = .true.
        else if (tycha(1:4) .eq. 'TENS') then
            tens = .true.
        end if
    end if
!
! ----------------------------------------------------------------------
!          IMPRESSION D'UN CHAMP TENSORIEL
! ----------------------------------------------------------------------
    if (tens) then
!
!        ECRITURE DE L'ENTETE DE View
!        ****************************
        nocmp = 'TENSEUR '
        call irgmpv(ifi, lresu, nomcon, chamsy, nbord2, &
                    para, nocmp, nbel, .false._1, .false._1, &
                    tens, versio)
!
! ---    BOUCLE SUR LES TYPES D'ELEMENTS SI NBEL>0
!        ON A RECUPERE L'ORDRE D'IMPRESSION PAR IRGMOR
        do ine = 1, neletr
            i = tord(ine)
            if (nbel(i) .ne. 0) then
                call irgnte(ifi, nbord2, coord, connx, point, &
                            nobj(i), nbel(i), cnsv, partie, jtype, &
                            cnsd)
            end if
        end do
!
!        FIN D'ECRITURE DE View
!        **********************
        write (ifi, 1000) '$EndView'
!
    end if
!
!
! ----------------------------------------------------------------------
!          IMPRESSION D'UN CHAMP VECTORIEL ( CMP = DX, DY, DZ )
! ----------------------------------------------------------------------
    if (vect) then
!
!        ECRITURE DE L'ENTETE DE View
!        ****************************
!
        nocmp = 'VECTEUR '
        call irgmpv(ifi, lresu, nomcon, chamsy, nbord2, &
                    para, nocmp, nbel, .false._1, vect, &
                    tens, versio)
!
!        LISTE DES COMPOSANTES
        if (versio .eq. 1) then
            tbcmp(1) = 'DX      '
            tbcmp(2) = 'DY      '
            tbcmp(3) = 'DZ      '
        else if (versio .eq. 2) then
            tbcmp(3) = '        '
            do i = 1, nbcmpi
                tbcmp(i) = vnocmp(i)
            end do
        end if
!
! ---    BOUCLE SUR LES TYPES D'ELEMENTS SI NBEL>0
!        ON A RECUPERE L'ORDRE D'IMPRESSION PAR IRGMOR
        do ine = 1, neletr
            i = tord(ine)
            if (nbel(i) .ne. 0) then
                call irgnal(ifi, nbord2, coord, connx, point, &
                            tbcmp, 3, i, nobj(i), nbel(i), &
                            cnsc, cnsl, cnsv, partie, jtype, &
                            cnsd)
            end if
        end do
!
!        FIN D'ECRITURE DE View
!        **********************
!
        write (ifi, 1000) '$EndView'
!
    end if
!
! ----------------------------------------------------------------------
!           IMPRESSION D'UN CHAMP SCALAIRE ( AUTRE CMP )
! ----------------------------------------------------------------------
!
    if (scal) then
        do k = 1, ncmpu
            nocmp = vnocmp(k)
!
!        ECRITURE DE L'ENTETE DE View
!        ****************************
!
            call irgmpv(ifi, lresu, nomcon, chamsy, nbord2, &
                        para, nocmp, nbel, scal, .false._1, &
                        tens, versio)
!
!        LISTE DES COMPOSANTES
            tbcmp(1) = nocmp
!
! ---    BOUCLE SUR LES TYPES D'ELEMENTS SI NBEL>0
!        ON A RECUPERE L'ORDRE D'IMPRESSION PAR IRGMOR
            do ine = 1, neletr
                i = tord(ine)
                if (nbel(i) .ne. 0) then
                    call irgnal(ifi, nbord2, coord, connx, point, &
                                tbcmp, 1, i, nobj(i), nbel(i), &
                                cnsc, cnsl, cnsv, partie, jtype, &
                                cnsd)
                end if
            end do
!
!        FIN D'ECRITURE DE View
!        **********************
!
            write (ifi, 1000) '$EndView'
!
        end do
    end if
!
    AS_DEALLOCATE(vi=cnsd)
    AS_DEALLOCATE(vi=cnsc)
    AS_DEALLOCATE(vi=cnsv)
    AS_DEALLOCATE(vi=cnsl)
    AS_DEALLOCATE(vk8=vnocmp)
    call jedetr('&&IRGMCN.TYPE')
    call jedema()
!
1000 format(a8)
!
end subroutine
