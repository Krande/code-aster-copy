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

subroutine moco99(nomres, resul, nbmod, lrang, iorne,&
                  seul)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
!
    integer :: nbmod, lrang(nbmod), iorne
    character(len=8) :: nomres, resul
    aster_logical :: seul
!
!     BUT:
!       POINTER LES PREMIERS MODES PROPRES D'UNE STRUCTURE RESULTAT
!       DE TYPE MODE_MECA DANS UNE AUTRE STRUCTURE DE TYPE BASE_MODALE
!       DEJA EXISTANTE A PARTIR D'UN NUMERO D'ORDRE
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NOMRES    : NOM UTILISATEUR DE LA STRUCTURE RESULTAT A REMPLIR
! IN   RESUL     : NOM DE LA STRUCTURE RESULTAT A POINTER
! IN   NBMOD     : NOMBRE DE MODES A POINTER
! IN   LRANG     : LISTE DES ANCIENS NUMEROS DE RANGEMENT A POINTER
! IN   SEUL      : INDICATEUR POUR UNE DEUXIEME OCCURENCE.
!                  MODE_MECA OU MODE_STAT ETC. SOUS LE MOT CLE "RITZ"
!
!      SORTIE :
!-------------
! OUT  IORNE     : NUMERO D'ORDRE DU PREMIER CHAMPS 'DEPL' A POINTER
!
! ......................................................................
!
!
!
!
    integer :: nbpabm
    parameter   (nbpabm=10)
    integer :: ldpar(nbpabm), ldpa2(nbpabm)
    integer :: nbcham, nbold(1), nbtrou, vali
    integer :: i, ii, jtyp, ier, iorol, ire, ibid, lpain(3),lpaou(3)
    integer :: llkge, llmge, llncp, llom2, lltmo, llval2, llvalo
!
    real(kind=8) :: genek, genem, omeg2, rbid, epsi
!
    character(len=8) :: k8bid, interf, typi, param(3)
    character(len=16) :: typres, bmpara(nbpabm), chmeca, typmo
    character(len=19) :: chamol, chamne
    character(len=24) :: type, typeba
!
    complex(kind=8) :: cbid
    integer, pointer :: ordr(:) => null()
!
!-----------------------------------------------------------------------
!
    data  bmpara /'NUME_MODE','FREQ','NORME','NOEUD_CMP','TYPE_DEFO',&
     &              'OMEGA2','MASS_GENE','RIGI_GENE','TYPE_MODE',&
     &              'AMOR_REDUIT'/
!
!-----------------------------------------------------------------------
    param(1) = 'MODELE'
    param(2) = 'CHAMPMAT'
    param(3) = 'CARAELEM'

!
! --- CAS DE L'ABSENCE D'UN MODE_MECA
!
    if (resul .eq. '          ' .or. nbmod .eq. 0) goto 999
!
!
! --- DETERMINATION DU NOMBRE DE MODES DANS LA STRUCTURE A POINTER
!
    call rsorac(resul, 'LONUTI', ibid, rbid, k8bid,&
                cbid, epsi, 'ABSOLU', nbold, 1,&
                nbtrou)
!
    if (nbmod .gt. nbold(1)) then
        vali = nbold(1)
        call utmess('I', 'ALGORITH13_48', si=vali)
        nbmod=nbold(1)
    endif
    nbmod=min(nbmod,nbold(1))
!
    if (nbmod .eq. 0) goto 999
!     --- ON RECUPERE LE TYPE D'INTERFACE ---
!
    call dismoi('REF_INTD_PREM', nomres, 'RESU_DYNA', repk=interf, arret='C',&
                ier=ire)
    if (interf .ne. ' ') then
        type = interf//'.IDC_TYPE'
        call jeveuo(type, 'L', jtyp)
        typi = zk8(jtyp)
    endif
!
! --- RECHERCHE DE L'ADRESSE DES ANCIENNES VALEURS PROPRES
!
    call dismoi('TYPE_RESU', resul, 'RESULTAT', repk=typres)
!
    typeba=' '
    if (typres .ne. 'MULT_ELAS') then
        call dismoi('TYPE_BASE', resul, 'RESU_DYNA', repk=typeba, arret='C',&
                    ier=ire)
    endif
!
    typmo=' '
!
!
! ----- RECUPERATION DU NOMBRE DE CHAMP POSSIBLE DE LA SD
    call jelira(resul//'           .DESC', 'NOMMAX', nbcham)
!
    call jeveuo(resul//'           .ORDR', 'L', vi=ordr)
    do i = 1, nbmod
!       LRANG : CONTIENT DES NUMEROS DE RANGEMENT : 1,2, ..., NBMOD
        ASSERT(lrang(i).eq.i)
        iorol=ordr(lrang(i))
!
! ------BOUCLE SUR LA LISTE DE CHAMPS AFIN D'IDENTIFIER CEUX QUI
!       SONT PRESENTES DANS L'ANCIEN RESULT
!
        do ii = 1, nbcham
! ----- RECUPERATION DU NOM DU CHAMP POSSIBLE DE LA SD
            call jenuno(jexnum(resul//'           .DESC', ii), chmeca)
!
! ----- REQUETE NOM ET ADRESSE ANCIEN CHAMNO
            call rsexch(' ', resul, chmeca, iorol, chamol,&
                        ier)
            if (ier .eq. 0) then
                call rsexch(' ', nomres, chmeca, iorne, chamne,&
                            ier)
                if (chamol .ne. chamne) then
                    call copisd('CHAMP', 'G', chamol, chamne)
                endif
                call rsnoch(nomres, chmeca, iorne)
            endif
        end do
!
        if (typres .ne. 'MODE_MECA') goto 11
!
! ----- RECUPERATION DES VALEURS GENERALISEES ET PULSATION CARREE
!
        call rsadpa(resul, 'L', 1, 'RIGI_GENE', iorol,&
                    0, sjv=llkge, styp=k8bid, istop=0)
        genek=zr(llkge)
        call rsadpa(resul, 'L', 1, 'MASS_GENE', iorol,&
                    0, sjv=llmge, styp=k8bid, istop=0)
        genem=zr(llmge)
        call rsadpa(resul, 'L', 1, 'OMEGA2', iorol,&
                    0, sjv=llom2, styp=k8bid, istop=0)
        omeg2=zr(llom2)
!
        call rsadpa(resul, 'L', 1, 'TYPE_MODE', iorol,&
                    0, sjv=lltmo, styp=k8bid, istop=0)
        typmo=zk16(lltmo)

!       RECUPERATION ET ECRITURE DES MODELES, CHAMPMAT ET CARA_ELEM

        call rsadpa(resul, 'L', 3, param, iorol,&
                    0, tjv=lpain, styp=k8bid, istop=0)

        call rsadpa(nomres, 'E', 3, param, iorol,&
                    0, tjv=lpaou, styp=k8bid, istop=0)
        zk8(lpaou(1)) = zk8(lpain(1))
        zk8(lpaou(2)) = zk8(lpain(2))
        zk8(lpaou(3)) = zk8(lpain(3))
!
 11     continue
!
! ----- ECRITURE DES NOUVEAUX PARAMETRES
!
        call rsadpa(nomres, 'E', nbpabm, bmpara, iorne,&
                    0, tjv=ldpar, styp=k8bid)
        zi(ldpar(1))=iorne
        if (typmo(1:8) .ne. 'MODE_DYN') then
            if (typeba(1:1) .ne. ' ') then
                call rsadpa(resul, 'L', nbpabm, bmpara, iorol,&
                            0, tjv=ldpa2, styp=k8bid, istop=0)
                zr(ldpar(2)) = zr(ldpa2(2))
                zk24(ldpar(3)) = zk24(ldpa2(3))
                zk16(ldpar(4)) = zk16(ldpa2(4))
                zk16(ldpar(5)) = zk16(ldpa2(5))
                zr(ldpar(6)) = zr(ldpa2(6))
                zr(ldpar(7)) = zr(ldpa2(7))
                zr(ldpar(8)) = zr(ldpa2(8))
                if (.not.seul) then
                    zk16(ldpar(9)) ='MODE_STA'
                else if (seul) then
                    zk16(ldpar(9)) = zk16(ldpa2(9))
                endif
                goto 12
            endif
            zr(ldpar(2))=0.d0
            zk24(ldpar(3))=' '
            zk16(ldpar(4))=' '
            zk16(ldpar(5))=' '
            if (typmo(1:8) .eq. 'MODE_STA') then
                call rsadpa(resul, 'L', 1, 'NOEUD_CMP', iorol,&
                            0, sjv=llncp, styp=k8bid, istop=0)
                zk16(ldpar(4)) = zk16(llncp)
                zk16(ldpar(5)) = 'STATIQUE'
                if (interf .ne. ' ') then
                    if (typi .eq. 'CRAIGB') zk16(ldpar(5)) = 'CONTRAINT'
                    if (typi .eq. 'MNEAL') zk16(ldpar(5)) = 'ATTACHE'
                endif
            endif
            zr(ldpar(6))=0.d0
            zr(ldpar(7))=0.d0
            zr(ldpar(8))=0.d0
! ON ETEND LA DECLARATION MODE_STA A TOUS LES MODES NON DYNA PAR
! EXEMPLE LES MULT_ELAS POUR ETRE COMPTABILISES COMME MODES STATIQUES
! DANS LES BASES DE RITZ PAR APPEL A DISMOI
            zk16(ldpar(9))='MODE_STA'
            goto 12
        endif
        call rsadpa(resul, 'L', 1, 'FREQ', iorol,&
                    0, sjv=llvalo, styp=k8bid, istop=0)
        call rsadpa(resul, 'L', 1, 'AMOR_REDUIT', iorol,&
                    0, sjv=llval2, styp=k8bid, istop=0)
        zr(ldpar(2))=zr(llvalo)
        zk24(ldpar(3))=' '
        zk16(ldpar(4))=' '
        zk16(ldpar(5))='PROPRE '
        zr(ldpar(6))=omeg2
        zr(ldpar(7))=genem
        zr(ldpar(8))=genek
        if (.not.seul) then
            zk16(ldpar(9)) ='MODE_STA'
        else if (seul) then
            zk16(ldpar(9)) = 'MODE_DYN'
        endif
        zr(ldpar(10))=zr(llval2)
!
 12     continue
!
! ----- INCREMENTATION DU NUMERO ORDRE
!
        iorne=iorne+1
!
    end do
!
!
999 continue
end subroutine
