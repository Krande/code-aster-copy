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

subroutine crlidd(nomres, mailla)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 19/02/91
!-----------------------------------------------------------------------
!  BUT:  CREER LE TABLEAU DESCRIPTEUR DES DEFORMEES A CALCULER
!
!  REMPLIR LA PREMIERE LIGNE PAR LES NUMERO (MAILLAGE) DES NOEUDS
!
!
!  REMPLIR LA DEUXIEME LIGNE PAR UN CODE DE TYPE D'INTERFACE
!   (-1 MAC NEAL) (-2 CRAIG BAMPTON)(-3  CB_HARMO)(-4 AUCUN)
!  CETTE COLONNE CONTIENDRA PLUS TARD LE NUMERO DE LA PREMIERE
!    DEFORMEE STATIQUE DU NOEUDS
!
! REMPLIR LA TROISIEME COLONNE  PAR LE  DU CUMUL DES
!       MASQUES AU NOEUDS (DDL SUR LESQUELS SERONT CALCULEES
!      LES DEFORMES A PRIORI)
!  CET ENTIER SERA REACTUALISE PLUS TARD EN FONCTION DES DDL
!    REELLEMENT ACTIFS AUX NOEUDS
!
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTLISATEUR DU CONCEPT RESULTAT DE L'OPERATEUR
! MAILLA /I/: NOM UTLISATEUR DU MAILLAGE DE LA SOUS-STRUCTURE
!
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cheris.h"
#include "asterfort/copvis.h"
#include "asterfort/dismoi.h"
#include "asterfort/isgeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nodoub.h"
#include "asterfort/utmess.h"
#include "asterfort/uttrii.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, mailla, cb, mn, ha, au, nomtyp
    character(len=24) :: temmn, temcb, temha, temau, temmas
    character(len=24) :: typint, noeint, desdef
    aster_logical :: doubok
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, ino, iran, j, lldes
    integer(kind=8) :: llnin, lltyp, ltau, ltcb, ltha, ltmas, ltmn
    integer(kind=8) :: nbau, nbcb, nbcmp, nbec, nbha, nbint, nbmn
    integer(kind=8) :: nbno, nbtem, nbto
!-----------------------------------------------------------------------
    data cb, mn, ha, au/'CRAIGB', 'MNEAL', 'CB_HARMO', 'AUCUN'/
    data doubok/.false./
!-----------------------------------------------------------------------
!
!-----------------INITIALISATION DES NOM LES PLUS UTILISES--------------
!
    call jemarq()
    typint = nomres//'.IDC_TYPE'
    noeint = nomres//'.IDC_LINO'
!
!----------------RECUPERATION DU NOMBRE D'ENTIERS CODES DE LA GRANDEUR-
!
    call dismoi('NB_EC', nomres, 'INTERF_DYNA', repi=nbec)
!
!----------------------RECUPERATION TYPE INTERFACE----------------------
!
    call jeveuo(typint, 'L', lltyp)
!
!--------------------RECHERCHE DU NOMBRE D'INTERFACES-------------------
!
    call jelira(noeint, 'NMAXOC', nbint)
!
!------- BOUCLE DE COMPTAGE DES NOEUDS PAR TYPE INTERFACES--------------
!
    nbcb = 0
    nbmn = 0
    nbau = 0
    nbha = 0
    do i = 1, nbint
        nomtyp = zk8(lltyp+i-1)
!
!    COMPTEUR CRAIG BAMPTON
!
        if (nomtyp .eq. cb) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            nbcb = nbcb+nbno
        end if
!
!    COMPTEUR MAC NEAL
!
        if (nomtyp .eq. mn) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            nbmn = nbmn+nbno
        end if
!
!    COMPTEUR CRAIG-BAMPTON-HARMONIQUE
!
        if (nomtyp .eq. ha) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            nbha = nbha+nbno
        end if
!
!    COMPTEUR AUCUN
!
        if (nomtyp .eq. au) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            nbau = nbau+nbno
        end if
!
    end do
!
!---------ALLOCATION DU VECTEUR DES NOEUDS CRAIG BAMPTON----------------
!
    temcb = '&&CRLIDD.NOE.CB'
    if (nbcb .gt. 0) then
        call wkvect(temcb, 'V V I', nbcb, ltcb)
    else
        ltcb = 1
    end if
!
!---------ALLOCATION DU VECTEUR DES NOEUDS MAC NEAL---------------------
!
    temmn = '&&CRLIDD.NOE.MN'
    if (nbmn .gt. 0) then
        call wkvect(temmn, 'V V I', nbmn, ltmn)
    else
        ltmn = 1
    end if
!
!---------ALLOCATION DU VECTEUR DES NOEUDS CB_HARMO---------------------
!
    temha = '&&CRLIDD.NOE.HA'
    if (nbha .gt. 0) then
        call wkvect(temha, 'V V I', nbha, ltha)
    else
        ltha = 1
    end if
!
!------------ALLOCATION DU VECTEUR DES NOEUDS AUCUN---------------------
!
    temau = '&&CRLIDD.NOE.AU'
    if (nbau .gt. 0) then
        call wkvect(temau, 'V V I', nbau, ltau)
    else
        ltau = 1
    end if
!
!---------BOUCLE DE REMPLISSAGE DES 3 VECTEURS TEMCB TEMMN TEMAU--------
!
    nbcb = 0
    nbmn = 0
    nbau = 0
    nbha = 0
!
    do i = 1, nbint
        nomtyp = zk8(lltyp+i-1)
!
!    NOEUDS DE CRAIG BAMPTON
!
        if (nomtyp .eq. cb) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            call jeveuo(jexnum(noeint, i), 'L', llnin)
            do j = 1, nbno
                zi(ltcb+nbcb) = zi(llnin+j-1)
                nbcb = nbcb+1
            end do
        end if
!
!    NOEUD DE MAC NEAL
!
        if (nomtyp .eq. mn) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            call jeveuo(jexnum(noeint, i), 'L', llnin)
            do j = 1, nbno
                zi(ltmn+nbmn) = zi(llnin+j-1)
                nbmn = nbmn+1
            end do
        end if
!
!    NOEUD DE CB_HARMO
!
        if (nomtyp .eq. ha) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            call jeveuo(jexnum(noeint, i), 'L', llnin)
            do j = 1, nbno
                zi(ltha+nbha) = zi(llnin+j-1)
                nbha = nbha+1
            end do
        end if
!
!    NOEUD DE AUCUN
!
        if (nomtyp .eq. au) then
            call jelira(jexnum(noeint, i), 'LONMAX', nbno)
            call jeveuo(jexnum(noeint, i), 'L', llnin)
            do j = 1, nbno
                zi(ltau+nbau) = zi(llnin+j-1)
                nbau = nbau+1
            end do
        end if
!
    end do
!
!-------SUPPRESSION DES DOUBLES ET ORDRE DES LISTES TROUVEES------------
!
    if (nbcb .ne. 0) call uttrii(zi(ltcb), nbcb)
    if (nbmn .ne. 0) call uttrii(zi(ltmn), nbmn)
    if (nbha .ne. 0) call uttrii(zi(ltha), nbha)
    if (nbau .ne. 0) call uttrii(zi(ltau), nbau)
!
!-----------COMPARAISON LISTES MN ET CB POUR DETECTION INTERSECTION-----
!
    call nodoub(nbmn, nbcb, zi(ltmn), zi(ltcb), mn, &
                cb, mailla, doubok)
    call nodoub(nbmn, nbha, zi(ltmn), zi(ltha), mn, &
                ha, mailla, doubok)
    call nodoub(nbha, nbcb, zi(ltha), zi(ltcb), ha, &
                cb, mailla, doubok)
    if (doubok) then
        call utmess('F', 'ALGORITH12_67')
    end if
!
!-----------ALLOCATION TABLEAU DESCRIPTION DEFORMEES CALCULEES----------
!
    nbto = nbcb+nbmn+nbha+nbau
!
    nbtem = (2+nbec)*nbto
!
    desdef = nomres//'.IDC_DEFO'
    call wkvect(desdef, 'G V I', nbtem, lldes)
!
!-----REMPLISSAGE DU TABLEAU PAR ORDRE DES NOEUDS (MN CB AU)------------
!
    nbtem = 0
!
!    NOEUD MAC NEAL
!
    call copvis(nbmn, zi(ltmn), zi(lldes+nbtem))
    nbtem = nbtem+nbmn
!
!
!    NOEUD CRAIG BAMPTON
!
    call copvis(nbcb, zi(ltcb), zi(lldes+nbtem))
    nbtem = nbtem+nbcb
!
!    NOEUD CRAIG-BAMPTON-HARMONIQUE
!
    call copvis(nbha, zi(ltha), zi(lldes+nbtem))
    nbtem = nbtem+nbha
!
!    NOEUD AUCUN
!
    call copvis(nbau, zi(ltau), zi(lldes+nbtem))
    nbtem = nbtem+nbau
!
!--RECUPERATION DU NOMBRE DE COMPOSANTES DE LA GRANDEUR SOUS-JACENTE----
!
    call dismoi('NB_CMP_MAX', nomres, 'INTERF_DYNA', repi=nbcmp)
!
!-------------DEFINITION DU NOM  OBJET MASQUE AUX INTERFACES------------
!
    temmas = '&&DEFINT'//'.MASQUE'
!
!---------MODIFICATION NUMEROTATION DANS DEFINITION INTERFACES----------
! LE NUMERO DU NOEUD DANS LE MAILLAGE DEVIENT LE NUMERO DANS LA
!        LISTE DES NOEUDS D'INTERFACE
!
! REMPLISSAGE COMME INDIQUE DANS L'ENTETE DES COLONNES 2 ET 3 DU
! DESCRIPTEUR DES DEFORMEES
!
    do i = 1, nbint
        nomtyp = zk8(lltyp+i-1)
        call jelira(jexnum(noeint, i), 'LONMAX', nbno)
        call jeveuo(jexnum(noeint, i), 'E', llnin)
        call jeveuo(jexnum(temmas, i), 'L', ltmas)
!
!    NOEUD DE MAC NEAL
!
        if (nomtyp .eq. mn) then
!
            nbtem = 0
!
            do j = 1, nbno
                ino = zi(llnin+j-1)
                call cheris(nbmn, zi(lldes+nbtem), ino, iran)
                zi(llnin+j-1) = iran+nbtem
!
                iad = lldes+nbto+iran+nbtem-1
                zi(iad) = -1
!
                iad = lldes+nbto*2+(iran-1)*nbec+nbtem-1+1
                call isgeco(zi(ltmas+(j-1)*nbec), zi(iad), nbcmp, 1, zi(iad))
            end do
        end if
!
!    NOEUDS DE CRAIG BAMPTON
!
        if (nomtyp .eq. cb) then
!
            nbtem = nbmn
!
            do j = 1, nbno
                ino = zi(llnin+j-1)
                call cheris(nbcb, zi(lldes+nbtem), ino, iran)
                zi(llnin+j-1) = iran+nbtem
!
                iad = lldes+nbto+iran+nbtem-1
                zi(iad) = -2
!
                iad = lldes+nbto*2+(iran-1)*nbec+nbtem-1+1
                call isgeco(zi(ltmas+(j-1)*nbec), zi(iad), nbcmp, 1, zi(iad))
!
            end do
        end if
!
!    NOEUDS DE CRAIG-BAMPTON-HARMONIQUE
!
        if (nomtyp .eq. ha) then
!
            nbtem = nbmn+nbcb
!
            do j = 1, nbno
                ino = zi(llnin+j-1)
                call cheris(nbha, zi(lldes+nbtem), ino, iran)
                zi(llnin+j-1) = iran+nbtem
!
                iad = lldes+nbto+iran+nbtem-1
                zi(iad) = -3
!
                iad = lldes+nbto*2+(iran-1)*nbec+nbtem-1+1
                call isgeco(zi(ltmas+(j-1)*nbec), zi(iad), nbcmp, 1, zi(iad))
!
            end do
        end if
!
!    NOEUD DE AUCUN
!
        if (nomtyp .eq. au) then
!
            nbtem = nbmn+nbcb+nbha
!
            do j = 1, nbno
                ino = zi(llnin+j-1)
                call cheris(nbau, zi(lldes+nbtem), ino, iran)
                zi(llnin+j-1) = iran+nbtem
!
                iad = lldes+nbto+iran+nbtem-1
                zi(iad) = -4
!
                iad = lldes+nbto*2+(iran-1)*nbec+nbtem-1+1
                call isgeco(zi(ltmas+(j-1)*nbec), zi(iad), nbcmp, 1, zi(iad))
!
            end do
        end if
!
    end do
!
!---------SAUVEGARDE DE LA DEFINITION DES INTERFACES--------------------
!
    call jedetr('&&DEFINT.MASQUE')
    call jedetr('&&DEFINT.NOM.INTF')
    if (nbmn .gt. 0) call jedetr(temmn)
    if (nbcb .gt. 0) call jedetr(temcb)
    if (nbau .gt. 0) call jedetr(temau)
    if (nbha .gt. 0) call jedetr(temha)
!
    call jedema()
end subroutine
