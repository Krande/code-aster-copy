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

subroutine crelil(kstop, nbmat, vlimat, lili, base, &
                  nomma, pref, gd, mailla, nec, &
                  ncmp, ilimo, nlili, nbelm, nume_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbcmp.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmat, gd, nec, ilimo, nlili, iconx1, iconx2, nbelm
    integer(kind=8) :: iadnem, iadlie
    character(len=*) :: lili, nomma, pref, mailla
    character(len=1) :: base, kstop
    character(len=*) :: vlimat(nbmat)
    character(len=14), optional :: nume_
!
!
!
! ---- OBJET : CREATION DU CHAMP .LILI A PARTIR D'UNE LISTE DE MATR_ELEM
!
! ---- DESCRIPTION DES PARAMETRES
! IN  K1   KSTOP  : (VOIR L'ARGUMENT "OUT" NLILI)
! IN  I    NBMAT  : NBRE DE MATR_ELEM
! IN  K* VLIMAT : LISTE DES MAT_ELE
! IN  K*24 LILI   : NOM DE L OBJET LILI QUI SERA CREE
! IN  K1   BASE   : ' G ' POUR CREER LILI SUR BASE GLOBALE
!                   ' V ' POUR CREER LILI SUR BASE VOLATILE
! IN  K*24 NOMMA  : NOM FORFAITAIRE DU LIGREL &MAILLA
! IN  K*19 PREF   : PREFIXE DES OBJETS TEMPORAIRES CREES
!
! OUT I    GD     : GRANDEUR SIMPLE ASSOCIEE AU NUME_DDL
! OUT K*8  MAILLA : NOM DU MAILLAGE
! OUT I    NEC    : NBRE D ENTIERS CODES POUR GD
! OUT I    NCMP   : NBRE DE CMP POUR GD
! OUT I    ILIMO  : NUMERO DANS LE REPERTOIRE LILI DU NOM DU 1ER LIGREL
!                   APPARTENANT A UNE S.D. DE TYPE MODELE , =0 SINON
! OUT I    NLILI  : DIMENSION DE L OBJET CREE K24LIL
!                   SI NLILI=1, C'EST QU'ON N'A TROUVE AUCUN RESUELEM
!                   => ERREUR FATALE SI KSTOP='F'
!                   => ON CONTINUE   SI KSTOP='C'
! OUT I    NBELM  : NBRE D ELEMENTS DU MAILLAGE
!    --- DESCRIPTION DES OBJETS ADNE ET ADLI ---
!     ADNE (1          ) = NBRE DE MAILLES DU MAILLAGE
!     ADNE (2          ) = 0
!     ADNE (3          ) = 0
!     ADLI (1          ) = NBRE DE MAILLES DU MAILLAGE
!     ADLI (2          ) = 0
!     ADLI (3          ) = 0
!     POUR 2<=ILI<=NLILI
!     ADNE (3*(ILI-1)+1) = NBRE MAX D'OBJETS DE LA COLLECTION
!                            LILI(ILI).NEMA
!     ADNE (3*(ILI-1)+2) = ADRESSE DE L'OBJET LILI(ILI).NEMA
!     ADNE (3*(ILI-1)+3) = ADRESSE DU VECTEUR DES LONG. CUMULEES DE
!                            LILI(ILI).NEMA
!     ADLI (3*(ILI-1)+1) = NBRE MAX D'OBJETS DE LA COLLECTION
!                            LILI(ILI).LIEL
!     ADLI (3*(ILI-1)+2) = ADRESSE DE L'OBJET LILI(ILI).LIEL
!     ADLI (3*(ILI-1)+3) = ADRESSE DU VECTEUR DES LONG. CUMULEES DE
!                            LILI(ILI).LIEL
!-----------------------------------------------------------------------
!     FONCTIONS JEVEUX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    character(len=8) :: modele, models, exiss1, exiss2
    character(len=16) :: pheno
!----------------------------------------------------------------------
!     VARIABLES LOCALES
!----------------------------------------------------------------------
    character(len=19) :: prefix, matel
    character(len=24) :: resu, nomli, k24lil, kmaill
!-----------------------------------------------------------------------
!     FONCTIONS LOCALES D'ACCES AUX DIFFERENTS CHAMPS DES
!     S.D. MANIPULEES DANS LE SOUS PROGRAMME
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!                DEBUT DES INSTRUCTIONS
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iad, icomp, idimli, idlres
    integer(kind=8) :: ili, imat, iresu, iret, iret1, n1
    integer(kind=8) :: nbgr, nbmo, nbresu, nbsup, ncmp
    character(len=24), pointer :: rerr(:) => null()
    aster_logical :: l_parallel_mesh
!
!-----------------------------------------------------------------------
    k24lil = lili
    kmaill = nomma
    prefix = pref
!
!---- CALCUL DU NBRE DE LIGRELS REFERENCES
!
    idimli = 1
    models = ' '
    exiss2 = 'NON'
!
!     -- VERIFICATION DES MATR_ELEM :
!     -------------------------------
    l_parallel_mesh = isParallelMesh(mailla)
    do imat = 1, nbmat
        matel = vlimat(imat)
        call jeexin(matel//'.RERR', iret1)
        ASSERT(iret1 .gt. 0)
        call jeveuo(matel//'.RERR', 'L', vk24=rerr)
        modele = rerr(1) (1:8)
        if ((modele .ne. models) .and. (models .ne. ' ')) then
            call utmess('F', 'ASSEMBLA_18')
        end if
        models = modele
!
        call dismoi('NB_SS_ACTI', matel, 'MATR_ELEM', repi=n1)
        if (n1 .gt. 0) then
            exiss1 = 'OUI'
            exiss2 = 'OUI'
        else
            exiss1 = 'NON'
        end if
!
!       testing the presence of .RELR is not mandatory in parallel
        call jeexin(matel//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(matel//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(matel//'.RELR', 'L', idlres)
            idimli = idimli+nbresu
        else
            if (.not. l_parallel_mesh) then
                if (exiss1(1:3) .eq. 'NON') then
                    call utmess('F', 'ASSEMBLA_19')
                end if
            end if
        end if
    end do
!
!     --SI IL EXISTE DES SOUS-STRUCTURES, ON COMPTE 1 LIGREL DE PLUS:
!
    if (exiss2(1:3) .eq. 'OUI') idimli = idimli+1
!
    if (present(nume_)) then
        call dismoi('NOM_MAILLA', nume_, 'NUME_DDL', repk=mailla(1:8))
        call dismoi('NUM_GD_SI', nume_, 'NUME_DDL', repi=gd)
    else
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mailla(1:8))
        call dismoi('PHENOMENE', modele, 'MODELE', repk=pheno)
        call dismoi('NUM_GD', pheno, 'PHENOMENE', repi=gd)
    end if
!
!---- CALCUL DE NEC ET NCMP
!
    nec = nbec(gd)
    ncmp = nbcmp(gd)
!
!
!
!---- CREATION DU REPERTOIRE .LILI DE TOUS LES NOMS DE LIGRELS /=
!     TROUVES DANS LES RESUELEMS ATTEINT A PARTIR DE LA LISTE DES
!     MATR_ELEM + EN 1ER LE MOT '&MAILLA'
!
    call jeexin(k24lil, iret)
    if (iret .gt. 0) call jedetr(k24lil)
    call jecreo(k24lil, base//' N  K24 ')
    call jeecra(k24lil, 'NOMMAX', idimli)
!---- LILI(1)= '&MAILLA'
    call jecroc(jexnom(k24lil, kmaill))
!---- SI LES RESUELEM SONT VIDES .LILI NE CONTIENT QUE &MAILLA
    nlili = 1
    if (idimli .eq. 1) goto 101
!
!---- CALCUL DE LILI
!
    do imat = 1, nbmat
        matel = vlimat(imat)
        call jeexin(matel//'.RELR', iret)
        if (iret .eq. 0) goto 110
        call jelira(matel//'.RELR', 'LONUTI', nbresu)
        if (nbresu .gt. 0) call jeveuo(matel//'.RELR', 'L', idlres)
        do iresu = 1, nbresu
            resu = zk24(idlres+iresu-1)
!
            call jeexin(resu(1:19)//'.NOLI', iret)
            if (iret .eq. 0) goto 120
            call jeveuo(resu(1:19)//'.NOLI', 'L', iad)
            nomli = zk24(iad)
            call jenonu(jexnom(k24lil, nomli), ili)
            if (ili .eq. 0) then
!
!           ---- SI CE LIGREL N EST PAS DANS LILI ON LE MET
!
                call jeexin(nomli(1:19)//'.NBNO', iret)
                if (iret .ne. 0) then
                    call jeveuo(nomli(1:19)//'.NBNO', 'L', iad)
                else
                end if
                nlili = nlili+1
                call jecroc(jexnom(k24lil, nomli))
            end if
120         continue
        end do
110     continue
    end do
!
!     -- ON REGARDE SI ON DOIT AJOUTER LE LIGREL DE MODELE POUR LES
!     -- SOUS-STRUCTURES:
    if (exiss2(1:3) .eq. 'OUI') then
        icomp = 0
        do ili = 1, nlili
            call jenuno(jexnum(k24lil, ili), nomli)
            if (nomli(1:8) .eq. modele) icomp = 1
        end do
        if (icomp .eq. 0) then
            nlili = nlili+1
            call jecroc(jexnom(k24lil, modele//'.MODELE'))
        end if
    end if
!
!
!
    if (nlili .eq. 1) then
        if (kstop .eq. 'C') then
            goto 999
        else
            ASSERT(kstop .eq. 'F')
            call utmess('F', 'ASSEMBLA_20')
        end if
    end if
!
101 continue
!
!---- RECUPERATION ADRESSES DE CONNEX: ICONX1, ICONX2  ET NBELM
!
    call dismoi('NB_MA_MAILLA', mailla(1:8), 'MAILLAGE', repi=nbelm)
    if (nbelm .gt. 0) then
        call jeveuo(mailla(1:8)//'.CONNEX', 'L', iconx1)
        call jeveuo(jexatr(mailla(1:8)//'.CONNEX', 'LONCUM'), 'L', iconx2)
    end if
!
!---- CREATION DES OBJETS ADNE ET ADLI SUR 'V'
!
    call jeexin(prefix//'.ADNE', iret)
    if (iret .gt. 0) call jedetr(prefix//'.ADNE')
    if (iret .gt. 0) call jedetr(prefix//'.ADLI')
    call wkvect(prefix//'.ADNE', ' V V I', 3*nlili, iadnem)
    call wkvect(prefix//'.ADLI', ' V V I', 3*nlili, iadlie)
!---- ADNE(1)= NBELM
    zi(iadnem) = nbelm
    nbmo = 0
    ilimo = 0
    do ili = 2, nlili
        call jenuno(jexnum(k24lil, ili), nomli)
!
!---- CALCUL DU NBRE DE LIGRELS DE MODELE : NBMO ET DE ILIMO
!
        if (nomli(9:15) .eq. '.MODELE') then
            nbmo = nbmo+1
            if (nbmo .eq. 1) ilimo = ili
            if (nbmo .gt. 1) then
                call utmess('F', 'ASSEMBLA_21')
            end if
        end if
        call jeexin(nomli(1:19)//'.NEMA', iret)
        if (iret .ne. 0) then
!
!---- ADNE(3*(ILI-1)+1)=NBRE DE MAILLES SUP DU LIGREL NOMLI
!
            call jelira(nomli(1:19)//'.NEMA', 'NUTIOC', nbsup)
            zi(iadnem+3*(ili-1)) = nbsup
            call jeveut(nomli(1:19)//'.NEMA', 'L', iad)
            zi(iadnem+3*(ili-1)+1) = iad
            call jeveut(jexatr(nomli(1:19)//'.NEMA', 'LONCUM'), 'L', iad)
            zi(iadnem+3*(ili-1)+2) = iad
        else
            zi(iadnem+3*(ili-1)) = 0
            zi(iadnem+3*(ili-1)+1) = 2**30
            zi(iadnem+3*(ili-1)+2) = 2**30
        end if
!
!---- ADLI(3*(ILI-1)+1)=NBRE DE MAILLES DU LIGREL NOMLI
!
        call jeexin(nomli(1:19)//'.LIEL', iret)
        if (iret .gt. 0) then
            call jelira(nomli(1:19)//'.LIEL', 'NUTIOC', nbgr)
            zi(iadlie+3*(ili-1)) = nbgr
            call jeveut(nomli(1:19)//'.LIEL', 'L', iad)
            zi(iadlie+3*(ili-1)+1) = iad
            call jeveut(jexatr(nomli(1:19)//'.LIEL', 'LONCUM'), 'L', iad)
            zi(iadlie+3*(ili-1)+2) = iad
        else
            zi(iadlie+3*(ili-1)) = 0
            zi(iadlie+3*(ili-1)+1) = 2**30
            zi(iadlie+3*(ili-1)+2) = 2**30
        end if
    end do
999 continue
end subroutine
