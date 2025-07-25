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

subroutine lgtlgr(basez, ligrey, ligrez)
    implicit none
!     BUT : CREATION ET AFFECTATION DE L'OBJET DE TYPE LIGREL
!           ET DE NOM LIGREZ A PARTIR DE L'OBJET DE TYPE LIGRET
!           ET DE NOM LIGREY SUR LA BASE BASEZ.
!
!
!  ARGUMENT       E/S    TYPE          ROLE
!
!  BASEZ          IN      K1      NOM DE LA BASE
!  LIGREY         IN      K19     NOM DU LIGRET SERVANT A CREER
!                                 LE LIGREL LIGREZ
!  LIGREZ         IN      K19     NOM DU LIGREL A CREER ET AFFECTER
!                 JXVAR
!-------------------------------------------------------------
!
! ====================== DEBUT DES DECLARATIONS ========================
#include "jeveux.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/initel.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
!
! ----- ARGUMENTS
    character(len=*) :: basez, ligrey, ligrez
! ----- VARIABLES LOCALES -------------------------------
    character(len=1) :: base
    character(len=8) :: moloc
    character(len=19) :: ligret, ligrel
! ====================== DEBUT DU CODE EXECUTABLE ======================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, idligi
    integer(kind=8) ::   idnema, idphen
    integer(kind=8) :: ij, imodl, iret, j, jdnbno, jdpm, k
    integer(kind=8) :: k1, lonlie, nbapma, nbapno, nbmato, nbmaty, nbno
    integer(kind=8) :: nbno2, nbnoto, ntypoi, nutyp1, nutype
    integer(kind=8), pointer :: pono(:) => null()
    integer(kind=8), pointer :: lity(:) => null()
    integer(kind=8), pointer :: poma(:) => null()
    integer(kind=8), pointer :: lima(:) => null()
    integer(kind=8), pointer :: apno(:) => null()
    character(len=16), pointer :: mode(:) => null()
    integer(kind=8), pointer :: lino(:) => null()
    integer(kind=8), pointer :: apma(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    base = basez
    ligret = ligrey
    ligrel = ligrez
!
! --- ON VERIFIE SI LE LIGREL EXISTE ,S'IL N'EXISTE PAS, ON
! --- S'ARRETE EN ERREUR FATALE :
!     -------------------------
    call jeexin(ligret//'.LGRF', iret)
    ASSERT(iret .ne. 0)
!
! --- NUMERO DU TYPE ASSOCIE A DES MAILLES TARDIVES :
!     ---------------------------------------------
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ntypoi)
!
! --- RECUPERATION DES ATTRIBUTS DU LIGRET :
!     ====================================
!
! --- VECTEUR DE LA LISTE DES MAILLES CUMULEES DU LIGRET :
!     --------------------------------------------------
    call jeveuo(ligret//'.LIMA', 'L', vi=lima)
!
! --- VECTEUR DES TYPES DES MAILLES CUMULEES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LITY', 'L', vi=lity)
!
! --- NOM DE LA MODELISATION :
!     ----------------------
    call jeveuo(ligret//'.MODE', 'L', vk16=mode)
!
! --- NOM DU PHENOMENE :
!     ----------------
    call jeveuo(ligret//'.PHEN', 'L', idphen)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES MAILLES :
!     ----------------------------------------------
    call jeveuo(ligret//'.POMA', 'L', vi=poma)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES NOEUDS :
!     ---------------------------------------------
    call jeveuo(ligret//'.PONO', 'L', vi=pono)
!
! --- VECTEUR DE LA LISTE DES NOEUDS CUMULES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LINO', 'L', vi=lino)
!
! --- NOMBRE D'AFFECTATIONS DE MAILLES AU LIGRET :
!     ------------------------------------------
    call jeveuo(ligret//'.APMA', 'L', vi=apma)
!
! --- NOMBRE D'AFFECTATIONS DE NOEUDS AU LIGRET :
!     -----------------------------------------
    call jeveuo(ligret//'.APNO', 'L', vi=apno)
!
    nbapma = apma(1)
    nbapno = apno(1)
!
! --- NOMBRE TOTAL DE MAILLES :
!     -----------------------
    nbmato = poma(nbapma+1)
!
! --- NOMBRE TOTAL DE NOEUDS :
!     ----------------------
    nbnoto = pono(nbapno+1)
!
    k1 = 1
    nutyp1 = lity(1)
    do i = 1, nbmato
        nutype = lity(i)
        if (nutype .ne. nutyp1) then
            nutyp1 = nutype
            k1 = k1+1
        end if
    end do
!
! --- ON CREE LE .LIEL SI LE NOMBRE DE MAILLES EST NON NUL :
!     ----------------------------------------------------
    if (nbmato+nbnoto .gt. 0) then
!
        call jecrec(ligrel//'.LIEL', base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                    k1+nbnoto)
!
! ---   LONGUEUR DU LIGREL.LIEL :
!       -----------------------
        lonlie = nbmato+k1+nbnoto+nbapno
!
! ---   AFFECTATION DE LA LONGUEUR :
!       --------------------------
        call jeecra(ligrel//'.LIEL', 'LONT', lonlie)
!
! ---   CREATION DE L'OBJET LIGREL.NBNO :
!       -------------------------------
        call jeexin(ligrel//'.NBNO', iret)
        if (iret .eq. 0) then
            call wkvect(ligrel//'.NBNO', base//' V I', 1, jdnbno)
            zi(jdnbno) = nbnoto
        end if
!
! ---   CREATION DE L'OBJET LIGREL.LGRF :
!       -------------------------------
        call jedupo(ligret//'.LGRF', 'V', ligrel//'.LGRF', .false._1)
        call jeecra(ligrel//'.LGRF', 'DOCU', ibid, 'MECA')
!
        nutyp1 = lity(1)
        k = 0
        nbmaty = 0
        ij = 0
!
! ---   BOUCLE SUR LE NOMBRE D'AFFECTATIONS DU LIGRET PAR DES MAILLES :
!       -------------------------------------------------------------
        do i = 1, nbmato
!
            nutype = lity(i)
!
            if (nutype .ne. nutyp1) then
!
                k = k+1
!
                if (k .eq. k1) goto 32
!
! ---       CREATION DU IEME OBJET DE COLLECTION :
!           ------------------------------------
                call jecroc(jexnum(ligrel//'.LIEL', k))
!
! ---       LONGUEUR DU IEME OBJET DE COLLECTION :
!           ------------------------------------
                call jeecra(jexnum(ligrel//'.LIEL', k), 'LONMAX', nbmaty+1)
!
! ---       AFFECTATION DU IEME OBJET DE COLLECTION :
!           ---------------------------------------
                call jeveuo(jexnum(ligrel//'.LIEL', k), 'E', idligi)
!
                do j = 1, nbmaty
                    ij = ij+1
                    zi(idligi+j-1) = lima(ij)
                end do
!
                zi(idligi+nbmaty) = nutyp1
                nutyp1 = nutype
                nbmaty = 1
!
            else
!
                nbmaty = nbmaty+1
!
            end if
!
!
32          continue
!
        end do
!
        if (nutyp1 .ne. 0) then
            k = k+1
!
! ---       CREATION DU IEME OBJET DE COLLECTION :
!           ------------------------------------
            call jecroc(jexnum(ligrel//'.LIEL', k))
!
! ---       LONGUEUR DU IEME OBJET DE COLLECTION :
!           ------------------------------------
            call jeecra(jexnum(ligrel//'.LIEL', k), 'LONMAX', nbmaty+1)
!
! ---       AFFECTATION DU IEME OBJET DE COLLECTION :
!           ---------------------------------------
            call jeveuo(jexnum(ligrel//'.LIEL', k), 'E', idligi)
!
            do j = 1, nbmaty
                ij = ij+1
                zi(idligi+j-1) = lima(ij)
            end do
!
            zi(idligi+nbmaty) = nutyp1
        end if
!
! --- RECHERCHE DU TYPE DES POI1 :
!     --------------------------
        call jenonu(jexnom('&CATA.'//zk16(idphen) (1:13)//'.MODL', mode(1)), imodl)
        call jeveuo(jexnum('&CATA.'//zk16(idphen), imodl), 'L', jdpm)
!
        nutype = zi(jdpm+ntypoi-1)
!
! ---   BOUCLE SUR LE NOMBRE D'AFFECTATIONS DU LIGRET PAR DES NOEUDS :
!       ------------------------------------------------------------
        do i = 1, nbapno
!
! ---     NOMBRE DE NOEUDS POUR LA IEME OCCURENCE :
!         ---------------------------------------
            nbno = pono(1+i)-pono(i)
!
! ---     CREATION DU IEME OBJET DE COLLECTION :
!         ------------------------------------
            call jecroc(jexnum(ligrel//'.LIEL', i+k))
!
! ---     LONGUEUR DU IEME OBJET DE COLLECTION :
!         ------------------------------------
            call jeecra(jexnum(ligrel//'.LIEL', i+k), 'LONMAX', nbno+1)
!
! ---     AFFECTATION DU IEME OBJET DE COLLECTION :
!         ---------------------------------------
            call jeveuo(jexnum(ligrel//'.LIEL', i+k), 'E', idligi)
!
            do j = 1, nbno
                zi(idligi+j-1) = -j
            end do
!
            zi(idligi+nbno) = nutype
!
        end do
!
    end if
!
! --- ON CREE LE .NEMA SI LE NOMBRE DE NOEUDS EST NON NUL :
!     ---------------------------------------------------
    if (nbnoto .gt. 0) then
!
        call jecrec(ligrel//'.NEMA', base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                    nbnoto)
!
! ---   AFFECTATION DE LA LONGUEUR DU LIGREL.NEMA :
!       -----------------------------------------
        call jeecra(ligrel//'.NEMA', 'LONT', 2*nbnoto)
!
! ---   BOUCLE SUR LE NOMBRE D'AFFECTATIONS DU LIGRET PAR DES NOEUDS :
!       ------------------------------------------------------------
        do i = 1, nbnoto
!
! ---     NOMBRE DE NOEUDS POUR LA IEME OCCURENCE :
!         ---------------------------------------
            nbno = pono(1+i)-pono(i)
            if (nbno .gt. 0) then
                nbno2 = nbno
            else
                nbno2 = -nbno
            end if
!
! ---     CREATION DU IEME OBJET DE COLLECTION :
!         ------------------------------------
            call jecroc(jexnum(ligrel//'.NEMA', i))
!
! ---     LONGUEUR DU IEME OBJET DE COLLECTION :
!         ------------------------------------
            call jeecra(jexnum(ligrel//'.NEMA', i), 'LONMAX', 2*nbno2)
!
! ---     AFFECTATION DU IEME OBJET DE COLLECTION :
!         ---------------------------------------
            call jeveuo(jexnum(ligrel//'.NEMA', i), 'E', idnema)
!
            zi(idnema+1-1) = lino(1+pono(i))
            zi(idnema+2-1) = 1
!
        end do
!
    end if
!
! --- RECUPERATION DU MODE LOCAL ASSOCIE AU PHENOMENE :
!     -----------------------------------------------
    call dismoi('NOM_MOLOC', zk16(idphen), 'PHENOMENE', repk=moloc)
!
! --- ADAPTATION DE LA TAILLE DES GRELS :
!     ---------------------------------
    call adalig(ligrel)
    call cormgi(base, ligrel)
    call initel(ligrel)
!
    call jedema()
!
end subroutine
