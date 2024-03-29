! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine ajellt(ligrez, nomaz, nbma, limaz, typelz, &
                  phenoz, modelz, nbno, linoz)
    implicit none
!     BUT : AFFECTATION DE L'OBJET DE TYPE LIGRET
!           ET DE NOM LIGREZ
!           ON STOCKE LA LISTE DE MAILLES LIMA DANS LE VECTEUR
!           LIGRET//'.LIMA'
!           QUAND IL S'AGIT DE MAILLES TARDIVES, ON A LIMA(1) = 0,
!           ALORS ON STOCKE LES NOEUDS DE LA MAILLE TARDIVE DANS
!           LE VECTEUR LIGRET//'.LINO'
!           QUAND TYPEL N'EST PAS DEFINI, ON RECUPERE LE TYPE DES
!           MAILLES VIA LA MODELISATION ET LE PHENOMENE .
!
!
!  ARGUMENT       E/S    TYPE          ROLE
!
!  LIGREZ         IN      K19     NOM DU LIGRET A AFFECTER
!                 JXVAR
!  NOMAZ          IN      K8      NOM DU MAILLAGE SUR LEQUEL S'APPUIE
!                                 LE LIGRET
!  NBMA           IN      I       NOMBRE DE MAILLES A AFFECTER
!  LIMAZ          IN      K24     NOM DU VECTEUR JEVEUX CONTENANT
!                                 LA LISTE DES NUMEROS DE MAILLES
!  TYPELZ         IN      K16     TYPE DES MAILLES A AFFECTER
!  PHENOZ         IN      K16     PHENOMENE ASSOCIE AU MODELE
!  MODELZ         IN      K16     MODELISATION ASSOCIEE AU MODELE
!  NBNO           IN      I       NOMBRE DE NOEUDS DE LA MAILLE
!                                 TARDIVE A AFFECTER
!  LINOZ          IN      K24     NOM DU VECTEUR JEVEUX CONTENANT
!                                 LA LISTE DES NUMEROS DE NOEUDS
!-------------------------------------------------------------
!
! ====================== DEBUT DES DECLARATIONS ========================
#include "jeveux.h"
#include "asterfort/crelgt.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/wkvect.h"
!
!
! ----- ARGUMENTS
    character(len=*) :: ligrez, nomaz, typelz, phenoz, modelz, limaz, linoz
! ----- VARIABLES LOCALES -------------------------------
!-----------------------------------------------------------------------
    integer :: i, idlima, idlino, idlity
    integer :: idpoma, idpono, imodl
    integer :: iret, iret1, iret2, ityp, jdlima, jdlino, jdpm
    integer :: jdtm, lolima, lolimx, lolino, lolinx, lopomx, loponx
    integer :: matard, nbma, nbmadi, nbmail, nbmax, nbno, nbnodi
    integer :: nbnox, nlolim, nlolin, numail, nutypm
!-----------------------------------------------------------------------
    parameter(nbmail=10000)
!
    character(len=8) :: noma
    character(len=16) :: pheno, modeli, typel
    character(len=19) :: ligret
    character(len=24) :: lima, lino, typmai
    integer, pointer :: apma(:) => null()
    character(len=16), pointer :: phen(:) => null()
    integer, pointer :: apno(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer, pointer :: vnbma(:) => null()
    integer, pointer :: mata(:) => null()
    character(len=16), pointer :: mode(:) => null()
! ====================== DEBUT DU CODE EXECUTABLE ======================
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    noma = nomaz
    pheno = phenoz
    modeli = modelz
    ligret = ligrez
    typel = typelz
    lima = limaz
    lino = linoz
    typmai = noma//'.TYPMAIL'
!
    matard = 0
!
! --- ON VERIFIE SI LE LIGREL EXISTE ,S'IL N'EXISTE PAS, ON LE CREE :
!     -------------------------------------------------------------
    call jeexin(ligret//'.LGRF', iret)
!
    if (iret .eq. 0) then
        call crelgt('V', ligret)
    end if
!
! --- VECTEUR DE LA LISTE DES MAILLES CUMULEES DU LIGRET :
!     --------------------------------------------------
    call jeveuo(ligret//'.LIMA', 'E', idlima)
!
! --- VECTEUR DES TYPES DES MAILLES CUMULEES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LITY', 'E', idlity)
!
! --- NOM DE LA MODELISATION :
!     ----------------------
    call jeveuo(ligret//'.MODE', 'E', vk16=mode)
!
! --- NOM DU PHENOMENE :
!     ----------------
    call jeveuo(ligret//'.PHEN', 'E', vk16=phen)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES MAILLES :
!     ----------------------------------------------
    call jeveuo(ligret//'.POMA', 'E', idpoma)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES NOEUDS :
!     ---------------------------------------------
    call jeveuo(ligret//'.PONO', 'E', idpono)
!
! --- NOM DU MAILLAGE :
!     ---------------
    call jeveuo(ligret//'.LGRF', 'E', vk8=lgrf)
!
! --- NOMBRE DE MAILLES TARDIVES :
!     --------------------------
    call jeveuo(ligret//'.MATA', 'E', vi=mata)
!
! --- VECTEUR DE LA LISTE DES NOEUDS CUMULES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LINO', 'E', idlino)
!
! --- NOMBRE D'AFFECTATIONS DE MAILLES :
!     --------------------------------
    call jeveuo(ligret//'.APMA', 'E', vi=apma)
!
! --- NOMBRE D'AFFECTATIONS DE NOEUDS :
!     -------------------------------
    call jeveuo(ligret//'.APNO', 'E', vi=apno)
!
! --- NOMBRE D'AFFECTATIONS DE MAILLES :
!     --------------------------------
    call jeveuo(ligret//'.NBMA', 'E', vi=vnbma)
!
    vnbma(1) = vnbma(1)+nbma
!
! --- ON AFFECTE UNE FOIS POUR TOUTES LE NOM DU MAILLAGE :
!     --------------------------------------------------
    if (iret .eq. 0) then
        lgrf(1) = noma
    end if
!
! --- RECUPERATION DE LA LISTE DES MAILLES A AFFECTER :
!     ===============================================
    if (lima .ne. ' ') then
        call jeexin(lima, iret1)
        if (iret1 .eq. 0) then
            call wkvect(lima, 'V V I', 1, jdlima)
        else
            call jeveuo(lima, 'L', jdlima)
        end if
    end if
!
! --- RECUPERATION DE LA LISTE DES NOEUDS A AFFECTER :
!     ===============================================
    if (lino .ne. ' ') then
        call jeexin(lino, iret2)
        if (iret2 .eq. 0) then
            call wkvect(lino, 'V V I', 1, jdlino)
        else
            call jeveuo(lino, 'L', jdlino)
        end if
    end if
!
! --- VERIFICATION DE L'ADEQUATION DE L'AFFECTATION DES MAILLES
! --- A LA LISTE DES MAILLES CUMULEES :
!     ===============================
    if (zi(jdlima) .gt. 0 .and. nbma .ge. 1) then
!
! ---   NOMBRE DE MAILLES DEJA AFFECTEES :
!       --------------------------------
        call jelira(ligret//'.LIMA', 'LONUTI', lolima)
!
! ---   LONGUEUR DU VECTEUR LIGRET.LIMA :
!       -------------------------------
        call jelira(ligret//'.LIMA', 'LONMAX', lolimx)
!
! ---   NOMBRE DE MAILLES DISPONIBLES :
!       -----------------------------
        nbmadi = lolimx-lolima
!
! ---   REAJUSTEMENT EVENTUEL DES VECTEURS LIMA ET LITY :
!       -----------------------------------------------
        if (nbma .gt. nbmadi) then
            nlolim = nbma-nbmadi
            nbmax = lolimx+max(nlolim, nbmail)
            call juveca(ligret//'.LIMA', nbmax)
            call jeveuo(ligret//'.LIMA', 'E', idlima)
            call juveca(ligret//'.LITY', nbmax)
            call jeveuo(ligret//'.LITY', 'E', idlity)
        end if
!
! ---   VERIFICATION DE L'ADEQUATION DE LA TAILLE DU VECTEUR
! ---   DES POINTEURS DANS LA LISTE DE MAILLES :
!       --------------------------------------
!
! ---   NOMBRE D'AFFECTATIONS DE MAILLES :
!       --------------------------------
        apma(1) = apma(1)+1
!
! ---   LONGUEUR DU VECTEUR LIGRET.POMA :
!       -------------------------------
        call jelira(ligret//'.POMA', 'LONMAX', lopomx)
!
! ---   REAJUSTEMENT EVENTUEL DU VECTEUR POMA :
!       -------------------------------------
        if (apma(1) .ge. lopomx) then
            call juveca(ligret//'.POMA', 2*lopomx)
            call jeveuo(ligret//'.POMA', 'E', idpoma)
        end if
!
    end if
!
! --- VERIFICATION DE L'ADEQUATION DE L'AFFECTATION DES NOEUDS
! --- A LA LISTE DES NOEUDS CUMULES :
!     =============================
    if (zi(jdlima) .eq. 0 .and. nbma .eq. 1) then
!
! ---   NOMBRE DE NOEUDS DEJA AFFECTES :
!       ------------------------------
        call jelira(ligret//'.LINO', 'LONUTI', lolino)
!
! ---   LONGUEUR DU VECTEUR LIGRET.LINO :
!       -------------------------------
        call jelira(ligret//'.LINO', 'LONMAX', lolinx)
!
! ---   NOMBRE DE NOEUDS DISPONIBLES :
!       ----------------------------
        nbnodi = lolinx-lolino
!
! ---   REAJUSTEMENT EVENTUEL DU VECTEUR LINO :
!       -------------------------------------
        if (nbno .gt. nbnodi) then
            nlolin = nbno-nbnodi
            nbnox = lolinx+max(nlolin, nbmail)
            call juveca(ligret//'.LINO', nbnox)
            call jeveuo(ligret//'.LINO', 'E', idlino)
        end if
!
! ---   VERIFICATION DE L'ADEQUATION DE LA TAILLE DU VECTEUR
! ---   DES POINTEURS DANS LA LISTE DE NOEUDS :
!       -------------------------------------
!
! ---   NOMBRE D'AFFECTATIONS DE NOEUDS :
!       -------------------------------
        apno(1) = apno(1)+1
!
! ---   LONGUEUR DU VECTEUR LIGRET.PONO :
!       -------------------------------
        call jelira(ligret//'.PONO', 'LONMAX', loponx)
!
! ---   REAJUSTEMENT EVENTUEL DU VECTEUR PONO :
!       -------------------------------------
        if (apno(1) .gt. loponx) then
            call juveca(ligret//'.PONO', 2*loponx)
            call jeveuo(ligret//'.PONO', 'E', idpono)
        end if
!
    end if
!
! --- AFFECTATION DES MAILLES TARDIVES :
!     ================================
    if (zi(jdlima) .eq. 0 .and. nbma .eq. 1) then
!
! ---   ON INCREMENTE LE NOMBRE DE MAILLES TARDIVES :
!       -------------------------------------------
        mata(1) = mata(1)+1
        matard = matard+1
!
! ---   AFFECTATION DU VECTEUR DES NOEUDS DU LIGRET :
!       -------------------------------------------
        do i = 1, nbno
            zi(idlino+zi(idpono+matard-1)+i-1) = zi(jdlino+i-1)
        end do
!
        zi(idpono+matard) = zi(idpono+matard-1)+nbno
!
        call jeecra(ligret//'.LINO', 'LONUTI', zi(idpono+matard))
!
! --- AFFECTATION DES MAILLES PHYSIQUES :
!     =================================
    else
!
! ---   AFFECTATION DU TYPE DES MAILLES :
!       -------------------------------
        if (typel .eq. ' ') then
            call jenonu(jexnom('&CATA.'//pheno(1:13)//'.MODL', modeli), imodl)
            call jeveuo(jexnum('&CATA.'//pheno, imodl), 'L', jdpm)
            call jeveuo(typmai, 'L', jdtm)
        else
            call jenonu(jexnom('&CATA.TE.NOMTE', typel), ityp)
        end if
!
!
! ---   AFFECTATION DE LA LISTE DES MAILLES CUMULEES :
!       --------------------------------------------
        do i = 1, nbma
            zi(idlima+zi(idpoma+apma(1)-1)+i-1) = zi(jdlima+i-1)
            if (typel .eq. ' ') then
                numail = zi(jdlima+i-1)
                nutypm = zi(jdtm+numail-1)
                ityp = zi(jdpm+nutypm-1)
            else
                call jenonu(jexnom('&CATA.TE.NOMTE', typel), ityp)
            end if
!
            zi(idlity+zi(idpoma+apma(1)-1)+i-1) = ityp
        end do
!
        vnbma(1) = vnbma(1)+nbma
!
! ---   VECTEUR DE POINTEURS DANS LE VECTEUR DES MAILLES :
!       ------------------------------------------------
        zi(idpoma+apma(1)) = zi(idpoma+apma(1)-1)+nbma
!
        call jeecra(ligret//'.LIMA', 'LONUTI', zi(idpoma+apma(1)))
!
    end if
!
! --- AFFECTATION DE LA MODELISATION AU LIGRET :
!     ----------------------------------------
    mode(1) = modeli
!
! --- AFFECTATION DU PHENOMENE AU LIGRET :
!     ----------------------------------
    phen(1) = pheno
!
    call jedema()
!
end subroutine
