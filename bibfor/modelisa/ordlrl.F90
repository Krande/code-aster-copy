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

subroutine ordlrl(charge, lisrel, nomgd)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/ordrel.h"
#include "asterfort/utmess.h"
!
! person_in_charge: jacques.pellet at edf.fr
!
    character(len=19), intent(in) :: lisrel
    character(len=8), intent(in) :: charge
    character(len=8), intent(in) :: nomgd
!
! -----------------------------------------------------------------
!     MISE A JOUR DE L'OBJET DE TYPE  LISTE_RELA ET DE NOM
!     LISREL  :
!               LES RELATIONS SONT REORDONNEES  PAR ORDRE DE NOEUD
!               CROISSANT ET POUR UN NOEUD DONNE PAR DDL CROISSANT
!
!               LES RELATIONS STRICTEMENT EGALES SONT ELIMINEES
!               (I.E. ON NE GARDE QUE LA RELATION DE PLUS GRAND
!                     INDICE DANS LISREL)
! -----------------------------------------------------------------
!     L'OBJET LISREL DOIT OBLIGATOIREMENT EXISTER
! -----------------------------------------------------------------
!  CHARGE        - IN/JXIN    - K8  - : NOM DE LA SD_CHARGE
!  LISREL        - IN/JXVAR   - K19  - : NOM DE LA LISTE_RELA
! -------------------------------------------------------
!
!
! --------- VARIABLES LOCALES ---------------------------
    character(len=24) :: valk(2)
    integer(kind=8) :: nmocl
    parameter(nmocl=320)
    complex(kind=8) :: coproc, rapcoc, dcmplx
    character(len=4) :: typcoe
    character(len=8) :: nomnoe
    character(len=8) :: noma, mod, cmp, nomcmp(nmocl)
    character(len=16) :: kidrel
    character(len=19) :: ligrmo
! --------- FIN  DECLARATIONS  VARIABLES LOCALES --------
    real(kind=8) :: copror, difrel, eps1, eps2, epsrel, rapcoe, coemax
    integer(kind=8) :: i, icmp, icomp, iddl, iddl1, iddl2, ideca1, ideca2
    integer(kind=8) :: idecal, in, indmax, ino, ier
    integer(kind=8) ::  inom, ipntr1, ipntr2, ipntrl, irela, irela1
    integer(kind=8) :: irela2
    integer(kind=8) ::  jprnm, jrlco, jrlco1, jrlco2, jrlcof
    integer(kind=8) :: jrldd
    integer(kind=8) :: jrlno, idnoe1, idnoe2, idnoeu
    integer(kind=8) :: nbcmp, nbec, nbrela, nbtema, nbter1, nbter2, nbterm
    integer(kind=8) :: nddla, nidrel
    integer(kind=8), pointer :: rlnt(:) => null()
    character(len=8), pointer :: rltc(:) => null()
    integer(kind=8), pointer :: rlsu(:) => null()
    integer(kind=8), pointer :: rlpo(:) => null()
    integer(kind=8), pointer :: rlnr(:) => null()
    complex(kind=8), pointer :: coef_c(:) => null()
    integer(kind=8), pointer :: coefmax(:) => null()
    real(kind=8), pointer :: coef_r(:) => null()
    integer(kind=8), pointer :: noeud_occ(:) => null()
    integer(kind=8), pointer :: noeud_rela(:) => null()
    aster_logical :: lcolle
!
    call jemarq()
!
    eps1 = 1.d4*r8prem()
    eps2 = 1.d0/r8gaem()
!
! - Mesh and model
!
    call dismoi('NOM_MODELE', charge, 'CHARGE', repk=mod)
    call dismoi('NOM_LIGREL', mod, 'MODELE', repk=ligrmo)
    call dismoi('NOM_MAILLA', ligrmo, 'LIGREL', repk=noma)
    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', inom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', nbcmp)
    nddla = nbcmp-1
    ASSERT(nddla .le. nmocl)
!
    do i = 1, nbcmp
        nomcmp(i) = zk8(inom-1+i)
    end do
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
    call jeveuo(ligrmo//'.PRNM', 'L', jprnm)
!
! --- ACCES AUX COMPOSANTES DE LA LISTE_RELA
!
    call jeveuo(lisrel//'.RLCO', 'E', jrlco)
    call jeveuo(lisrel//'.RLDD', 'E', jrldd)
    call jeveuo(lisrel//'.RLNO', 'E', jrlno)
    call jeveuo(lisrel//'.RLNT', 'E', vi=rlnt)
    call jeveuo(lisrel//'.RLPO', 'E', vi=rlpo)
    call jeveuo(lisrel//'.RLSU', 'E', vi=rlsu)
    call jeveuo(lisrel//'.RLTC', 'L', vk8=rltc)
!
! --- TYPE DE VALEUR DES COEFFICIENTS DES RELATIONS ---
!
    typcoe = rltc(1) (1:4)
!
! --- NOMBRE DE RELATIONS DE LA LISTE_RELA
!
    call jeveuo(lisrel//'.RLNR', 'L', vi=rlnr)
    nbrela = rlnr(1)
!
! --- NOMBRE DE TERMES  MAX IMPLIQUES DANS UNE RELATION
!
    nbtema = 0
    do irela = 1, nbrela
        if (nbtema .lt. rlnt(irela)) nbtema = rlnt(irela)
    end do
!
! --- CREATION D'UN VECTEUR DE TRAVAIL DESTINE A CONTENIR
! --- L'INDICE DU PLUS GRAND COEFFICIENT EN VALEUR ABSOLUE
! --- (MODULE) D'UNE RELATION
!
    AS_ALLOCATE(vi=coefmax, size=nbrela)
!
! --- CREATION D'UN VECTEUR DE TRAVAIL DESTINE A CONTENIR
! --- LES NUMEROS DES NOEUDS D'UNE RELATION
!
    AS_ALLOCATE(vi=noeud_rela, size=nbtema)
!
! --- CREATION D'UN VECTEUR DE TRAVAIL DESTINE A CONTENIR
! --- LE NOMBRE D'OCCURENCES DE CHAQUE NOEUD APPRAISSANT
! --- DANS UNE RELATION
!
    AS_ALLOCATE(vi=noeud_occ, size=nbtema)
!
! --- CREATION D'UN VECTEUR DE TRAVAIL DESTINE A CONTENIR
! --- LES COEFFICIENTS REELS D'UNE RELATION
!
    AS_ALLOCATE(vr=coef_r, size=nbtema)
!
! --- CREATION D'UN VECTEUR DE TRAVAIL DESTINE A CONTENIR
! --- LES COEFFICIENTS COMPLEXES D'UNE RELATION
!
    AS_ALLOCATE(vc=coef_c, size=nbtema)
!
!
!
!     0. ON ORDONNE LES TERMES DE CHAQUE RELATION POUR POUVOIR
!        LES COMPARER PLUS FACILEMENT ET DETECTER LES DOUBLONS
!     ----------------------------------------------------------
    do irela = 1, nbrela
        ipntrl = rlpo(irela)
        nbterm = rlnt(irela)
        idecal = ipntrl-nbterm
        jrlcof = jrlco+idecal
        idnoeu = jrlno+idecal
        iddl = jrldd+idecal
!
        if (typcoe .eq. 'COMP') then
            do ino = 1, nbterm
                coef_c(ino) = zc(jrlcof+ino-1)
            end do
        else if (typcoe .eq. 'REEL') then
            do ino = 1, nbterm
                coef_r(ino) = zr(jrlcof+ino-1)
            end do
        else
            ASSERT(.false.)
        end if
!
        do ino = 1, nbterm
            nomnoe = zk8(idnoeu+ino-1)
            in = char8_to_int(nomnoe, lcolle, noma, "NOEUD")
            noeud_rela(ino) = in
            cmp = zk8(iddl+ino-1)
            icmp = indik8(nomcmp, cmp, 1, nbcmp)
            if (.not. exisdg(zi(jprnm-1+(in-1)*nbec+1), icmp)) then
                valk(1) = cmp
                valk(2) = nomnoe
                call utmess('F', 'CHARGES2_31', nk=2, valk=valk)
            end if
        end do
!
! ----- Rearrangement of linear relation tables in ascending order of nodes and dof for given node
!
        call ordrel(noeud_rela, zk8(idnoeu), zk8(iddl), coef_r, coef_c, &
                    noeud_occ, nbterm, zk8(inom), nddla)
!
!       -- REAFFECTATION DU TABLEAU DES COEFFICIENTS
        if (typcoe .eq. 'COMP') then
            do ino = 1, nbterm
                zc(jrlcof+ino-1) = coef_c(ino)
            end do
        else if (typcoe .eq. 'REEL') then
            do ino = 1, nbterm
                zr(jrlcof+ino-1) = coef_r(ino)
            end do
        else
            ASSERT(.false.)
        end if
!
        coemax = 0.0d0
        if (typcoe .eq. 'COMP') then
            do ino = 1, nbterm
                if (abs(coef_c(ino)) .gt. coemax) then
                    coemax = abs(coef_c(ino))
                    indmax = ino
                end if
            end do
        else if (typcoe .eq. 'REEL') then
            do ino = 1, nbterm
                if (abs(coef_r(ino)) .gt. coemax) then
                    coemax = abs(coef_r(ino))
                    indmax = ino
                end if
            end do
        else
            ASSERT(.false.)
        end if
        coefmax(irela) = indmax
    end do
!
!
!
!   1. IDENTIFICATION DES RELATIONS REDONDANTES A 1 TERME
!   ----------------------------------------------------------------
    call jecreo('&&ORDLRL.KIDREL', 'V N K16')
    call jeecra('&&ORDLRL.KIDREL', 'NOMMAX', nbrela)
    do irela1 = nbrela, 1, -1
        nbter1 = rlnt(irela1)
        if (nbter1 .le. 1) then
            ipntr1 = rlpo(irela1)
            ideca1 = ipntr1-nbter1
            idnoe1 = jrlno+ideca1
            iddl1 = jrldd+ideca1
            kidrel = zk8(idnoe1)//zk8(iddl1)
            call jenonu(jexnom('&&ORDLRL.KIDREL', kidrel), nidrel)
            if (nidrel .eq. 0) then
                call jecroc(jexnom('&&ORDLRL.KIDREL', kidrel))
            else
                rlsu(irela1) = 1
            end if
        end if
    end do
    call jedetr('&&ORDLRL.KIDREL')
!
!
!
!   2. SETTING TO ZERO TERMS BELOW NUMERICAL PRECISION
!   --------------------------------------------------
    do irela1 = 1, nbrela
        nbter1 = rlnt(irela1)
        if (nbter1 .eq. 1) cycle
        ipntr1 = rlpo(irela1)
        ideca1 = ipntr1-nbter1
        jrlco1 = jrlco+ideca1
!
        indmax = coefmax(irela1)
!
!       set to 0 terms below numerical precision
        do ino = 1, nbter1
            if (typcoe .eq. 'COMP') then
                if (abs(zc(jrlco1+ino-1)) < eps1*abs(zc(jrlco1+indmax-1))) then
                    zc(jrlco1+ino-1) = dcmplx(0.d0, 0.d0)
                end if
            else if (typcoe .eq. 'REEL') then
                if (abs(zr(jrlco1+ino-1)) < eps1*abs(zr(jrlco1+indmax-1))) then
                    zr(jrlco1+ino-1) = 0.d0
                end if
            else
                ASSERT(.false.)
            end if
        end do
    end do

!   3. IDENTIFICATION DES RELATIONS REDONDANTES A PLUSIEURS TERMES
!   ----------------------------------------------------------------
    do irela1 = nbrela, 2, -1
        nbter1 = rlnt(irela1)
        if (nbter1 .eq. 1) goto 170
        ipntr1 = rlpo(irela1)
        ideca1 = ipntr1-nbter1
        jrlco1 = jrlco+ideca1
        idnoe1 = jrlno+ideca1
        iddl1 = jrldd+ideca1
!
        indmax = coefmax(irela1)
!
        if (typcoe .eq. 'COMP') then
            if (abs(zc(jrlco1+indmax-1)) .lt. eps2) then
                call utmess('F', 'CHARGES2_32')
            end if
        else if (typcoe .eq. 'REEL') then
            if (abs(zr(jrlco1+indmax-1)) .lt. eps2) then
                call utmess('F', 'CHARGES2_32')
            end if
        else
            ASSERT(.false.)
        end if
!
!
!       --  CAS DES COEF. COMPLEXES
!       -----------------------------------
        if (typcoe .eq. 'COMP') then
            do irela2 = 1, irela1-1
                nbter2 = rlnt(irela2)
                ipntr2 = rlpo(irela2)
                ideca2 = ipntr2-nbter2
                jrlco2 = jrlco+ideca2
                idnoe2 = jrlno+ideca2
                iddl2 = jrldd+ideca2
                coproc = zc(jrlco2+indmax-1)/zc(jrlco1+indmax-1)
!
                if (nbter1 .eq. nbter2) then
                    icomp = 0
                    do ino = 1, nbter1
                        if (zk8(idnoe1+ino-1) .eq. zk8(idnoe2+ino-1)) then
                            if (zk8(iddl1+ino-1) .eq. zk8(iddl2+ino-1)) then
                                rapcoc = coproc*zc(jrlco1+ino-1)
                                epsrel = eps1*abs(zc(jrlco1+indmax-1))
                                difrel = abs(zc(jrlco2+ino-1)-rapcoc)
                                if (difrel .le. epsrel) goto 110
                                icomp = 1
                                goto 120
                            else
                                icomp = 1
                                goto 120
                            end if
                        else
                            icomp = 1
                            goto 120
                        end if
110                     continue
                    end do
120                 continue
                    if (icomp .eq. 0) rlsu(irela2) = 1
                end if
            end do
!
!
!       --  CAS DES COEF. REEL
!       -----------------------------------
        else if (typcoe .eq. 'REEL') then
            do irela2 = 1, irela1-1
                nbter2 = rlnt(irela2)
                ipntr2 = rlpo(irela2)
                ideca2 = ipntr2-nbter2
                jrlco2 = jrlco+ideca2
                idnoe2 = jrlno+ideca2
                iddl2 = jrldd+ideca2
                copror = zr(jrlco2+indmax-1)/zr(jrlco1+indmax-1)
!
                if (nbter1 .eq. nbter2) then
                    icomp = 0
                    do ino = 1, nbter1
                        if (zk8(idnoe1+ino-1) .eq. zk8(idnoe2+ino-1)) then
                            if (zk8(iddl1+ino-1) .eq. zk8(iddl2+ino-1)) then
                                rapcoe = copror*zr(jrlco1+ino-1)
                                epsrel = eps1*abs(zr(jrlco1+indmax-1))
                                difrel = abs(zr(jrlco2+ino-1)-rapcoe)
                                if (difrel .le. epsrel) goto 140
                                icomp = 1
                                goto 150
                            else
                                icomp = 1
                                goto 150
                            end if
                        else
                            icomp = 1
                            goto 150
                        end if
140                     continue
                    end do
150                 continue
                    if (icomp .eq. 0) rlsu(irela2) = 1
                end if
            end do
        else
            ASSERT(.false.)
        end if
170     continue
    end do
!
!
! ---  MENAGE  ---
    AS_DEALLOCATE(vi=noeud_rela)
    AS_DEALLOCATE(vi=noeud_occ)
    AS_DEALLOCATE(vi=coefmax)
    AS_DEALLOCATE(vr=coef_r)
    AS_DEALLOCATE(vc=coef_c)
!
    call jedema()
end subroutine
