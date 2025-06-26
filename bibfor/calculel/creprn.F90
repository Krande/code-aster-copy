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

subroutine creprn(ligrez, molocz, basez, prnmz, prnsz)
! person_in_charge: jacques.pellet at edf.fr
!
! aslint: disable=
    implicit none
!-----------------------------------------------------------------
!  BUT : CREATION :
!      .DU VECTEUR PRNM DES ENTIERS CODES DECRIVANT
!       LA NATURE DES DDLS DES NOEUDS PHYSIQUES DU LIGREL.
!
!      .DU VECTEUR PRNS DES ENTIERS CODES DECRIVANT
!       LES NOEUDS LAGRANGE DU LIGREL
!       (S'IL Y EN A)
!-----------------------------------------------------------------
!   ARGUMENT        E/S  TYPE         ROLE
!    LIGREZ          IN    K*     NOM DU LIGREL
!    MOLOCZ          IN    K*   / NOM DU MODE_LOCAL PERMETTANT
!                                 DE CHOISIR LES DDLS.
!                               / ' '
!    BASEZ           IN    K*     NOM DE LA BASE
!    PRNMZ      IN/JXOUT    K24   NOM DU VECTEUR DES ENTIERS CODES
!                                 DECRIVANT LA NATURE DES DDLS DES
!                                 NOEUDS PHYSIQUES
!   (CE VECTEUR EST TOUJOURS CREE)
!    PRNSZ      IN/JXOUT    K24   NOM DU VECTEUR DES ENTIERS CODES
!                                 DECRIVANT LES NOEUDS LAGRANGE
!   (CE VECTEUR N'EST CREE QUI SI LIGREL CONTIENT DES NOEUDS TARDIFS)
!-----------------------------------------------------------------
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/entcod.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/nbno.h"
#include "asterfort/typele.h"
#include "asterfort/wkvect.h"
!
!
! -----  VARIABLES LOCALES
    character(len=*) :: ligrez, molocz, basez, prnmz, prnsz
    integer(kind=8) :: gd, i, iamaco, iamail, iamsco, jmoloc, iancmp
    integer(kind=8) ::  iaprnm, iaprno, iaprns, icmp
    integer(kind=8) :: icodla, iec, igr, illiel, ilmaco, ilmsco
    integer(kind=8) :: ima, imode, ino, inold, iret, ite, j, k, l, lgncmp, nbnm
    integer(kind=8) :: nbnoms, nbsma, nbssa, nec, nel, nl, nm, nnoe, numa, nunoel
    integer(kind=8) :: admodl, lcmodl
    integer(kind=8) :: lshift
    character(len=1) :: base
    character(len=8) :: noma, nomgd, exiel, nomacr, moloc
    character(len=16) :: nomte
    character(len=14) :: num2
    character(len=16) :: phenom
    character(len=19) :: ligrel
    character(len=24) :: prnm, prns
    integer(kind=8), pointer :: liel(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
!
! -----  FONCTIONS FORMULES
!     NUMAIL(IGR,IEL)=NUMERO DE LA MAILLE ASSOCIEE A L'ELEMENT IEL
#define numail(igr,iel) liel(zi(illiel+igr-1)+iel-1)
!     NUMGLM(IMA,INO)=NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
!                     IMA ETANT UNE MAILLE DU MAILLAGE.
#define numglm(ima,ino) zi(iamaco-1+zi(ilmaco+ima-1)+ino-1)
!     NUMGLS(IMA,INO)=NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
!                     IMA ETANT UNE MAILLE SUPPLEMENTAIRE DU LIGREL
#define numgls(ima,ino) zi(iamsco-1+zi(ilmsco+ima-1)+ino-1)
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    base = basez
    moloc = molocz
    ligrel = ligrez
    prnm = prnmz
    prns = prnsz
    nm = 0
    nl = 0
    nbnoms = 0
    nbsma = 0
    nbssa = 0
!
!
    call dismoi('EXI_ELEM', ligrel, 'LIGREL', repk=exiel)
    call dismoi('NB_SS_ACTI', ligrel, 'LIGREL', repi=nbssa)
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=noma)
!
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nm)
    call dismoi('NB_NL_MAILLA', noma, 'MAILLAGE', repi=nl)
    call dismoi('NB_SM_MAILLA', noma, 'MAILLAGE', repi=nbsma)
!
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', lcmodl)
    call jeveuo(jexnum('&CATA.TE.MODELOC', 1), 'L', admodl)
!
    if (moloc .eq. ' ') then
        call dismoi('PHENOMENE', ligrel, 'LIGREL', repk=phenom)
        call dismoi('NOM_MOLOC', phenom, 'PHENOMENE', repk=moloc)
    end if
!
!
!     -- DETERMINATION DE LA GRANDEUR NOMGD A PARTIR DE MOLOC :
!     ---------------------------------------------------------
    if (exiel(1:3) .eq. 'OUI') then
        do igr = 1, nbgrel(ligrel)
            ite = typele(ligrel, igr)
            call jenuno(jexnum('&CATA.TE.NOMTE', ite), nomte)
            call jenonu(jexnom('&CATA.TE.NOMMOLOC', nomte//moloc), imode)
            if (imode .gt. 0) then
                call jeveuo(jexnum('&CATA.TE.MODELOC', imode), 'L', jmoloc)
                call jenuno(jexnum('&CATA.GD.NOMGD', zi(jmoloc-1+2)), nomgd)
                goto 20
            end if
        end do
!       -- IL PEUT ARRIVER QUE NBGREL=0. ON S'EN SORT AVEC MOLOC :
        if (moloc .eq. 'DDL_MECA') then
            nomgd = 'DEPL_R'
        else if (moloc .eq. 'DDL_THER') then
            nomgd = 'TEMP_R'
        else
            ASSERT(.false.)
        end if
20      continue
!
    else
!       -- SI IL N'Y A PAS D'ELEMENTS FINIS
!          ON EST EN SOUS-STRUCTURATION STATIQUE => MECANIQUE.
        call dismoi('NOM_GD', 'MECANIQUE', 'PHENOMENE', repk=nomgd)
!
!
    end if
!
! --- SI QUE l'ON RENTRE AVEC MODE_LOCAL
    if (molocz .ne. ' ') then
! --- Pour certains RESU_ELEM, il faut changer le nom
! Il en manque sûrement (voir nulili.F90)
        if (nomgd(1:4) == "MDEP" .or. nomgd(1:4) == "VDEP" &
            .or. nomgd(1:4) == "MDNS") then
            nomgd = "DEPL_"//nomgd(6:6)
        elseif (nomgd(1:4) == "MTEM" .or. nomgd(1:4) == "VTEM" &
                .or. nomgd(1:4) == "MTNS") then
            nomgd = "TEMP_"//nomgd(6:6)
        elseif (nomgd(1:4) == "MPRE" .or. nomgd(1:4) == "VPRE") then
            nomgd = "PRES_"//nomgd(6:6)
        elseif (nomgd(1:4) == "MSIZ" .or. nomgd(1:4) == "VSIZ") then
            nomgd = "SIZZ_"//nomgd(6:6)
        elseif (nomgd(1:4) == "MZNS" .or. nomgd(1:4) == "VNEU") then
            nomgd = "NEUT_"//nomgd(6:6)
        end if
    end if
!
!     -- CALCUL DE GD ET NEC :
!     ---------------------------------------------------------
    call dismoi('NUM_GD_SI', nomgd, 'GRANDEUR', repi=gd)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
!
!
    call jeexin(noma//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(noma//'.CONNEX', 'L', iamaco)
        call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', ilmaco)
    else
        iamaco = 1
        ilmaco = 1
    end if
!
    call jeexin(ligrel(1:19)//'.NEMA', iret)
    if (iret .gt. 0) then
        call jeveuo(ligrel(1:19)//'.NEMA', 'L', iamsco)
        call jeveuo(jexatr(ligrel(1:19)//'.NEMA', 'LONCUM'), 'L', ilmsco)
    else
        iamsco = 1
        ilmsco = 1
    end if
!
!
! --- ALLOCATION DE PRNM :
!     ------------------
    call jeexin(prnm, iret)
    if (iret .eq. 0) then
        call wkvect(prnm, base//' V I', (nm+nl)*nec, iaprnm)
    else
        call jeveuo(prnm, 'L', iaprnm)
    end if

! - allocation de prns (pour un ligrel contenant des noeuds tardifs):
!   ----------------------------------------------------------------
    call dismoi('NB_NO_SUP', ligrel, 'LIGREL', repi=nbnoms)
    call jeexin(prns, iret)
    if (iret .eq. 0) then
        if (nbnoms .gt. 0) call wkvect(prns, base//' V I', nbnoms*nec, iaprns)
    else
        call jeveuo(prns, 'L', iaprns)
    end if

! - traitement des elements finis classiques :
!   ----------------------------------------
    if (exiel(1:3) .eq. 'NON') goto 90
    call jeveuo(ligrel(1:19)//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel(1:19)//'.LIEL', 'LONCUM'), 'L', illiel)

    do igr = 1, nbgrel(ligrel)
!
! ---   calcul de imode (mode_local) :
!       ------------------------------
        ite = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', ite), nomte)
        call jenonu(jexnom('&CATA.TE.NOMMOLOC', nomte//moloc), imode)

        if (imode .gt. 0) then
            nnoe = nbno(imode)
            nel = nbelem(ligrel, igr)
            do j = 1, nel
                numa = numail(igr, j)
                if (numa .gt. 0) then

!                   -- il s'agit d'une maille physique du maillage :
!                   ------------------------------------------------
                    do k = 1, nnoe
                        nunoel = numglm(numa, k)
                        do l = 1, nec
                            iec = entcod(admodl, lcmodl, nec, imode, k, l)
                            zi(iaprnm-1+nec*(nunoel-1)+l) = ior( &
                                                            zi(iaprnm-1+nec*(nunoel-1)+l), &
                                                            iec &
                                                            )
                        end do
                    end do
                else

!                   -- il s'agit d'une maille tardive :
!                   -----------------------------------
                    numa = -numa
                    do k = 1, nnoe
                        nunoel = numgls(numa, k)
                        do l = 1, nec
                            iec = entcod(admodl, lcmodl, nec, imode, k, l)
                            if (nunoel .gt. 0) then
                                zi(iaprnm-1+nec*(nunoel-1)+l) = &
                                    ior(zi(iaprnm-1+nec*(nunoel-1)+l), &
                                        iec)
                            else
                                zi(iaprns-1+nec*(-nunoel-1)+l) = &
                                    ior(zi(iaprns-1+nec*(-nunoel-1)+l), &
                                        iec)
                            end if
                        end do
                    end do
                end if
            end do
        end if
    end do

90  continue

! --- BOUCLE SUR LES SUPERELEMENTS :
!     ----------------------------
    if (nbssa .gt. 0) then
!
        call jeveuo(ligrel//'.SSSA', 'L', vi=sssa)
!
! ---   le seul ddl porte par un noeud de lagrange est 'lagr' :
!       -------------------------------------------------------
        call jeveuo(noma//'.NOMACR', 'L', vk8=vnomacr)
!
        call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iancmp)
        call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', lgncmp)
        icmp = indik8(zk8(iancmp), 'LAGR', 1, lgncmp)
! on ne trouve pas la composante "LAGR" dans la grandeur
        ASSERT(icmp .ne. 0)
! il est imprévu d avoir la composante "LAGR" au delà de 30
        ASSERT(icmp .le. 30)
!
        icodla = lshift(1, icmp)
!
        do ima = 1, nbsma
            nomacr = vnomacr(ima)
            call dismoi('NOM_NUME_DDL', nomacr, 'MACR_ELEM_STAT', repk=num2)
            call jeveuo(nomacr//'.CONX', 'L', vi=conx)
            call jeveuo(jexnum(num2//'.NUME.PRNO', 1), 'L', iaprno)
            if (sssa(ima) .eq. 1) then
                call jeveuo(jexnum(noma//'.SUPMAIL', ima), 'L', iamail)
                call jelira(jexnum(noma//'.SUPMAIL', ima), 'LONMAX', nbnm)
!
                do i = 1, nbnm
                    ino = zi(iamail-1+i)
                    inold = conx(3*(i-1)+2)
                    if (ino .gt. nm) then
!
! ---                   CAS D'UN NOEUD DE LAGRANGE :
!                       --------------------------
                        zi(iaprnm-1+nec*(ino-1)+1) = ior(zi(iaprnm-1+nec*(ino-1)+1), icodla)

                    else if (inold .gt. 0) then
!
! ---                   CAS D'UN NOEUD PHYSIQUE DU MAILLAGE :
!                       -----------------------------------
                        do iec = 1, nec
                            zi(iaprnm-1+nec*(ino-1)+iec) = ior( &
                                                           zi(iaprnm-1+nec*(ino-1)+iec), &
                                                           zi(iaprno-1+(nec+2)*(inold-1)+2+iec))
                        end do
                    else
! on traite un super-élément  et le noeud courant n'est ni un noeud Lagrange,
! ni un noeud physique du maillage.
                        ASSERT(.false.)
                    end if
                end do
            end if
        end do
    end if
!
    call jedema()
!
end subroutine
