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

subroutine memaxm(typmx, champ, nocmp, nbcmp, lcmp, &
                  vr, nbmail, numail)
! aslint: disable=W1306
    implicit none
#include "asterc/r8maem.h"
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"
    character(len=*) :: typmx
    character(len=*) :: champ, nocmp, lcmp(*)
    integer(kind=8) :: nbcmp, nbmail, numail(*)
    real(kind=8) :: vr(*)
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
! BUT :  EXTRAIRE LE "MIN/MAX" DE COMPOSANTES
!        D'UN CHAMP (CHAM_ELEM OU CARTE) SUIVANT LA COMPOSANTE NOCMP
!
! EXEMPLE :
!   CONSIDERONS 2 MAILLES M1 ET M2 AVEC LES VALEURS D'UN CHAMP
!   (DX,DY,DZ):
!           SUR LA MAILLE M1 ---> (5,3,1)
!           SUR LA MAILLE M2 ---> (1,4,3)
!   SI LA COMPOSANTE CRITERE EST NOMCMP='DZ',ET QUE L'UTILISA-
!   TEUR DEMANDE LES VALEURS DES COMPOSANTES DX ET DY DU CHAMP
!   SUR L'ELEMENT OU LE MAX EST ATTEINT,LA FONCTION RETOURNERA:
!           VR=(1,4)
!   SI LA COMPOSANTE CRITERE EST NOMCMP='DZ',ET QUE L'UTILISA-
!   TEUR DEMANDE LA VALEUR DE LA COMPOSANTE DX DU CHAMP SUR
!   L'ELEMENT OU LE MIN EST ATTEINT,LA FONCTION RETOURNERA:
!           VR=5
!   SI LA COMPOSANTE CRITERE EST NOMCMP='DY',ET QUE L'UTILISA-
!   TEUR DEMANDE LA VALEUR DES COMPOSANTES DX,DY ET DZ DU CHAMP
!   SUR L'ELEMENT OU LE MAX EST ATTEINT,LA FONCTION RETOURNERA:
!           VR=(1,4,3)
!
! IN  : TYPMX  :  'MIN'/'MAX'/'MIN_ABS'/'MAX_ABS'
! IN  : CHAMP  :  NOM DU CHAMP A SCRUTER (VALEURS REELLES OU ENTIERES)
! IN  : NOCMP  :  NOM DE LA COMPOSANTE SUR LAQUELLE ON FAIT LE TEST
! IN  : NBCMP  :  NOMBRE DE COMPOSANTES DEMANDEES (=LONGUEUR DE VR)
! IN  : LICMP  :  NOM DES COMPOSANTES DEMANDEES PAR L'UTILISATEUR
! IN  : NBMAIL :  = 0   , COMPARAISON SUR TOUT LE MAILLAGE
!                 SINON , COMPARAISON SUR UNE PARTIE DU MAILLAGE
! IN  : NUMAIL :  NUMEROS DES MAILLES SUR LESQUELLES ON EFFECTUE LES
!                 COMPARAISONS (SI NBMAIL>0)
! OUT : VR     :  VECTEUR CONTENANT LES VALEURS DES COMPOSANTES DU CHAMP
!                 SUR L'ELEMENT (OU NOEUD OU POINT DE GAUSS) OU LE
!                 'MIN'/'MAX' EST ATTEINT
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: iret
    integer(kind=8) :: longt
    character(len=8) :: kmpic, typ1, nomgd, tsca, tych, mesh
    integer(kind=8) :: jcesd, jcesl, jcesv, nel, iel, nbpt, nbsspt, ncmp
    integer(kind=8) :: ipt, isp, icmp, ncp, iicmp, iadr1
    integer(kind=8) :: iadr2, iel1
    real(kind=8) :: valr, vmima
    character(len=19) :: chams, cham19
    integer(kind=8) :: tncomp(nbcmp)
    aster_logical :: copi, lmax, labs, lreel
    character(len=8), pointer :: cesc(:) => null()
    character(len=8), pointer :: cesk(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    cham19 = champ
!
!     -- TRANSFORMATION DU CHAMP EN CHAMP SIMPLE :
!     --------------------------------------------
    chams = '&&MEMAXM.CES'
    call dismoi('TYPE_CHAMP', cham19, 'CHAMP', repk=tych)
    if (tych(1:2) .eq. 'EL') then
        call celces(cham19, 'V', chams)
    else if (tych .eq. 'CART') then
        call carces(cham19, 'ELEM', ' ', 'V', chams, &
                    ' ', iret)
        ASSERT(iret .eq. 0)
    else
        ASSERT(.false.)
    end if
    call jelira(chams//'.CESV', 'TYPE', cval=typ1)
    ! ASSERT(typ1 .eq. 'R')
!
!
!
    call jeveuo(chams//'.CESD', 'L', jcesd)
    call jeveuo(chams//'.CESL', 'L', jcesl)
    call jeveuo(chams//'.CESC', 'L', vk8=cesc)
    call jeveuo(chams//'.CESK', 'L', vk8=cesk)
    call jeveuo(chams//'.CESV', 'L', jcesv)
!
    mesh = cesk(1)
    nomgd = cesk(2)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R' .or. tsca .eq. 'I')
    lreel = tsca .eq. 'R'
!
!
!
!     -- INITIALISATION DE VMIMA :
    ASSERT(typmx(1:3) .eq. 'MIN' .or. typmx(1:3) .eq. 'MAX')
    lmax = typmx(1:3) .eq. 'MAX'
    if (.not. lmax) vmima = +r8maem()
    if (len(typmx) .gt. 3) then
        ASSERT(len(typmx) .eq. 7)
        ASSERT(typmx(4:7) .eq. '_ABS')
        labs = .true.
        if (lmax) vmima = 0.d0
    else
        labs = .false.
        if (lmax) vmima = -r8maem()
    end if
!
!
!     INITIALISATION DE TNCOMP CONTENANT LES INDICES
!     DES CMP
!     ----------------------------------
    do icmp = 1, nbcmp
        tncomp(icmp) = 0
    end do
!
    ncmp = zi(jcesd-1+2)
    do icmp = 1, ncmp
        do iicmp = 1, nbcmp
            if (lcmp(iicmp) .eq. cesc(icmp)) then
                tncomp(iicmp) = icmp
            end if
        end do
    end do
!
!
!     COMPARAISON NOCMP AVEC TTES LES
!     AUTRES AFIN DE RECUPERER LE NUM DE LA COMPOSANTE
!     RECUPERE L'INDEX DE LA COMPOSANTE A TESTER DANS LE CHAMP
    ncp = 0
    do icmp = 1, ncmp
        if (cesc(icmp) .eq. nocmp) ncp = icmp
    end do
!
!     -- CAS : TOUTES LES MAILLES :
!     -----------------------------
    if (nbmail .le. 0) then
!       NOMBRE D'ELEMENTS DU MAILLAGE
        nel = zi(jcesd-1+1)
!     -- CAS : LISTE DE MAILLES :
!     ---------------------------
    else
        nel = nbmail
    end if
!
!
    do iel = 1, nel
!
        if (nbmail .le. 0) then
            iel1 = iel
        else
            iel1 = numail(iel)
        end if
!
!       NOMBRE DE PTS ET SSPTS POUR CHAQUE ELEMENT
        nbpt = zi(jcesd-1+5+4*(iel1-1)+1)
        nbsspt = zi(jcesd-1+5+4*(iel1-1)+2)
        ncmp = zi(jcesd-1+5+4*(iel1-1)+3)
!
!
        do ipt = 1, nbpt
            do isp = 1, nbsspt
                call cesexi('C', jcesd, jcesl, iel1, ipt, &
                            isp, ncp, iadr1)
                if (iadr1 .gt. 0) then
                    if (lreel) then
                        valr = zr(jcesv-1+iadr1)
                    else
                        valr = zi(jcesv-1+iadr1)
                    end if
                    if (labs) valr = abs(valr)
                    copi = .false.
                    if ((lmax) .and. (valr .gt. vmima)) copi = .true.
                    if ((.not. lmax) .and. (valr .lt. vmima)) copi = .true.
                    if (copi) then
                        vmima = valr
                        do iicmp = 1, nbcmp
                            call cesexi('C', jcesd, jcesl, iel1, ipt, &
                                        isp, tncomp(iicmp), iadr2)
                            if (iadr2 .eq. 0) then
                                vr(iicmp) = r8nnem()
                            else
                                if (lreel) then
                                    vr(iicmp) = zr(jcesv-1+iadr2)
                                else
                                    vr(iicmp) = zi(jcesv-1+iadr2)
                                end if
                            end if
                        end do
                    end if
                end if
            end do
        end do
    end do
!
!
    call detrsd('CHAMP', chams)
!
!     -- IL FAUT PARFOIS COMMUNIQUER LE RESULTAT ENTRE LES PROCS :
    call dismoi('MPI_COMPLET', champ, 'CHAMP', repk=kmpic)
    if (kmpic .eq. 'NON' .or. isParallelMesh(mesh)) then
        if (lmax) then
            call asmpi_comm_vect('MPI_MAX', 'R', nbval=longt, vr=vr)
        else
            call asmpi_comm_vect('MPI_MIN', 'R', nbval=longt, vr=vr)
        end if
    end if
!
    call jedema()
end subroutine
