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
subroutine mecact(base, nomcar, moclez, nomco, nomgdz, &
                  ncmp, nomcmp, si, sr, sc, &
                  sk, lnomcmp, vi, vr, vc, &
                  vk)
! aslint: disable=W1306
!-----------------------------------------------------------------------
    implicit none
!     CREER 1 CARTE CONSTANTE SUR 1 MODELE.
!-----------------------------------------------------------------------
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
    character(len=*), intent(in) :: base, nomcar, nomco
    character(len=*), intent(in) :: nomgdz, moclez
    character(len=8) :: nomgd
    character(len=6) :: mocle
    integer(kind=8), intent(in) :: ncmp
    character(len=*), intent(in), optional :: nomcmp, sk, lnomcmp(ncmp), vk(ncmp)
    integer(kind=8), intent(in), optional :: si, vi(ncmp)
    real(kind=8), intent(in), optional :: sr, vr(ncmp)
    complex(kind=8), intent(in), optional :: sc, vc(ncmp)
!
!   Liste des variables locales utilisées pour récuperer les arguments
    character(len=24) :: licmp(ncmp), kcmp(ncmp)
    integer(kind=8) :: icmp(ncmp)
    real(kind=8) :: rcmp(ncmp)
    complex(kind=8) :: ccmp(ncmp)
! ----------------------------------------------------------------------
!     ENTREES:
!      BASE       : BASE DE CREATION POUR LA CARTE (G/V/L)
!      NOMCAR     : NOM DE LA CARTE A CREER (SI ELLE EXISTE ON LA DETRUIT).
!      MOCLEZ     : 'MAILLA' , 'MODELE' OU 'LIGREL'
!      NOMCO      : NOM DU MAILLAGE SUPPORT DE LA CARTE (SI MOCLE='MAILLA')
!                 : NOM DU MODELE SUPPORT DE LA CARTE (SI MOCLE='MODELE')
!                 : NOM DU LIGREL SUPPORT DE LA CARTE.(SI MOCLE= 'LIGREL')
!                 : K19 SI MOCLEZ = LIGREL, K8 SINON
!      NOMGDZ     : NOM DE LA GRANDEUR ASSOCIEE A LA CARTE.
!      NCMP       : NOMBRE DE CMP A EDITER SUR LA CARTE.
! NOMCMP/LNOMCMP  : LISTE DES NOMS DE CMP A EDITER.
!      SI/VI      : LISTE DES VALEURS ENTIERES DES CMP A EDITER.(EVENTUEL).
!      SR/VR      : LISTE DES VALEURS REELLES  DES CMP A EDITER.(EVENTUEL).
!      SC/VC      : LISTE DES VALEURS COMPLEX  DES CMP A EDITER.(EVENTUEL).
!      SK/VK      : LISTE DES VALEURS CHAR*8,16,24  DES CMP A EDITER.(EVT).
!
!     SORTIES:
!       NOMCAR : EST REMPLI.
!
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
    character(len=24) :: nommo2, nomca2
!
!     VARIABLES LOCALES:
!     ------------------
    character(len=1) :: type, bas2
    character(len=8) :: noma
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ianoma, iret, jncmp, jvalv, ltyp, j
!
!-----------------------------------------------------------------------
    call jemarq()
!
!   DETERMINATION DU TYPE DE CARTE A AFFECTER
    ASSERT(EXCLUS2(si, vi))
    ASSERT(EXCLUS2(sr, vr))
    ASSERT(EXCLUS2(sc, vc))
    ASSERT(EXCLUS2(sk, vk))
    if (UN_PARMI2(si, vi)) then
        j = 1
    else if (UN_PARMI2(sr, vr)) then
        j = 2
    else if (UN_PARMI2(sc, vc)) then
        j = 3
    else
        ASSERT(UN_PARMI2(sk, vk))
        j = 4
    end if
!
!   AFFECTATION DES VARIABLES LOCALES EN FONCTION DES ARGUMENTS
    if (ncmp .gt. 1) then
        ASSERT(UN_PARMI4(vi, vr, vc, vk))
        do i = 1, ncmp
            licmp(i) = lnomcmp(i)
        end do
        select case (j)
        case (1)
            do i = 1, ncmp
                icmp(i) = vi(i)
            end do
        case (2)
            do i = 1, ncmp
                rcmp(i) = vr(i)
            end do
        case (3)
            do i = 1, ncmp
                ccmp(i) = vc(i)
            end do
        case (4)
            do i = 1, ncmp
                kcmp(i) = vk(i)
            end do
        end select
    else
        ASSERT(UN_PARMI4(si, sr, sc, sk))
        licmp(1) = nomcmp
        select case (j)
        case (1)
            icmp(1) = si
        case (2)
            rcmp(1) = sr
        case (3)
            ccmp(1) = sc
        case (4)
            kcmp(1) = sk
        end select
    end if
!
!
    mocle = moclez
    nomgd = nomgdz
!
    bas2 = base
!
!     -- RECUPERATION DU NOM DU MAILLAGE:
    nommo2 = nomco
    nomca2 = nomcar
    if (mocle(1:6) .eq. 'MAILLA') then
        noma = nommo2(1:8)
    else if (mocle(1:6) .eq. 'MODELE') then
        call jeveuo(nommo2(1:8)//'.MODELE    .LGRF', 'L', ianoma)
        noma = zk8(ianoma-1+1)
    else if (mocle(1:6) .eq. 'LIGREL') then
        call jeveuo(nommo2(1:19)//'.LGRF', 'L', ianoma)
        noma = zk8(ianoma-1+1)
    else
        ASSERT(.false.)
    end if
!
!     -- SI LA CARTE EXISTE DEJA , ON LA DETRUIT COMPLETEMENT:
!
    call jeexin(nomca2(1:19)//'.NOMA', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.NOMA')
    call jeexin(nomca2(1:19)//'.NOLI', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.NOLI')
    call jeexin(nomca2(1:19)//'.DESC', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.DESC')
    call jeexin(nomca2(1:19)//'.LIMA', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.LIMA')
    call jeexin(nomca2(1:19)//'.VALE', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.VALE')
    call jeexin(nomca2(1:19)//'.NCMP', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.NCMP')
    call jeexin(nomca2(1:19)//'.VALV', iret)
    if (iret .gt. 0) call jedetr(nomca2(1:19)//'.VALV')
!
!     -- ON ALLOUE LA CARTE:
!
    call alcart(bas2, nomca2, noma, nomgd)
    call jeveuo(nomca2(1:19)//'.NCMP', 'E', jncmp)
    call jeveuo(nomca2(1:19)//'.VALV', 'E', jvalv)
    call jelira(nomca2(1:19)//'.VALV', 'TYPE', cval=type)
    call jelira(nomca2(1:19)//'.VALV', 'LTYP', ltyp)
    do i = 1, ncmp
        zk8(jncmp-1+i) = licmp(i)
        if (type(1:1) .eq. 'R') then
            zr(jvalv-1+i) = rcmp(i)
        end if
        if (type(1:1) .eq. 'C') then
            zc(jvalv-1+i) = ccmp(i)
        end if
        if (type(1:1) .eq. 'I') then
            zi(jvalv-1+i) = icmp(i)
        end if
        if (type(1:1) .eq. 'K') then
            if (ltyp .eq. 8) then
                zk8(jvalv-1+i) = kcmp(i)
            else if (ltyp .eq. 16) then
                zk16(jvalv-1+i) = kcmp(i)
            else if (ltyp .eq. 24) then
                zk24(jvalv-1+i) = kcmp(i)
            else
                ASSERT(.false.)
            end if
        end if
    end do
!
!     -- ON NOTE DANS LA CARTE LES VALEURS VOULUES :
!
    call nocart(nomca2, 1, ncmp)
!
    call jedetr(nomca2(1:19)//'.VALV')
    call jedetr(nomca2(1:19)//'.NCMP')
    call jedema()
end subroutine
