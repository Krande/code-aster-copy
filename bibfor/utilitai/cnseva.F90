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

subroutine cnseva(cnsf, npara, lpara, cnsr)
! person_in_charge: jacques.pellet at edf.fr
! A_UTIL
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/cnscre.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: npara
    character(len=*) :: cnsf, lpara(npara), cnsr
! ---------------------------------------------------------------------
! BUT: EVALUER LE CHAM_NO_S DE FONCTIONS CNSF EN UTILISANT
!      LES PARAMETRES TROUVES DANS LES CHAM_NO_S LPARA
! ---------------------------------------------------------------------
! ARGUMENTS:
! CNSF  IN/JXIN  K19 : SD CHAM_NO_S A EVALUER
! NPARA IN       I   : NOMBRE DE CHAM_NO_S PARAMETRES (LPARA)
! LPARA IN/JXIN  V(K19) : LISTE DES CHAM_NO_S PARAMETRES
! CNSR  IN/JXOUT K19  : SD CHAM_NO RESULTAT DE L'EVALUATION
!
! REMARQUES :
!  EN CHAQUE POINT DE DISCRETISATION DE CNSF, ON FERA "CALL FOINTE"
!  POUR EVALUER LES FONCTIONS AFFECTEES A CE POINT.
!  ON PASSERA EN PARAMETRES DES FONCTIONS, LES VALEURS DES CHAMPS
!  DE LPARA AVEC COMME NOM DE PARAMETRE LE NOM DE LA CMP
!
!  ON NE TRAITE QUE LES FONCTIONS REELLES F : R * R(* R,...) -> R
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) ::    jfl
    integer(kind=8) :: jpd, jpc, jpv, jpl
    integer(kind=8) :: jrd, jrc, jrl, jrk
    integer(kind=8) :: nbno, k, ino, ncmp, nbpu, ier, nbpumx
    integer(kind=8) :: k2, ncmp2, ipara, ibid
    parameter(nbpumx=50)
    character(len=8) :: ma, nomgdf, nomgdr, fo, nompu(nbpumx)
    character(len=8) :: ma2, nomgd2
    character(len=3) :: tsca
    character(len=19) :: f, p, r
    character(len=24) :: valk
    real(kind=8) :: x, valpu(nbpumx)
    integer(kind=8), pointer :: fd(:) => null()
    character(len=8), pointer :: fc(:) => null()
    character(len=8), pointer :: fk(:) => null()
    character(len=8), pointer :: pk(:) => null()
    character(len=8), pointer :: fv(:) => null()
    real(kind=8), pointer :: rv(:) => null()
    integer(kind=8), pointer :: vjad1(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     1- RECUPERATIONS D'INFOS DANS LE CHAMP DE FONCTIONS :
!     ------------------------------------------------------------
    f = cnsf
    call jeveuo(f//'.CNSK', 'L', vk8=fk)
    call jeveuo(f//'.CNSD', 'L', vi=fd)
    call jeveuo(f//'.CNSC', 'L', vk8=fc)
    call jeveuo(f//'.CNSV', 'L', vk8=fv)
    call jeveuo(f//'.CNSL', 'L', jfl)
!
    ma = fk(1)
    nomgdf = fk(2)
    nbno = fd(1)
    ncmp = fd(2)
!
    call dismoi('TYPE_SCA', nomgdf, 'GRANDEUR', repk=tsca)
    if (tsca .ne. 'K8') then
        call utmess('F', 'UTILITAI_16')
    end if
!
!
!     2- ALLOCATION DU CHAM_NO_S RESULTAT ET RECUPERATION
!        DES ADRESSES DE SES OBJETS   :
!     ------------------------------------------------------------
    r = cnsr
    nomgdr = nomgdf(1:4)//'_R'
    call cnscre(ma, nomgdr, ncmp, fc, 'V', &
                r)
    call jeveuo(r//'.CNSK', 'L', jrk)
    call jeveuo(r//'.CNSD', 'L', jrd)
    call jeveuo(r//'.CNSC', 'L', jrc)
    call jeveuo(r//'.CNSV', 'E', vr=rv)
    call jeveuo(r//'.CNSL', 'E', jrl)
!
!
!     3- ON MET EN MEMOIRE LES OBJETS UTILES DES CHAMPS PARAMETRES :
!     --------------------------------------------------------------
    AS_ALLOCATE(vi=vjad1, size=4*npara)
    do ipara = 1, npara
        p = lpara(ipara)
        call jeveuo(p//'.CNSK', 'L', vk8=pk)
        call jeveuo(p//'.CNSD', 'L', jpd)
        call jeveuo(p//'.CNSC', 'L', jpc)
        call jeveuo(p//'.CNSV', 'L', jpv)
        call jeveuo(p//'.CNSL', 'L', jpl)
        ma2 = pk(1)
        nomgd2 = pk(2)
!
        call dismoi('TYPE_SCA', nomgd2, 'GRANDEUR', repk=tsca)
        if (tsca .ne. 'R') then
            call utmess('F', 'UTILITAI_17')
        end if
        if (ma2 .ne. ma) then
            call utmess('F', 'UTILITAI_18')
        end if
        vjad1(4*(ipara-1)+1) = jpc
        vjad1(4*(ipara-1)+2) = jpd
        vjad1(4*(ipara-1)+3) = jpl
        vjad1(4*(ipara-1)+4) = jpv
    end do
!
!
!
!     4- EVALUATION DES FONCTIONS :
!     ---------------------------------------
!     ON BOUCLE D'ABORD SUR LES CMPS POUR AVOIR PLUS DE CHANCES
!     DE FAIRE PLUSIEURS FOINTE SUCCESSIFS AVEC LA MEME FONCTION.
!
    do k = 1, ncmp
        do ino = 1, nbno
            if (zl(jfl-1+(ino-1)*ncmp+k)) then
                zl(jrl-1+(ino-1)*ncmp+k) = .true.
                fo = fv((ino-1)*ncmp+k)
                if (fo .eq. ' ') goto 40
!
!           4.1 FABRICATION DE LA LISTE DES PARAMETRES POUR FOINTE:
!           -------------------------------------------------------
                nbpu = 0
                do ipara = 1, npara
                    jpc = vjad1(4*(ipara-1)+1)
                    jpd = vjad1(4*(ipara-1)+2)
                    jpl = vjad1(4*(ipara-1)+3)
                    jpv = vjad1(4*(ipara-1)+4)
                    ncmp2 = zi(jpd-1+2)
                    do k2 = 1, ncmp2
                        if (zl(jpl-1+(ino-1)*ncmp2+k2)) then
                            nbpu = nbpu+1
                            if (nbpu .gt. nbpumx) then
                                call utmess('F', 'CALCULEL2_66')
                            end if
!
!                 -- ON VERIFIE QU'UN MEME PARAMETRE N'EST PAS AJOUTE
!                    PLUSIEURS FOIS:
                            ibid = indik8(nompu, zk8(jpc-1+k2), 1, nbpu-1)
                            if (ibid .gt. 0) then
                                call utmess('F', 'CALCULEL2_78', sk=zk8(jpc-1+k2))
                            end if
!
                            nompu(nbpu) = zk8(jpc-1+k2)
                            valpu(nbpu) = zr(jpv-1+(ino-1)*ncmp2+k2)
                        end if
                    end do
                end do
!
!
!           4.2 APPEL A FOINTE :
!           --------------------
                call fointe('E', fo, nbpu, nompu, valpu, &
                            x, ier)
                if (ier .ne. 0) then
                    call utmess('F+', 'FONCT0_9', sk=fo)
                    valk = int_to_char8(ino)
                    call utmess('F', 'FONCT0_53', sk=valk)
                end if
!
!           4.3 STOCKAGE DU RESULTAT :
!           --------------------------
                rv((ino-1)*ncmp+k) = x
!
            end if
40          continue
        end do
    end do
!
!
!     5- MENAGE :
!     ---------------------------------------
    AS_DEALLOCATE(vi=vjad1)
!
    call jedema()
end subroutine
