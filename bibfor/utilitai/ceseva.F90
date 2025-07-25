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

subroutine ceseva(cesf, npara, lpara, cesr)
! person_in_charge: jacques.pellet at edf.fr
! A_UTIL
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/cestas.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: npara
    character(len=*) :: cesf, lpara(npara), cesr
! ---------------------------------------------------------------------
! BUT: EVALUER LE CHAM_ELEM_S DE FONCTIONS CESF EN UTILISANT
!      LES PARAMETRES TROUVES DANS LES CHAM_ELEM_S LPARA
! ---------------------------------------------------------------------
! ARGUMENTS:
! CESF  IN/JXIN  K19 : SD CHAM_ELEM_S A EVALUER
! NPARA IN       I   : NOMBRE DE CHAM_ELEM_S PARAMETRES (LPARA)
! LPARA IN/JXIN  V(K19) : LISTE DES CHAM_ELEM_S PARAMETRES
! CESR  IN/JXOUT K19  : SD CHAM_ELEM RESULTAT DE L'EVALUATION
!
! REMARQUES :
!  EN CHAQUE POINT DE DISCRETISATION DE CESF, ON FERA "CALL FOINTE"
!  POUR EVALUER LES FONCTIONS AFFECTEES A CE POINT.
!  ON PASSERA EN PARAMETRES DES FONCTIONS, LES VALEURS DES CHAMPS
!  DE LPARA AVEC COMME NOM DE PARAMETRE LE NOM DE LA CMP
!  ON NE TRAITE QUE LES FONCTIONS REELLES F : R * R(* R,...) -> R
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: jfd, jfl
    integer(kind=8) :: jpd, jpc, jpv, jpl
    integer(kind=8) :: jrd, jrc, jrl, jrk
    integer(kind=8) :: nbma, k, ima, ncmp, nbpu, ier, nbpumx, ibid
    integer(kind=8) :: ncmp2, ipara, ncmpmx, nspmx, nptmx
    integer(kind=8) :: k2, iadf, iadr, iadp, nbpt, nbsp, ipt, isp, jnompu, jvalpu
    character(len=8) :: ma, nomgdf, nomgdr, fo
    character(len=8) :: ma2, nomgd2, typces
    character(len=3) :: tsca
    character(len=19) :: f, p, r
    character(len=24) :: valk
    real(kind=8) :: x
    character(len=8), pointer :: fv(:) => null()
    real(kind=8), pointer :: rv(:) => null()
    character(len=8), pointer :: fc(:) => null()
    character(len=8), pointer :: fk(:) => null()
    character(len=8), pointer :: pk(:) => null()
    integer(kind=8), pointer :: vjad1(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
!     -- 2 VECTEURS POUR STOCKER LE NOM ET LES VALEURS DES PARAMETRES
!        DES FONCTIONS :
    nbpumx = 10
    call wkvect('&&CESEVA.NOMPU', 'V V K8', nbpumx, jnompu)
    call wkvect('&&CESEVA.VALPU', 'V V R', nbpumx, jvalpu)
!
!     1- RECUPERATIONS D'INFOS DANS LE CHAMP DE FONCTIONS :
!     ------------------------------------------------------------
    f = cesf
    call jeveuo(f//'.CESK', 'L', vk8=fk)
    call jeveuo(f//'.CESD', 'L', jfd)
    call jeveuo(f//'.CESC', 'L', vk8=fc)
    call jeveuo(f//'.CESV', 'L', vk8=fv)
    call jeveuo(f//'.CESL', 'L', jfl)
!
    ma = fk(1)
    nomgdf = fk(2)
    typces = fk(3)
    nbma = zi(jfd-1+1)
    ncmp = zi(jfd-1+2)
    nptmx = zi(jfd-1+3)
    nspmx = zi(jfd-1+4)
    ncmpmx = zi(jfd-1+5)
!
    call dismoi('TYPE_SCA', nomgdf, 'GRANDEUR', repk=tsca)
    if (tsca .ne. 'K8') then
        call utmess('F', 'UTILITAI_16')
    end if
!
!     2- ALLOCATION DU CHAM_ELEM_S RESULTAT ET RECUPERATION
!        DES ADRESSES DE SES OBJETS   :
!     ------------------------------------------------------------
    r = cesr
    nomgdr = nomgdf(1:4)//'_R'
    call cescre('V', r, typces, ma, nomgdr, &
                ncmp, fc, [-nptmx], [-nspmx], [-ncmpmx])
    call jeveuo(r//'.CESK', 'L', jrk)
    call jeveuo(r//'.CESD', 'L', jrd)
    call jeveuo(r//'.CESC', 'L', jrc)
    call jeveuo(r//'.CESV', 'E', vr=rv)
    call jeveuo(r//'.CESL', 'E', jrl)
!
!
!     3- ON MET EN MEMOIRE LES OBJETS UTILES DES CHAMPS PARAMETRES :
!     --------------------------------------------------------------
    AS_ALLOCATE(vi=vjad1, size=4*npara)
    do ipara = 1, npara
        p = lpara(ipara)
        call jeveuo(p//'.CESK', 'L', vk8=pk)
        call jeveuo(p//'.CESD', 'L', jpd)
        call jeveuo(p//'.CESC', 'L', jpc)
        call jeveuo(p//'.CESV', 'L', jpv)
        call jeveuo(p//'.CESL', 'L', jpl)
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
!     4- EVALUATION DES FONCTIONS :
!     ---------------------------------------
!     ON BOUCLE D'ABORD SUR LES CMPS POUR AVOIR PLUS DE CHANCES
!     DE FAIRE PLUSIEURS FOINTE SUCCESSIFS AVEC LA MEME FONCTION.
!
    do k = 1, ncmp
        do ima = 1, nbma
            nbpt = zi(jfd-1+5+4*(ima-1)+1)
            nbsp = zi(jfd-1+5+4*(ima-1)+2)
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    call cesexi('C', jfd, jfl, ima, ipt, &
                                isp, k, iadf)
                    if (iadf .le. 0) goto 40
!
                    fo = fv(iadf)
!
                    call cesexi('C', jrd, jrl, ima, ipt, &
                                isp, k, iadr)
                    ASSERT(iadr .lt. 0)
                    zl(jrl-1-iadr) = .true.
!
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
                            call cesexi('C', jpd, jpl, ima, ipt, &
                                        isp, k2, iadp)
                            if (iadp .le. 0) goto 20
!
                            nbpu = nbpu+1
                            if (nbpu .gt. nbpumx) then
!                    -- ON AGRANDIT .NOMPU ET .VALPU :
                                nbpumx = 2*nbpumx
                                call juveca('&&CESEVA.NOMPU', nbpumx)
                                call juveca('&&CESEVA.VALPU', nbpumx)
                                call jeveuo('&&CESEVA.NOMPU', 'E', jnompu)
                                call jeveuo('&&CESEVA.VALPU', 'E', jvalpu)
                            end if
!
!                 -- ON VERIFIE QU'UN MEME PARAMETRE N'EST PAS AJOUTE
!                    PLUSIEURS FOIS:
                            ibid = indik8(zk8(jnompu), zk8(jpc-1+k2), 1, &
                                          nbpu-1)
                            if (ibid .gt. 0) then
                                call utmess('F', 'CALCULEL2_78', sk=zk8(jpc-1+k2))
                            end if
!
                            zk8(jnompu-1+nbpu) = zk8(jpc-1+k2)
                            zr(jvalpu-1+nbpu) = zr(jpv-1+iadp)
20                          continue
                        end do
                    end do
!
!
!           4.2 APPEL A FOINTE :
!           --------------------
                    call fointe('E', fo, nbpu, zk8(jnompu), zr(jvalpu), &
                                x, ier)
                    if (ier .ne. 0) then
                        call utmess('F+', 'FONCT0_9', sk=fo)
                        valk = int_to_char8(ima)
                        call utmess('F', 'FONCT0_10', sk=valk)
                    end if
!
!           4.3 STOCKAGE DU RESULTAT :
!           --------------------------
                    rv(1-1-iadr) = x
!
40                  continue
                end do
            end do
        end do
    end do
!
    call cestas(cesr)
!
!     5- MENAGE :
!     ---------------------------------------
    AS_DEALLOCATE(vi=vjad1)
    call jedetr('&&CESEVA.NOMPU')
    call jedetr('&&CESEVA.VALPU')
!
    call jedema()
end subroutine
