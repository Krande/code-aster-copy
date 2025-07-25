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
subroutine vprecu(modes, nomsy, nbvect, lposi, nomvec, &
                  nbpara, nopara, nomvai, nomvar, nomvak, &
                  neq, nbmode, typmod, nbpari, nbparr, &
                  nbpark)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/irparb.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: modes, nomsy, nomvec, nopara
    character(len=*), intent(out) :: typmod
    character(len=*), intent(in) :: nomvai, nomvar, nomvak
    integer(kind=8), intent(in) :: nbvect, lposi(*), nbpara
    integer(kind=8), intent(out) :: neq, nbmode, nbpari, nbparr, nbpark
!     RECUPERATION DES VALEURS ET VECTEURS PROPRES
!     ------------------------------------------------------------------
! IN  MODES  : K8 : NOM DE LA STRUCTURE MODE_MECA
! IN  NOMSY  : K16: NOM DU CHAMP A RECUPERER
! IN  NBVECT : IS : NOMBRE DE MODE A RECUPERER
!                   =0   AUCUN CHAMP A RECUPERER
!                   >0   RECUPERATION DE LA LISTE PAR NUMERO D'ORDRE
!                   <0   RECUPERATION DE TOUS LES MODES
! IN  LPOSI  : IS : LISTE DES NUMEROS D'ORDRE
! IN  NOMVEC : K24: NOM DU TABLEAU CONTENANT LES VECTEURS PROPRES
! IN  NBPARA : IS : NOMBRE DE PARAMETRES A RECUPERER
!                   =0   AUCUN PARAMETRE A RECUPERER
!                   >0   RECUPERATION A PARTIR D'UNE LISTE
!                   <0   RECUPERATION DE TOUS LES PARAMETRES
! VAR NOPARA : K16: LISTE DES NOMS DES PARAMETRES
! IN  NOMVAI : K24: NOM DU TABLEAU CONTENANT LES PARAMETRES ENTIERS
! IN  NOMVAR : K24: NOM DU TABLEAU CONTENANT LES PARAMETRES REELS
! IN  NOMVAK : K24: NOM DU TABLEAU CONTENANT LES PARAMETRES CHARACTERS
! OUT NEQ    : IS : NOMBRE D'EQUATIONS
! OUT NBMODE : IS : NOMBRE DE MODE(S) RECUPERE(S)
! OUT NBPARI : IS : NOMBRE DE PARAMETRE(S) ENTIER(S) RECUPERE(S)
! OUT NBPARR : IS : NOMBRE DE PARAMETRE(S) REEL(S) RECUPERE(S)
! OUT NBPARK : IS : NOMBRE DE PARAMETRE(S) CHARACTER(S) RECUPERE(S)
! OUT TYPMOD : K  : TYPE DU CHAMP RECUPERE
!     ------------------------------------------------------------------
!
!     *** ATTENTION ***
!
!     LES VECTEURS PROPRES SONT LES UNS DERRIERE LES AUTRES
!         ACCES PAR   ZR(IDDL,IMODE)
!
!     LES PARAMETRES MODAUX SONT PAR CATEGORIES
!         ACCES PAR   ZR(IMODE,IPARA)
!
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    real(kind=8) :: rbid(2)
    character(len=4) :: type
    character(len=8) :: k8b
    character(len=24) :: vale, nomjv
    character(len=24) :: valk(2)
    complex(kind=8) :: c16b
    aster_logical :: recunp
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ieq, ii, ik, imode, ir, tmod(1)
    integer(kind=8) :: iret, j, jpara, lmode, lnopar, lnume, lnumor
    integer(kind=8) :: lresui, lresuk, lresur, lvale, nbmodt, nbout, nbtpar
    integer(kind=8) :: nbtrou, neq1, nordr
!-----------------------------------------------------------------------
    data vale/'                   .VALE'/
!     ------------------------------------------------------------------
    call jemarq()
!
    if (nbvect .lt. 0) then
!        --- TOUS LES MODES ---
        call rsorac(modes, 'LONUTI', ibid, rbid(1), k8b, &
                    c16b, 0.0d0, k8b, tmod, 1, &
                    nbtrou)
        nbmode = tmod(1)
        call wkvect('&&VPRECU.NUMERO.ORDRE', 'V V I', nbmode, lnumor)
        call rsorac(modes, 'TOUT_ORDRE', ibid, rbid(1), k8b, &
                    c16b, 0.0d0, k8b, zi(lnumor), nbmode, &
                    nbtrou)
    else if (nbvect .gt. 0) then
!        --- A PARTIR D'UNE LISTE DE NUMEROS D'ORDRE ---
        nbmode = nbvect
        call wkvect('&&VPRECU.NUMERO.ORDRE', 'V V I', nbmode, lnumor)
        do i = 1, nbvect
            zi(lnumor+i-1) = lposi(i)
        end do
!
    else
!        --- RIEN ---
        nbmode = 0
        typmod = '?'
        call rsorac(modes, 'LONUTI', ibid, rbid(1), k8b, &
                    c16b, 0.0d0, k8b, tmod, 1, &
                    nbtrou)
        nbmodt = tmod(1)
        call wkvect('&&VPRECU.NUMERO.ORDRE', 'V V I', nbmodt, lnumor)
        call rsorac(modes, 'TOUT_ORDRE', ibid, rbid(1), k8b, &
                    c16b, 0.0d0, k8b, zi(lnumor), nbmodt, &
                    nbtrou)
        goto 100
    end if
!     ------------------------------------------------------------------
!
!          *************** ON RECUPERE LES CHAMPS ***************
!
!     --- RECUPERATION DE NEQ ---
    call rsexch('F', modes, nomsy, zi(lnumor), vale(1:19), &
                iret)
    call jeexin(vale(1:19)//'.VALE', ibid)
    if (ibid .gt. 0) then
        vale(20:24) = '.VALE'
    else
        vale(20:24) = '.CELV'
    end if
!
    call jelira(vale, 'LONMAX', neq)
    call jelira(vale, 'TYPE', cval=typmod)
!
!     --- CREATION DES OBJETS DE NOM NOMVEC ET NOMVAL ---
    if (typmod(1:1) .eq. 'R') then
        call wkvect(nomvec, 'V V R', neq*nbmode, lmode)
    else if (typmod(1:1) .eq. 'C') then
        call wkvect(nomvec, 'V V C', neq*nbmode, lmode)
    else
        valk(1) = typmod(1:1)
        call utmess('F', 'ALGELINE4_80', sk=valk(1))
    end if
!
!        --- VECTEUR PROPRE ---
    do imode = 1, nbmode
        nordr = zi(lnumor-1+imode)
        call rsexch('F', modes, nomsy, nordr, vale(1:19), &
                    iret)
        call jeexin(vale(1:19)//'.VALE', ibid)
        if (ibid .gt. 0) then
            vale(20:24) = '.VALE'
        else
            vale(20:24) = '.CELV'
        end if
!
        call jeveuo(vale, 'L', lvale)
        call jelira(vale, 'LONMAX', neq1)
        call jelira(vale, 'TYPE', cval=k8b)
        if (typmod(1:1) .ne. k8b(1:1)) then
            call utmess('F', 'ALGELINE3_70')
        else if (neq .eq. neq1) then
            if (typmod(1:1) .eq. 'R') then
                do ieq = 0, neq-1
                    zr(lmode+neq*(imode-1)+ieq) = zr(lvale+ieq)
                end do
            else if (typmod(1:1) .eq. 'C') then
                do ieq = 0, neq-1
                    zc(lmode+neq*(imode-1)+ieq) = zc(lvale+ieq)
                end do
            end if
            call jelibe(vale)
        else
            call utmess('F', 'ALGELINE3_71')
        end if
    end do
100 continue
!     ------------------------------------------------------------------
!
!        *************** ON RECUPERE LES PARAMETRES ***************
!
    recunp = .false.
    if (nbpara .lt. 0) then
!        --- TOUS LES PARAMETRES ---
        nomjv = '&&VPRECU.NOM_PARA'
        call irparb(modes, nbpara, nopara, nomjv, nbout)
        call jeveuo(nomjv, 'L', jpara)
        recunp = .true.
    else if (nbpara .gt. 0) then
!        --- A PARTIR D'UNE LISTE DE PARAMETRES ---
        nomjv = '&&VPRECU.NOM_PARA'
        call irparb(modes, nbpara, nopara, nomjv, nbout)
        call jeveuo(nomjv, 'L', jpara)
!
    else
!        --- RIEN ---
        nbpari = 0
        nbparr = 0
        nbpark = 0
        goto 200
    end if
!
    nbpari = 0
    nbparr = 0
    nbpark = 0
    do i = 1, nbout
        call rsadpa(modes, 'L', 1, zk16(jpara+i-1), zi(lnumor), &
                    i, sjv=lnume, styp=type, istop=0)
        if (type(1:1) .eq. 'I') then
            nbpari = nbpari+1
        else if (type(1:1) .eq. 'R') then
            nbparr = nbparr+1
        else if (type(1:1) .eq. 'K') then
            nbpark = nbpark+1
        else
        end if
    end do
    if (recunp) then
        nbtpar = nbpari+nbparr+nbpark
        call wkvect(nopara, 'V V K16', nbtpar, lnopar)
    end if
!
    if (nbpari .ne. 0) call wkvect(nomvai, 'V V I', nbpari*nbmode, lresui)
    if (nbparr .ne. 0) call wkvect(nomvar, 'V V R', nbparr*nbmode, lresur)
    if (nbpark .ne. 0) call wkvect(nomvak, 'V V K24', nbpark*nbmode, lresuk)
!
    ii = 0
    ir = 0
    ik = 0
    do i = 1, nbout
        do j = 1, nbmode
            nordr = zi(lnumor-1+j)
            call rsadpa(modes, 'L', 1, zk16(jpara+i-1), nordr, &
                        i, sjv=lnume, styp=type, istop=0)
            if (type(1:1) .eq. 'I') then
                ii = ii+1
                zi(lresui+ii-1) = zi(lnume)
                if (recunp .and. j .eq. 1) then
                    zk16(lnopar+ii-1) = zk16(jpara+i-1)
                end if
            else if (type(1:1) .eq. 'R') then
                ir = ir+1
                zr(lresur+ir-1) = zr(lnume)
                if (recunp .and. j .eq. 1) then
                    zk16(lnopar+nbpari+ir-1) = zk16(jpara+i-1)
                end if
            else if (type(1:1) .eq. 'K') then
                ik = ik+1
                if (type(2:3) .eq. '24') then
                    zk24(lresuk+ik-1) = zk24(lnume)
                else
                    zk24(lresuk+ik-1) = zk16(lnume)//'        '
                end if
                if (recunp .and. j .eq. 1) then
                    zk16(lnopar+nbpari+nbparr+ik-1) = zk16(jpara+i-1)
                end if
            end if
        end do
    end do
200 continue
!
!     --- DESTRUCTION DES OBJET DE TRAVAIL ---
    call jeexin('&&VPRECU.NOM_PARA', iret)
    if (iret .ne. 0) call jedetr('&&VPRECU.NOM_PARA')
    call jeexin('&&VPRECU.NUMERO.ORDRE', iret)
    if (iret .ne. 0) call jedetr('&&VPRECU.NUMERO.ORDRE')
!
    call jedema()
end subroutine
