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
subroutine vpstor(ineg, typ, modes, nbmode, neq, &
                  vecpr8, vecpc8, mxresf, nbpari, nbparr, &
                  nbpark, nopara, mod45, resufi, resufr, &
                  resufk, iprec, modelz, matez, caraz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/indk24.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jenonu.h"
#include "asterfort/juveca.h"
#include "asterfort/refdaj.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsexis.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/vtcreb.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
!
    integer(kind=8) :: ineg, nbmode, neq, mxresf, nbpari, nbparr, nbpark
    integer(kind=8) :: iprec, resufi(mxresf, *)
    character(len=4) :: mod45
    character(len=*) :: typ, modes, resufk(mxresf, *), nopara(*)
    real(kind=8) :: vecpr8(neq, *), resufr(mxresf, *)
    complex(kind=8) :: vecpc8(neq, *)
    character(len=8), intent(in), optional :: modelz, matez, caraz
!     STOCKAGE DES VALEURS PROPRES
!
!     REMARQUE:
!        DANS NOPARA, ON A LES NOMS DE PARAMETRES DE TYP ENTIER
!                     ENSUITE LES NOMS DE PARAMETRES DE TYP CHARACTER
!                     ENSUITE LES NOMS DE PARAMETRES DE TYP REEL
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: imode, jmode, ier, nmin, imin, nmax, imax
    integer(kind=8) :: vali(3), jpara
    integer(kind=8) :: nmin1, kmode, nordr, iarg, i, ladpa, lmode, lvale
    integer(kind=8) :: nbpast, irang, iret, jmodg
    integer(kind=8) :: jmod2, igd, jrefe
    parameter(nbpast=19)
    character(len=8) :: res, k8b, modele, chmat, carael, basemo, mesh, nomgd
    character(len=16) :: typcon, nomcmd, nosy, typmod
    character(len=19) :: chamno, sd2
    character(len=24) :: nume, nopast(nbpast)
    character(len=24) :: valk, typeba, raide, raide2, k24b
    aster_logical :: lrefd, lbasm, lstock
    integer(kind=8), pointer :: p_desc(:) => null()
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
! --- PARAMETRES STOCKES DANS LA SD RESULTAT DYNAMIQUE
    data nopast/'NUME_MODE',&
     &  'NORME', 'TYPE_MODE', 'NOEUD_CMP',&
     &  'FREQ', 'OMEGA2', 'AMOR_REDUIT',&
     &  'MASS_GENE', 'RIGI_GENE', 'AMOR_GENE',&
     &  'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &  'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &  'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY', 'MASS_EFFE_UN_DZ'/
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call getres(res, typcon, nomcmd)
!
!     POUR POUVOIR UTILISER VPSTOR DANS STAT_NON_LINE VIA NMOP45
    if (typcon .eq. 'EVOL_NOLI') then
        typcon = 'MODE_FLAMB'
        if (mod45 .eq. 'VIBR') typcon = 'MODE_MECA'
        if (mod45 .eq. 'STAB') typcon = 'MODE_STAB'
    end if
!
    if (typcon .eq. 'MODE_ACOU') then
        nosy = 'PRES'
    else
        nosy = 'DEPL'
    end if
!
    lbasm = .false.
    lrefd = .true.
    lstock = .false.
!
    call jeexin(modes(1:8)//'           .REFD', iret)
    if (iret .eq. 0) then
        call refdaj(' ', modes, -1, ' ', 'INIT', &
                    ' ', iret)
        lrefd = .false.
    end if
!
    typeba = ' '
    call dismoi('TYPE_BASE', modes, 'RESU_DYNA', repk=typeba, arret='C', &
                ier=iret)
    if (typeba(1:1) .ne. ' ') lbasm = .true.
    if (lbasm) then
        call getvid(' ', 'RAIDE', scal=raide, nbret=ier)
        if (ier .ne. 0) then
            call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=nume)
        else
            call dismoi('REF_RIGI_PREM', modes, 'RESU_DYNA', repk=raide, arret='C')
            call dismoi('NUME_DDL', modes, 'RESU_DYNA', repk=nume, arret='C')
        end if
    else
        call dismoi('REF_RIGI_PREM', modes, 'RESU_DYNA', repk=k24b, arret='C')
        raide = k24b(1:8)
        call exisd('MATR_ASSE', raide, iret)
        if (iret .ne. 0) then
            call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=nume)
            lstock = .true.
        else
            call dismoi('NUME_DDL', modes, 'RESU_DYNA', repk=nume, arret='C', &
                        ier=iret)
        end if
    end if
!
!     --- CONTROLE PREALABLE ---
    do imode = 1, nbmode
        jmode = resufi(imode, 1)
        if (jmode .lt. 1 .and. ineg .gt. 0) then
            call utmess('A', 'ALGELINE3_79')
        end if
    end do
!
!     --- STOCKAGE DES MODES ---
    call rsexis(modes, ier)
    if (ier .eq. 0) then
        call utmess('F', 'ALGELINE3_80')
    end if
!
    nmin = resufi(1, 1)
    imin = 1
    nmax = resufi(1, 1)
    imax = 1
    do imode = 2, nbmode
        if (resufi(imode, 1) .lt. nmin) then
            nmin = resufi(imode, 1)
            imin = imode
        end if
        if (resufi(imode, 1) .gt. nmax) then
            nmax = resufi(imode, 1)
            imax = imode
        end if
    end do
    nmin1 = nmax
!
!     ON RECUPERE LE NOM DE LA MATRICE DE RAIDEUR AFIN DE
!     DETERMINER LE NOM DU MODELE, DU MATERIAU ET DES
!     CARACTERISTIQUES ELEMENTAIRES
    modele = '        '
    chmat = '        '
    carael = '        '
    if (lstock) then
        if (present(modelz)) then
            modele = modelz
            chmat = matez
            carael = caraz
            go to 39
        end if
        if (typcon(1:9) .eq. 'MODE_MECA' .or. typcon(1:9) .eq. 'MODE_ACOU' .or. &
            typcon(1:10) .eq. 'MODE_FLAMB' .or. typcon(1:9) .eq. 'MODE_STAB') then
            call dismoi('NOM_MODELE', raide, 'MATR_ASSE', repk=modele)
        else if (typcon(1:9) .eq. 'MODE_GENE') then
!            ON EST PASSE PAR UN PROJ_MATR_BASE
            call jeveuo(raide(1:19)//'.REFA', 'L', jmodg)
            basemo = zk24(jmodg) (1:8)
            if (basemo(1:8) .ne. '        ') then
                call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=raide2, arret='C')
                if (raide2 .eq. ' ') then
                    call jeveuo(jexnum(basemo//'           .TACH', 1), 'L', jmod2)
                    sd2 = zk24(jmod2) (1:8)
                    call rsadpa(sd2, 'L', 1, 'MODELE', 1, &
                                0, sjv=jpara, styp=k8b)
                    modele = zk8(jpara)
                    call rsadpa(sd2, 'L', 1, 'CHAMPMAT', 1, &
                                0, sjv=jpara, styp=k8b)
                    chmat = zk8(jpara)
                    call rsadpa(sd2, 'L', 1, 'CARAELEM', 1, &
                                0, sjv=jpara, styp=k8b)
                    carael = zk8(jpara)
                    goto 39
                end if
            end if
        end if
        if (chmat .eq. '        ') then
            call getvid(' ', 'CHAM_MATER', scal=chmat, nbret=ier)
            if (ier .eq. 0) chmat = '        '
        end if
        if (carael .eq. '        ') then
            call getvid(' ', 'CARA_ELEM', scal=carael, nbret=ier)
            if (ier .eq. 0) carael = '        '
        end if
    end if
!
39  continue
!
    do imode = 1, nbmode
!
        if ((typcon .ne. 'MODE_FLAMB') .or. (nomcmd .eq. 'NORM_MODE')) then
!           STOCKAGE DES FREQUENCES PAR ORDRE CROISSANT DE NUMERO
            if (imode .eq. 1) then
                kmode = imin
            else if (imode .eq. nbmode) then
                kmode = imax
            else
                do lmode = 1, nbmode
                    if (resufi(lmode, 1) .gt. nmin .and. resufi(lmode, 1) .lt. nmin1) then
                        nmin1 = resufi(lmode, 1)
                        kmode = lmode
                    end if
                end do
                nmin = nmin1
                nmin1 = nmax
            end if
!
            jmode = resufi(kmode, 1)
        else
!           stockage par ordre croissant de charge critique
!           cad par ordre decroissant de frequence
!           Rq : dans ce cas resufi(imode,1) = imode (ce qui simplifie les choses)
            kmode = nbmode-imode+1
            jmode = kmode
        end if
!
        nordr = iprec+imode
!
!        --- VECTEUR PROPRE ---
        call rsexch(' ', modes, nosy, nordr, chamno, &
                    ier)
        if (ier .eq. 0) then
            continue
        else if (ier .eq. 100 .and. lrefd) then
            call vtcreb(chamno, 'G', typ(1:1), nume_ddlz=nume, nb_equa_outz=neq)
        else
            vali(1) = kmode
            vali(2) = jmode
            vali(3) = ier
            valk = chamno
            call utmess('F', 'ALGELINE4_85', sk=valk, ni=3, vali=vali)
        end if
        if (typcon .eq. 'MODE_GENE' .or. typcon .eq. 'HARM_GENE') then
! GLUTE CAR ON A UTILISE VTCRE[ABM] POUR UN CHAM_GENE QUI A UN .REFE
! DE TAILLE 2 ET NON 4 COMME UN CHAM_NO ET PAS DE .DESC
            call wkvect(chamno//'.DESC', 'G V I', 2, vi=p_desc)
            call dismoi("NOM_GD", nume, "NUME_DDL", repk=nomgd)
            call jenonu(jexnom('&CATA.GD.NOMGD', nomgd(1:5)//typ(1:1)), igd)
            p_desc(1) = igd
            p_desc(2) = 1
            call jeecra(chamno//'.DESC', 'DOCU', iarg, 'VGEN')
            call jeecra(chamno//'.REFE', 'DOCU', iarg, 'VGEN')
            call juveca(chamno//'.REFE', 2)
!
            call dismoi("NOM_MAILLA", nume, "NUME_DDL", repk=mesh)
            call jeveuo(chamno//'.REFE', 'E', jrefe)
            zk24(jrefe) = mesh
            zk24(jrefe+1) = nume
        end if
        call jeveuo(chamno//'.VALE', 'E', lvale)
        if (typ(1:1) .eq. 'R') then
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vecpr8(1, kmode), b_incx, zr(lvale), b_incy)
        else if (typ(1:1) .eq. 'C') then
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, vecpc8(1, kmode), b_incx, zc(lvale), b_incy)
        end if
!       SI LE CHAMP A DEJA ETE NOTE PAR SEMOCO, ON NE LE REFAIT PAS
        if (ier .ne. 0) call rsnoch(modes, nosy, nordr)
!
! ----- ON STOCKE 'NUME_MODE'
!
        irang = indk24(nopara, nopast(1), 1, nbpari)
        if (irang .gt. 0) then
            call rsadpa(modes, 'E', 1, nopast(1), nordr, &
                        0, sjv=ladpa, styp=k8b)
            if ((typcon .ne. 'MODE_FLAMB') .or. (nomcmd .eq. 'NORM_MODE')) then
                zi(ladpa) = resufi(kmode, irang)
            else
                zi(ladpa) = nordr
            end if
        end if
!
! ----- ON STOCKE 'NORME'
!
        irang = indk24(nopara(nbpari+1), nopast(2), 1, nbpark)
        if (irang .gt. 0) then
            call rsadpa(modes, 'E', 1, nopast(2), nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk24(ladpa) = resufk(kmode, irang)
        end if
!
! ----- ON STOCKE 'TYPE_MODE' POUR LES MODES PROPRES 'MODE_MECA'
!
        if (typcon(1:9) .eq. 'MODE_MECA' .or. typcon(1:9) .eq. 'MODE_GENE') then
!
            irang = indk24(nopara(nbpari+1), nopast(3), 1, nbpark)
            if (irang .gt. 0) then
                typmod = resufk(kmode, irang)
                if (typmod(1:8) .eq. '        ') then
                    typmod = 'MODE_DYN'
                end if
            else
                typmod = "MODE_DYN"
            end if
!
            call rsadpa(modes, 'E', 1, nopast(3), nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk16(ladpa) = typmod
!
        end if
!
! ----- ON STOCKE 'NOEUD_CMP'
!
        irang = indk24(nopara(nbpari+1), nopast(4), 1, nbpark)
        if (irang .gt. 0) then
            call rsadpa(modes, 'E', 1, nopast(4), nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk16(ladpa) = resufk(kmode, irang)
        end if
!
! ----- ON STOCKE : MODELE, CARA_ELEM, CHAM_MATER
!
        if (lstock .and. (nomcmd .ne. 'NORM_MODE')) then
            call rsadpa(modes, 'E', 1, 'MODELE', nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk8(ladpa) = modele
            call rsadpa(modes, 'E', 1, 'CHAMPMAT', nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk8(ladpa) = chmat
            call rsadpa(modes, 'E', 1, 'CARAELEM', nordr, &
                        0, sjv=ladpa, styp=k8b)
            zk8(ladpa) = carael
        end if
!
!
! ----- ON STOCKE LES PARAMETRES REELS
!
        if (typcon .eq. 'MODE_FLAMB') then
            call rsadpa(modes, 'E', 1, 'CHAR_CRIT', nordr, &
                        0, sjv=ladpa, styp=k8b)
            if (nomcmd .eq. 'NORM_MODE') then
                zr(ladpa) = resufr(kmode, 1)
            else
                zr(ladpa) = -resufr(kmode, 2)
            end if
        else if (typcon .eq. 'MODE_STAB') then
            call rsadpa(modes, 'E', 1, 'CHAR_STAB', 1, &
                        0, sjv=ladpa, styp=k8b)
            zr(ladpa) = resufr(kmode, 1)
        else
            do i = 5, nbpast
                irang = indk24(nopara(nbpari+nbpark+1), nopast(i), 1, nbparr)
                if (irang .gt. 0) then
                    call rsadpa(modes, 'E', 1, nopast(i), nordr, &
                                0, sjv=ladpa, styp=k8b)
                    zr(ladpa) = resufr(kmode, irang)
                end if
            end do
        end if
!
    end do
!
    call jedema()
end subroutine
