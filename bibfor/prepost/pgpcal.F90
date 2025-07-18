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

subroutine pgpcal(sd_pgp)

    use DynaGene_module

    implicit none
! Calculate the physical fields of interest based on the user requests
! and preprocessed data in the pgp data structure
!
! Saves the results into a standard Code_Aster table structure
!
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pgpget.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!   ====================================================================
!   = 0 =   Variable declarations and initialization
!   ====================================================================
!   -0.1- Input/output arguments
    character(len=8), intent(in):: sd_pgp
!   -0.2- Local variables
    real(kind=8)      :: physvalr
    complex(kind=8)   :: physvalc
    integer(kind=8)           :: nbobs, iobs, physlen, length, nbmodes
    integer(kind=8)           :: i, j, iord, dec1, nord
    integer(kind=8)           :: jevol, jvecr, jvecc, jtblp, lc, shift, i_chreco
    character(len=4)  :: chreco, typcha, typsc, typres
    character(len=8)  :: resin, result
    character(len=12) :: bl11pt
    character(len=16) :: champ
    character(len=24) :: nomjv

    integer(kind=8), pointer :: lordr(:) => null()
    real(kind=8), pointer :: vectr(:) => null()
    complex(kind=8), pointer :: vectc(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: v_resu(:) => null()
    type(DynaGene) :: dyna_gene

!   ------------------------------------------------------------------------------------
!   Definition of statement functions giving the appropriate (i,j) term in the basis
!   vector
#define vr(m,n) vectr((n-1)*physlen+m)
#define evolr(n,p) v_resu((p-1)*nbmodes+n)
#define vc(m,n) vectc((n-1)*physlen+m)
#define evolc(n,p) zc(jevol+(p-1)*nbmodes+n-1)

!   -0.3- Initialization

    bl11pt = '           .'

    call jemarq()

    call pgpget(sd_pgp, 'RESU_OUT', kscal=result)
    call pgpget(sd_pgp, 'NB_OBSER', iscal=nbobs)

    call jeveuo(result//'           .TBLP', 'L', jtblp)
    nomjv = zk24(jtblp+4*(9-1)+2)
    call jeveuo(nomjv, 'E', jvecr)
    nomjv = zk24(jtblp+4*(10-1)+2)
    call jeveuo(nomjv, 'E', jvecc)

    call pgpget(sd_pgp, 'RESU_IN ', kscal=resin)
    call pgpget(sd_pgp, 'TYP_RESU ', kscal=typres)

    if (typres(1:4) .eq. 'TRAN') then
        call dyna_gene%init(resin)
    end if

    call jeveuo(resin//'           .DESC', 'L', vi=desc)
    nbmodes = desc(2)

!   Line counter, across the whole table (for different observations)
    lc = 0
    do iobs = 1, nbobs
        call pgpget(sd_pgp, 'NOM_CHAM ', iobs=iobs, kscal=champ)
        call pgpget(sd_pgp, 'TYP_CHAM ', iobs=iobs, kscal=typcha)
!
        chreco = 'DEPL'
        i_chreco = dyna_gene%depl
        if (champ(1:4) .eq. 'VITE') then
            chreco = 'VITE'
            i_chreco = dyna_gene%vite
        else if (champ(1:4) .eq. 'ACCE') then
            chreco = 'ACCE'
            i_chreco = dyna_gene%acce
        end if
!
        call pgpget(sd_pgp, 'NUM_ORDR', iobs=iobs, lonvec=nord)
        AS_ALLOCATE(vi=lordr, size=nord)
        call pgpget(sd_pgp, 'NUM_ORDR', iobs=iobs, ivect=lordr)
!

        call pgpget(sd_pgp, 'REF_COMP', iobs=iobs, lonvec=physlen)
        call pgpget(sd_pgp, 'TYP_SCAL', iobs=iobs, kscal=typsc)

        if (typsc(1:1) .eq. 'R') then

            call pgpget(sd_pgp, 'VEC_PR_R ', iobs=iobs, lonvec=length)
            AS_ALLOCATE(vr=vectr, size=length)
            call pgpget(sd_pgp, 'VEC_PR_R ', iobs=iobs, rvect=vectr)

            if (typres(1:4) .eq. 'TRAN') then
                do iord = 1, nord
                    call dyna_gene%get_values_by_index(i_chreco, lordr(iord), shift, vr=v_resu)

                    dec1 = lc+(iord-1)*physlen

                    do i = 1, physlen
                        physvalr = 0.d0
                        do j = 1, nbmodes
                            physvalr = physvalr+vr(i, j)*evolr(j, lordr(iord)-shift)
                        end do
                        zr(jvecr+dec1+i-1) = physvalr
                    end do
                end do

            else if (typres(1:4) .eq. 'HARM') then
                call jeveuo(resin//bl11pt//chreco, 'L', jevol)
                do iord = 1, nord
                    dec1 = lc+(iord-1)*physlen
                    do i = 1, physlen
                        physvalc = dcmplx(0.d0, 0.d0)
                        do j = 1, nbmodes
                            physvalc = physvalc+dcmplx(vr(i, j), 0.d0)*evolc(j, lordr(iord))
                        end do
                        zc(jvecc+dec1+i-1) = physvalc
                    end do
                end do

            end if
            AS_DEALLOCATE(vr=vectr)

        else if (typsc(1:1) .eq. 'C') then

            call pgpget(sd_pgp, 'VEC_PR_C ', iobs=iobs, lonvec=length)
            AS_ALLOCATE(vc=vectc, size=length)
            call pgpget(sd_pgp, 'VEC_PR_C ', iobs=iobs, cvect=vectc)

            if (typres(1:4) .eq. 'TRAN') then
                do iord = 1, nord
                    call dyna_gene%get_values_by_index(i_chreco, lordr(iord), shift, vr=v_resu)
                    dec1 = lc+(iord-1)*physlen
                    do i = 1, physlen
                        physvalc = dcmplx(0.d0, 0.d0)
                        do j = 1, nbmodes
                            physvalc = physvalc+vc(i, j)*dcmplx(evolr(j, lordr(iord)-shift), 0.d0)
                        end do
                        zc(jvecc+dec1+i-1) = physvalc
                    end do
                end do

            else if (typres(1:4) .eq. 'HARM') then
                call jeveuo(resin//bl11pt//chreco, 'L', jevol)
                do iord = 1, nord
                    dec1 = lc+(iord-1)*physlen
                    do i = 1, physlen
                        physvalc = dcmplx(0.d0, 0.d0)
                        do j = 1, nbmodes
                            physvalc = physvalc+vc(i, j)*evolc(j, lordr(iord))
                        end do
                        zc(jvecc+dec1+i-1) = physvalc
                    end do
                end do
            end if
            AS_DEALLOCATE(vc=vectc)
        end if
!
        lc = lc+nord*physlen
        AS_DEALLOCATE(vi=lordr)
!
    end do

    if (typres(1:4) .eq. 'TRAN') then
        call dyna_gene%free
    else
        call jelibe(resin//bl11pt//chreco)
    end if

    call jedema()

end subroutine
