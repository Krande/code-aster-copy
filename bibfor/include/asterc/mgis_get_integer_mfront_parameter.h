interface
subroutine mgis_get_integer_mfront_parameter(extern_addr, nopara, value)
    character(len=16), intent(in) :: extern_addr
    character(len=*), intent(in) :: nopara
    integer(kind=8), intent(out) :: value
end subroutine mgis_get_integer_mfront_parameter
end interface
