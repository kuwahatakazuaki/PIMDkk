module element_table
  use iso_c_binding, only: c_long
  use utility, only: lowerchr, program_abort
  implicit none
  private
  public :: symbol_to_atomic_number

contains

  integer(c_long) function symbol_to_atomic_number(symbol) result(z)
    character(len=*), intent(in) :: symbol
    character(len=len(symbol)) :: lower_symbol

    lower_symbol = lowerchr(trim(symbol))

    select case(trim(lower_symbol))
      case('h', 'd', 'mu')
        z = 1_c_long
      case('he')
        z = 2_c_long
      case('li')
        z = 3_c_long
      case('be')
        z = 4_c_long
      case('b')
        z = 5_c_long
      case('c')
        z = 6_c_long
      case('n')
        z = 7_c_long
      case('o')
        z = 8_c_long
      case('f')
        z = 9_c_long
      case('ne')
        z = 10_c_long
      case('na')
        z = 11_c_long
      case('mg')
        z = 12_c_long
      case('al')
        z = 13_c_long
      case('si')
        z = 14_c_long
      case('p')
        z = 15_c_long
      case('s')
        z = 16_c_long
      case('cl')
        z = 17_c_long
      case('ar')
        z = 18_c_long
      case('k')
        z = 19_c_long
      case('ca')
        z = 20_c_long
      case('sc')
        z = 21_c_long
      case('ti')
        z = 22_c_long
      case('v')
        z = 23_c_long
      case('cr')
        z = 24_c_long
      case('mn')
        z = 25_c_long
      case('fe')
        z = 26_c_long
      case('co')
        z = 27_c_long
      case('ni')
        z = 28_c_long
      case('cu')
        z = 29_c_long
      case('zn')
        z = 30_c_long
      case('ga')
        z = 31_c_long
      case('ge')
        z = 32_c_long
      case('as')
        z = 33_c_long
      case('se')
        z = 34_c_long
      case('br')
        z = 35_c_long
      case('kr')
        z = 36_c_long
      case('rb')
        z = 37_c_long
      case('sr')
        z = 38_c_long
      case('y')
        z = 39_c_long
      case('zr')
        z = 40_c_long
      case('nb')
        z = 41_c_long
      case('mo')
        z = 42_c_long
      case('tc')
        z = 43_c_long
      case('ru')
        z = 44_c_long
      case('rh')
        z = 45_c_long
      case('pd')
        z = 46_c_long
      case('ag')
        z = 47_c_long
      case('cd')
        z = 48_c_long
      case('in')
        z = 49_c_long
      case('sn')
        z = 50_c_long
      case('sb')
        z = 51_c_long
      case('te')
        z = 52_c_long
      case('i')
        z = 53_c_long
      case('xe')
        z = 54_c_long
      case('cs')
        z = 55_c_long
      case('ba')
        z = 56_c_long
      case('la')
        z = 57_c_long
      case('ce')
        z = 58_c_long
      case('pt')
        z = 78_c_long
      case('au')
        z = 79_c_long
      case('hg')
        z = 80_c_long
      case default
        call program_abort('ERROR!!! Unknown element: '//trim(symbol))
    end select
  end function symbol_to_atomic_number

end module element_table
