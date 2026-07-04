module python_pfp_interface
  use iso_c_binding
  implicit none
  private
  public :: pfp_initialize
  public :: pfp_finalize
  public :: pfp_calculate_energy_and_forces

  interface
     subroutine Py_Initialize() bind(C, name="Py_Initialize")
     end subroutine Py_Initialize

     integer(c_int) function Py_IsInitialized() bind(C, name="Py_IsInitialized")
       import :: c_int
     end function Py_IsInitialized

     integer(c_int) function Py_FinalizeEx() bind(C, name="Py_FinalizeEx")
       import :: c_int
     end function Py_FinalizeEx

     type(c_ptr) function PyImport_ImportModule(name) bind(C, name="PyImport_ImportModule")
       import :: c_char, c_ptr
       character(kind=c_char), dimension(*), intent(in) :: name
     end function PyImport_ImportModule

     type(c_ptr) function PyObject_GetAttrString(obj, name) bind(C, name="PyObject_GetAttrString")
       import :: c_ptr, c_char
       type(c_ptr), value :: obj
       character(kind=c_char), dimension(*), intent(in) :: name
     end function PyObject_GetAttrString

     integer(c_int) function PyCallable_Check(obj) bind(C, name="PyCallable_Check")
       import :: c_ptr, c_int
       type(c_ptr), value :: obj
     end function PyCallable_Check

     type(c_ptr) function PyList_New(size) bind(C, name="PyList_New")
       import :: c_ptr, c_size_t
       integer(c_size_t), value :: size
     end function PyList_New

     integer(c_int) function PyList_SetItem(list, index, item) bind(C, name="PyList_SetItem")
       import :: c_ptr, c_ptrdiff_t, c_int
       type(c_ptr), value :: list
       integer(c_ptrdiff_t), value :: index
       type(c_ptr), value :: item
     end function PyList_SetItem

     type(c_ptr) function PyList_GetItem(list, index) bind(C, name="PyList_GetItem")
       import :: c_ptr, c_ptrdiff_t
       type(c_ptr), value :: list
       integer(c_ptrdiff_t), value :: index
     end function PyList_GetItem

     type(c_ptr) function PyBool_FromLong(value) bind(C, name="PyBool_FromLong")
       import :: c_long, c_ptr
       integer(c_long), value :: value
     end function PyBool_FromLong

     integer(c_ptrdiff_t) function PyList_Size(list) bind(C, name="PyList_Size")
       import :: c_ptr, c_ptrdiff_t
       type(c_ptr), value :: list
     end function PyList_Size

     type(c_ptr) function PyLong_FromLong(value) bind(C, name="PyLong_FromLong")
       import :: c_long, c_ptr
       integer(c_long), value :: value
     end function PyLong_FromLong

     type(c_ptr) function PyFloat_FromDouble(value) bind(C, name="PyFloat_FromDouble")
       import :: c_double, c_ptr
       real(c_double), value :: value
     end function PyFloat_FromDouble

     type(c_ptr) function PyTuple_New(length) bind(C, name="PyTuple_New")
       import :: c_size_t, c_ptr
       integer(c_size_t), value :: length
     end function PyTuple_New

     integer(c_int) function PyTuple_SetItem(tup, index, item) bind(C, name="PyTuple_SetItem")
       import :: c_ptr, c_ptrdiff_t, c_int
       type(c_ptr), value :: tup
       integer(c_ptrdiff_t), value :: index
       type(c_ptr), value :: item
     end function PyTuple_SetItem

     type(c_ptr) function PyTuple_GetItem(tup, index) bind(C, name="PyTuple_GetItem")
       import :: c_ptr, c_ptrdiff_t
       type(c_ptr), value :: tup
       integer(c_ptrdiff_t), value :: index
     end function PyTuple_GetItem

     type(c_ptr) function PyObject_CallObject(callable, args) bind(C, name="PyObject_CallObject")
       import :: c_ptr
       type(c_ptr), value :: callable
       type(c_ptr), value :: args
     end function PyObject_CallObject

     type(c_ptr) function PyUnicode_FromString(string) bind(C, name="PyUnicode_FromString")
       import :: c_char, c_ptr
       character(kind=c_char), dimension(*), intent(in) :: string
     end function PyUnicode_FromString

     real(c_double) function PyFloat_AsDouble(obj) bind(C, name="PyFloat_AsDouble")
       import :: c_ptr, c_double
       type(c_ptr), value :: obj
     end function PyFloat_AsDouble

     integer(c_int) function PyRun_SimpleString(command) bind(C, name="PyRun_SimpleString")
       import :: c_char, c_int
       character(kind=c_char), dimension(*), intent(in) :: command
     end function PyRun_SimpleString

     subroutine PyErr_Print() bind(C, name="PyErr_Print")
     end subroutine PyErr_Print
  end interface

  logical, save :: python_ready = .false.
  logical, save :: owns_interpreter = .false.
  type(c_ptr), save :: cached_module_ptr = c_null_ptr
  type(c_ptr), save :: cached_func_ptr = c_null_ptr

  character(len=*), parameter :: module_name = 'pfp_function'
  character(len=*), parameter :: function_name = 'calculate_energy_and_forces'
  character(len=*), parameter :: tolist_name = 'tolist'

contains

  pure function to_c_chars(str) result(chars)
    character(len=*), intent(in) :: str
    character(kind=c_char), dimension(len_trim(str)+1) :: chars
    integer :: i, n
    n = len_trim(str)
    do i = 1, n
       chars(i) = str(i:i)
    end do
    chars(n+1) = c_null_char
  end function to_c_chars

  subroutine ensure_ptr(ptr, context)
    type(c_ptr), intent(in) :: ptr
    character(len=*), intent(in) :: context
    if (.not. c_associated(ptr)) then
       call PyErr_Print()
       stop 'Python call failed during '//trim(context)
    end if
  end subroutine ensure_ptr

  subroutine ensure_ok(status, context)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: context
    if (status /= 0_c_int) then
       call PyErr_Print()
       stop 'Python error during '//trim(context)
    end if
  end subroutine ensure_ok

  subroutine pfp_initialize(extra_sys_path)
    implicit none
    character(len=*), intent(in), optional :: extra_sys_path
    character(len=:), allocatable :: insert_path
    character(len=:), allocatable :: path_setup

    if (python_ready) return

    if (Py_IsInitialized() == 0_c_int) then
      call Py_Initialize()
      owns_interpreter = .true.
    end if

    if (present(extra_sys_path)) then
      insert_path = trim(extra_sys_path)
    else
      insert_path = '.'
    end if

    path_setup = 'import sys'//char(10)//'sys.path.insert(0, ' // &
      char(34)//trim(insert_path)//char(34)//')'
    call ensure_ok(PyRun_SimpleString(to_c_chars(path_setup)), 'extend sys.path')

    if (.not. c_associated(cached_module_ptr)) then
      cached_module_ptr = PyImport_ImportModule(to_c_chars(module_name))
      call ensure_ptr(cached_module_ptr, 'import module')
    end if

    if (.not. c_associated(cached_func_ptr)) then
      cached_func_ptr = PyObject_GetAttrString(cached_module_ptr, to_c_chars(function_name))
      call ensure_ptr(cached_func_ptr, 'get function')
      if (PyCallable_Check(cached_func_ptr) == 0_c_int) then
        stop 'Python object is not callable'
      end if
    end if

    python_ready = .true.
  end subroutine pfp_initialize

  subroutine pfp_finalize()
    implicit none

    if (.not. python_ready) return

    if (owns_interpreter) then
      if (Py_FinalizeEx() /= 0) then
        stop 'Python finalization failed'
      end if
      owns_interpreter = .false.
    end if

    python_ready = .false.
    cached_module_ptr = c_null_ptr
    cached_func_ptr = c_null_ptr
  end subroutine pfp_finalize

  subroutine pfp_calculate_energy_and_forces( &
       n_atoms, atomic_numbers, positions, cell, energy, &
       forces, extra_sys_path, pbc, calc_mode)
    implicit none
    integer, intent(in) :: n_atoms
    integer(c_long), intent(in) :: atomic_numbers(n_atoms)
    real(c_double), intent(in) :: positions(n_atoms, 3)
    real(c_double), intent(in) :: cell(3, 3)
    real(c_double), intent(out) :: energy
    real(c_double), intent(out) :: forces(n_atoms, 3)
    character(len=*), intent(in), optional :: extra_sys_path
    logical, intent(in) :: pbc
    character(len=*), intent(in) :: calc_mode

    type(c_ptr) :: atomic_list_ptr, positions_list_ptr
    type(c_ptr) :: cell_list_ptr
    type(c_ptr) :: row_ptr, value_ptr
    type(c_ptr) :: args_ptr, result_ptr
    type(c_ptr) :: energy_obj_ptr, forces_obj_ptr
    type(c_ptr) :: tolist_method_ptr, tolist_result_ptr
    type(c_ptr) :: empty_tuple_ptr
    type(c_ptr) :: py_pbc
    type(c_ptr) :: py_calc_mode
    integer :: i, j
    integer(c_ptrdiff_t) :: list_size

    if (.not. python_ready) then
      if (present(extra_sys_path)) then
        call pfp_initialize(extra_sys_path)
      else
        call pfp_initialize()
      end if
    end if

    call ensure_ptr(cached_func_ptr, 'cached Python function')

    atomic_list_ptr = PyList_New(int(n_atoms, kind=c_size_t))
    call ensure_ptr(atomic_list_ptr, 'create atomic number list')

    do i = 1, n_atoms
      value_ptr = PyLong_FromLong(int(atomic_numbers(i), kind=c_long))
      call ensure_ptr(value_ptr, 'convert atomic number')
      call ensure_ok( &
           PyList_SetItem(atomic_list_ptr, int(i-1, kind=c_ptrdiff_t), value_ptr), &
           'set atomic list item')
    end do

    positions_list_ptr = PyList_New(int(n_atoms, kind=c_size_t))
    call ensure_ptr(positions_list_ptr, 'create positions list')

    do i = 1, n_atoms
      row_ptr = PyList_New(int(3, kind=c_size_t))
      call ensure_ptr(row_ptr, 'create position row')
      do j = 1, 3
        value_ptr = PyFloat_FromDouble(positions(i, j))
        call ensure_ptr(value_ptr, 'convert position value')
        call ensure_ok( &
             PyList_SetItem(row_ptr, int(j-1, kind=c_ptrdiff_t), value_ptr), &
             'set position component')
      end do
      call ensure_ok( &
           PyList_SetItem(positions_list_ptr, int(i-1, kind=c_ptrdiff_t), row_ptr), &
           'set positions row')
    end do

    cell_list_ptr = PyList_New(int(3, kind=c_size_t))
    call ensure_ptr(cell_list_ptr, 'create cell list')

    do i = 1, 3
      row_ptr = PyList_New(int(3, kind=c_size_t))
      call ensure_ptr(row_ptr, 'create cell row')
      do j = 1, 3
        value_ptr = PyFloat_FromDouble(cell(i, j))
        call ensure_ptr(value_ptr, 'convert cell value')
        call ensure_ok( &
             PyList_SetItem(row_ptr, int(j-1, kind=c_ptrdiff_t), value_ptr), &
             'set cell component')
      end do
      call ensure_ok( &
           PyList_SetItem(cell_list_ptr, int(i-1, kind=c_ptrdiff_t), row_ptr), &
           'set cell row')
    end do

    py_pbc = PyBool_FromLong(merge(1_c_long, 0_c_long, pbc))
    call ensure_ptr(py_pbc, 'convert pbc value')

    py_calc_mode = PyUnicode_FromString(to_c_chars(calc_mode))
    call ensure_ptr(py_calc_mode, 'convert calc_mode value')

    ! Fixed 5-element positional tuple: (atomic_numbers, positions, cell, pbc, calc_mode)
    args_ptr = PyTuple_New(int(5, kind=c_size_t))
    call ensure_ptr(args_ptr, 'create argument tuple')

    call ensure_ok( &
         PyTuple_SetItem(args_ptr, 0_c_ptrdiff_t, atomic_list_ptr), &
         'set atomic_numbers argument')
    call ensure_ok( &
         PyTuple_SetItem(args_ptr, 1_c_ptrdiff_t, positions_list_ptr), &
         'set positions argument')
    call ensure_ok( &
         PyTuple_SetItem(args_ptr, 2_c_ptrdiff_t, cell_list_ptr), &
         'set cell argument')
    call ensure_ok( &
         PyTuple_SetItem(args_ptr, 3_c_ptrdiff_t, py_pbc), &
         'set pbc argument')
    call ensure_ok( &
         PyTuple_SetItem(args_ptr, 4_c_ptrdiff_t, py_calc_mode), &
         'set calc_mode argument')

    result_ptr = PyObject_CallObject(cached_func_ptr, args_ptr)
    if (.not. c_associated(result_ptr)) then
      call PyErr_Print()
      stop 'Python function call failed'
    end if

    energy_obj_ptr = PyTuple_GetItem(result_ptr, 0_c_ptrdiff_t)
    forces_obj_ptr = PyTuple_GetItem(result_ptr, 1_c_ptrdiff_t)

    energy = PyFloat_AsDouble(energy_obj_ptr)

    tolist_method_ptr = PyObject_GetAttrString(forces_obj_ptr, to_c_chars(tolist_name))
    call ensure_ptr(tolist_method_ptr, 'get tolist method')

    empty_tuple_ptr = PyTuple_New(int(0, kind=c_size_t))
    call ensure_ptr(empty_tuple_ptr, 'create empty tuple')

    tolist_result_ptr = PyObject_CallObject(tolist_method_ptr, empty_tuple_ptr)
    call ensure_ptr(tolist_result_ptr, 'call tolist')

    list_size = PyList_Size(tolist_result_ptr)
    if (list_size /= n_atoms) then
      stop 'Unexpected forces list size'
    end if

    do i = 1, n_atoms
      row_ptr = PyList_GetItem(tolist_result_ptr, int(i-1, kind=c_ptrdiff_t))
      do j = 1, 3
        value_ptr = PyList_GetItem(row_ptr, int(j-1, kind=c_ptrdiff_t))
        forces(i, j) = PyFloat_AsDouble(value_ptr)
      end do
    end do
  end subroutine pfp_calculate_energy_and_forces

end module python_pfp_interface
