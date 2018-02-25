module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs
