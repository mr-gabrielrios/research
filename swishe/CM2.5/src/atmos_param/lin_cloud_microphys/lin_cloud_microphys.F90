module lin_cld_microphys_mod

 use time_manager_mod,  only: time_type
 use fms_mod,           only: error_mesg, FATAL

 implicit none
 private

 public  lin_cld_microphys_driver, lin_cld_microphys_init, lin_cld_microphys_end, sg_conv
 public  qsmith_init, qsmith, es2_table1d, es3_table1d, esw_table1d, g_sum, wqsat_moist, wqsat2_moist, sat_adj2

 contains

!-----------------------------------------------------------------------------------------------------------------------
  subroutine lin_cld_microphys_driver(qv,    ql,    qr,    qi,    qs,    qg,    qa,  &
                               qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
                               pt_dt, pt, p3, dz,  delp, area, dt_in,                &
                               land,  rain, snow, ice, graupel,                      &
                               hydrostatic, phys_hydrostatic,                        &
                               iis,iie, jjs,jje, kks,kke, ktop, kbot, time)

  type(time_type), intent(in):: time
  logical,         intent(in):: hydrostatic, phys_hydrostatic
  integer,         intent(in):: iis,iie, jjs,jje
  integer,         intent(in):: kks,kke
  integer,         intent(in):: ktop, kbot
  real,            intent(in):: dt_in

  real, intent(in   ), dimension(:,:)  :: area
  real, intent(in   ), dimension(:,:)  :: land
  real, intent(out  ), dimension(:,:)  :: rain, snow, ice, graupel
  real, intent(in   ), dimension(:,:,:):: p3, delp, dz
  real, intent(in   ), dimension(:,:,:):: pt, qv, ql, qr, qi, qs, qg, qa
  real, intent(inout), dimension(:,:,:):: pt_dt,  qa_dt
  real, intent(inout), dimension(:,:,:):: qv_dt, ql_dt, qr_dt, qi_dt,  &
                                          qs_dt, qg_dt

  call error_mesg('lin_cld_microphys_driver','The null version of lin_cld_microphys_driver should not be called.',FATAL)

 end subroutine lin_cld_microphys_driver
!-----------------------------------------------------------------------------------------------------------------------

 subroutine sat_adj2(mdt, is, ie, js, je, ng, km, k, hydrostatic, consv_te, &
                     te0, qv, ql, qi, qr, qs, qa, area, peln, delz, pt, dp, last_step)
 real, intent(in):: mdt
 integer, intent(in):: is, ie, js, je, km, ng, k
 logical, intent(in):: hydrostatic, last_step
 logical, intent(in):: consv_te
 real, intent(in), dimension(is-ng:ie+ng,js-ng:je+ng):: dp, area
 real, intent(in):: delz(is:ie,js:je)
 real, intent(in):: peln(is:ie,km+1,js:je)
 real, intent(inout), dimension(is-ng:ie+ng,js-ng:je+ng):: pt, qv, ql, qi, qr, qs, qa
 real, intent(inout):: te0(is:ie,js:je)

 call error_mesg('sat_adj2','The null version of sat_adj2 should not be called.',FATAL)

 end subroutine sat_adj2
!-----------------------------------------------------------------------------------------------------------------------

 subroutine lin_cld_microphys_init(id, jd, kd, axes, time)

    integer,         intent(in) :: id, jd, kd
    integer,         intent(in) :: axes(4)
    type(time_type), intent(in) :: time

 call error_mesg('lin_cld_microphys_init','The null version of lin_cld_microphys_init should not be called.',FATAL)

 end subroutine lin_cld_microphys_init
!-----------------------------------------------------------------------------------------------------------------------

 subroutine lin_cld_microphys_end

 call error_mesg('lin_cld_microphys_end','The null version of lin_cld_microphys_end should not be called.',FATAL)

 end subroutine lin_cld_microphys_end
!-----------------------------------------------------------------------------------------------------------------------

 subroutine qsmith_init

 call error_mesg('qsmith_init','The null version of qsmith_init should not be called.',FATAL)

 end subroutine qsmith_init
!-----------------------------------------------------------------------------------------------------------------------

 real function wqsat2_moist(ta, qv, pa, dqdt)
 real, intent(in):: ta, pa, qv
 real, intent(out):: dqdt

 call error_mesg('wqsat2_moist','The null version of wqsat2_moist should not be called.',FATAL)

 end function wqsat2_moist
!-----------------------------------------------------------------------------------------------------------------------

 real function wqsat_moist(ta, qv, pa)
 real, intent(in):: ta, pa, qv

 call error_mesg('wqsat_moist','The null version of wqsat_moist should not be called.',FATAL)

 end function wqsat_moist
!-----------------------------------------------------------------------------------------------------------------------

 subroutine esw_table1d(ta, es, n)
 integer, intent(in):: n
 real, intent(in)::  ta(n)
 real, intent(out):: es(n)

 call error_mesg('esw_table1d','The null version of esw_table1d should not be called.',FATAL)

 end subroutine esw_table1d
!-----------------------------------------------------------------------------------------------------------------------

 subroutine es2_table1d(ta, es, n)
 integer, intent(in):: n
 real, intent(in)::  ta(n)
 real, intent(out):: es(n)

 call error_mesg('es2_table1d','The null version of es2_table1d should not be called.',FATAL)

 end subroutine es2_table1d
!-----------------------------------------------------------------------------------------------------------------------

 subroutine es3_table1d(ta, es, n)
 integer, intent(in):: n
 real, intent(in)::  ta(n)
 real, intent(out):: es(n)

 call error_mesg('es3_table1d','The null version of es3_table1d should not be called.',FATAL)

 end subroutine es3_table1d
!-----------------------------------------------------------------------------------------------------------------------

 subroutine qsmith(im, km, ks, t, p, q, qs, dqdt)

 integer, intent(in):: im, km, ks
 real, intent(in),dimension(im,km):: t, p, q
 real, intent(out),dimension(im,km):: qs
 real, intent(out), optional:: dqdt(im,km)

 call error_mesg('qsmith','The null version of qsmith should not be called.',FATAL)

 end subroutine qsmith
!-----------------------------------------------------------------------------------------------------------------------

 subroutine sg_conv(is, ie, js, je, isd, ied, jsd, jed,               &
                    km, nq, dt, tau,             &
                    delp, phalf, pm, zfull, zhalf, ta, qa, ua, va, w, &
                    u_dt, v_dt, t_dt, q_dt, nqv, nql, nqi, nqr, nqs, nqg, &
                    hydrostatic, phys_hydrostatic)

      logical, intent(in):: hydrostatic, phys_hydrostatic
      integer, intent(in):: is, ie, js, je, km, nq
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: tau
      integer, intent(in):: nqv, nql, nqi
      integer, intent(in):: nqr, nqs, nqg
      real, intent(in):: dt
      real, intent(in):: phalf(is:ie,js:je,km+1)
      real, intent(in):: pm(is:ie,js:je,km)
      real, intent(in):: zfull(is:ie,js:je,km)
      real, intent(in):: zhalf(is:ie,js:je,km+1)
      real, intent(in):: delp(isd:ied,jsd:jed,km)

      real, intent(inout), dimension(is:ie,js:je,km)::  ta, ua, va
      real, intent(inout):: qa(is:ie,js:je,km,nq)
      real, intent(inout):: w(isd:ied,jsd:jed,km)
      real, intent(inout):: u_dt(isd:ied,jsd:jed,km)
      real, intent(inout):: v_dt(isd:ied,jsd:jed,km)
      real, intent(inout):: t_dt(is:ie,js:je,km)
      real, intent(inout):: q_dt(is:ie,js:je,km,nq)

  call error_mesg('sg_conv','The null version of sg_conv should not be called.',FATAL)

 end subroutine sg_conv
!-----------------------------------------------------------------------------------------------------------------------

 real function g_sum(p, ifirst, ilast, jfirst, jlast, area, mode)
 use mpp_mod,           only: mpp_sum
 integer, intent(IN) :: ifirst, ilast
 integer, intent(IN) :: jfirst, jlast
 integer, intent(IN) :: mode
 real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)
 real, intent(IN) :: area(ifirst:ilast,jfirst:jlast)
 real gsum

 call error_mesg('g_sum','The null version of g_sum should not be called.',FATAL)

 end function g_sum
!-----------------------------------------------------------------------------------------------------------------------

end module lin_cld_microphys_mod
