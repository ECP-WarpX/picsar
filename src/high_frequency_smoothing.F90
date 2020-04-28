
SUBROUTINE high_frequency_smoothing
  USE fields , ONLY: jm_h, jp_h, jl_h, rho_h, rhoold_h
  USE gpstd_solver , ONLY : tau_filter
  IMPLICIT NONE

  ! Apply filter
  rhoold_h(:,:,:) = rhoold_h(:,:,:) * tau_filter(:,:,:)
  rho_h(:,:,:) = rho_h(:,:,:) * tau_filter(:,:,:)
  jp_h (:,:,:) = jp_h (:,:,:) * tau_filter(:,:,:)
  jm_h (:,:,:) = jm_h (:,:,:) * tau_filter(:,:,:)
  jl_h (:,:,:) = jl_h (:,:,:) * tau_filter(:,:,:)

END SUBROUTINE high_frequency_smoothing

