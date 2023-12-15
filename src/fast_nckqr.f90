SUBROUTINE fast_nckqr(delta, Kmat, nobs, y, nlam1, ulam1, nlam2, &
  & ulam2, ntau, utau, eps, maxit, gamma, is_exact, anlam1, anlam2, npass, &
  & jerr, alpmat)
! --------------------------------------------------
IMPLICIT NONE
! - - - arg types - - -
INTEGER :: nobs, nlam1, anlam1, jerr(2), maxit, npass(nlam1, nlam2), is_exact
INTEGER :: nlam2, anlam2, ntau
DOUBLE PRECISION :: ulam1 (nlam1), ulam2 (nlam2)
DOUBLE PRECISION :: delta, eps, Kmat (nobs, nobs), y (nobs)
DOUBLE PRECISION :: gamma, alpmat ((nobs+1), ntau, nlam1, nlam2), utau(ntau)
! - - - local declarations - - -
INTEGER :: j, l1, l2, nlmax, lwork, tt, ierr, uu, info, npass_proj
INTEGER :: delta_id, delta_save, nn, is_exit, elbowid (nobs), elbchk
INTEGER, PARAMETER :: delta_len = 7, mproj = 5
DOUBLE PRECISION :: lam = 1.0D-3, del = 3.051758e-05, bigdel
DOUBLE PRECISION :: zvec(nobs), tau1, vareps, KKTeps, told, tnew, mul
DOUBLE PRECISION :: r1(nobs, ntau), dif((nobs + 1), ntau), yvec(nobs)
DOUBLE PRECISION :: alpvec ((nobs+1), ntau), dvec ((nobs + 1), ntau)
DOUBLE PRECISION :: d1vec ((nobs + 1), ntau), r2(nobs, (ntau - 1))
DOUBLE PRECISION :: d2vec ((nobs + 1), ntau), a1, cont(delta_len)
DOUBLE PRECISION :: d2dif (nobs), d2dif_sum, Kd (nobs), Kinv (nobs, nobs)
DOUBLE PRECISION :: Umat (nobs, nobs), eigens (nobs), Usum (nobs), rnobs
DOUBLE PRECISION :: lpinv (nobs, delta_len), lpUsum (nobs, delta_len)
DOUBLE PRECISION :: eigens1 (nobs), v1vec (nobs, delta_len)
DOUBLE PRECISION :: svec (nobs, delta_len), vvec (nobs, delta_len)
DOUBLE PRECISION :: gval(delta_len), hval, lpinv1 (nobs, delta_len)
DOUBLE PRECISION :: kktvec1 ((nobs + 1), ntau), kktvec2 ((nobs + 1), ntau)
DOUBLE PRECISION :: KKT ((nobs + 1), ntau), rmg, alptmp((nobs+1), ntau)
DOUBLE PRECISION :: alp_old((nobs+1), ntau), mdd, theta(nobs), ka(nobs)
DOUBLE PRECISION, ALLOCATABLE :: works (:)
! - - - begin - - -
anlam1 = 0
anlam2 = 0
npass = 0
jerr = 0
rnobs = Real(nobs)
nlmax = 8 * nobs
zvec = 0.0D0
yvec = 0.0D0
svec = 0.0D0
vvec = 0.0D0
gval = 0.0D0
v1vec = 0.0D0
vareps = 1.0D-8
dvec = 0.0D0
d1vec = 0.0D0
d2vec = 0.0D0
kktvec1 = 0.0D0
kktvec2 = 0.0D0
KKT = 0.0D0
alpmat = 0.0D0
alpvec = 0.0D0
r1 = 0.0D0
r2 = 0.0D0
dif = 0.0D0
a1 = 0.0D0
Kd = 0.0D0
rmg = 0.0D0
KKTeps = 1.0D-3
Kinv = 0.0D0
mdd = 0.0D0
alptmp = 0.0D0
alp_old = 0.0D0
theta = 0.0D0
elbowid = 0
rmg = 0.0D0
mul = 1.0D0

Umat = Kmat
ALLOCATE(works(nlmax))
lwork = -1
CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
lwork = MIN(nlmax, INT(works(1)))
CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
DEALLOCATE(works)
eigens = eigens+gamma
Usum = Sum(Umat, dim=1)

r_loop: DO tt=1, ntau
  r1(1:nobs, tt) = y
ENDDO r_loop

lambda2_loop: DO l2 = 1, nlam2
  lambda1_loop: DO l1 = 1, nlam1
    delta = 0.125D0
    delta_id = 0
    delta_save = 0
    a1 = 1 + 4*ulam1(l1)*rnobs
    eigens1 = a1*eigens
    update_delta: do
      delta_id = delta_id + 1
      IF (delta_id > delta_save) THEN
        cont(delta_id) = 2*rnobs*delta*vareps + rnobs + &
                                 & 4*ulam1(l1)*rnobs**2 + lam*ulam1(l1)*rnobs
        lpinv(:,delta_id) = 1/(a1*eigens + rnobs*delta*ulam2(l2) &
                                & + ulam1(l1)*lam*rnobs/eigens)
        lpinv1(:,delta_id) = a1 * lpinv(:,delta_id)
        lpUsum(:,delta_id) = lpinv1(:,delta_id) * Usum
        v1vec(:,delta_id) = Matmul(Umat, eigens1 * lpUsum(:,delta_id))
        vvec(:,delta_id) = Matmul(Umat, eigens * lpUsum(:,delta_id))
        svec(:,delta_id) = Matmul(Umat, lpUsum(:,delta_id))
        gval(delta_id) = 1.0D0 / (cont(delta_id) - Sum(v1vec(:,delta_id)))
        delta_save = delta_id
      ENDIF
      dif = 0
      told = 1.0D0
      update_alpha: DO
        z_loop: DO tt = 1, ntau
          tau1 = utau(tt)-1.0D0
          DO j = 1, nobs
            IF (r1(j,tt) <= - delta) THEN
              zvec(j) = - tau1
            ELSE IF (r1(j,tt) >= delta) THEN
              zvec(j) = -utau(tt)
            ELSE
              zvec(j) = -r1(j,tt)/(2.0D0 * delta) - utau(tt) + 0.50D0
            END IF
          ENDDO
          d1vec(1,tt) = sum(zvec) + 2.0D0 * rnobs * vareps * alpvec(1,tt)
          d1vec(2:(nobs+1),tt) = zvec + rnobs * ulam2(l2) * alpvec(2:nobs+1, tt)
        ENDDO z_loop

        d2vec=0.0D0
        bigdel=max(delta, del)
        y_loop: DO uu=1, (ntau-1)
          DO j = 1, nobs
            IF (r2(j,uu) <= - bigdel) THEN
              yvec(j) =  0
            ELSE IF (r2(j,uu) >= bigdel) THEN
              yvec(j) = 1
            ELSE
              yvec(j) = 0.5 * r2(j,uu)/bigdel + 0.5
            END IF
          ENDDO
          d2vec(1,uu) = d2vec(1,uu) + ulam1(l1) * sum(yvec)*rnobs
          d2vec(2:(nobs+1), uu) = d2vec(2:(nobs+1),uu) + ulam1(l1)*yvec*rnobs
          d2vec(1, (uu+1)) = d2vec(1, (uu+1)) - ulam1(l1) * sum(yvec)*rnobs
          d2vec(2:(nobs+1), (uu+1)) = d2vec(2:(nobs+1), (uu+1)) - ulam1(l1)*yvec*rnobs
        ENDDO y_loop
        dvec = d1vec + d2vec

        tau_loop: DO tt = 1, ntau
          hval = dvec(1,tt) - Dot_product(vvec(:,delta_id), dvec(2:(nobs+1),tt))
          dif(1,tt) = - delta * gval(delta_id) * hval
          dif(2:(nobs+1),tt) = - dif(1,tt) * svec(:,delta_id) - delta * &
            & Matmul(Umat, Matmul(dvec(2:(nobs+1),tt), Umat) * lpinv(:,delta_id))
          dif(:,tt) = mul * dif(:,tt)
          alpvec(:,tt) = alpvec(:,tt) + dif(:,tt)
          r1(:,tt) = y - (alpvec(1, tt) + Matmul(Kmat, alpvec(2:(nobs + 1), tt)))
        ENDDO tau_loop

        DO uu=1, (ntau-1)
          r2(:, uu) = Matmul(Kmat, alpvec(2:(nobs + 1), uu)-&
              & alpvec(2:(nobs + 1), (uu+1))) + &
              & alpvec(1, uu)-alpvec(1, (uu+1))
        ENDDO

        npass(l1, l2) = npass(l1, l2) + 1
        IF (Maxval(dif * dif) < eps) EXIT
        IF (Sum(npass) > maxit) THEN
          EXIT
        ENDIF
      ENDDO update_alpha

      z_loop2: DO tt = 1, ntau
          tau1 = utau(tt)-1.0D0
          DO j = 1, nobs
            IF (r1(j,tt) <= - 1.0D-9) THEN
              zvec(j) = -tau1
            ELSE IF (r1(j,tt) >= 1.0D-9) THEN
              zvec(j) = -utau(tt)
            ELSE
              zvec(j) =  -r1(j,tt)/(2.0D0 * 1.0D-9) - utau(tt) + 0.50D0
            END IF
          ENDDO
          kktvec1(1,tt) = sum(zvec)/rnobs + 2*vareps*alpvec(1,tt)
          kktvec1(2:(nobs+1),tt) = Matmul(Kmat, zvec/rnobs + ulam2(l2) * &
                                      & alpvec(2:(nobs+1),tt))
        ENDDO z_loop2
        kktvec2 = 0.0D0
        y_loop2: DO uu = 1, (ntau - 1)
          DO j = 1, nobs
            IF (r2(j,uu) <= - del) THEN
              yvec(j) =  0
            ELSE IF (r2(j,uu) >= del) THEN
              yvec(j) = 1
            ELSE
              yvec(j) = 0.5 * r2(j,uu)/del + 0.5
            END IF
          ENDDO
          d2dif = ulam1(l1) * yvec * rnobs
          d2dif_sum = Sum(d2dif)
          Kd = Matmul(Kmat, d2dif)
          kktvec2(1,uu) = kktvec2(1,uu) + d2dif_sum
          kktvec2(2:(nobs+1), uu) = kktvec2(2:(nobs+1),uu) + Kd
          kktvec2(1, (uu+1)) = kktvec2(1, (uu+1)) - d2dif_sum
          kktvec2(2:(nobs+1), (uu+1)) = kktvec2(2:(nobs+1), (uu+1)) - Kd
        ENDDO y_loop2
        KKT = kktvec1 + kktvec2
        IF (maxval(KKT * KKT) < KKTeps) THEN
          IF (Maxval(dif * dif) < nobs * eps * mul * mul) THEN
            IF (is_exact .EQ. 0) THEN
              EXIT
            ELSE
              is_exit = 0
              npass_proj = 0
              ! projection step
              tau_loop2: DO tt=1, ntau
                alptmp(:, tt) = alpvec(:, tt)
                projection: DO nn = 1, mproj
                  elbowid = 0
                  elbchk = 1
                  DO j = 1, nobs
                    rmg = r1(j, tt)
                    IF (rmg < delta .AND. rmg > -delta ) THEN
                      IF (rmg > 1.0D-3) THEN
                        elbchk = 0
                      ENDIF
                      elbowid(j) = 1
                    ENDIF
                  ENDDO
                  IF (elbchk .EQ. 1) THEN
                    EXIT
                  ENDIF
                  told = 1.0D0
                  update_alpha_proj: DO
                    ka = Matmul(Kmat, alptmp(2:(nobs+1), tt))
                    r1(:, tt) = y - (alptmp(1, tt) + ka)
                    DO j = 1, nobs
                      IF (r1(j, tt) <= - delta) THEN
                        zvec(j) = -(utau(tt)-1.0D0)
                      ELSE IF (r1(j, tt) >= delta) THEN
                        zvec(j) = -utau(tt)
                      ELSE
                        zvec(j) = -r1(j, tt)/(2.0D0 * delta) - utau(tt) + 0.50D0
                      END IF
                    ENDDO
                    d1vec(1,tt) = sum(zvec) + 2.0D0 * rnobs * vareps * alptmp(1,tt)
                    d1vec(2:(nobs+1),tt) = zvec + rnobs * ulam2(l2) * alptmp(2:(nobs+1), tt)

                    d2vec(:, tt) = 0.0D0
                    DO uu = 1, (ntau-1)
                      r2(:, uu) = Matmul(Kmat, alptmp(2:(nobs + 1), uu)-&
                          & alptmp(2:(nobs + 1), (uu+1))) + &
                          & alptmp(1, uu) - alptmp(1, (uu+1))
                    ENDDO
                    IF (tt==1) THEN
                      DO j = 1, nobs
                        IF (r2(j,tt) <= - bigdel) THEN
                          yvec(j) =  0
                        ELSE IF (r2(j,tt) >= bigdel) THEN
                          yvec(j) = 1
                        ELSE
                          yvec(j) = 0.5 * r2(j,tt)/bigdel + 0.5
                        END IF
                      ENDDO
                      d2dif = ulam1(l1) * yvec * rnobs
                      d2dif_sum = Sum(d2dif)
                      d2vec(1, tt) = d2dif_sum
                      d2vec(2:(nobs+1), tt) = d2dif
                    ELSE IF (tt==ntau) THEN
                      DO j = 1, nobs
                        IF (r2(j,tt-1) <= - bigdel) THEN
                          yvec(j) =  0
                        ELSE IF (r2(j,tt-1) >= bigdel) THEN
                          yvec(j) = 1
                        ELSE
                          yvec(j) = 0.5 * r2(j,tt-1)/bigdel + 0.5
                        END IF
                      ENDDO
                      d2dif = ulam1(l1) * yvec * rnobs
                      d2dif_sum = Sum(d2dif)
                      d2vec(1, tt) =  - d2dif_sum
                      d2vec(2:(nobs+1), tt) = - d2dif
                    ELSE
                      DO j = 1, nobs
                        IF (r2(j,tt-1) <= - bigdel) THEN
                          yvec(j) =  0
                        ELSE IF (r2(j,tt-1) >= bigdel) THEN
                          yvec(j) = 1
                        ELSE
                          yvec(j) = 0.5 * r2(j,tt-1)/bigdel + 0.5
                        END IF
                      ENDDO
                      d2dif = ulam1(l1) * yvec * rnobs
                      d2dif_sum = Sum(d2dif)
                      d2vec(1,tt) = - d2dif_sum
                      d2vec(2:(nobs+1), tt) = - d2dif
                      DO j = 1, nobs
                        IF (r2(j,tt) <= - bigdel) THEN
                          yvec(j) =  0
                        ELSE IF (r2(j,tt) >= bigdel) THEN
                          yvec(j) = 1
                        ELSE
                          yvec(j) = 0.5 * r2(j,tt)/bigdel + 0.5
                        END IF
                      ENDDO
                      d2dif = ulam1(l1) * yvec * rnobs
                      d2dif_sum = Sum(d2dif)
                      d2vec(1,tt) = d2vec(1,tt) + d2dif_sum
                      d2vec(2:(nobs+1), tt) = d2vec(2:(nobs+1), tt) + d2dif
                    ENDIF
                    dvec(:, tt) = d1vec(:, tt) + d2vec(:, tt)
                    tnew = 0.5 + 0.5 * Sqrt(1.0 + 4.0 * told * told)!
                    mul = 1.0 + (told - 1.0) / tnew
                    told = tnew

                    hval = dvec(1,tt) - Dot_product(vvec(:,delta_id), &
                      & dvec(2:(nobs+1), tt))
                    dif(1,tt) = - mul * delta * gval(delta_id) * hval / rnobs / ulam1(l1)
                    dif(2:(nobs+1),tt) = - mul *dif(1,tt) * svec(:,delta_id) - mul * &
                      & delta * Matmul(Umat, Matmul(dvec(2:(nobs+1),tt), Umat) &
                      & * lpinv(:,delta_id))
                    alptmp(:,tt) = alptmp(:,tt) + dif(:,tt)
                    ka = Matmul(Kmat, alptmp(2:(nobs+1), tt))
                    r1(:, tt) = y - (alptmp(1, tt) + ka)
                    alp_old(:,tt) = alptmp(:,tt)
                    IF (Sum(elbowid) > 1) THEN
                      ! CALL INTPR("sum_elbowid", -1, elbowid, 1)
                      theta = ka
                      DO j = 1, nobs
                        IF (elbowid(j) .EQ. 1) THEN
                          theta(j) = theta(j) +  r1(j, tt)!
                        ENDIF
                      ENDDO
                      alptmp(2:nobs+1, tt) = Matmul(Kinv, theta)
                    ENDIF
                    dif(:, tt) = dif(:, tt) + alptmp(:,tt) - alp_old(:,tt)
                    ! r1(:, tt) = y - (alptmp(1, tt) + Matmul(Kmat, alptmp(2:(nobs+1), tt)))
                    npass(l1, l2) = npass(l1, l2) + 1
                    npass_proj = npass_proj + 1
                    mdd = Maxval(dif * dif)

                    IF (mdd < 1e-02 .OR. npass_proj > 1e+02) THEN
                      EXIT
                    ELSEIF (mdd > nobs .AND. npass(l1,l2) > 2) THEN
                      is_exit = 1
                      EXIT
                    ENDIF
                    IF (Sum(npass) > maxit) THEN
                      is_exit = 1
                      EXIT
                    ENDIF
                  ENDDO update_alpha_proj
                ENDDO projection
              ENDDO tau_loop2

              IF (is_exit .EQ. 1) THEN
                EXIT
              ENDIF
              !KKT
              DO tt=1, ntau
                r1(:, tt) = y - (alptmp(1, tt) + Matmul(Kmat, alptmp(2:(nobs+1), tt)))
              ENDDO
              DO uu = 1, (ntau-1)
                r2(:, uu) = Matmul(Kmat, alptmp(2:(nobs + 1), uu)-&
                    & alptmp(2:(nobs + 1), (uu+1))) + &
                    & alptmp(1, uu) - alptmp(1, (uu+1))
              ENDDO

              DO tt = 1, ntau
                tau1 = utau(tt)-1.0D0
                DO j = 1, nobs
                  IF (r1(j,tt) <= - 1.0D-9) THEN
                    zvec(j) = -tau1
                  ELSE IF (r1(j,tt) >= 1.0D-9) THEN
                    zvec(j) = -utau(tt)
                  ELSE
                    zvec(j) =  -r1(j,tt)/(2.0D0 * 1.0D-9) - utau(tt) + 0.50D0
                  END IF
                ENDDO
                kktvec1(1,tt) = sum(zvec)/rnobs + 2*vareps*alptmp(1,tt)
                kktvec1(2:(nobs+1),tt) = Matmul(Kmat, zvec/rnobs + ulam2(l2) * &
                                            & alptmp(2:(nobs+1),tt))
              ENDDO
              kktvec2 = 0.0D0
              DO uu = 1, (ntau - 1)
                DO j = 1, nobs
                  IF (r2(j,uu) <= - del) THEN
                    yvec(j) =  0
                  ELSE IF (r2(j,uu) >= del) THEN
                    yvec(j) = 1
                  ELSE
                    yvec(j) = 0.5 * r2(j,uu)/del + 0.5
                  END IF
                ENDDO
                d2dif = ulam1(l1) * yvec * rnobs
                d2dif_sum = Sum(d2dif)
                Kd = Matmul(Kmat, d2dif)
                kktvec2(1,uu) = kktvec2(1,uu) + d2dif_sum
                kktvec2(2:(nobs+1), uu) = kktvec2(2:(nobs+1),uu) + Kd
                kktvec2(1, (uu+1)) = kktvec2(1, (uu+1)) - d2dif_sum
                kktvec2(2:(nobs+1), (uu+1)) = kktvec2(2:(nobs+1), (uu+1)) - Kd
              ENDDO
              KKT = kktvec1 + kktvec2

              IF (maxval(KKT * KKT) < KKTeps) THEN
                alpvec = alptmp
                EXIT
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF (delta_id .EQ. delta_len) EXIT
        delta = 0.125D0 * delta
      ENDDO update_delta
      alpmat(:, :, l1, l2) = alpvec
      anlam1 = l1
      anlam2 = l2
      IF (Sum(npass(:, :)) > maxit) THEN
        jerr(1) = -l1
        jerr(2) = -l2
        EXIT
      ENDIF
    ENDDO lambda1_loop
  ENDDO lambda2_loop
END SUBROUTINE fast_nckqr
