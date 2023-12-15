SUBROUTINE fast_kqr(delta, Kmat, nobs, y, nlam, ulam, tau, &
  & eps, maxit, gamma, is_exact, anlam, npass, jerr, alpmat)
! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, nlam, anlam, jerr, maxit, npass(nlam), is_exact
  DOUBLE PRECISION :: delta, eps, Kmat (nobs, nobs), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gamma, alpmat (nobs+1, nlam), tau
  ! - - - local declarations - - -
  INTEGER :: i, j, l, nn, info, delta_id, delta_save, nlmax, lwork
  INTEGER, PARAMETER :: delta_len = 4, mproj = 5
  DOUBLE PRECISION :: zvec (nobs), rds, vareps, told, tnew, mul = 1.0D0
  DOUBLE PRECISION :: r (nobs), dif (nobs+1), alpvec (nobs+1) , uo
  DOUBLE PRECISION :: KKTeps, cvec(nobs+1), dvec(nobs+1)
  DOUBLE PRECISION :: oldalpvec(nobs+1), KKT(nobs+1)
  DOUBLE PRECISION :: aka, Ka (nobs), obj0, obj1, alp0, mdd
  DOUBLE PRECISION :: Kinv (nobs, nobs), eU(nobs, nobs), einv (nobs)
  DOUBLE PRECISION :: theta (nobs), rmg, alp_old (nobs+1), alptmp (nobs+1)
  DOUBLE PRECISION :: gamvec(nobs), gval(delta_len), hval
  DOUBLE PRECISION :: Umat (nobs, nobs), eigens (nobs), Usum (nobs), rnobs
  DOUBLE PRECISION :: lpinv (nobs, delta_len), lpUsum (nobs, delta_len)
  DOUBLE PRECISION :: svec (nobs, delta_len), vvec (nobs,delta_len)
  INTEGER :: elbowid (nobs), elbchk, is_exit
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  !usestdlib_linalg, only: diag
  ! - - - begin - - -
  npass = 0
  nlmax = 8 * nobs
  rnobs = Real(nobs)
  r = 0.0D0
  alpmat = 0.0D0
  alpvec = 0.0D0
  alptmp = 0.0D0
  alp_old = 0.0D0
  theta = 0.0D0
  elbowid = 0
  rmg = 0.0D0
  zvec = 0.0D0
  svec = 0.0D0
  vvec = 0.0D0
  rds = 0.0D0
  vareps = 1.0D-8
  KKTeps = 1.0D-4
  Kinv = 0.0D0
  mdd = 0.0D0
  aka = 0.0D0
  Ka = 0.0D0
  obj0 = 0.0D0
  obj1 = 0.0D0
  alp0 = 0.0D0
  jerr = 0
  Umat = Kmat
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens+gamma
  Usum = Sum(Umat, dim=1)
  einv = 1/eigens
  DO j = 1, nobs
    eU(j,:) = einv(j) * Umat(:,j)
  ENDDO
  DO i = 1, nobs
    Kinv(:,i) = Matmul(Umat, eU(:,i))
  ENDDO
  r = y
  lambda_loop: DO l = 1,nlam
    delta = 0.125D0
    delta_id = 0
    delta_save = 0
    oldalpvec = 0.0D0
    update_delta: do
      ! CALL INTPR("delta_id", -1, delta_id, 1)
      delta_id = delta_id + 1
      IF (delta_id > delta_save) THEN
        lpinv(:,delta_id) = 1.0D0 / (eigens + rnobs * delta * ulam(l))
        lpUsum(:,delta_id) = lpinv(:,delta_id) * Usum
        vvec(:,delta_id) = Matmul(Umat, eigens * lpUsum(:,delta_id))
        svec(:,delta_id) = Matmul(Umat, lpUsum(:,delta_id))
        gval(delta_id) = 1.0D0 / (rnobs + 2.0D0 * rnobs * delta * vareps &
          & - Sum(vvec(:,delta_id))) !
        delta_save = delta_id
      ENDIF

      dif = 0.0D0
      told = 1.0D0
      ! - - - update alpha - - -
      update_alpha: DO
        DO j = 1, nobs
          IF (r(j) <= - delta) THEN
            zvec(j) = -(tau-1.0D0)
          ELSE IF (r(j) >= delta) THEN
            zvec(j) = -tau
          ELSE
            zvec(j) = -r(j)/(2.0D0 * delta) - tau + 0.50D0
          END IF
        ENDDO
        gamvec = zvec + nobs * ulam(l) * alpvec(2:(nobs + 1))
        rds = sum(zvec) + 2.0D0 * nobs * vareps * alpvec(1)
        hval = rds - Dot_product(vvec(:,delta_id), gamvec) !
        tnew = 0.5 + 0.5 * Sqrt(1.0 + 4.0 * told * told)
        mul = 1.0 + (told - 1.0) / tnew
        told = tnew
        dif(1) = - mul * delta * gval(delta_id) * hval !
        dif(2:(nobs+1)) = - dif(1) * svec(:,delta_id) - mul * delta * &
          & Matmul(Umat, Matmul(gamvec, Umat) * lpinv(:,delta_id))
        alpvec = alpvec + dif
        r = y - (alpvec(1) + Matmul(Kmat, alpvec(2:(nobs + 1))))
        npass(l) = npass(l) + 1
        IF (Maxval(dif * dif) < eps * mul * mul) EXIT
        IF (Sum(npass) > maxit) EXIT
      ENDDO update_alpha

      IF (Sum(npass) > maxit) EXIT

      dif = oldalpvec - alpvec
      ka = Matmul(Kmat, alpvec(2:(nobs+1)))
      aKa = Dot_product(ka, alpvec(2:(nobs+1)))
      CALL objfun(alpvec(1), aka, ka, y, ulam(l), nobs, tau, 1.0D-9, obj0) !object function
      CALL opt_int(-1.0D2, 1.0D2, nobs, ka, aka, y, ulam(l), tau, 1.0D-9, obj1, alp0)

      IF (obj1 < obj0) THEN
        dif(1) = dif(1) + alp0 - alpvec(1)
        r = r - (alp0 - alpvec(1))
        alpvec(1) = alp0
      ENDIF
      oldalpvec = alpvec
      !calculate kkt
      DO j = 1, nobs
        IF (r(j) <= - 1.0D-9) THEN
          zvec(j) = -(tau-1.0D0)
        ELSE IF (r(j) >= 1.0D-9) THEN
          zvec(j) = -tau
        ELSE
          zvec(j) =  -r(j)/(2.0D0 * 1.0D-9) - tau + 0.50D0
        END IF
      ENDDO
      cvec(1) = sum(zvec)
      cvec(2:(nobs+1)) = Matmul(Kmat,zvec)
      dvec(1) = 2*vareps*alpvec(1)
      dvec(2:(nobs+1)) = ulam(l) * Matmul(Kmat,alpvec(2:(nobs+1)))
      KKT = cvec/nobs + dvec
      uo = Max(ulam(l), 1.0D0)
      IF (Sum(KKT * KKT) / uo / uo < KKTeps) THEN
        IF (Maxval(dif * dif) < nobs * eps * mul * mul) THEN
          IF (is_exact .EQ. 0) THEN
            EXIT
          ELSE
            is_exit = 0
            ! projection step
            alptmp = alpvec
            projection: DO nn = 1, mproj
              elbowid = 0
              elbchk = 1
              DO j = 1, nobs
                rmg = r(j)
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
                !CALL INTPR("eigen_elbowid", -1, elbowid, 1)
                ka = Matmul(Kmat, alptmp(2:(nobs+1)))
                aKa = Dot_product(ka, alptmp(2:(nobs+1)))
                CALL objfun(alptmp(1), aka, ka, y, ulam(l), nobs, tau, 1.0D-9, obj0)
                CALL opt_int(-1.0D2, 1.0D2, nobs, ka, aka, y, ulam(l), tau, 1.0D-9, obj1, alp0)
                IF (obj1 < obj0) THEN
                  dif(1) = dif(1) + alp0 - alptmp(1)
                  alptmp(1) = alp0
                ENDIF
                r = y - (alptmp(1) + ka)
                DO j = 1, nobs
                  IF (r(j) <= - delta) THEN
                    zvec(j) = -(tau-1.0D0)
                  ELSE IF (r(j) >= delta) THEN
                    zvec(j) = -tau
                  ELSE
                    zvec(j) = -r(j)/(2.0D0 * delta) - tau + 0.50D0
                  END IF
                ENDDO

                gamvec = zvec + nobs * ulam(l) * alptmp(2:(nobs + 1))
                rds = sum(zvec) + 2.0D0 * nobs * vareps * alptmp(1)
                hval = rds - Dot_product(vvec(:,delta_id), gamvec) !
                tnew = 0.5 + 0.5 * Sqrt(1.0 + 4.0 * told * told)
                mul = 1.0 + (told - 1.0) / tnew
                told = tnew
                dif(1) = - mul * delta * gval(delta_id) * hval
                dif(2:(nobs+1)) = - dif(1) * svec(:,delta_id) - mul * delta * &
                  & Matmul(Umat, Matmul(gamvec, Umat) * lpinv(:,delta_id))
                alptmp = alptmp + dif
                ka = Matmul(Kmat, alptmp(2:(nobs+1)))
                r = y - (alptmp(1) + ka)
                alp_old = alptmp
                IF (Sum(elbowid) > 1) THEN
                  ! CALL INTPR("sum_elbowid", -1, elbowid, 1)
                  theta = Matmul(Kmat, alptmp(2:(nobs+1)))
                  DO j = 1, nobs
                    IF (elbowid(j) .EQ. 1) THEN
                      theta(j) = theta(j) +  r(j)!
                    ENDIF
                  ENDDO
                  alptmp(2:nobs+1) = Matmul(Kinv, theta)
                ENDIF
                dif = dif + alptmp - alp_old
                r = y - (alptmp(1) + Matmul(Kmat, alptmp(2:(nobs+1))))
                npass(l) = npass(l) + 1
                mdd = Maxval(dif * dif)
                IF (mdd < eps * mul * mul) THEN
                  EXIT
                ELSEIF (mdd > nobs .AND. npass(l) > 2) THEN
                  is_exit = 1
                  EXIT
                ENDIF
                IF (Sum(npass) > maxit) THEN
                  is_exit = 1
                  EXIT
                ENDIF
              ENDDO update_alpha_proj
            ENDDO projection
            IF (is_exit .EQ. 1) THEN
              EXIT
            ENDIF
            DO j = 1, nobs
              IF (r(j) <= - 1.0D-9) THEN
                zvec(j) = -(tau-1.0D0)
              ELSE IF (r(j) >= 1.0D-9) THEN
                zvec(j) = -tau
              ELSE
                zvec(j) =  -r(j)/(2.0D0 * 1.0D-9) - tau + 0.50D0
              END IF
            ENDDO
            cvec(1) = sum(zvec)
            cvec(2:(nobs+1)) = Matmul(Kmat,zvec)
            dvec(1) = 2*vareps*alptmp(1)
            dvec(2:(nobs+1)) = ulam(l) * Matmul(Kmat,alptmp(2:(nobs+1)))
            KKT = cvec/nobs + dvec
            uo = Max(ulam(l), 1.0D0)
            ! If KKT fails, reduce delta
            IF (Sum(KKT * KKT) / uo / uo < KKTeps) THEN
              alpvec = alptmp
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF (delta_id .EQ. delta_len) EXIT
      delta = 0.125D0 * delta
    ENDDO update_delta
    ! CALL INTPR("delta_id", -1, delta_id, 1)
    alpmat(:, l) = alpvec
    anlam = l
    IF (Sum(npass) > maxit) THEN
      jerr = -l
      EXIT
    ENDIF
  ENDDO lambda_loop
END SUBROUTINE fast_kqr
