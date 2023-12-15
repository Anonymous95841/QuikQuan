SUBROUTINE objfun (intcpt, aka, ka, y, lam, nobs, tau, del, objval) !change
  IMPLICIT NONE
  INTEGER :: nobs, j
  DOUBLE PRECISION :: intcpt, aka, ka (nobs), y (nobs), lam, objval
  DOUBLE PRECISION :: fh (nobs), xi (nobs), xi_tmp, tau, del, ttau
  xi = 0.0D0
  !check loss
  DO j = 1, nobs
    fh(j) = ka(j) + intcpt
    xi_tmp = y(j) - fh(j)
    ttau = tau-1.0D0
    IF (xi_tmp <= -del) THEN
      xi_tmp = xi_tmp * ttau
    ELSE IF (xi_tmp >= del) THEN
      xi_tmp = xi_tmp * tau
    ELSE
      xi_tmp = xi_tmp * xi_tmp/(4.0D0*del) + (tau-0.5D0) * xi_tmp + del/4.0D0
    END IF
    xi(j) = xi_tmp
  ENDDO
  objval = (lam/2.0D0) * aka + Sum(xi)/(nobs) + 1.0D-8*intcpt*intcpt
END SUBROUTINE objfun


! This code is partially modified from fmim.c from in the source code of R
SUBROUTINE opt_int (lmin, lmax, nobs, ka, aka, y, lam, tau, del, objval, lhat)
  IMPLICIT NONE
  INTEGER :: nobs
  DOUBLE PRECISION :: lmin, lmax, lhat, objval, ka (nobs), aka, y (nobs), lam, tau
  DOUBLE PRECISION :: a, b, d, e, p, q, r, u, v, w, x, del
  DOUBLE PRECISION :: t2, fu, fv, fw, fx, xm, tol, tol1, tol3
  REAL(KIND = SELECTED_REAL_KIND(10, 99)) :: eps
  DOUBLE PRECISION, PARAMETER :: gold = (3.0D0 - Sqrt(5.0D0)) * 0.5D0

  eps = 3.14 !any double precision value
  eps = Epsilon(eps)
  tol = eps ** 0.25D0
  tol1 = eps + 1.0D0
  eps = Sqrt(eps)

  a = lmin
  b = lmax
  v = a + gold * (b - a)
  w = v
  x = v

  d = 0.0D0
  e = 0.0D0
  objval = 0.0D0
  CALL objfun(x, aka, ka, y, lam, nobs, tau, del, objval)
  fx = objval
  fv = fx
  fw = fx
  tol3 = tol / 3.0D0

  main_loop: DO
    xm = (a + b) * 0.5D0
    tol1 = eps * Abs(x) + tol3
    t2 = tol1 * 2.0D0

    IF (Abs(x - xm) .LE. t2 - (b - a) * 0.5D0) EXIT
    p = 0.0D0
    q = 0.0D0
    r = 0.0D0
    IF (Abs(e) > tol1) THEN
      r = (x - w) * (fx - fv)
      q = (x - v) * (fx - fw)
      p = (x - v) * q - (x - w) * r
      q = (q - r) * 2.0D0
      IF (q > 0.0D0) THEN
        p = -p
      ELSE
        q = -q
      ENDIF
      r = e
      e = d
    ENDIF

    IF ((Abs(p) .GE. Abs(q * 0.5D0 * r)) .OR. (p .LE. q * (a - x)) .OR. &
      &  (p .GE. q * (b - x))) THEN
      IF (x < xm) THEN
        e = b - x
      ELSE
        e = a - x
      ENDIF
      d = gold * e
    ELSE
      d = p / q
      u = x + d
      IF (u - a < t2 .OR. b - u < t2) THEN
        d = tol1
        IF (x .GE. xm) d = -d
      ENDIF
    ENDIF

    IF (Abs(d) .GE. tol1) THEN
      u = x + d
    ELSEIF (d > 0.0D0) THEN
      u = x + tol1
    ELSE
      u = x - tol1
    ENDIF

    objval = 0.0D0
    CALL objfun(u, aka, ka, y, lam, nobs, tau, del, objval)
    fu = objval

    IF (fu .LE. fx) THEN
      IF (u < x) THEN
        b = x
      ELSE
        a = x
      ENDIF
      v = w
      w = x
      x = u
      fv = fw
      fw = fx
      fx = fu
    ELSE
      IF (u < x) THEN
        a = u
      ELSE
        b = u
      ENDIF

      IF ((fu .LE. fw) .OR. (w .EQ. x)) THEN
        v = w
        fv = fw
        w = u
        fw = fu
      ELSEIF ((fu .LE. fv) .OR. (v .EQ. x) .OR. (v .EQ. w)) THEN
        v = u
        fv = fu
      ENDIF
    ENDIF

  ENDDO main_loop
  lhat = x
  objval = 0.0D0
  CALL objfun(x, aka, ka, y, lam, nobs, tau, del, objval)

END SUBROUTINE opt_int


! --------------------------------------------------
SUBROUTINE rbfdot(X1, X2, nobs1, nobs2, p1, p2, sigma, Kmat)
! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs1, nobs2, p1, p2
  DOUBLE PRECISION :: sigma, X1(nobs1, p1), X2(nobs2, p2), Kmat(nobs1, nobs2)
  ! - - - local declarations - - -
  INTEGER :: i, j

  Do i = 1, nobs1
    Do j = 1, nobs2
      If (i > j) THEN
        Kmat(i,j) = Kmat(j,i)
      ELSE
        Kmat(i,j) = exp(sigma * (-sum((X1(i,:) - X2(j,:))**2)))
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE rbfdot
