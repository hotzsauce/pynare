!****************************************************************************
!
!  MODULE: globalvariablesfunctions
!
!****************************************************************************

MODULE globalvariablesfunctions

IMPLICIT NONE

!	MODEL PARAMETERS
DOUBLE PRECISION, PARAMETER		::	alpha				= 0.333D0
DOUBLE PRECISION, PARAMETER		::	beta 				= 0.994D0
DOUBLE PRECISION, PARAMETER		::	delta0 				= 0.025D0
DOUBLE PRECISION, PARAMETER		::	delta1 				= beta**(-1D0) - 1D0 + delta0
DOUBLE PRECISION, PARAMETER		::	delta2 				= 0.00031D0
DOUBLE PRECISION, PARAMETER		::	eta 				= 0.346674415887757D0
DOUBLE PRECISION, PARAMETER		::	phik				= 10D0
DOUBLE PRECISION, PARAMETER		::	phip				= 100.0D0
DOUBLE PRECISION, PARAMETER		::	sigma				= 15D0
DOUBLE PRECISION, PARAMETER		::	ies					= 0.95D0
DOUBLE PRECISION, PARAMETER		::	thetavf				= (1D0-sigma)/(1D0-1D0/ies)
DOUBLE PRECISION, PARAMETER		::	thetamu				= 6.0D0
DOUBLE PRECISION, PARAMETER		::	piess				= 1.00D0
DOUBLE PRECISION, PARAMETER		::	rss					= piess/beta
DOUBLE PRECISION, PARAMETER		::	rhopie				= 1.5D0
DOUBLE PRECISION, PARAMETER		::	rhoy				= 0.2D0
DOUBLE PRECISION, PARAMETER		::	ass					= 1.0D0
DOUBLE PRECISION, PARAMETER		::	yss					= 1.0D0
DOUBLE PRECISION, PARAMETER		::	rhoa				= 0.93564D0
DOUBLE PRECISION, PARAMETER		::	rhovola				= 0.74227D0
DOUBLE PRECISION, PARAMETER		::	volass				= 0.0026251D0
DOUBLE PRECISION, PARAMETER		::	volvola				= 0.0015D0
DOUBLE PRECISION, PARAMETER		::	fixedcost			= (thetamu/(thetamu-1D0)-1D0)*yss
DOUBLE PRECISION, PARAMETER		::	productionconstant	= 1.14997551663338D0
DOUBLE PRECISION, PARAMETER		::	utilityconstant		= 0.00588552318290974D0

!	PROGRAM AND GRID PARAMETERS
INTEGER, PARAMETER				::	nstatevars			= 4
INTEGER, PARAMETER				::	ndecisionvars		= 5
INTEGER, PARAMETER				::	nshocks				= 2

INTEGER, PARAMETER				::	na					= 14
DOUBLE PRECISION, PARAMETER		::	mina				= 0.90D0
DOUBLE PRECISION, PARAMETER		::	maxa				= 1.03D0
DOUBLE PRECISION, PARAMETER		::	stepa				= (maxa - mina) / (na - 1D0)
DOUBLE PRECISION				::	grida(na)			

INTEGER, PARAMETER				::	nk					= 17
DOUBLE PRECISION, PARAMETER		::	mink				= 10.4D0
DOUBLE PRECISION, PARAMETER		::	maxk				= 11.2D0
DOUBLE PRECISION, PARAMETER		::	stepk				= (maxk - mink) / (nk - 1D0)
DOUBLE PRECISION				::	gridk(nk)

INTEGER, PARAMETER				::	nvola				= 10
DOUBLE PRECISION, PARAMETER		::	minvola				= 0.001D0
DOUBLE PRECISION, PARAMETER		::	maxvola				= 0.00775D0
DOUBLE PRECISION, PARAMETER		::	stepvola			= (maxvola - minvola) / (nvola - 1D0)
DOUBLE PRECISION				::	gridvola(nvola)

INTEGER, PARAMETER				::	nlagy				= 13
DOUBLE PRECISION, PARAMETER		::	minlagy				= 0.97D0
DOUBLE PRECISION, PARAMETER		::	maxlagy				= 1.03D0
DOUBLE PRECISION, PARAMETER		::	steplagy			= (maxlagy - minlagy) / (nlagy - 1D0)
DOUBLE PRECISION				::	gridlagy(nlagy)

INTEGER, PARAMETER				::	nepsilona			= 17
DOUBLE PRECISION, PARAMETER		::	rangestdepsilona	= 4.0D0
DOUBLE PRECISION, PARAMETER		::	minepsilona			= -1D0*rangestdepsilona
DOUBLE PRECISION, PARAMETER		::	maxepsilona			= 1D0*rangestdepsilona
DOUBLE PRECISION, PARAMETER		::	stepepsilona		= (maxepsilona - minepsilona) / (nepsilona - 1D0)
DOUBLE PRECISION				::	gridepsilona(nepsilona)	

INTEGER, PARAMETER				::	nepsilonvola		= 17
DOUBLE PRECISION, PARAMETER		::	rangestdepsilonvola	= 4.0D0
DOUBLE PRECISION, PARAMETER		::	minepsilonvola		= -1D0*rangestdepsilonvola
DOUBLE PRECISION, PARAMETER		::	maxepsilonvola		= 1D0*rangestdepsilonvola
DOUBLE PRECISION, PARAMETER		::	stepepsilonvola		= (maxepsilonvola - minepsilonvola) / (nepsilonvola - 1D0)
DOUBLE PRECISION				::	gridepsilonvola(nepsilonvola)	

INTEGER							:: state(nstatevars)
DOUBLE PRECISION				:: rule(na,nk,nvola,nlagy,ndecisionvars)

!$OMP THREADPRIVATE(state)

CONTAINS

!****************************************************************************
!
!  FUNCTION: transformguess
!
!****************************************************************************

FUNCTION transformguess(guess,nequations)

IMPLICIT NONE

INTEGER, INTENT(IN)					::	nequations
DOUBLE PRECISION, INTENT(IN)		::	guess(nequations)
DOUBLE PRECISION					::	transformguess(nequations)

transformguess(1) = guess(1)/0.4D0 - 0.2D0/0.4D0
transformguess(1) = DLOG(transformguess(1)/(1D0-transformguess(1)))

transformguess(2) = guess(2)/0.45D0 - 0.05D0/0.45D0
transformguess(2) = DLOG(transformguess(2)/(1D0-transformguess(2)))

transformguess(3) = guess(3)/0.2D0 - 0.90D0/0.2D0
transformguess(3) = DLOG(transformguess(3)/(1D0-transformguess(3)))

transformguess(4) = guess(4)/1.0D0 - 0.50D0/1.0D0
transformguess(4) = DLOG(transformguess(4)/(1D0-transformguess(4)))

transformguess(5) = (guess(5) - 1.0D-6)**(0.5D0)


END FUNCTION transformguess


!****************************************************************************
!
!  FUNCTION: untransformguess
!
!****************************************************************************

FUNCTION untransformguess(guess,nequations)

IMPLICIT NONE

INTEGER, INTENT(IN)					::	nequations
DOUBLE PRECISION, INTENT(IN)		::	guess(nequations)
DOUBLE PRECISION					::	untransformguess(nequations)

untransformguess(1)	= 0.2D0 + 0.4D0*DEXP(guess(1))/(1D0+DEXP(guess(1)))
untransformguess(2)	= 0.05D0 + 0.45D0*DEXP(guess(2))/(1D0+DEXP(guess(2)))
untransformguess(3)	= 0.9D0 + 0.2D0*DEXP(guess(3))/(1D0+DEXP(guess(3)))
untransformguess(4)	= 0.5D0 + 1.0D0*DEXP(guess(4))/(1D0+DEXP(guess(4)))
untransformguess(5)	= 1.0D-6 + guess(5)**(2D0)

END FUNCTION untransformguess


END MODULE globalvariablesfunctions



!****************************************************************************
!
!  PROGRAM: basubundickeffectivedemand
!
!****************************************************************************


PROGRAM basubundickeffectivedemand

USE globalvariablesfunctions

IMPLICIT NONE
INTEGER				:: ia, ik, ivola, ilagy, iepsilona, iepsilonvola, converge, idecisionvars
DOUBLE PRECISION	:: convergetol, nonlineareqtol, transformedguess(ndecisionvars), fnorm, lastupdatedrule(ndecisionvars)
DOUBLE PRECISION	:: updatedrule(na,nk,nvola,nlagy,ndecisionvars), percentdiffrule(ndecisionvars)

EXTERNAL fcn

DO ia = 1, na
	grida(ia) = mina + (DBLE(ia) - 1D0)*stepa
END DO

DO ik = 1, nk
	gridk(ik) = mink + (DBLE(ik) - 1D0)*stepk
END DO

DO ivola = 1, nvola
	gridvola(ivola) = minvola + (DBLE(ivola) - 1D0)*stepvola
END DO

DO ilagy = 1, nlagy
	gridlagy(ilagy) = minlagy + (DBLE(ilagy) - 1D0)*steplagy
END DO

DO iepsilona = 1, nepsilona
	gridepsilona(iepsilona) = minepsilona + (DBLE(iepsilona) - 1D0)*stepepsilona
END DO

DO iepsilonvola = 1, nepsilonvola
	gridepsilonvola(iepsilonvola) = minepsilonvola + (DBLE(iepsilonvola) - 1D0)*stepepsilonvola
END DO

OPEN(UNIT = 10, FILE = "decisionrulen.dat", STATUS = "OLD", ACTION = "READ")
OPEN(UNIT = 11, FILE = "decisionrulei.dat", STATUS = "OLD", ACTION = "READ")
OPEN(UNIT = 12, FILE = "decisionrulepie.dat", STATUS = "OLD", ACTION = "READ")
OPEN(UNIT = 13, FILE = "decisionruleu.dat", STATUS = "OLD", ACTION = "READ")
OPEN(UNIT = 14, FILE = "decisionruleexpvf1sigma.dat", STATUS = "OLD", ACTION = "READ")

DO ia = 1, na
	DO ik = 1, nk
		DO ivola = 1,nvola	
			DO ilagy = 1,nlagy	
			
			READ(10,'(F18.16)') rule(ia,ik,ivola,ilagy,1)
			READ(11,'(F18.16)') rule(ia,ik,ivola,ilagy,2)
			READ(12,'(F18.16)') rule(ia,ik,ivola,ilagy,3)
			READ(13,'(F18.16)') rule(ia,ik,ivola,ilagy,4)
			READ(14,'(F18.16)') rule(ia,ik,ivola,ilagy,5)
			
			END DO
		END DO
	END DO
END DO


CLOSE(10)
CLOSE(11)
CLOSE(12)
CLOSE(13)
CLOSE(14)

converge = 0
convergetol = 1D-6
nonlineareqtol = 1D-12
updatedrule = rule 



DO WHILE (converge .EQ. 0)	

	converge = 1

	!$OMP PARALLEL PRIVATE(ia,ik,ivola,ilagy,transformedguess,fnorm,lastupdatedrule) 

	CALL ERSET(0, 0, 0) 

	!$OMP DO COLLAPSE(nstatevars) SCHEDULE(DYNAMIC)
	DO ia = 1, na
		DO ik = 1, nk
			DO ivola = 1, nvola
				DO ilagy = 1, nlagy

				state = (/ia,ik,ivola,ilagy/)
				
				
				CALL DNEQNF(fcn, nonlineareqtol, ndecisionvars, 100, &
				transformguess(rule(ia,ik,ivola,ilagy,:),ndecisionvars), transformedguess, fnorm)
				
				IF (fnorm .GT. ndecisionvars*nonlineareqtol**2) THEN
	
					CALL DNEQNF(fcn, nonlineareqtol, ndecisionvars, 100, &
					transformguess(lastupdatedrule,ndecisionvars), transformedguess, fnorm)
			
					IF (fnorm .GT. ndecisionvars*nonlineareqtol**2) THEN
						WRITE(*,*) "Failed To Find Solution To System Of Non-Linear Equations"
						WRITE(*,*) "State: ",state
						WRITE(*,*) "norm(Residuals): ",fnorm
	
					END IF
										
				END IF
			
				updatedrule(ia,ik,ivola,ilagy,:) = untransformguess(transformedguess,ndecisionvars)
				
				lastupdatedrule = updatedrule(ia,ik,ivola,ilagy,:)
				
				END DO			
			END DO
		END DO
	END DO	
	!$OMP END DO 
	!$OMP END PARALLEL

	DO ia = 1, na
		DO ik = 1, nk
			DO ivola = 1, nvola
				DO ilagy = 1, nlagy
				
					IF (converge .EQ. 1) THEN
						DO idecisionvars = 1,ndecisionvars
						percentdiffrule(idecisionvars) = &
						ABS(LOG(updatedrule(ia,ik,ivola,ilagy,idecisionvars)/rule(ia,ik,ivola,ilagy,idecisionvars)))
					IF (percentdiffrule(idecisionvars) .GT. convergetol) THEN
								converge = 0
							END IF
						END DO
					END IF
					
				END DO						
			END DO						
		END DO			
	END DO

	
	
	WRITE(*,*) "--------------------------------------------------------------"
	WRITE(*,*) " "
	WRITE(*,*) "Minimum N: ",minval(minval(minval(minval(updatedrule(:,:,:,:,1),4),3),2),1) 
	WRITE(*,*) "Maximum N: ",maxval(maxval(maxval(maxval(updatedrule(:,:,:,:,1),4),3),2),1)
	WRITE(*,*) "  "
	WRITE(*,*) "Minimum I: ",minval(minval(minval(minval(updatedrule(:,:,:,:,2),4),3),2),1) 
	WRITE(*,*) "Maximum I: ",maxval(maxval(maxval(maxval(updatedrule(:,:,:,:,2),4),3),2),1)
	WRITE(*,*) "  "
	WRITE(*,*) "Minimum PIE: ",minval(minval(minval(minval(updatedrule(:,:,:,:,3),4),3),2),1) 
	WRITE(*,*) "Maximum PIE: ",maxval(maxval(maxval(maxval(updatedrule(:,:,:,:,3),4),3),2),1)
	WRITE(*,*) "  "
	WRITE(*,*) "Minimum U: ",minval(minval(minval(minval(updatedrule(:,:,:,:,4),4),3),2),1) 
	WRITE(*,*) "Maximum U: ",maxval(maxval(maxval(maxval(updatedrule(:,:,:,:,4),4),3),2),1)
	WRITE(*,*) "  "	
	WRITE(*,*) "Minimum EXPVF1SIGMA: ",minval(minval(minval(minval(updatedrule(:,:,:,:,5),4),3),2),1) 
	WRITE(*,*) "Maximum EXPVF1SIGMA: ",maxval(maxval(maxval(maxval(updatedrule(:,:,:,:,5),4),3),2),1)
	WRITE(*,*) "  "
	
	WRITE(*,*) "Maximum N Difference:",maxval(maxval(maxval(maxval(abs(log(updatedrule(:,:,:,:,1)/rule(:,:,:,:,1))),4),3),2),1)
	WRITE(*,*) "Maximum I Difference:",maxval(maxval(maxval(maxval(abs(log(updatedrule(:,:,:,:,2)/rule(:,:,:,:,2))),4),3),2),1)
	WRITE(*,*) "Maximum PIE Difference:",maxval(maxval(maxval(maxval(abs(log(updatedrule(:,:,:,:,3)/rule(:,:,:,:,3))),4),3),2),1)
	WRITE(*,*) "Maximum U Difference:",maxval(maxval(maxval(maxval(abs(log(updatedrule(:,:,:,:,4)/rule(:,:,:,:,4))),4),3),2),1)
	WRITE(*,*) "Maximum EXPVF1SIGMA Difference:",&
				maxval(maxval(maxval(maxval(abs(log(updatedrule(:,:,:,:,5)/rule(:,:,:,:,5))),4),3),2),1)
	WRITE(*,*) "  "
	
	rule = updatedrule
	
END DO

OPEN(UNIT = 10, FILE = "decisionrulen.dat", STATUS = "REPLACE", ACTION = "WRITE")
OPEN(UNIT = 11, FILE = "decisionrulei.dat", STATUS = "REPLACE", ACTION = "WRITE")
OPEN(UNIT = 12, FILE = "decisionrulepie.dat", STATUS = "REPLACE", ACTION = "WRITE")
OPEN(UNIT = 13, FILE = "decisionruleu.dat", STATUS = "REPLACE", ACTION = "WRITE")
OPEN(UNIT = 14, FILE = "decisionruleexpvf1sigma.dat", STATUS = "REPLACE", ACTION = "WRITE")



DO ia = 1, na
	DO ik = 1, nk
		DO ivola = 1, nvola
			DO ilagy = 1, nlagy

				WRITE(10,'(F18.16)') rule(ia,ik,ivola,ilagy,1)
				WRITE(11,'(F18.16)') rule(ia,ik,ivola,ilagy,2)
				WRITE(12,'(F18.16)') rule(ia,ik,ivola,ilagy,3)
				WRITE(13,'(F18.16)') rule(ia,ik,ivola,ilagy,4)
				WRITE(14,'(F18.16)') rule(ia,ik,ivola,ilagy,5)
							
			END DO
		END DO
	END DO
END DO

CLOSE(10)
CLOSE(11)
CLOSE(12)
CLOSE(13)
CLOSE(14)



END PROGRAM basubundickeffectivedemand



!****************************************************************************
!
!  SUBROUTINE: fcn
!
!****************************************************************************

SUBROUTINE fcn(guess,residuals,nequations)

USE globalvariablesfunctions

IMPLICIT NONE

INTEGER, INTENT(IN)					::	nequations
DOUBLE PRECISION, INTENT(IN)		::	guess(nequations)
DOUBLE PRECISION, INTENT(OUT)		::	residuals(nequations)

DOUBLE PRECISION					::	untransformedguess(nequations), interpolateddecisionvariables(nequations)


DOUBLE PRECISION					::	n, i, pie, u, expvf1sigma, a, k, vola, lagy 
DOUBLE PRECISION					::	y, deltau, deltauprime, k1, q, c, w, mu, rrk, rd, r, vf, epsilonvola1, epsilona1
DOUBLE PRECISION					::	n1, i1, pie1, u1, expvf1sigma1, a1, vola1
DOUBLE PRECISION					::	y1, deltau1, deltauprime1, k2, q1, c1, w1, mu1, rrk1, rd1, r1, vf1, sdf1

DOUBLE PRECISION					::	normdensepsilona1, normdensepsilonvola1
DOUBLE PRECISION					::	sdf1pie1epsilonvolaepsilona(nepsilona), sdf1pie1epsilonvola(nepsilonvola)
DOUBLE PRECISION					::	sdf1rrk1q1epsilonvolaepsilona(nepsilona), sdf1rrk1q1epsilonvola(nepsilonvola)
DOUBLE PRECISION					::	sdf1pie1y1epsilonvolaepsilona(nepsilona), sdf1pie1y1epsilonvola(nepsilonvola) 
DOUBLE PRECISION					::	vf1epsilonvolaepsilona(nepsilona), vf1epsilonvola(nepsilonvola)

DOUBLE PRECISION					::	expsdf1pie1, expsdf1rrk1q1, expsdf1pie1y1, expvf1

INTEGER								::	iepsilona1, iepsilonvola1
DOUBLE PRECISION, PARAMETER			::	pi = 3.14159265358979D0

EXTERNAL linearinterp4

untransformedguess = untransformguess(guess,nequations)

n 			= untransformedguess(1)
i 			= untransformedguess(2)
pie			= untransformedguess(3)
u			= untransformedguess(4)
expvf1sigma	= untransformedguess(5)

a			= grida(state(1))
k			= gridk(state(2))
vola		= gridvola(state(3))
lagy		= gridlagy(state(4))

y			= productionconstant*(u*k)**(alpha)*(n)**(1D0-alpha) - fixedcost
deltau		= delta0 + delta1*(u-1D0) + (delta2/2D0)*(u-1D0)**(2D0);
deltauprime	= delta1 + delta2*(u-1D0);
k1			= ( (1D0-deltau) - (phik/2D0) * (i/k-delta0)**(2D0) ) *k + i
q			= (1D0 - phik * (i/k - delta0))**(-1D0)
c			= y - i - (phip/2D0) * (pie/piess - 1D0)**(2D0) * y
w			= ( (1D0-eta) / eta ) * c / (1D0-n)
mu			= (1D0-alpha) * (y + fixedcost) / (w * n) 
rrk			= alpha * (y + fixedcost) / (u * k * mu)
rd			= DEXP(DLOG(rss) + rhopie*DLOG(pie/piess) + rhoy*DLOG(y/lagy) )
vf			= (utilityconstant*a*(c**(eta)*(1D0-n)**(1D0-eta))**((1D0-sigma)/thetavf) & 
				+ beta * expvf1sigma**(1D0/thetavf))**(thetavf/(1D0-sigma))

IF (rd .LT. 1.0D0) THEN

	r = 1.00D0

ELSE

	r = rd

END IF

DO iepsilonvola1 = 1, nepsilonvola

		epsilonvola1 = gridepsilonvola(iepsilonvola1)
		vola1 = rhovola*vola + (1D0-rhovola)*volass + volvola*epsilonvola1
		normdensepsilonvola1 = DEXP(-0.5D0 * ((epsilonvola1)**2D0))	/ (DSQRT(2D0 * pi))

		IF (vola1 .LT. 1D-4) THEN
			vola1 = 1D-4
		END IF

	DO iepsilona1 = 1, nepsilona

		epsilona1 = gridepsilona(iepsilona1)
		a1 = rhoa*a + (1D0-rhoa) + vola*epsilona1
		
		normdensepsilona1 = DEXP(-0.5D0 * ((epsilona1)**2D0))	/ (DSQRT(2D0 * pi))

		CALL linearinterp4(interpolateddecisionvariables,nequations,na,grida,stepa,a1,nk,gridk,stepk,k1,&
							nvola,gridvola,stepvola,vola1,nlagy,gridlagy,steplagy,rd,rule(:,:,:,:,:))
							
		n1 				= interpolateddecisionvariables(1)
		i1 				= interpolateddecisionvariables(2)
		pie1 			= interpolateddecisionvariables(3)
		u1				= interpolateddecisionvariables(4)
		expvf1sigma1	= interpolateddecisionvariables(5)
			
		y1				= productionconstant * (u1*k1)**(alpha) * (n1)**(1D0-alpha) - fixedcost
		deltau1			= delta0 + delta1*(u1-1D0) + (delta2/2D0)*(u1-1D0)**(2D0);
		deltauprime1	= delta1 + delta2*(u1-1D0);
		k2			= ( (1D0-deltau1) - (phik/2D0) * (i1/k1-delta0)**(2D0) ) *k1 + i1
		q1			= (1D0 - phik * (i1/k1 - delta0))**(-1D0)
		c1			= y1 - i1 - (phip/2D0) * (pie1/piess - 1D0)**(2D0) * y1
		w1			= ( (1D0-eta) / eta ) * c1 / (1D0-n1)
		mu1			= (1D0-alpha) * (y1 + fixedcost) / (w1 * n1) 
		rrk1		= alpha * (y1 + fixedcost) / (u1 * k1 * mu1)
		rd1			= DEXP(DLOG(rss) + rhopie*DLOG(pie1/piess) + rhoy*DLOG(y1/y) )
		vf1			= (utilityconstant * a1 * (c1**(eta)*(1D0-n1)**(1D0-eta))**((1D0-sigma)/thetavf) & 
						+ beta * expvf1sigma1**(1D0/thetavf))**(thetavf/(1D0-sigma))
							
		IF (rd1 .LT. 1.00D0) THEN
		
			r1 = 1.00D0
		
		ELSE
		
			r1 = rd1
		
		END IF

	
		sdf1 = beta*(a1/a)*((c1**(eta)*(1D0-n1)**(1D0-eta))/(c**(eta)*(1D0-n)**(1D0-eta)))**((1D0-sigma)/thetavf)*(c/c1)&
			*(vf1**(1D0-sigma)/expvf1sigma)**(1D0-1D0/thetavf)
		
		sdf1pie1epsilonvolaepsilona(iepsilona1) = normdensepsilonvola1 * normdensepsilona1 * sdf1 / pie1
		
		sdf1rrk1q1epsilonvolaepsilona(iepsilona1) = normdensepsilonvola1 * normdensepsilona1 * sdf1 * (u1*rrk1 + q1*((1D0-deltau1) &
										- (phik/2D0)*(i1/k1-delta0)**(2D0) + phik*(i1/k1-delta0)*(i1/k1) ) )

		sdf1pie1y1epsilonvolaepsilona(iepsilona1) = normdensepsilonvola1 * normdensepsilona1 * sdf1 * &
												phip * (pie1/piess - 1D0) * (y1/y) * (pie1/piess)
		
		vf1epsilonvolaepsilona(iepsilona1) = normdensepsilonvola1 * normdensepsilona1 * vf1**(1D0-sigma)
		
	
	END DO

	sdf1pie1epsilonvola(iepsilonvola1) = (stepepsilona / 2D0) * & 
		(2D0 * SUM(sdf1pie1epsilonvolaepsilona) - sdf1pie1epsilonvolaepsilona(1) - sdf1pie1epsilonvolaepsilona(nepsilona))
	
	sdf1rrk1q1epsilonvola(iepsilonvola1) = (stepepsilona / 2D0) * &
		(2D0 * SUM(sdf1rrk1q1epsilonvolaepsilona) - sdf1rrk1q1epsilonvolaepsilona(1) - sdf1rrk1q1epsilonvolaepsilona(nepsilona))
	
	sdf1pie1y1epsilonvola(iepsilonvola1) = (stepepsilona / 2D0) * &
		(2D0 * SUM(sdf1pie1y1epsilonvolaepsilona) - sdf1pie1y1epsilonvolaepsilona(1) - sdf1pie1y1epsilonvolaepsilona(nepsilona))
		
	vf1epsilonvola(iepsilonvola1) = (stepepsilona / 2D0) * &
		(2D0 * SUM(vf1epsilonvolaepsilona) - vf1epsilonvolaepsilona(1) - vf1epsilonvolaepsilona(nepsilona))

		
END DO




expsdf1pie1 = (stepepsilonvola / 2D0) * &
		(2D0 * SUM(sdf1pie1epsilonvola) - sdf1pie1epsilonvola(1) - sdf1pie1epsilonvola(nepsilonvola))

expsdf1rrk1q1 = (stepepsilonvola / 2D0) * &
		(2D0 * SUM(sdf1rrk1q1epsilonvola) - sdf1rrk1q1epsilonvola(1) - sdf1rrk1q1epsilonvola(nepsilonvola))

expsdf1pie1y1 = (stepepsilonvola / 2D0) * &
		(2D0 * SUM(sdf1pie1y1epsilonvola) - sdf1pie1y1epsilonvola(1) - sdf1pie1y1epsilonvola(nepsilonvola))

expvf1 = (stepepsilonvola / 2D0) * &
		(2D0 * SUM(vf1epsilonvola) - vf1epsilonvola(1) - vf1epsilonvola(nepsilonvola))


residuals(1)	= 1D0 - r * expsdf1pie1
residuals(2)	= q - expsdf1rrk1q1
residuals(3)	= phip * (pie/piess - 1D0) * (pie/piess) - (1D0 - thetamu) - thetamu / mu - expsdf1pie1y1 
residuals(4)	= alpha * (y + fixedcost) * mu**(-1D0) - q * deltauprime * u * k
residuals(5)	= expvf1sigma - expvf1


!WRITE(*,*) n, i, pie, u, expvf1sigma
!WRITE(*,*) a, k, vola, lagy 
!WRITE(*,*) state 
!WRITE(*,*) residuals 

!PAUSE 

RETURN
END SUBROUTINE fcn


!****************************************************************************
!
!  SUBROUTINE: linearinterp4
!
!****************************************************************************

SUBROUTINE linearinterp4(zi,ndecisionrules,nx1,x1grid,x1step,x1i,nx2,x2grid,x2step,x2i,nx3,x3grid,x3step,x3i,nx4,x4grid,x4step,x4i,z)

IMPLICIT NONE

INTEGER, INTENT(IN)					::	ndecisionrules, nx1, nx2, nx3, nx4
DOUBLE PRECISION, INTENT(IN)		::	x1i, x1grid(nx1), x1step
DOUBLE PRECISION, INTENT(IN)		::	x2i, x2grid(nx2), x2step
DOUBLE PRECISION, INTENT(IN)		::	x3i, x3grid(nx3), x3step
DOUBLE PRECISION, INTENT(IN)		::	x4i, x4grid(nx4), x4step

DOUBLE PRECISION, INTENT(IN)		::	z(nx1,nx2,nx3,nx4,ndecisionrules)
DOUBLE PRECISION, INTENT(OUT)		::	zi(ndecisionrules)

INTEGER								::	location1, location2, location3, location4
INTEGER								::	mz, m1, m2, m3, m4
DOUBLE PRECISION					::	w1(2), w2(2), w3(2), w4(2)

location1 = MIN(nx1-1,MAX(1,FLOOR((x1i-x1grid(1))/x1step)+1))
w1(2) = (x1i-x1grid(location1)) / x1step
w1(1) = 1D0-w1(2)    

location2 = MIN(nx2-1,MAX(1,FLOOR((x2i-x2grid(1))/x2step)+1))
w2(2) = (x2i-x2grid(location2)) / x2step
w2(1) = 1D0-w2(2)    

location3 = MIN(nx3-1,MAX(1,FLOOR((x3i-x3grid(1))/x3step)+1))
w3(2) = (x3i-x3grid(location3)) / x3step
w3(1) = 1D0-w3(2)    

location4 = MIN(nx4-1,MAX(1,FLOOR((x4i-x4grid(1))/x4step)+1))
w4(2) = (x4i-x4grid(location4)) / x4step
w4(1) = 1D0-w4(2)    

DO mz = 1, ndecisionrules

	zi(mz) = 0D0

		DO m1 = 0, 1
			DO m2 = 0, 1
				DO m3 = 0, 1
					DO m4 = 0, 1		
					zi(mz) = zi(mz) + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*z(location1+m1,location2+m2,location3+m3,location4+m4,mz)
					END DO
				END DO
			END DO
		END DO

END DO

RETURN
END SUBROUTINE linearinterp4