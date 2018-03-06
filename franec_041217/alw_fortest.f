      PROGRAM FORTEST
c      REAL*8 TESTSUM
C      TESTSUM = 0
c      OPEN(UNIT=46,FILE='outtest.txt')
c      DO 10 I = 1,10
c	TESTSUM = TESTSUM + I
c	WRITE(46,*) 'Sum = ',TESTSUM
c      CONTINUE
      implicit none

      write ( *, '(a)' ) '  Hello, world!'
c      WRITE(*,*)'END OF TEST PROGRAM'
      STOP
      END
