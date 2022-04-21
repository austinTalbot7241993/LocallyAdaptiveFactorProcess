matrix_utils.o: matrix_utils.c
	gcc -c matrix_utils.c -I/usr/local/include

clean:
	rm *.o

testMatrixUtils.o: testMatrixUtils.c
	gcc -c testMatrixUtils.c -I/usr/local/include

testMatrixUtils: testmatrixUtils.o gmlib.o matrix_utils.o testMatrixUtils.c
	gcc -o testMatrixUtils testMatrixUtils.o gmlib.o matrix_utils.o -lgsl -lgslcblas -lmarray 

mdarray.o: mdarray.c
	gcc -c mdarray.c -I/usr/local/include

testMdarray.o: testMdarray.c
	gcc -c testMdarray.c

testMdarray: testMdarray.o gmlib.o matrix_utils.o mdarray.o
	gcc -o testMdarray testMdarray.o gmlib.o matrix_utils.o mdarray.o -lgsl -lgslcblas -lmarray

mvnorm.o: mvnorm.c
	gcc -c mvnorm.c -I/usr/local/include

marraytools.o: marraytools.c
	gcc -c marraytools.c -I/usr/local/include

sstools2.o: sstools2.c
	gcc -c sstools2.c -I/usr/local/include

sstools2Test: sstools2.o marraytools.o mvnorm.o matrix_utils.o
	gcc -o sstools2Test sstools2.o matrix_utils.o marraytools.o mvnorm.o -lgsl -lgslcblas -lmarray -lm

gmlib.o: gmlib.c
	gcc -c gmlib.c -I/usr/local/include

KalmanFilter.o: KalmanFilter.c
	gcc -c KalmanFilter.c -I/usr/local/include

KalmanFilter2.o:KalmanFilter2.c
	gcc -c KalmanFilter2.c

testFile: testKalman.o marraytools.o matrix_utils.o KalmanFilter.o gmlib.o mvnorm.o
	gcc -o testFile testKalman.o KalmanFilter.o marraytools.o mvnorm.o matrix_utils.o gmlib.o -lgsl -lgslcblas -lmarray -lm

testKalman.o: testKalman.c
	gcc -c testKalman.c

testKalman2.o: testKalman2.c
	gcc -c testKalman2.c

testKalman2: testKalman2.o mdarray.o matrix_utils.o gmlib.o KalmanFilter2.o
	gcc -o testKalman2 testKalman2.o KalmanFilter2.o gmlib.o matrix_utils.o mdarray.o -lgsl -lgslcblas -lmarray

fastStateSmoother.o: fastStateSmoother.c
	gcc -c fastStateSmoother.c -I/usr/local/include

testFSS.o: testFSS.c
	gcc -c testFSS.c

testFSS: testFSS.o marraytools.o matrix_utils.o fastStateSmoother.o gmlib.o
	gcc -o testFSS testFSS.o fastStateSmoother.o marraytools.o mvnorm.o matrix_utils.o gmlib.o -lgsl -lgslcblas -lmarray -lm

testFSS2.o: testFSS2.c
	gcc -c testFSS2.c

testFSS2: testFSS2.o mdarray.o matrix_utils.o gmlib.o fastStateSmoother2.o
	gcc -o testFSS2 testFSS2.o fastStateSmoother2.o gmlib.o matrix_utils.o mdarray.o -lgsl -lgslcblas -lmarray
testNGPtools.o: testNGPtools.c marraytools.o matrix_utils.o NGPtools.o mvnorm.o KalmanFilter.o fastStateSmoother.o gmlib.o
	gcc -c testNGPtools.c

testNGPtools: testNGPtools.o NGPtools.o matrix_utils.o marraytools.o mvnorm.o KalmanFilter.o fastStateSmoother.o gmlib.o
	gcc -o testNGPtools testNGPtools.o matrix_utils.o marraytools.o mvnorm.o KalmanFilter.o fastStateSmoother.o gmlib.o NGPtools.o -lmarray -lgsl -lgslcblas

testSSsimulate.o: testSSsimulate.c
	gcc -c testSSsimulate.c

testSSsimulate: testSSsimulate.o KalmanFilter.o SSsimulate.o fastStateSmoother.o gmlib.o marraytools.o mvnorm.o matrix_utils.o 
	gcc -o testSSsimulate testSSsimulate.o SSsimulate.o matrix_utils.o marraytools.o mvnorm.o KalmanFilter.o fastStateSmoother.o gmlib.o -lmarray -lgsl -lgslcblas -lm

SSsimulate2.o: SSsimulate2.c
	gcc -c SSsimulate2.c

testSSsimulate2.o: testSSsimulate2.c
	gcc -c testSSsimulate2.c

testSSsimulate2: testSSsimulate2.o gmlib.o matrix_utils.o mdarray.o KalmanFilter2.o fastStateSmoother2.o mvnorm.o SSsimulate2.o
	gcc -o testSSsimulate2 testSSsimulate2.o gmlib.o matrix_utils.o mdarray.o mvnorm.o KalmanFilter2.o fastStateSmoother2.o SSsimulate2.o -lmarray -lgsl -lgslcblas


testNGPsample.o: testNGPsample.c
	gcc -c testNGPsample.c

testNGPsample: testNGPsample.o NGPsample.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o matrix_utils.o marraytools.o mvnorm.o gmlib.o 
	gcc -o testNGPsample testNGPsample.o NGPsample.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o matrix_utils.o marraytools.o mvnorm.o gmlib.o -lmarray -lgsl -lgslcblas -lm

testMVNorm.o: testMVNorm.c
	gcc -c testMVNorm.c

testMVNorm: testMVNorm.o mvnorm.o matrix_utils.o 
	gcc -o testMVNorm testMVNorm.o mvnorm.o matrix_utils.o -lgsl -lgslcblas 

testLAFsamplingRoutines.o: testLAFsamplingRoutines.c
	gcc -c testLAFsamplingRoutines.c

testLAFsamplingRoutines: testLAFsamplingRoutines.o matrix_utils.o gmlib.o marraytools.o mvnorm.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o NGPsample.o
	gcc -o testLAFsamplingRoutines testLAFsamplingRoutines.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o NGPsample.o matrix_utils.o mvnorm.o gmlib.o marraytools.o -lgsl -lmarray -lgslcblas -lm


NGPlaf.o: NGPlaf.c
	gcc -c NGPlaf.c -I/usr/local/include

NGPlafTest: NGPlaf.o gmlib.o matrix_utils.o marraytools.o mvnorm.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o NGPsample.o
	gcc -O3 -o NGPlafTest NGPlaf.o gmlib.o mvnorm.o matrix_utils.o marraytools.o KalmanFilter.o fastStateSmoother.o SSsimulate.o NGPtools.o NGPsample.o -lgsl -lmarray -lgslcblas

NGPlaf2.o: NGPlaf2.c
	gcc -c NGPlaf2.c -I/usr/local/include

NGPlafTest2: NGPlaf2.o gmlib.o matrix_utils.o mdarray.o mvnorm.o KalmanFilter2.o fastStateSmoother2.o SSsimulate2.o NGPtools2.o NGPsample2.o
	gcc -o NGPlafTest2 NGPlaf2.o gmlib.o mvnorm.o matrix_utils.o mdarray.o KalmanFilter2.o fastStateSmoother2.o SSsimulate2.o NGPtools2.o NGPsample2.o -lgsl -lmarray -lgslcblas

testNGPmcmc.o: testNGPmcmc.c
	gcc -c testNGPmcmc.c

testNGPmcmc: testNGPmcmc.o gmlib.o matrix_utils.o mdarray.o SSsimulate2.o NGPmcmc.o NGPtools2.o
	gcc -o testNGPmcmc testNGPmcmc.o gmlib.o matrix_utils.o mvnorm.o mdarray.o KalmanFilter2.o NGPtools2.o fastStateSmoother2.o SSsimulate2.o NGPmcmc.o -lgsl -lgslcblas -lpthread

runNGPmcmc.o: runNGPmcmc.c
	gcc -c runNGPmcmc.c

runNGPmcmc: runNGPmcmc.o gmlib.o matrix_utils.o mdarray.o SSsimulate2.o NGPmcmc.o NGPtools2.o
	gcc -o runNGPmcmc runNGPmcmc.o gmlib.o matrix_utils.o mvnorm.o mdarray.o KalmanFilter2.o NGPtools2.o fastStateSmoother2.o SSsimulate2.o NGPmcmc.o -lgsl -lgslcblas -lpthread



