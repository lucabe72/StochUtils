This set of utilities and test programs implement the Gamma-bound presented
in the paper "Efficient and Robust Probabilistic Guarantees for Real-Time
Tasks". The utility functions in the "src" directory are implemented in fairly
portable C code, and can be used for on-line admission tests. The test programs
from the "Tools" directory show how to use the programming interface of such
utility functions.

The test programs are:
avg:	Gets as input a PMF for the execution times and a PMF for z (number of
	server periods in an inter-arrival time). Then, computes the average
	execution time, the average z, and the average utilization of the task.
	These values are used to compute the minimum (based on the average
	utilization) and maximum (based on the worst-case utilization) values
	for the reservation's budget Q^s.
	Then, the PMF of Y is computed for the values of Q^s between minimum
	and maximum, and the sum in Equation (6) of the paper is computed for
	various values of Gamma.
pmf-yt:	Gets as input a PMF for the execution times and a PMF for the
	inter-arrival times. Inter-arrival times are not assumed to be multiple
	of the server period, and their PMF is resampled to obtain the PMF of z
	(computing a pessimistic bound with inter-arrival times multiple of
	the server period).
	Then, the PMF of Y is computed, and a value for Gamma is found by using
	dichotomic search. Such a value is used to compute a bound for the PMF
	of v, and then a pessimistic bound for the probabilistic deadlines is
	computed based on v.
	Parameters:
		-t <T^s>: set the reservation's period
		-q <Q^s>: set the reservation's maximum budget
		-T <maxy>: set the maximum value of Y which is interesting in the analysis
		-s <samples>: set the size of the execution times PMF (in samples)
pmf-y:	Similar to pmf-yt, but directly reads the PMF of z
confint:	Compute average and confidence interval of a set of numbers.

