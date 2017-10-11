package gridGenerator;

/**
 * Curve fits a set of three points with a parabola at an angle theta such that the middle point
 * (between (x1,y1) and (x2,y2)) is the vertex.
 * 
 * @author Chaitanya Varier
 * @version 08/05/2016
 */

public class TiltedParabolaFitter {
	
	public static final double PI = 3.141592653589793238;
	
	// Solves trig equation for theta using Bisection method
	public static double solveTheta(double x1, double y1, double x2, double y2){
		
		double leftPrev = -1, rightPrev = 1, leftNext, rightNext, m, aInitial = 0.0, bInitial = PI;
		
		// Find theta using bisection method
		
		leftNext = aInitial;
		rightNext = bInitial;
		
		while (leftNext != leftPrev || rightNext != rightPrev) {
			leftPrev = leftNext;
			rightPrev = rightNext;
			
			m = (leftNext+rightNext)/2;
			
			if (f(leftNext, x1, y1, x2, y2) * f(m, x1, y1, x2, y2) < 0 )
				rightNext = m;
			else
				leftNext = m;
		}
		
		return (leftNext + rightNext)/2; 
		
	}
	
	// The trig equation in functional form to be solved for
	public static double f (double theta, double x1, double y1, double x2, double y2) {
		
		return (-x1*Math.sin(theta) + y1*Math.cos(theta)) 
				*(x2*Math.cos(theta) + y2*Math.sin(theta))*(x2*Math.cos(theta) + y2*Math.sin(theta)) 
				- (-x2*Math.sin(theta) + y2*Math.cos(theta)) 
				*(x1*Math.cos(theta) + y1*Math.sin(theta))*(x1*Math.cos(theta) + y1*Math.sin(theta));
		
	}

}