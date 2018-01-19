package gridGenerator;

/**
 * Class containing methods which perform the main algorithmic computations throughout the mesh generation process.
 * 
 * @author Chaitanya Varier
 * @version 08/05/2016
 */

public class MeshSolver {

    public static final double e = 2.7182818284590452;
    private double[][] x_new, y_new, x_old, y_old;
    private double deltaZi, deltaEta;
    private int length, height;
    
    public MeshSolver(double[][] x_new, double[][] y_new, int length, int height, double deltaZi, double deltaEta) {
        this.x_new = x_new;
        this.y_new = y_new;
        this.length = length;
        this.height = height;
        this.deltaZi = deltaZi;
        this.deltaEta = deltaEta;
    }
    
    public MeshSolver(double[][] x_old, double[][] y_old, int length, int height) {
        this.x_old = x_old;
        this.y_old = y_old;
        this.length = length;
        this.height = height;
    }
    
    // Make an initial guess for the grid using Transfinite Interpolation
    public void initializeGrid() {

        for (int j = 1; j < height - 1; j++) {
            double dy = j / (height - 1.0);

            for (int i = 1; i < length - 1; i++) {
                double dx = i / (length - 1.0);
                x_old[j][i] = (1.0 - dy)
                        * x_old[0][i]
                        + dy
                        * x_old[height - 1][i]
                        + (1.0 - dx)
                        * x_old[j][0]
                        + dx
                        * x_old[j][length - 1]
                        - (dx * dy * x_old[height - 1][length - 1] + (1.0 - dx)
                                * dy * x_old[height - 1][0] + dx * (1.0 - dy)
                                * x_old[0][length - 1] + (1.0 - dx)
                                * (1.0 - dy) * x_old[0][0]);

                y_old[j][i] = (1.0 - dy)
                        * y_old[0][i]
                        + dy
                        * y_old[height - 1][i]
                        + (1.0 - dx)
                        * y_old[j][0]
                        + dx
                        * y_old[j][length - 1]
                        - (dx * dy * y_old[height - 1][length - 1] + (1.0 - dx)
                                * dy * y_old[height - 1][0] + dx * (1.0 - dy)
                                * y_old[0][length - 1] + (1.0 - dx)
                                * (1.0 - dy) * y_old[0][0]);
            }
        }

    }

    // Assemble the coefficients for the tridiagonal matrices if stretching is
    // disabled
    public void assembleCoeffNoStretch(double[][] b, double[][] a, double[][] deTerm, double[][] dTerm, 
                    double[][] eTerm) {

        for (int i = 1; i < length - 1; i++) {
            for (int j = 1; j < height - 1; j++) {

                double xiNext = x_new[j][i + 1];
                double xiPrev = x_new[j][i - 1];
                double xjNext = x_new[j + 1][i];
                double xjPrev = x_new[j - 1][i];

                double yiNext = y_new[j][i + 1];
                double yiPrev = y_new[j][i - 1];
                double yjNext = y_new[j + 1][i];
                double yjPrev = y_new[j - 1][i];

                double xijNext = x_new[j + 1][i + 1];
                double xijPrev = x_new[j - 1][i - 1];
                double xiPrevjNext = x_new[j + 1][i - 1];
                double xiNextjPrev = x_new[j - 1][i + 1];

                double yijNext = y_new[j + 1][i + 1];
                double yijPrev = y_new[j - 1][i - 1];
                double yiPrevjNext = y_new[j + 1][i - 1];
                double yiNextjPrev = y_new[j - 1][i + 1];

                double x1 = 0.5 * (xiNext - xiPrev) / deltaZi;
                double x2 = 0.5 * (xjNext - xjPrev) / deltaEta;
                double y1 = 0.5 * (yiNext - yiPrev) / deltaZi;
                double y2 = 0.5 * (yjNext - yjPrev) / deltaEta;

                double g11 = x1 * x1 + y1 * y1;
                double g22 = x2 * x2 + y2 * y2;
                double g12 = x1 * x2 + y1 * y2;

                b[j][i] = 2.0 * (g11 / (deltaEta * deltaEta) + g22
                        / (deltaZi * deltaZi));
                a[j][i] = g11 / (deltaEta * deltaEta);
                deTerm[j][i] = g22 / (deltaZi * deltaZi);
                dTerm[j][i] = -0.5 * g12
                        * (xijNext + xijPrev - xiNextjPrev - xiPrevjNext)
                        / (deltaZi * deltaEta);
                eTerm[j][i] = -0.5 * g12
                        * (yijNext + yijPrev - yiNextjPrev - yiPrevjNext)
                        / (deltaZi * deltaEta);
            }
        }

    }

    // Assemble the coefficients for the tridiagonal matrices if stretching is
    // enabled
    public void assembleCoeffStretch(double endX, double alpha, double beta,
                    double[][] b, double[][] a, double[][] c, double[][] dTerm, 
                    double[][] eTerm) {

        for (int i = 1; i < length - 1; i++) {
            for (int j = 1; j < height - 1; j++) {

                double xiNext = x_new[j][i + 1];
                double xiPrev = x_new[j][i - 1];
                double xjNext = x_new[j + 1][i];
                double xjPrev = x_new[j - 1][i];

                double yiNext = y_new[j][i + 1];
                double yiPrev = y_new[j][i - 1];
                double yjNext = y_new[j + 1][i];
                double yjPrev = y_new[j - 1][i];

                double x1 = 0.5 * (xiNext - xiPrev) / deltaZi;
                double x2 = 0.5 * (xjNext - xjPrev) / deltaEta;
                double y1 = 0.5 * (yiNext - yiPrev) / deltaZi;
                double y2 = 0.5 * (yjNext - yjPrev) / deltaEta;

                double g11 = x1 * x1 + y1 * y1;
                double g22 = x2 * x2 + y2 * y2;

                double f1Prime = 0, f1DoublePrime = 0;
                if (i <= (length - 2)/2) {
                    alpha = Math.abs(alpha);
                    f1Prime = alpha * (Math.pow(e, alpha * i))
                            / (Math.pow(e, alpha) - 1);
                    f1DoublePrime = alpha * alpha * (Math.pow(e, alpha * i))
                            / (Math.pow(e, alpha) - 1);
                } else {
                    alpha = -Math.abs(alpha);
                    double iReverse = (length -2)/2 - i;
                    f1Prime = alpha * (Math.pow(e, alpha * iReverse))
                            / (Math.pow(e, alpha) - 1);
                    f1DoublePrime = alpha * alpha * (Math.pow(e, alpha * iReverse))
                            / (Math.pow(e, alpha) - 1);
                }

                double f2Prime = beta * (Math.pow(e, beta * j))
                        / (Math.pow(e, beta) - 1);
                double f2DoublePrime = beta * beta * (Math.pow(e, beta * j))
                        / (Math.pow(e, beta) - 1);

                b[j][i] = 2.0 * (g11 / (deltaEta * deltaEta) + g22
                        / (deltaZi * deltaZi));
                a[j][i] = (f2DoublePrime / f2Prime * g11 / (2 * deltaEta) + g11
                        / (deltaEta * deltaEta));
                c[j][i] = -f2DoublePrime / f2Prime * g11 / (2 * deltaEta) + g11
                        / (deltaEta * deltaEta);
                dTerm[j][i] = g22
                        / (deltaZi)
                        * (xiPrev / deltaZi + f1DoublePrime / f1Prime * xiPrev
                                / 2 + xiNext / deltaZi - f1DoublePrime
                                / f1Prime * xiNext / 2);
                eTerm[j][i] = g22
                        / (deltaZi)
                        * (yiPrev / deltaZi + f1DoublePrime / f1Prime * yiPrev
                                / 2 + yiNext / deltaZi - f1DoublePrime
                                / f1Prime * yiNext / 2);
            }
        }

    }

    // Solve the tridiagonal matrices using the TDMA algorithm if stretching is
    // disabled
    public void solveTDMANoStretch(double[][] phi, double[][] a, double[][] b, double[][] deTerm, double[][] dTerm) {

        int imax = length - 2, jmax = height - 2;
        int imin = 0, jmin = 0;
        double[] P = new double[length], Q = new double[length], bArr = new double[length];

        // Set P(1) to 0 since a(1) = c(1) = 0
        P[jmin] = 0.0;

        // Start West-East sweep
        for (int i = imin + 1; i <= imax; i++) {

            // Set Q(1) to x(1) since x(i) = P(i)x(i+1) + Q(i) and P(1) = 0
            Q[jmin] = phi[jmin][i];

            // Start South-North traverse
            for (int j = jmin + 1; j <= jmax; j++) {

                // Assemble TDMA coefficients, rename North, South, East and
                // West as follows

                // Store a's = c's in P
                P[j] = a[j][i];
                // Store d's/e's in Qx/Qy
                Q[j] = deTerm[j][i] * (phi[j][i + 1] + phi[j][i - 1])
                        + dTerm[j][i];
                // Store b's in bArr
                bArr[j] = b[j][i];

                // Calculate coefficients of recursive formuli
                double term = 1.0 / (bArr[j] - P[j] * P[j - 1]);
                Q[j] = (Q[j] + P[j] * Q[j - 1]) * term;
                P[j] = P[j] * term;

            }

            // Obtain new values of phi (either x or y)
            for (int j = jmax - 1; j > jmin; j--)
                phi[j][i] = P[j] * phi[j + 1][i] + Q[j];

        }

    }

    // Solve the tridiagonal matrices using the TDMA algorithm if stretching is
    // enabled
    public void solveTDMAStretch(double[][] phi, double[][] a, double[][] b, double[][] c, double[][] dTerm) {

        int imax = length - 2, jmax = height - 2;
        int imin = 0, jmin = 0;
        double[] P = new double[length], Q = new double[length], bArr = new double[length], aArr = new double[length];

        // Set P(1) to 0 since a(1) = c(1) = 0
        P[jmin] = 0.0;

        // Start West-East sweep
        for (int i = imin + 1; i <= imax; i++) {

            // Set Q(1) to x(1) since x(i) = P(i)x(i+1) + Q(i) and P(1) = 0
            Q[jmin] = phi[jmin][i];

            // Start South-North traverse
            for (int j = jmin + 1; j <= jmax; j++) {

                // Assemble TDMA coefficients, rename North, South, East and
                // West as follows

                // Store c's in P's
                P[j] = c[j][i];
                // Store d's/e's in Qx/Qy
                Q[j] = dTerm[j][i];
                // Store b's in bArr
                bArr[j] = b[j][i];
                // Store a's in aArr
                aArr[j] = a[j][i];

                // Calculate coefficients of recursive formuli
                double term = 1.0 / (bArr[j] - aArr[j] * P[j - 1]);
                Q[j] = (Q[j] + aArr[j] * Q[j - 1]) * term;
                P[j] = P[j] * term;

            }

            // Obtain new values of phi (either x or y)
            for (int j = jmax - 1; j > jmin; j--)
                phi[j][i] = P[j] * phi[j + 1][i] + Q[j];

        }

    }
    
 // Adjust boundary nodes for orthogonality
    public void orthogonalizeBoundary(double[][] phi1, double[][] phi2, 
                    boolean orthogonalizeSouth, boolean orthogonalizeWest, 
                    boolean orthogonalizeNorth, boolean orthogonalizeEast) {

        // Adjust nodes on following boundaries in order: South, West, North,
        // East
        // 1. Get slope at node on boundary curve
        // 2. Fit tilted parabolas such that current node is always vertex
        // 3. Solve for the current node's new x and y

        double phi2_bPrime = 0;
        double phi1_b, phi2_b, phi1_p, phi2_p;
        double phi1_bNext, phi2_bNext, phi1_bPrev, phi2_bPrev;

        if (orthogonalizeSouth) {
            for (int i = length - 2; i > 0; i--) {
                phi1_b = phi1[0][i];
                phi2_b = phi2[0][i];

                phi1_p = phi1[1][i];
                phi2_p = phi2[1][i];

                phi1_bNext = phi1[0][i - 1] - phi1_b;
                phi2_bNext = phi2[0][i - 1] - phi2_b;
                phi1_bPrev = phi1[0][i + 1] - phi1_b;
                phi2_bPrev = phi2[0][i + 1] - phi2_b;

                if (Math.abs(phi2_bNext - phi2_bPrev) < 1e-12) {
                    phi1[0][i] = phi1_p;
                    phi2[0][i] = phi2_p;
                } else {
                    phi2_bPrime = Math.tan(TiltedParabolaFitter.solveTheta(
                            phi1_bPrev, phi2_bPrev, phi1_bNext, phi2_bNext));
                    phi1[0][i] = (phi2_bPrime * phi1_b + (1.0 / phi2_bPrime)
                            * phi1_p + phi2_p - phi2_b)
                            / (phi2_bPrime + 1.0 / phi2_bPrime);
                    phi2[0][i] = -(1.0 / phi2_bPrime) * (phi1[0][i] - phi1_p)
                            + phi2_p;
                }
            }
        }

        if (orthogonalizeWest) {
            for (int j = 1; j < height - 1; j++) {
                phi1_b = phi1[j][0];
                phi2_b = phi2[j][0];

                phi1_p = phi1[j][1];
                phi2_p = phi2[j][1];

                phi1_bNext = phi1[j + 1][0] - phi1_b;
                phi2_bNext = phi2[j + 1][0] - phi2_b;
                phi1_bPrev = phi1[j - 1][0] - phi1_b;
                phi2_bPrev = phi2[j - 1][0] - phi2_b;

                if (Math.abs(phi2_bNext - phi2_bPrev) < 1e-12) {
                    phi1[j][0] = phi1_p;
                    // phi2[j][0] = phi2_p;
                } else {
                    phi2_bPrime = Math.tan(TiltedParabolaFitter.solveTheta(
                            phi1_bPrev, phi2_bPrev, phi1_bNext, phi2_bNext));
                    phi1[j][0] = (phi2_bPrime * phi1_b + (1.0 / phi2_bPrime)
                            * phi1_p + phi2_p - phi2_b)
                            / (phi2_bPrime + 1.0 / phi2_bPrime);
                    phi2[j][0] = -(1.0 / phi2_bPrime) * (phi1[j][0] - phi1_p)
                            + phi2_p;
                }

            }
        }

        if (orthogonalizeNorth) {
            for (int i = 1; i < length - 1; i++) {
                phi1_b = phi1[height - 1][i];
                phi2_b = phi2[height - 1][i];

                phi1_p = phi1[height - 2][i];
                phi2_p = phi2[height - 2][i];

                phi1_bNext = phi1[height - 1][i + 1] - phi1_b;
                phi2_bNext = phi2[height - 1][i + 1] - phi2_b;
                phi1_bPrev = phi1[height - 1][i - 1] - phi1_b;
                phi2_bPrev = phi2[height - 1][i - 1] - phi2_b;

                if (Math.abs(phi2_bNext - phi2_bPrev) < 1e-12) {
                    phi1[height - 1][i] = phi1_p;
                    // phi2[height-1][i] = phi2_p;
                } else {
                    phi2_bPrime = Math.tan(TiltedParabolaFitter.solveTheta(
                            phi1_bPrev, phi2_bPrev, phi1_bNext, phi2_bNext));
                    phi1[height - 1][i] = (phi2_bPrime * phi1_b
                            + (1.0 / phi2_bPrime) * phi1_p + phi2_p - phi2_b)
                            / (phi2_bPrime + 1.0 / phi2_bPrime);
                    phi2[height - 1][i] = -(1.0 / phi2_bPrime)
                            * (phi1[height - 1][i] - phi1_p) + phi2_p;
                }

            }
        }

        if (orthogonalizeEast) {
            for (int j = height - 2; j > 0; j--) {
                phi1_b = phi1[j][length - 1];
                phi2_b = phi2[j][length - 1];

                phi1_p = phi1[j][length - 2];
                phi2_p = phi2[j][length - 2];

                phi1_bNext = phi1[j - 1][length - 1] - phi1_b;
                phi2_bNext = phi2[j - 1][length - 1] - phi2_b;
                phi1_bPrev = phi1[j + 1][length - 1] - phi1_b;
                phi2_bPrev = phi2[j + 1][length - 1] - phi2_b;

                if (Math.abs(phi2_bNext - phi2_bPrev) < 1e-12) {
                    phi1[j][length - 1] = phi1_p;
                    // phi2[j][length-1] = phi2_p;
                } else {
                    phi2_bPrime = Math.tan(TiltedParabolaFitter.solveTheta(
                            phi1_bPrev, phi2_bPrev, phi1_bNext, phi2_bNext));
                    phi1[j][length - 1] = (phi2_bPrime * phi1_b
                            + (1.0 / phi2_bPrime) * phi1_p + phi2_p - phi2_b)
                            / (phi2_bPrime + 1.0 / phi2_bPrime);
                    phi2[j][length - 1] = -(1.0 / phi2_bPrime)
                            * (phi1[j][length - 1] - phi1_p) + phi2_p;
                }
            }
        }

    }

    // Adjust interior nodes for orthogonality
    public void orthogonalizeInterior(double[][] phi1, double[][] phi2) {

        double phi2_bPrime = 0;
        double phi1_p1, phi2_p1, phi1_p2, phi2_p2;
        double phi1_p1Next, phi2_p1Next, phi1_p1Prev, phi2_p1Prev;

        // Adjust interior grid nodes in a line-by line fashion
        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < length - 1; i++) {
                phi1_p1 = phi1[j][i];
                phi2_p1 = phi2[j][i];

                phi1_p2 = phi1[j + 1][i];
                phi2_p2 = phi2[j + 1][i];

                phi1_p1Next = phi1[j][i - 1] - phi1_p1;
                phi2_p1Next = phi2[j][i - 1] - phi2_p1;
                phi1_p1Prev = phi1[j][i + 1] - phi1_p1;
                phi2_p1Prev = phi2[j][i + 1] - phi2_p1;

                if (Math.abs(phi2_p1Next - phi2_p1Prev) < 1e-12) {
                    phi1[j][i] = phi1_p2;
                    // phi2[j][i] = phi2_p;
                } else {
                    phi2_bPrime = Math
                            .tan(TiltedParabolaFitter.solveTheta(phi1_p1Prev,
                                    phi2_p1Prev, phi1_p1Next, phi2_p1Next));
                    phi1[j][i] = (phi2_bPrime * phi1_p1 + (1.0 / phi2_bPrime)
                            * phi1_p2 + phi2_p2 - phi2_p1)
                            / (phi2_bPrime + 1.0 / phi2_bPrime);
                    phi2[j][i] = -(1.0 / phi2_bPrime) * (phi1[j][i] - phi1_p2)
                            + phi2_p2;
                }
            }
        }

    }
    
}
