package gridGenerator;

import java.awt.Color;
import java.awt.Font;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.math.plot.Plot2DPanel;
import org.math.plot.plotObjects.BaseLabel;

/**
 * 2D orthogonal elliptic grid generator using Laplace equations, univariate
 * stretching functions and tilted parabola tangent fitter.
 * 
 * Direction of solution: West to East Sweep with South to North traverse.
 * 
 * @author Chaitanya Varier
 * @version 09/19/2016
 */

public class OrthogonalStretchedWinslow2D {

    public static final double PI = 3.141592653589793238;
    public static final double e = 2.7182818284590452;

    private static double deltaX = 1.0, deltaY = 1.0;
    private static int boundaryType = 1;
    private static String boundaryName = "";
    private static boolean orthogonalizeBoundary = false,
            orthogonalizeInterior = false, doStretch = false;
    private static boolean orthogonalizeSouth = false,
            orthogonalizeWest = false, orthogonalizeNorth = false,
            orthogonalizeEast = false;
    private static double amplitude, stepSize, radius;
    private static double startX = 0.0, startY = 0.0;
    private static double endX = 100.0, endY = 100.0;
    private static double[][] x_old, y_old, x_new, y_new;
    private static int length;
    private static int height;
    private static double b[][], a[][], c[][], deTerm[][], dTerm[][],
            eTerm[][];
    private static double deltaZi = 1.0, deltaEta = 1.0;
    private static double xdiffmax = 100.0, ydiffmax = 100.0;
    private static double alpha = 1E-10, beta = 1E-10;
    private static int count = 0;
    private static final double diffThreshold = 1E-30, omegaSOR = 1000000.0;
    private static final double checkAngleThreshold = 10.0;
    private static Font plotFont;
    private static Plot2DPanel plot;
    private static Scanner in;
    private static PrintWriter outputInfo, outputInitial, outputFinal;

    public static void main(String args[]) {

        String initialFileName = "Initial_Grid.txt";
        String infoFileName = "Grid_Info.txt";
        String finalFileName = "Final_Grid.txt";

        try {
            outputInitial = new PrintWriter(new FileOutputStream(
                    initialFileName, false));
            outputFinal = new PrintWriter(new FileOutputStream(finalFileName,
                    false));
            outputInfo = new PrintWriter(new FileOutputStream(infoFileName,
                    false));
        } catch (FileNotFoundException e) {
            System.out.println("File error. Program aborted.");
            System.exit(0);
        }

        in = new Scanner(System.in);
        System.out.println("Enter a starting X value");
        startX = in.nextDouble();
        System.out.println("Enter an ending X value");
        endX = in.nextDouble();
        System.out.println("Enter a starting Y value");
        startY = in.nextDouble();
        System.out.println("Enter an ending Y value");
        endY = in.nextDouble();
        outputInfo.println("X interval: [" + startX + ", " + endX + "]");
        outputInfo.println("Y interval: [" + startY + ", " + endY + "]");
        System.out
                .println("Enter a resolution value (distance between adjacent points)");
        deltaX = in.nextDouble();
        deltaY = deltaX;
        length = (int) ((endX - startX) / deltaX) + 1;
        height = (int) ((endY - startY) / deltaY) + 1;
        outputInfo.println("Resolution value: " + deltaX);

        System.out
                .println("\nChoose a boundary type from the following options:\n");
        System.out.println("1. Rectangular");
        System.out.println("2. Gaussian hill");
        System.out.println("3. Absolute value");
        System.out.println("4. Greatest integer (tessellated)");
        System.out.println("5. Forwards step");
        System.out.println("6. Semi-ellipse");
        boundaryType = in.nextInt();

        switch (boundaryType) {

            case 1:
                boundaryName = "rectangular";
                break;
            case 2:
                boundaryName = "Gaussian hill";
                System.out.println("\nEnter an amplitude scale factor [0-1]");
                amplitude = in.nextDouble();
                break;
            case 3:
                boundaryName = "absolute value";
                System.out.println("\nEnter an amplitude");
                amplitude = in.nextDouble();
                break;
            case 4:
                boundaryName = "greatest integer (tessellated)";
                System.out.println("\nEnter an amplitude");
                amplitude = in.nextDouble();
                System.out
                        .println("\nEnter a stepsize (the greater, the wider the step)");
                stepSize = in.nextDouble();
                break;
            case 5:
                boundaryName = "forward step";
                System.out.println("\nEnter a step height");
                amplitude = in.nextDouble();
                break;
            case 6:
                boundaryName = "semi-elliptical";
                System.out.println("\nEnter an amplitude");
                amplitude = in.nextDouble();
                System.out.println("\nEnter a radius");
                radius = in.nextDouble();
                break;

        }

        outputInfo.println("\nBoundary type: " + boundaryName);
        outputInfo.println();
        in.nextLine();

        System.out
                .println("\nDo you wish for the solution grid to have its boundary orthogonalized? (y/n)");

        String oBInput = in.nextLine();

        if (oBInput.equals("y"))
            orthogonalizeBoundary = true;
        else
            orthogonalizeBoundary = false;

        if (orthogonalizeBoundary) {
            System.out
                    .println("\nDo you wish for the solution grid to have its South boundary orthogonalized? (y/n)");

            String sBInput = in.nextLine();

            if (sBInput.equals("y"))
                orthogonalizeSouth = true;
            else
                orthogonalizeSouth = false;

            System.out
                    .println("\nDo you wish for the solution grid to have its West boundary orthogonalized? (y/n)");

            String wBInput = in.nextLine();

            if (wBInput.equals("y"))
                orthogonalizeWest = true;
            else
                orthogonalizeWest = false;

            System.out
                    .println("\nDo you wish for the solution grid to have its North boundary orthogonalized? (y/n)");

            String nBInput = in.nextLine();

            if (nBInput.equals("y"))
                orthogonalizeNorth = true;
            else
                orthogonalizeNorth = false;

            System.out
                    .println("\nDo you wish for the solution grid to have its East boundary orthogonalized? (y/n)");

            String eBInput = in.nextLine();

            if (eBInput.equals("y"))
                orthogonalizeEast = true;
            else
                orthogonalizeEast = false;
        }

        System.out
                .println("\nDo you wish for the solution grid to have its interior orthogonalized? (y/n)");

        String oIInput = in.nextLine();

        if (oIInput.equals("y"))
            orthogonalizeInterior = true;
        else
            orthogonalizeInterior = false;

        System.out
                .println("\nDo you wish for the solution grid to be stretched? (y/n)");

        String dSInput = in.nextLine();

        if (dSInput.equals("y"))
            doStretch = true;
        else
            doStretch = false;

        if (doStretch) {
            System.out
                    .println("\nPlease enter a value for the X-stretching parameter");
            alpha = in.nextDouble();

            System.out
                    .println("\nPlease enter a value for the Y-stretching parameter");
            beta = in.nextDouble();

            outputInfo.println("\nX-stretching parameter: " + alpha);
            outputInfo.println("\nY-Stretching parameter: " + beta);
            outputInfo.println();
        }

        x_old = new double[length][height];
        y_old = new double[length][height];
        x_new = new double[length][height];
        y_new = new double[length][height];

        // ----------DEFINE GRID BOUNDARIES----------

        setBoundaries();

        // ----------LINEAR INTERPOLATION FOR INTERIOR NODES----------

        initializeGrid();

        // Puncture grid

        // punctureGrid(x_old);
        // punctureGrid(y_old);

        /*
         * Display plot of the initial guess grid
         */

        System.out.println("\nInitial\n");
        outputInfo.println("\nInitial\n");
        outputInfo.println();
        checkOrthogonalityBoundary(x_old, y_old);
        System.out.println();
        outputInfo.println("");
        checkOrthogonalityInterior(x_old, y_old);
        System.out.println();
        outputInfo.println("");

        double avgARI = calcAvgAR(x_old, y_old);

        System.out.println("\nThe average aspect ratio of all cells is: "
                + avgARI);
        System.out.println("The standard deviation of all aspect ratios is: "
                + calcStdDevAR(x_old, y_old, avgARI));

        outputInfo.println("\nThe average aspect ratio of all cells is: "
                + avgARI);
        outputInfo.println("The standard deviation of all aspect ratios is: "
                + calcStdDevAR(x_old, y_old, avgARI));
        outputInfo.println();

        String initialName = "Initial grid (Transfinite Interpolation) with "
                + boundaryName + " boundary";

        outputInitial.println(initialName);
        outputInitial.println();
        outputInitial.println("Grid points");
        outputInitial.println();
        outputInitial.println("   X\t\t\t\tY");
        outputInitial.println();

        setUpGrid();
        plotInterior(initialName, Color.blue, x_old, y_old, true);
        plotHorizontalGridLines(initialName, Color.blue, x_old, y_old);
        plotVerticalGridLines(initialName, Color.blue, x_old, y_old);
        printGrid(initialName, x_old, y_old);

        /*
         * Set old matrices to new
         */

        copyMatricesOldToNew(x_new, x_old);
        copyMatricesOldToNew(y_new, y_old);

        while (xdiffmax > diffThreshold || ydiffmax > diffThreshold) {

            // ----------SET UP TRIDIAGONAL MATRICES----------

            b = new double[height][length];
            a = new double[height][length];
            c = new double[height][length];
            deTerm = new double[height][length];
            dTerm = new double[height][length];
            eTerm = new double[height][length];

            // Assemble cofficients
            if (doStretch)
                assembleCoeffStretch();
            else
                assembleCoeffNoStretch();

            // ----------SOLVE MATRICES USING THE THOMAS ALGORITHM----------

            // Calculate new values of x and y
            if (doStretch)
                solveTDMAStretch(x_new);
            else
                solveTDMANoStretch(x_new);

            // Put eTerms in dTerm array for calculating solution of y
            for (int i = 1; i < length - 1; i++) {
                for (int j = 1; j < height - 1; j++) {
                    dTerm[j][i] = eTerm[j][i];
                }
            }

            if (doStretch)
                solveTDMAStretch(y_new);
            else
                solveTDMANoStretch(y_new);

            // ----------ADJUST BOUNDARY NODES FOR ORTHOGONALITY----------

            if (orthogonalizeBoundary)
                orthogonalizeBoundary(x_new, y_new);
            if (orthogonalizeInterior)
                orthogonalizeInterior(x_new, y_new);

            // ----------FIND MAXIMAL DIFFERENCES BETWEEN OLD AND NEW
            // MATRICES----------

            xdiffmax = computeMaxDiff(x_new, x_old);
            ydiffmax = computeMaxDiff(y_new, y_old);

            /*
             * Set old matrices to new
             */

            copyMatricesNewToOld(x_new, x_old);
            copyMatricesNewToOld(y_new, y_old);

            count++;

        }

        /*
         * Display plot of the final transformed grid
         */
        Color darkGreen = new Color(71, 168, 54);
        String isOrthogonalized = (orthogonalizeBoundary || orthogonalizeInterior) ? "orthogonalized"
                : "non-orthogonalized";
        String isStretched = doStretch ? "stretched" : "non-stretched";

        String finalName = "Final " + isOrthogonalized + ", " + isStretched
                + " grid with " + boundaryName + " boundary";

        if (doStretch)
            finalName += " (X-stretching parameter of " + alpha
                    + " and Y-stretching parameter of " + beta + ")";

        outputFinal.println(finalName);
        outputFinal.println();
        outputFinal.println("Grid points");
        outputFinal.println();
        outputFinal.println("   X\t\t\t\tY");
        outputFinal.println();

        setUpGrid();
        plotInterior(finalName, darkGreen, x_new, y_new, false);
        plotHorizontalGridLines(finalName, darkGreen, x_new, y_new);
        plotVerticalGridLines(finalName, darkGreen, x_new, y_new);
        printGrid(finalName, x_new, y_new);

        /*
         * Assess the grid's quality by determining its orthogonality and aspect
         * ratio statistics
         */

        System.out.println("\nFinal" + "\n");
        outputInfo.println("\nFinal" + "\n");
        System.out.println("Completed in " + count + " iterations\n");
        outputInfo.println();
        outputInfo.println("Completed in " + count + " iterations\n");
        checkOrthogonalityBoundary(x_new, y_new);
        System.out.println();
        outputInfo.println();
        checkOrthogonalityInterior(x_new, y_new);

        double avgARF = calcAvgAR(x_new, y_new);

        System.out.println("\nThe average aspect ratio of all cells is: "
                + avgARF);
        System.out.println("The standard deviation of all aspect ratios is: "
                + calcStdDevAR(x_new, y_new, avgARF));

        outputInfo.println();
        outputInfo.println("\nThe average aspect ratio of all cells is: "
                + avgARF);
        outputInfo.println("The standard deviation of all aspect ratios is: "
                + calcStdDevAR(x_new, y_new, avgARF));

        outputInitial.close();
        outputFinal.close();
        outputInfo.close();

    }

    // Set the boundaries for the grid
    public static void setBoundaries() {

        /*
         * Set grid boundary nodes in following order: south, west, north, east
         */

        double x = startX;
        double y = startY;

        for (int i = 0; i < length; i++) {
            y = startY;
            
            switch (boundaryType) {
                case 2:
                    // Gaussian hill
                    y = amplitude
                            * (endY - startY)
                            * Math.exp(-(x - (endX + startX) / 2.0)
                                    * (x - (endX + startX) / 2.0)
                                    / (2 * (endX - startX) * (0.06)
                                            * (endX - startX) * (0.06))) + startY;
                    break;
                    
                case 3:
                    // absolute value
                    if (x <= (startX + endX) / 2.0) {
                        y = amplitude * (x - startX) + startY;
                    } else {
                        y = amplitude * (endX - x) + startY;
                    }
                    break;
                
                case 4:
                    // greatest integer (tessellated)
                    if (x <= (startX + endX) / 2.0) {
                        y = amplitude * Math.floor(1.0 / stepSize * (x - startX))
                                + startY;
                    } else {
                        y = -amplitude
                                * Math.floor(1.0 / stepSize * (x - endX + stepSize))
                                + startY;
                    }
                    break;
                    
                case 5:
                    // square wave
                    y = amplitude
                            * Math.floor(1.0 / ((endX - startX) / 2.0)
                                    * (x - startX)) + startY;
                    break;
                    
                case 6:
                    // semi-ellipse
                    if (x >= (startX + endX) / 2.0 - radius
                            && x <= (startX + endX) / 2.0 + radius)
                        y = amplitude
                                * Math.sqrt(radius * radius
                                        - (x - (startX + endX) / 2)
                                        * (x - (startX + endX) / 2)) + startY;
                    else
                        y = startY;
                    break;
            }
            
            x_old[0][i] = x;
            y_old[0][i] = y;
            x += deltaX;
        }

        // West, north and east boundaries defined as solid edges
        x = startX;
        y = startY;

        for (int j = 0; j < height; j++) {
            x_old[j][0] = x;
            y_old[j][0] = y;
            y += deltaY;
        }

        x = startX;
        y = endY;

        for (int i = 0; i < length; i++) {
            x_old[height - 1][i] = x;
            y_old[height - 1][i] = y;
            x += deltaX;
        }

        x = endX;
        y = startY;

        if (boundaryType == 5) {
            y = amplitude;
        }

        for (int j = 0; j < height; j++) {

            x_old[j][length - 1] = x;
            y_old[j][length - 1] = y;

            if (boundaryType == 5) {
                y += deltaY * (endY - amplitude) / (endY - startY);
            } else {
                y += deltaY;
            }

        }

    }

    // Make an initial guess for the grid using Transfinite Interpolation
    public static void initializeGrid() {

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

    // Displace nodes on the grid randomly
    public static void punctureGrid(double[][] phi) {

        for (int iter = 1; iter <= 9; iter++) {
            int i = (int) (Math.random() * length);
            int j = (int) (Math.random() * height);
            int p = (int) (Math.random() * length);

            phi[j][i] = p;
        }

    }

    // Assemble the coefficients for the tridiagonal matrices if stretching is
    // disabled
    public static void assembleCoeffNoStretch() {

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
    public static void assembleCoeffStretch() {

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

                // f1 = alpha*cos(1/50*pi*x)+k, f2 = (e^(beta*eta)-1)/(e^beta-1)
                double f1Prime = 0, f1DoublePrime = 0;
                if (i <= endX / 2) {
                    f1Prime = alpha * (Math.pow(e, alpha * i))
                            / (Math.pow(e, alpha) - 1);
                    f1DoublePrime = alpha * alpha * (Math.pow(e, alpha * i))
                            / (Math.pow(e, alpha) - 1);
                } else {
                    f1Prime = -alpha * (Math.pow(e, alpha * i))
                            / (Math.pow(e, alpha) - 1);
                    f1DoublePrime = -alpha * alpha * (Math.pow(e, alpha * i))
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
    public static void solveTDMANoStretch(double[][] phi) {

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
    public static void solveTDMAStretch(double[][] phi) {

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

    // Compute the maximum difference between corresponding elements from two 2D
    // grids
    public static double computeMaxDiff(double[][] phi_new, double[][] phi_old) {

        double phidiffmax = Math.abs(phi_new[1][1] - phi_old[1][1]);

        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < length - 1; i++) {
                double phidiff = Math.abs(phi_new[j][i] - phi_old[j][i]);

                if (phidiff > phidiffmax)
                    phidiffmax = phidiff;
            }
        }

        return phidiffmax;

    }

    // Copy the old matrix into the new one (during set-up phase)
    public static void copyMatricesOldToNew(double[][] phi_new,
            double[][] phi_old) {

        // Copy old matrices into new for difference comparison
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                phi_new[j][i] = phi_old[j][i];
            }
        }

    }

    // Copy the new matrix into the old one (during solving phase)
    public static void copyMatricesNewToOld(double[][] phi_new,
            double[][] phi_old) {

        // Copy new matrices into old for difference comparison with SOR
        // (successive over-relaxation)
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                phi_old[i][j] = phi_old[i][j] + omegaSOR
                        * (phi_new[i][j] - phi_old[i][j]);
            }
        }

    }

    // Create necessary grid objects
    public static void setUpGrid() {

        plotFont = new Font(Font.MONOSPACED, Font.PLAIN, 12);

        plot = new Plot2DPanel();

    }

    // Adjust boundary nodes for orthogonality
    public static void orthogonalizeBoundary(double[][] phi1, double[][] phi2) {

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
    public static void orthogonalizeInterior(double[][] phi1, double[][] phi2) {

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

    // Check the orthogonality of all boundary intersections (where vertical and
    // horizontal grid lines cross) by doing the following:
    // find the standard deviation of the angles, the mean angle, the angle that
    // is the farthest
    // away from 90 degrees and the percent of angles within one degree of 90
    // degrees.
    // Here "angle" references the angle formed between
    // intersecting vertical and horizontal grid lines. The intersections of
    // grid lines in this case
    // correspond to boundary grid nodes.

    public static void checkOrthogonalityBoundary(double[][] phi1,
            double[][] phi2) {

        // Find the angles at every grid point formed by intersecting horizontal
        // and vertical grid lines
        // on the boundaries and store them in an array (one for each side).
        double[] anglesS = new double[length], anglesW = new double[height], anglesN = new double[length], anglesE = new double[height];
        double thetaH = 0, thetaV = 0;

        // Obtain and store South boundary angles
        for (int i = 1; i < length - 1; i++) {

            double phi1_p = phi1[0][i], phi2_p = phi2[0][i];

            double phi1_HNext = phi1[0][i + 1] - phi1_p, phi2_HNext = phi2[0][i + 1]
                    - phi2_p, phi1_HPrev = phi1[0][i - 1] - phi1_p, phi2_HPrev = phi2[0][i - 1]
                    - phi2_p;

            double phi1_VNext = phi1[2][i] - phi1[1][i], phi2_VNext = phi2[2][i]
                    - phi2[1][i], phi1_VPrev = phi1_p - phi1[1][i], phi2_VPrev = phi2_p
                    - phi2[1][i];

            thetaH = TiltedParabolaFitter.solveTheta(phi1_HPrev, phi2_HPrev,
                    phi1_HNext, phi2_HNext);
            thetaV = TiltedParabolaFitter.solveTheta(phi1_VPrev, phi2_VPrev,
                    phi1_VNext, phi2_VNext);

            if (thetaH > PI / 2 && thetaV > PI / 2)
                thetaH = PI - thetaH;
            else if (thetaH < PI / 2 && thetaV < PI / 2)
                thetaV = PI - thetaV;

            // convert rad to deg
            anglesS[i] = Math.abs(thetaH - thetaV) * 180.0 / PI;

        }

        // Obtain and store West boundary angles
        for (int j = 1; j < height - 1; j++) {

            double phi1_p = phi1[j][0], phi2_p = phi2[j][0];

            double phi1_HNext = phi1[j][2] - phi1[j][1], phi2_HNext = phi2[j][2]
                    - phi2[j][1], phi1_HPrev = phi1_p - phi1[j][1], phi2_HPrev = phi2_p
                    - phi2[j][1];

            double phi1_VNext = phi1[j + 1][0] - phi1_p, phi2_VNext = phi2[j + 1][0]
                    - phi2_p, phi1_VPrev = phi1[j - 1][0] - phi1_p, phi2_VPrev = phi2[j - 1][0]
                    - phi2_p;

            thetaH = TiltedParabolaFitter.solveTheta(phi1_HPrev, phi2_HPrev,
                    phi1_HNext, phi2_HNext);
            thetaV = TiltedParabolaFitter.solveTheta(phi1_VPrev, phi2_VPrev,
                    phi1_VNext, phi2_VNext);

            if (thetaH > PI / 2 && thetaV > PI / 2)
                thetaH = PI - thetaH;
            else if (thetaH < PI / 2 && thetaV < PI / 2)
                thetaV = PI - thetaV;

            // convert rad to deg
            anglesW[j] = Math.abs(thetaH - thetaV) * 180.0 / PI;

        }

        // Obtain and store North boundary angles
        for (int i = 1; i < length - 1; i++) {

            double phi1_p = phi1[height - 1][i], phi2_p = phi2[height - 1][i];

            double phi1_HNext = phi1[height - 1][i + 1] - phi1_p, phi2_HNext = phi2[height - 1][i + 1]
                    - phi2_p, phi1_HPrev = phi1[height - 1][i - 1] - phi1_p, phi2_HPrev = phi2[height - 1][i - 1]
                    - phi2_p;

            double phi1_VNext = phi1_p - phi1[height - 2][i], phi2_VNext = phi2_p
                    - phi2[height - 2][i], phi1_VPrev = phi1[height - 3][i]
                    - phi1[height - 2][i], phi2_VPrev = phi2[height - 3][i]
                    - phi2[height - 2][i];

            thetaH = TiltedParabolaFitter.solveTheta(phi1_HPrev, phi2_HPrev,
                    phi1_HNext, phi2_HNext);
            thetaV = TiltedParabolaFitter.solveTheta(phi1_VPrev, phi2_VPrev,
                    phi1_VNext, phi2_VNext);

            if (thetaH > PI / 2 && thetaV > PI / 2)
                thetaH = PI - thetaH;
            else if (thetaH < PI / 2 && thetaV < PI / 2)
                thetaV = PI - thetaV;

            // convert rad to deg
            anglesN[i] = Math.abs(thetaH - thetaV) * 180.0 / PI;

        }

        // Obtain and store East boundary angles
        for (int j = 1; j < height - 1; j++) {

            double phi1_p = phi1[j][length - 1], phi2_p = phi2[j][length - 1];

            double phi1_HNext = phi1_p - phi1[j][length - 2], phi2_HNext = phi2_p
                    - phi2[j][length - 2], phi1_HPrev = phi1[j][length - 3]
                    - phi1[j][length - 2], phi2_HPrev = phi2[j][length - 3]
                    - phi2[j][length - 2];

            double phi1_VNext = phi1[j + 1][length - 1] - phi1_p, phi2_VNext = phi2[j + 1][length - 1]
                    - phi2_p, phi1_VPrev = phi1[j - 1][length - 1] - phi1_p, phi2_VPrev = phi2[j - 1][length - 1]
                    - phi2_p;

            thetaH = TiltedParabolaFitter.solveTheta(phi1_HPrev, phi2_HPrev,
                    phi1_HNext, phi2_HNext);
            thetaV = TiltedParabolaFitter.solveTheta(phi1_VPrev, phi2_VPrev,
                    phi1_VNext, phi2_VNext);

            if (thetaH > PI / 2 && thetaV > PI / 2)
                thetaH = PI - thetaH;
            else if (thetaH < PI / 2 && thetaV < PI / 2)
                thetaV = PI - thetaV;

            // convert rad to deg
            anglesE[j] = Math.abs(thetaH - thetaV) * 180.0 / PI;

        }

        // Find the mean of all angles and angle farthest away from 90 deg
        double mean = anglesS[1];
        double maxDev90 = Math.abs(anglesS[1] - 90.0);
        double angleMaxDev90 = anglesS[1];
        int iDev = 1, jDev = 0;

        // South
        for (int i = 2; i < length - 1; i++) {

            mean += anglesS[i];
            if (Math.abs(anglesS[i] - 90.0) > maxDev90) {
                maxDev90 = Math.abs(anglesS[i] - 90.0);
                angleMaxDev90 = anglesS[i];
                iDev = i;
                jDev = 0;
            }
        }

        // West
        for (int j = 1; j < height - 1; j++) {

            mean += anglesW[j];
            if (Math.abs(anglesW[j] - 90.0) > maxDev90) {
                maxDev90 = Math.abs(anglesW[j] - 90.0);
                angleMaxDev90 = anglesW[j];
                iDev = 0;
                jDev = j;
            }
        }

        // North
        for (int i = 1; i < length - 1; i++) {

            mean += anglesN[i];
            if (Math.abs(anglesN[i] - 90.0) > maxDev90) {
                maxDev90 = Math.abs(anglesN[i] - 90.0);
                angleMaxDev90 = anglesN[i];
                iDev = i;
                jDev = height - 1;
            }
        }

        // East
        for (int j = 1; j < height - 1; j++) {

            mean += anglesE[j];
            if (Math.abs(anglesE[j] - 90.0) > maxDev90) {
                maxDev90 = Math.abs(anglesE[j] - 90.0);
                angleMaxDev90 = anglesE[j];
                iDev = length - 1;
                jDev = j;
            }
        }

        mean /= ((length - 2) * 2 + (height - 2) * 2);

        // Find the standard deviation and the percent of angles that are within
        // 10 deg of 90 deg
        double var = 0;
        double stdDev = 0;
        double count90 = 0;
        double percent90 = 0;

        // South
        for (int i = 1; i < length - 1; i++) {

            var = (anglesS[i] - mean) * (anglesS[i] - mean);
            if (Math.abs(anglesS[i] - 90.0) <= checkAngleThreshold) {
                count90++;
            }
        }

        // West
        for (int j = 1; j < height - 1; j++) {

            var = (anglesW[j] - mean) * (anglesW[j] - mean);
            if (Math.abs(anglesW[j] - 90.0) <= checkAngleThreshold) {
                count90++;
            }
        }

        // North
        for (int i = 1; i < length - 1; i++) {

            var = (anglesN[i] - mean) * (anglesN[i] - mean);
            if (Math.abs(anglesN[i] - 90.0) <= checkAngleThreshold) {
                count90++;
            }
        }

        // East
        for (int j = 1; j < height - 1; j++) {

            var = (anglesE[j] - mean) * (anglesE[j] - mean);
            if (Math.abs(anglesE[j] - 90.0) <= checkAngleThreshold) {
                count90++;
            }
        }

        var /= ((length - 2) * 2 + (height - 2) * 2 - 1);
        stdDev = Math.sqrt(var);
        percent90 = count90 / ((length - 2) * 2 + (height - 2) * 2) * 100.0;

        // Print required information
        System.out.println("The standard deviation of all boundary angles is: "
                + stdDev + " deg");
        System.out.println("The mean boundary angle is: " + mean + " deg");
        System.out
                .println("The maximum deviation of all angles from 90 deg on the boundary is: "
                        + maxDev90 + " deg");
        System.out.println("The corresponding angle is: " + angleMaxDev90
                + " deg");
        System.out.println("and is located at position: [" + jDev + ", " + iDev
                + "]");
        System.out
                .println("The percent of angles within 10 deg of 90 deg on the boundary is: "
                        + percent90 + "%");

        outputInfo.println("The standard deviation of all boundary angles is: "
                + stdDev + " deg");
        outputInfo.println("The mean boundary angle is: " + mean + " deg");
        outputInfo
                .println("The maximum deviation of all angles from 90 deg on the boundary is: "
                        + maxDev90 + " deg");
        outputInfo.println("The corresponding angle is: " + angleMaxDev90
                + " deg");
        outputInfo.println("and is located at position: [" + jDev + ", " + iDev
                + "]");
        outputInfo
                .println("The percent of angles within 10 deg of 90 deg on the boundary is: "
                        + percent90 + "%");

    }

    // Check the orthogonality of all interior grid intersections (where
    // vertical and
    // horizontal grid lines cross) by doing the following:
    // find the standard deviation of the angles, the mean angle, the angle that
    // is the farthest
    // away from 90 degrees and the percent of . Here "angle" references the
    // angle formed between
    // intersecting vertical and horizontal grid lines. The intersections of
    // grid lines in this case
    // correspond to interior grid nodes.

    public static void checkOrthogonalityInterior(double[][] phi1,
            double[][] phi2) {

        // Find the angles at every grid point formed by intersecting horizontal
        // and vertical grid lines
        // and store them in a 2D array.
        double angles[][] = new double[height][length];
        double thetaH = 0, thetaV = 0;

        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < length - 1; i++) {

                double phi1_p = phi1[j][i], phi2_p = phi2[j][i];

                double phi1_HNext = phi1[j][i + 1] - phi1_p, phi2_HNext = phi2[j][i + 1]
                        - phi2_p, phi1_HPrev = phi1[j][i - 1] - phi1_p, phi2_HPrev = phi2[j][i - 1]
                        - phi2_p;

                double phi1_VNext = phi1[j + 1][i] - phi1_p, phi2_VNext = phi2[j + 1][i]
                        - phi2_p, phi1_VPrev = phi1[j - 1][i] - phi1_p, phi2_VPrev = phi2[j - 1][i]
                        - phi2_p;

                thetaH = TiltedParabolaFitter.solveTheta(phi1_HPrev,
                        phi2_HPrev, phi1_HNext, phi2_HNext);
                thetaV = TiltedParabolaFitter.solveTheta(phi1_VPrev,
                        phi2_VPrev, phi1_VNext, phi2_VNext);

                if (thetaH > PI / 2 && thetaV > PI / 2)
                    thetaH = PI - thetaH;
                else if (thetaH < PI / 2 && thetaV < PI / 2)
                    thetaV = PI - thetaV;

                // convert rad to deg
                angles[j][i] = Math.abs(thetaH - thetaV) * 180.0 / PI;
            }
        }

        // Find the mean of all angles and angle farthest away from 90 deg
        double mean = 0;
        double maxDev90 = Math.abs(angles[1][1] - 90.0);
        double angleMaxDev90 = angles[1][1];
        int iDev = 0, jDev = 0;

        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < length - 1; i++) {

                mean += angles[j][i];
                if (Math.abs(angles[j][i] - 90.0) > maxDev90) {
                    maxDev90 = Math.abs(angles[j][i] - 90.0);
                    angleMaxDev90 = angles[j][i];
                    iDev = i;
                    jDev = j;
                }
            }
        }

        mean /= ((height - 2) * (length - 2));

        // Find the standard deviation and the percent of angles that are within
        // 10 deg of 90 deg
        double var = 0;
        double stdDev = 0;
        double count90 = 0;
        double percent90 = 0;

        for (int j = 1; j < height - 1; j++) {
            for (int i = 1; i < length - 1; i++) {

                var = (angles[j][i] - mean) * (angles[j][i] - mean);
                if (Math.abs(angles[j][i] - 90.0) <= checkAngleThreshold) {
                    count90++;
                }
            }
        }

        var /= ((height - 2) * (length - 2) - 1);
        stdDev = Math.sqrt(var);
        percent90 = count90 / ((height - 2) * (length - 2)) * 100.0;

        // Print required information
        System.out.println("The standard deviation of all interior angles is: "
                + stdDev + " deg");
        System.out.println("The mean interior angle is: " + mean + " deg");
        System.out
                .println("The maximum deviation of all angles from 90 deg on the interior is: "
                        + maxDev90 + " deg");
        System.out.println("The corresponding angle is: " + angleMaxDev90
                + " deg");
        System.out.println("and is located at position: [" + jDev + ", " + iDev
                + "]");
        System.out
                .println("The percent of angles within 10 deg of 90 deg on the interior is: "
                        + percent90 + "%");

        outputInfo.println("The standard deviation of all interior angles is: "
                + stdDev + " deg");
        outputInfo.println("The mean interior angle is: " + mean + " deg");
        outputInfo
                .println("The maximum deviation of all angles from 90 deg on the interior is: "
                        + maxDev90 + " deg");
        outputInfo.println("The corresponding angle is: " + angleMaxDev90
                + " deg");
        outputInfo.println("and is located at position: [" + jDev + ", " + iDev
                + "]");
        outputInfo
                .println("The percent of angles within 10 deg of 90 deg on the interior is: "
                        + percent90 + "%");

    }

    // Calculate the average aspect ratio of all cells in the grid
    public static double calcAvgAR(double[][] phi1, double[][] phi2) {

        double AR = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, longest = s1, shortest = s2;

        for (int j = 0; j < height - 1; j++) {
            for (int i = 0; i < length - 1; i++) {
                s1 = Math.sqrt((phi1[j][i] - phi1[j][i + 1])
                        * (phi1[j][i] - phi1[j][i + 1])
                        + (phi2[j][i] - phi2[j][i + 1])
                        * (phi2[j][i] - phi2[j][i + 1]));
                s2 = Math.sqrt((phi1[j + 1][i] - phi1[j][i])
                        * (phi1[j + 1][i] - phi1[j][i])
                        + (phi2[j + 1][i] - phi2[j][i])
                        * (phi2[j + 1][i] - phi2[j][i]));
                s3 = Math.sqrt((phi1[j + 1][i + 1] - phi1[j + 1][i])
                        * (phi1[j + 1][i + 1] - phi1[j + 1][i])
                        + (phi2[j + 1][i + 1] - phi2[j + 1][i])
                        * (phi2[j + 1][i + 1] - phi2[j + 1][i]));
                s4 = Math.sqrt((phi1[j + 1][i + 1] - phi1[j][i + 1])
                        * (phi1[j + 1][i + 1] - phi1[j][i + 1])
                        + (phi2[j + 1][i + 1] - phi2[j][i + 1])
                        * (phi2[j + 1][i + 1] - phi2[j][i + 1]));

                longest = s1;
                shortest = s1;

                if (s2 > longest) {
                    longest = s2;
                }

                if (s3 > longest) {
                    longest = s3;
                }

                if (s4 > longest) {
                    longest = s4;
                }

                if (s2 < shortest) {
                    shortest = s2;
                }

                if (s3 < shortest) {
                    shortest = s3;
                }

                if (s4 < shortest) {
                    shortest = s4;
                }

                AR += longest / shortest;

            }
        }

        return AR / ((length - 1) * (height - 1));
    }

    public static double calcStdDevAR(double[][] phi1, double[][] phi2,
            double avgAR) {

        double AR = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, longest = s1, shortest = s2, var = 0, stdDev = 0;

        for (int j = 0; j < height - 1; j++) {
            for (int i = 0; i < length - 1; i++) {
                s1 = Math.sqrt((phi1[j][i] - phi1[j][i + 1])
                        * (phi1[j][i] - phi1[j][i + 1])
                        + (phi2[j][i] - phi2[j][i + 1])
                        * (phi2[j][i] - phi2[j][i + 1]));
                s2 = Math.sqrt((phi1[j + 1][i] - phi1[j][i])
                        * (phi1[j + 1][i] - phi1[j][i])
                        + (phi2[j + 1][i] - phi2[j][i])
                        * (phi2[j + 1][i] - phi2[j][i]));
                s3 = Math.sqrt((phi1[j + 1][i + 1] - phi1[j + 1][i])
                        * (phi1[j + 1][i + 1] - phi1[j + 1][i])
                        + (phi2[j + 1][i + 1] - phi2[j + 1][i])
                        * (phi2[j + 1][i + 1] - phi2[j + 1][i]));
                s4 = Math.sqrt((phi1[j + 1][i + 1] - phi1[j][i + 1])
                        * (phi1[j + 1][i + 1] - phi1[j][i + 1])
                        + (phi2[j + 1][i + 1] - phi2[j][i + 1])
                        * (phi2[j + 1][i + 1] - phi2[j][i + 1]));

                longest = s1;
                shortest = s1;

                if (s2 > longest) {
                    longest = s2;
                }

                if (s3 > longest) {
                    longest = s3;
                }

                if (s4 > longest) {
                    longest = s4;
                }

                if (s2 < shortest) {
                    shortest = s2;
                }

                if (s3 < shortest) {
                    shortest = s3;
                }

                if (s4 < shortest) {
                    shortest = s4;
                }

                AR = longest / shortest;
                var += (AR - avgAR) * (AR - avgAR);

            }
        }
        var /= ((length - 1) * (height - 1) - 1);
        stdDev = Math.sqrt(var);

        return stdDev;

    }

    // Plot interior grid nodes
    public static void plotInterior(String name, Color color, double[][] phi1,
            double[][] phi2, boolean printInitial) {

        // Get X and Y coordinates of all nodes into separate arrays
        double[] X = new double[length * height], Y = new double[length
                * height];
        int iter = 0;

        while (iter <= length * height - 1)
            for (int j = 0; j < height; j++) {
                for (int i = 0; i < length; i++) {
                    X[iter] = phi1[j][i];
                    Y[iter] = phi2[j][i];
                    if (printInitial) {
                        outputInitial.printf("%.10f \t\t %.10f", X[iter],
                                Y[iter]);
                        outputInitial.println();
                    } else {
                        outputFinal
                                .printf("%.10f \t\t %.10f", X[iter], Y[iter]);
                        outputFinal.println();
                    }
                    iter++;
                }
            }

        // Plot grid nodes
        // plot.addScatterPlot(name, color, X, Y);

    }

    // Draw horizontal grid lines
    public static void plotHorizontalGridLines(String name, Color color,
            double[][] phi1, double[][] phi2) {

        double[][] xHoriz = new double[height][length], yHoriz = new double[height][length];
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                xHoriz[j][i] = phi1[i][j];
                yHoriz[j][i] = phi2[i][j];
            }
        }

        for (int j = 0; j < height; j++) {
            plot.addLinePlot(name, color, xHoriz[j], yHoriz[j]);
        }

    }

    // Draw vertical grid lines
    public static void plotVerticalGridLines(String name, Color color,
            double[][] phi1, double[][] phi2) {

        double[][] xVert = new double[height][length], yVert = new double[height][length];
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                xVert[j][i] = phi1[j][i];
                yVert[j][i] = phi2[j][i];
            }
        }

        for (int j = 0; j < height; j++) {
            plot.addLinePlot(name, color, xVert[j], yVert[j]);
        }

    }

    // Display the grid
    public static void printGrid(String name, double[][] phi1, double[][] phi2) {

        // Find the farthest nodes away from the center on the boundaries
        double maxl = phi1[0][0], maxt = phi2[height - 1][0], maxr = phi1[0][length - 1], maxb = phi2[0][0];

        for (int j = 0; j < height; j++) {

            if (phi1[j][0] < maxl)
                maxl = phi1[j][0];
            if (phi1[j][length - 1] > maxr)
                maxr = phi1[j][length - 1];

        }

        for (int i = 0; i < length; i++) {

            if (phi2[0][i] < maxb)
                maxb = phi2[0][i];
            if (phi2[height - 1][i] > maxt)
                maxt = phi2[height - 1][i];

        }

        plot.setFixedBounds(0, maxl, maxr);
        plot.setFixedBounds(1, maxb, maxt);
        plot.setAxisLabels("X", "Y");
        plot.getAxis(0).setLabelPosition(0.5, -0.1);
        plot.getAxis(0).setLabelFont(plotFont);
        plot.getAxis(1).setLabelPosition(-0.15, 0.5);
        plot.getAxis(1).setLabelFont(plotFont);
        BaseLabel title1 = new BaseLabel(name, Color.BLACK, 0.5, 1.1);
        title1.setFont(plotFont);
        plot.addPlotable(title1);
        JFrame frame1 = new JFrame(name);
        frame1.setSize(1000, 1000);
        frame1.setContentPane(plot);
        frame1.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame1.setVisible(true);
        
    }

}