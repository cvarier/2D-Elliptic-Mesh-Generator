package gridGenerator;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import javax.imageio.stream.FileImageOutputStream;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.math.plot.Plot2DPanel;
import org.math.plot.plotObjects.BaseLabel;

/**
 * 2D orthogonal elliptic mesh generator using Winslow equations, univariate
 * stretching functions and a tilted parabola tangent fitter.
 * 
 * Direction of solution: West to East Sweep with South to North traverse.
 * 
 * @author Chaitanya Varier
 * @version 10/10/2017
 */

public class EllipticMeshGenerator2D {

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
    private static Font plotFont;
    private static Plot2DPanel plot;
    private static Scanner in;
    private static PrintWriter outputInfo, outputInitial, outputFinal;
    private static double animationFPS = 1;
    private static int animationFrameCount = 5;

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
        
        MeshHelper meshHelper = new MeshHelper(length, height, deltaX, deltaY);

        // ----------DEFINE GRID BOUNDARIES----------

        meshHelper.setBoundaries(startX, startY, endX, endY, boundaryType, amplitude, stepSize, radius, x_old, y_old);

        // ----------LINEAR INTERPOLATION FOR INTERIOR NODES----------

        MeshSolver meshSolverInitial = new MeshSolver(x_old, y_old, length, height);
        meshSolverInitial.initializeGrid();

        // Puncture grid

        /*punctureGrid(x_old);
          punctureGrid(y_old);*/

        /*
         * Display plot of the initial guess grid
         */
        
        MeshStatistics meshStatisticsInitial = new MeshStatistics(length, height, outputInfo);

        System.out.println("\nInitial\n");
        outputInfo.println("\nInitial\n");
        outputInfo.println();
        meshStatisticsInitial.checkOrthogonalityBoundary(x_old, y_old);
        System.out.println();
        outputInfo.println("");
        meshStatisticsInitial.checkOrthogonalityInterior(x_old, y_old);
        System.out.println();
        outputInfo.println("");

        double avgARI = meshStatisticsInitial.calcAvgAR(x_old, y_old);
        double arStdDevI = meshStatisticsInitial.calcStdDevAR(x_old, y_old, avgARI);

        System.out.println("\nThe average aspect ratio of all cells is: "
                + avgARI);
        System.out.println("The standard deviation of all aspect ratios is: "
                + arStdDevI);

        outputInfo.println("\nThe average aspect ratio of all cells is: "
                + avgARI);
        outputInfo.println("The standard deviation of all aspect ratios is: "
                + arStdDevI);
        outputInfo.println();

        String initialName = "Initial grid (Transfinite Interpolation) with "
                + boundaryName + " boundary";

        outputInitial.println(initialName);
        outputInitial.println();
        outputInitial.println("Grid points");
        outputInitial.println();
        outputInitial.println("   X\t\t\t\tY");
        outputInitial.println();
        
        // ----------- FOR ANIMATION ------------
        
        // Specifiy frames per second and slowdown factor
        double slowMo = 1.0;
        ImageOutputStream animationOutput = null;
        GifSequenceWriter animationWriter = null;

        try {
            // Initialize the gif writer
            animationOutput = new FileImageOutputStream(new File("animatedMesh.gif"));
            animationWriter = new GifSequenceWriter(animationOutput, 1, (int)(1000*slowMo/animationFPS), true);
        } catch (IOException e) {
            e.printStackTrace();
        }

        setUpGrid();
        plotInterior(initialName, Color.blue, x_old, y_old, true, true);
        plotHorizontalGridLines(initialName, Color.blue, x_old, y_old);
        plotVerticalGridLines(initialName, Color.blue, x_old, y_old);
        try {
            printGrid(initialName, x_old, y_old, true, animationWriter);
        } catch (IOException e2) {
            e2.printStackTrace();
        }

        /*
         * Set old matrices to new
         */

        meshHelper.copyMatricesOldToNew(x_new, x_old);
        meshHelper.copyMatricesOldToNew(y_new, y_old);
        
        // ----------- MAIN LOOP -----------

        while (xdiffmax > diffThreshold || ydiffmax > diffThreshold) {

            // ----------SET UP TRIDIAGONAL MATRICES----------

            b = new double[height][length];
            a = new double[height][length];
            c = new double[height][length];
            deTerm = new double[height][length];
            dTerm = new double[height][length];
            eTerm = new double[height][length];
            
            MeshSolver meshSolver = new MeshSolver(x_new, y_new, length, height, deltaZi, deltaEta);

            // Assemble cofficients
            if (doStretch)
                meshSolver.assembleCoeffStretch(endX, alpha, beta, b, a, c, dTerm, eTerm);
            else
                meshSolver.assembleCoeffNoStretch(b, a, deTerm, dTerm, eTerm);

            // ----------SOLVE MATRICES USING THE THOMAS ALGORITHM----------

            // Calculate new values of x and y
            if (doStretch)
                meshSolver.solveTDMAStretch(x_new, a, b, c, dTerm);
            else
                meshSolver.solveTDMANoStretch(x_new, a, b, deTerm, dTerm);

            // Put eTerms in dTerm array for calculating solution of y
            for (int i = 1; i < length - 1; i++) {
                for (int j = 1; j < height - 1; j++) {
                    dTerm[j][i] = eTerm[j][i];
                }
            }

            if (doStretch)
                meshSolver.solveTDMAStretch(y_new, a, b, c, dTerm);
            else
                meshSolver.solveTDMANoStretch(y_new, a, b, deTerm, dTerm);

            // ----------ADJUST BOUNDARY NODES FOR ORTHOGONALITY----------

            if (orthogonalizeBoundary)
                meshSolver.orthogonalizeBoundary(x_new, y_new, orthogonalizeSouth, orthogonalizeWest, 
                                orthogonalizeNorth, orthogonalizeEast);
            if (orthogonalizeInterior)
                meshSolver.orthogonalizeInterior(x_new, y_new);

            // ----------FIND MAXIMAL DIFFERENCES BETWEEN OLD AND NEW
            // MATRICES----------

            xdiffmax = meshHelper.computeMaxDiff(x_new, x_old);
            ydiffmax = meshHelper.computeMaxDiff(y_new, y_old);

            /*
             * Set old matrices to new
             */

            meshHelper.copyMatricesNewToOld(x_new, x_old, omegaSOR);
            meshHelper.copyMatricesNewToOld(y_new, y_old, omegaSOR);

            count++;
            
            // if count is a multiple of animationFrameCount, print the grid to produce an animated effect
            
            if (count % animationFrameCount == 0) {
                Color darkGreen = new Color(71, 168, 54);
                String isOrthogonalized = (orthogonalizeBoundary || orthogonalizeInterior) ? "orthogonalized"
                        : "non-orthogonalized";
                String isStretched = doStretch ? "stretched" : "non-stretched";

                String name = "Frame " + count + " " + isOrthogonalized + ", " + isStretched
                        + " grid with " + boundaryName + " boundary";
                
                setUpGrid();
                plotInterior(name, darkGreen, x_new, y_new, false, false);
                plotHorizontalGridLines(name, darkGreen, x_new, y_new);
                plotVerticalGridLines(name, darkGreen, x_new, y_new);
                try {
                    printGrid(name, x_new, y_new, true, animationWriter);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

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
        plotInterior(finalName, darkGreen, x_new, y_new, false, true);
        plotHorizontalGridLines(finalName, darkGreen, x_new, y_new);
        plotVerticalGridLines(finalName, darkGreen, x_new, y_new);
        try {
            printGrid(finalName, x_new, y_new, true, animationWriter);
        } catch (IOException e1) {
            e1.printStackTrace();
        }
        
        try {
            animationWriter.close();
            animationOutput.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        /*
         * Assess the grid's quality by determining its orthogonality and aspect
         * ratio statistics
         */
        
        MeshStatistics meshStatisticsFinal = new MeshStatistics(length, height, outputInfo);

        System.out.println("\nFinal" + "\n");
        outputInfo.println("\nFinal" + "\n");
        System.out.println("Completed in " + count + " iterations\n");
        outputInfo.println();
        outputInfo.println("Completed in " + count + " iterations\n");
        meshStatisticsFinal.checkOrthogonalityBoundary(x_new, y_new);
        System.out.println();
        outputInfo.println();
        meshStatisticsFinal.checkOrthogonalityInterior(x_new, y_new);

        double avgARF = meshStatisticsFinal.calcAvgAR(x_new, y_new);
        double arStdDevF = meshStatisticsFinal.calcStdDevAR(x_new, y_new, avgARF);

        System.out.println("\nThe average aspect ratio of all cells is: "
                + avgARF);
        System.out.println("The standard deviation of all aspect ratios is: "
                + arStdDevF);

        outputInfo.println();
        outputInfo.println("\nThe average aspect ratio of all cells is: "
                + avgARF);
        outputInfo.println("The standard deviation of all aspect ratios is: "
                + arStdDevF);

        outputInitial.close();
        outputFinal.close();
        outputInfo.close();
       

    }

    // Create necessary grid objects
    public static void setUpGrid() {

        plotFont = new Font(Font.MONOSPACED, Font.PLAIN, 12);
        plot = new Plot2DPanel();

    }

    // Plot interior grid nodes
    public static void plotInterior(String name, Color color, double[][] phi1,
            double[][] phi2, boolean printInitial, boolean writeToFile) {

        // Get X and Y coordinates of all nodes into separate arrays
        double[] X = new double[length * height], Y = new double[length
                * height];
        int iter = 0;

        while (iter <= length * height - 1)
            for (int j = 0; j < height; j++) {
                for (int i = 0; i < length; i++) {
                    X[iter] = phi1[j][i];
                    Y[iter] = phi2[j][i];
                    if (writeToFile) {
                        if (printInitial) {
                            outputInitial.printf("%.10f \t\t %.10f", X[iter],
                                    Y[iter]);
                            outputInitial.println();
                        } else {
                            outputFinal
                                    .printf("%.10f \t\t %.10f", X[iter], Y[iter]);
                            outputFinal.println();
                        }
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
    public static void printGrid(String name, double[][] phi1, double[][] phi2, boolean animate, 
                    GifSequenceWriter gifWriter) throws IOException {

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
        
        if (!animate) {
            frame1.setVisible(true);
        } else {
            //frame1.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
            frame1.setVisible(true);
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            // Create image from JFrame and write to the animation
            BufferedImage image = new BufferedImage(1000,1000, BufferedImage.TYPE_INT_RGB);
            Graphics2D graphics2D = image.createGraphics();
            frame1.paint(graphics2D);
            gifWriter.writeToSequence(image);
            //frame1.dispatchEvent(new WindowEvent(frame1, WindowEvent.WINDOW_CLOSING));
        }
        
    }

}