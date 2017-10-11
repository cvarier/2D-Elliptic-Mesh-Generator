package gridGenerator;

import java.io.PrintWriter;

public class MeshStatistics {
    
    private static final double PI = 3.141592653589793238;
    private static final double checkAngleThreshold = 10.0;
    
    
    private int length, height;
    private PrintWriter outputInfo;
    
    public MeshStatistics(int length, int height, PrintWriter outputInfo) {
        this.length = length;
        this.height = height;
        this.outputInfo = outputInfo;
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

    public void checkOrthogonalityBoundary(double[][] phi1,
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

    public void checkOrthogonalityInterior(double[][] phi1,
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
    public double calcAvgAR(double[][] phi1, double[][] phi2) {

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

    public double calcStdDevAR(double[][] phi1, double[][] phi2,
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

}
