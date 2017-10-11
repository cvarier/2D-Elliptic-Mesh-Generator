package gridGenerator;

public class MeshHelper {
    
    private int length, height;
    private double deltaX, deltaY;
    
    public MeshHelper(int length, int height, double deltaX, double deltaY) {
        this.length = length;
        this.height = height;
        this.deltaX = deltaX;
        this.deltaY = deltaY;
    }
    
    // Set the boundaries for the grid
    public void setBoundaries(double startX, double startY, double endX, double endY, 
                    int boundaryType, double amplitude, double stepSize, double radius, 
                    double[][] x_old, double[][] y_old) {

        /*
         * Set grid boundary nodes in following order: south, west, north, east
         */
        double x = startX;
        double y = startY;

        for (int i = 0; i < length; i++) {
            y = startY;
            
/*            // Crazy function
            y = 50*Math.exp(-((10*x/60)-3.5)*((10*x/60)-3.5))+50*Math.exp(-((6*x/60)-7.0)*((6*x/60)-7.0));*/
            
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
            x_old[j][0] = x; //10*Math.exp(-((6*y/60)-7.0)*((6*y/60)-7.0)); // <--crazy function; 
            y_old[j][0] = y;
            y += deltaY;
        }

        x = startX;
        y = endY;

        for (int i = 0; i < length; i++) {
            x_old[height - 1][i] = x;
            y_old[height - 1][i] = y; /*-5*Math.exp(-((20*x/60)-3.5)*((20*x/60)-3.5))-
                            20*Math.exp(-((10*x/60)-7.0)*((10*x/60)-7.0))-
                            15*Math.exp(-((26*x/100)-20.0)*((26*x/100)-20.0))+100;*///<-- crazy function; 
            x += deltaX;
        }

        x = endX;
        y = startY;

        if (boundaryType == 5) {
            y = amplitude;
        }

        for (int j = 0; j < height; j++) {

            x_old[j][length - 1] = x; //-10*Math.exp(-((10*y/60)-3.5)*((10*y/60)-3.5))-30*Math.exp(-((6*y/60)-7.0)*((6*y/60)-7.0))+100; 
                            // <-- crazy function 
            y_old[j][length - 1] = y;

            if (boundaryType == 5) {
                y += deltaY * (endY - amplitude) / (endY - startY);
            } else {
                y += deltaY;
            }

        }

    }

    // Compute the maximum difference between corresponding elements from two 2D
    // grids
    public double computeMaxDiff(double[][] phi_new, double[][] phi_old) {

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
    public void copyMatricesOldToNew(double[][] phi_new,
            double[][] phi_old) {

        // Copy old matrices into new for difference comparison
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                phi_new[j][i] = phi_old[j][i];
            }
        }

    }

    // Copy the new matrix into the old one (during solving phase)
    public void copyMatricesNewToOld(double[][] phi_new,
            double[][] phi_old, double omegaSOR) {

        // Copy new matrices into old for difference comparison with SOR
        // (successive over-relaxation)
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < length; i++) {
                phi_old[i][j] = phi_old[i][j] + omegaSOR
                        * (phi_new[i][j] - phi_old[i][j]);
            }
        }

    }
    
    // Displace nodes on the grid randomly
    public void punctureGrid(double[][] phi) {

        for (int iter = 1; iter <= 9; iter++) {
            int i = (int) (Math.random() * length);
            int j = (int) (Math.random() * height);
            int p = (int) (Math.random() * length);

            phi[j][i] = p;
        }

    }
    
}
