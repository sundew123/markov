import java.awt.*;
import java.util.Arrays;
public class entropy {
    public static Object[] tree = new Object[3];
    public static DrawingPanel drawingPanel = new DrawingPanel(11, 11);
    public static Graphics graphics = drawingPanel.getGraphics();
    public static boolean[][] pixels;
    public static Object[] horizontal;
    public static Object[] vertical;
    public static boolean a;
    public static double[][] A;
    public static double[][][][] lookup = new double[drawingPanel.getHeight()][drawingPanel.getWidth()][2][10];
    public static void main(String args[]) {
        a = true;
        //initiate image (press Shift+L to start)
        Image img = drawingPanel.loadImage("C:\\Users\\leomi\\Downloads\\loss.png");
        graphics.drawImage(img, 0, 0, drawingPanel);
        //alternative option to draw own image
        int[] pos = new int[2];
        drawingPanel.setAntiAlias(false);
        drawingPanel.onMouseDown((x, y) -> {if (a) { pos[0] = x; pos[1] = y; drawingPanel.setPixel(x, y, Color.BLACK); }});
        drawingPanel.onDrag((x, y) -> {if (a) { graphics.drawLine(pos[0], pos[1], x, y); pos[0] = x; pos[1] = y; }});
        drawingPanel.onKeyUp((e) -> {if (e == 'L' && a) { initiate(); }});
    }
    public static void initiate() {
        a = false;
        //convert pixels to boolean form
        int[][] pixels1 = drawingPanel.getPixelsRGB();
        pixels = new boolean[drawingPanel.getHeight()][drawingPanel.getWidth()];
        for (int i = 0; i < pixels1.length; i++) {
            for (int j = 0; j < pixels1[i].length; j++) {
                pixels[i][j] = (pixels1[i][j] == -1);
            }
        }
        //construct Markov chain tree
        addToTree(1);
        double[] temp = {((double)((int[])tree[2])[0] + 1) / ((double)((int[])tree[2])[1] + (double)((int[])tree[2])[0] + 2), ((double)((int[])tree[2])[1] + 1) / ((double)((int[])tree[2])[1] + (double)((int[])tree[2])[0] + 2)};
        tree[2] = temp;
        calculate(tree);
        horizontal = tree.clone();
        tree = new Object[3];
        transpose(pixels);
        addToTree(1);
        double[] temp2 = {((double)((int[])tree[2])[0] + 1) / ((double)((int[])tree[2])[1] + (double)((int[])tree[2])[0] + 2), ((double)((int[])tree[2])[1] + 1) / ((double)((int[])tree[2])[1] + (double)((int[])tree[2])[0] + 2)};
        tree[2] = temp2;
        calculate(tree);
        vertical = tree.clone();
        transpose(pixels);
        getLookup();
        A = getA();
        System.out.println(A);
        drawingPanel.clear();
        draw();
    }
    public static void draw() {
        //calculates and displays local entropy of each pixel
        double[][] buffer = new double[pixels.length][pixels[0].length];
        double max = 0;
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                double[] prob = get(i, j, pixels[i][j]);
                prob = flDivide(prob, flAdd(prob, get(i, j, !pixels[i][j])));
                buffer[i][j] = -Math.log(prob[0]) - prob[1] * Math.log(2);
                if (buffer[i][j] > max) {
                    max = buffer[i][j];
                }
            }
        }
        //prints total entropy
        double s = 0;
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                s += -Math.log(lookup[i][j][0][9] * A[Math.min(j, 9)][Math.min(i, 9)] + lookup[i][j][1][9] * (1 - A[Math.min(j, 9)][Math.min(i, 9)]));
            }
        }
        System.out.println(s);
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                Color cc = new Color((int)Math.floor(255 * buffer[i][j] / max), (int)Math.floor(255 * buffer[i][j] / max), (int)Math.floor(255 * buffer[i][j] / max));
                drawingPanel.setPixel(j, i, cc);
            }
        }
    }
    public static double[] get(int I, int J, boolean K) {
        //calculates value of local entropy given pixel
        double[] prod = {1, 0};
        boolean temp2 = pixels[I][J];
        pixels[I][J] = K;
        for (int i = I; i < pixels.length && i < I + 10; i++) {
            double hP;
            double vP;
            Object[] temp = horizontal;
            int k = -1;
            while (J + k >= 0 && temp[pixels[i][J + k] ? 1 : 0] != null) {
                temp = (Object[])temp[pixels[i][J + k] ? 1 : 0];
                k--;
            }
            hP = ((double[])temp[2])[pixels[i][J] ? 1 : 0];
            temp = vertical;
            k = -1;
            while (i + k >= 0 && temp[pixels[i + k][J] ? 1 : 0] != null) {
                temp = (Object[])temp[pixels[i + k][J] ? 1 : 0];
                k--;
            }
            vP = ((double[])temp[2])[pixels[i][J] ? 1 : 0];
            prod = flMultiply(prod, new double[]{hP * A[Math.min(J, 9)][Math.min(i, 9)] + vP * (1 - A[Math.min(J, 9)][Math.min(i, 9)]), 0});
        }
        for (int j = J + 1; j < pixels[I].length && j < J + 10; j++) {
            double hP;
            double vP;
            Object[] temp = horizontal;
            int k = -1;
            while (j + k >= 0 && temp[pixels[I][j + k] ? 1 : 0] != null) {
                temp = (Object[])temp[pixels[I][j + k] ? 1 : 0];
                k--;
            }
            hP = ((double[])temp[2])[pixels[I][j] ? 1 : 0];
            temp = vertical;
            k = -1;
            while (I + k >= 0 && temp[pixels[I + k][j] ? 1 : 0] != null) {
                temp = (Object[])temp[pixels[I + k][j] ? 1 : 0];
                k--;
            }
            vP = ((double[])temp[2])[pixels[I][j] ? 1 : 0];
            prod = flMultiply(prod, new double[]{hP * A[Math.min(j, 9)][Math.min(I, 9)] + vP * (1 - A[Math.min(j, 9)][Math.min(I, 9)]), 0});
        }
        pixels[I][J] = temp2;
        return prod;
    }
    /*
    public static double[] getProb(double[] seed) {
        double[] p = new double[2];
        p[0] = 1;
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels.length; j++) {
                Object[] temp = horizontal;
                int k = -1;
                while (j + k > 0 && temp[pixels[i][j + k] ? 1 : 0] != null) {
                    temp = (Object[])temp[pixels[i][j + k] ? 1 : 0];
                    k--;
                }
                double hP = ((double[])temp[2])[pixels[i][j] ? 1 : 0];
                temp = vertical;
                k = -1;
                while (i + k > 0 && temp[pixels[i + k][j] ? 1 : 0] != null) {
                    temp = (Object[])temp[pixels[i + k][j] ? 1 : 0];
                    k--;
                }
                double vP = ((double[])temp[2])[pixels[i][j] ? 1 : 0];
                p = flMultiply(p, flAdd(flMultiply(seed, flAdd(new double[]{hP, 0}, new double[]{-vP, 0})), new double[]{vP, 0}));
            }
        }
        return p;
    }
    public static double getA() {
        double[] wavg = new double[2];
        double[] avg = new double[2];
        for (int i = 0; i < 20000; i++) {
            double seed = Math.random();
            wavg = flAdd(wavg, flMultiply(new double[]{seed, 0}, getProb(new double[]{seed, 0})));
            avg = flAdd(avg, getProb(new double[]{seed, 0}));
        }
        //System.out.println(Arrays.toString(wavg));
        double[] sol = flDivide(wavg, avg);
        return sol[0] * Math.pow(2, sol[1]);
    }
    */
    public static void getLookup() {
        //store horizontal and vertical probabilities in table for slight speed increase
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                double[] hP = new double[10];
                double[] vP = new double[10];
                Object[] temp = horizontal;
                int k = 0;
                while (j - k >= 0 && temp[pixels[i][j - k] ? 1 : 0] != null) {
                    hP[k] = ((double[])temp[2])[pixels[i][j] ? 1 : 0];
                    if (j - k - 1 >= 0) {
                        temp = (Object[])temp[pixels[i][j - k - 1] ? 1 : 0];
                    }
                    k++;
                }
                for (; k < hP.length; k++) {
                    hP[k] = hP[k - 1];
                }
                temp = vertical;
                k = 0;
                while (i - k >= 0 && temp[pixels[i - k][j] ? 1 : 0] != null) {
                    vP[k] = ((double[])temp[2])[pixels[i][j] ? 1 : 0];
                    if (i - k - 1 >= 0) {
                        temp = (Object[])temp[pixels[i - k - 1][j] ? 1 : 0];
                    }
                    k++;
                }
                for (; k < vP.length; k++) {
                    vP[k] = vP[k - 1];
                }
                lookup[i][j][0] = hP;
                lookup[i][j][1] = vP;
            }
        }
    }
    public static double[][][] getProb(double[] seed) {
        //calculate P(a) for a given sample value of a, returns array for entire family of Markov chains
        double[][][] sol = new double[10][10][2];
        for (int i = 0; i < sol.length; i++) {
            for (int j = 0; j < sol[i].length; j++) {
                double[] temp = {1, 0};
                sol[i][j] = temp;
            }
        }
        for (int i = 0; i < lookup.length; i++) {
            for (int j = 0; j < lookup[i].length; j++) {
                for (int k = 0; k < lookup[i][j][0].length; k++) {
                    for (int l = 0; l < lookup[i][j][1].length; l++) {
                        sol[k][l] = flMultiply(sol[k][l], flAdd(flMultiply(new double[]{lookup[i][j][0][k], 0}, seed), flMultiply(new double[]{-lookup[i][j][1][l], 0}, flAdd(new double[]{-1, 0}, seed))));
                    }
                }
            }
        }
        return sol;
    }
    public static double[][] getA() {
        //calculate constant a, takes the longest here but can be parallelized by GPU
        double[][][] wavg = new double[10][10][2];
        double[][][] avg = new double[10][10][2];
        for (int i = 0; i < 20000; i++) {
            double seed = Math.random();
            double con[][][] = getProb(new double[]{seed, 0});
            for (int j = 0; j < wavg.length; j++) {
                for (int k = 0; k < wavg[j].length; k++) {
                    wavg[j][k] = flAdd(wavg[j][k], flMultiply(new double[]{seed, 0}, con[j][k]));
                    avg[j][k] = flAdd(avg[j][k], con[j][k]);
                }
            }
        }
        double[][] sol = new double[10][10];
        for (int j = 0; j < wavg.length; j++) {
            for (int k = 0; k < wavg[j].length; k++) {
                double[] temp = flDivide(wavg[j][k], avg[j][k]);
                sol[j][k] = temp[0] * Math.pow(2, temp[1]);
            }
        }
        return sol;
    }
    public static void transpose(boolean[][] x) {
        boolean[][] buffer = new boolean[x[0].length][x.length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[i].length; j++) {
                buffer[j][i] = x[i][j];
            }
        }
        x = buffer;
    }
    public static void addToTree(int len) {
        for (int i = 0; i < pixels.length; i++) {
            for (int j = len - 1; j < pixels[i].length; j++) {
                Object[] temp = tree;
                for (int k = -1; k > -len; k--) {
                    if (temp[pixels[i][j + k] ? 1 : 0] == null) {
                        temp[pixels[i][j + k] ? 1 : 0] = new Object[3];
                    }
                    temp = (Object[])temp[pixels[i][j + k] ? 1 : 0];
                }
                if (temp[2] == null) {
                    temp[2] = new int[2];
                }
                ((int[])temp[2])[pixels[i][j] ? 1 : 0] += 1;
            }
        }
        if (len + 1 <= 10) {
            addToTree(len + 1);
        }
    }
    public static double[] normalize(double[] f) {
        if ((Double.isNaN(f[0]) || Double.isNaN(f[1])) || (Double.isInfinite(f[0]) && f[1] == Double.NEGATIVE_INFINITY) || (f[0] == 0 && f[1] == Double.POSITIVE_INFINITY)) {
            return new double[]{Double.NaN, Double.NaN};
        } else if (f[0] == 0 || f[1] == Double.NEGATIVE_INFINITY) {
            return new double[]{0, Double.NEGATIVE_INFINITY};
        } else if (f[0] == Double.POSITIVE_INFINITY) {
            return new double[]{1, Double.POSITIVE_INFINITY};
        } else if (f[0] == Double.NEGATIVE_INFINITY) {
            return new double[]{-1, Double.POSITIVE_INFINITY};
        } else if (f[1] == Double.POSITIVE_INFINITY) {
            return new double[]{1, Double.POSITIVE_INFINITY};
        }
        double[] F = f.clone();
        F[0] *= Math.pow(2, F[1] % 1);
        F[1] = (int)F[1];
        double[] temp = {F[0] * Math.pow(2, -Math.floor(Math.log(Math.abs(F[0])) / Math.log(2))), F[1] + Math.floor(Math.log(Math.abs(F[0])) / Math.log(2))};
        return temp;
    }
    public static double[] flAdd(double[] f1, double[] f2) {
        double[] F1 = normalize(f1);
        double[] F2 = normalize(f2);
        if (F1[1] > F2[1]) {
            return normalize(new double[]{F1[0] + F2[0] * Math.pow(2, F2[1] - F1[1]), F1[1]});
        } else {
            return normalize(new double[]{F2[0] + F1[0] * Math.pow(2, F1[1] - F2[1]), F2[1]});
        }
    }
    public static double[] flMultiply(double[] f1, double[] f2) {
        double[] fl1 = normalize(f1.clone());
        double[] fl2 = normalize(f2.clone());
        double[] prodf = {fl1[0] * fl2[0], fl1[1] + fl2[1]};
        return normalize(prodf);
    }
    public static double[] flDivide(double[] f1, double[] f2) {
        double[] k = f2.clone();
        k[0] = 1 / k[0];
        k[1] *= -1;
        return flMultiply(f1, k);
    }
    public static double[] flexponent(double[] a, double[] b) {
        double[] A = normalize(a);
        double[] B = normalize(b);
        double[] first = normalize(new double[]{1, Math.log(A[0]) / Math.log(2) * B[0] * Math.pow(2, B[1])});
        double[] second = normalize(new double[]{1, A[1] * B[0] * Math.pow(2, B[1])});
        return flMultiply(first, second);
    }
    public static double getC(double m, double e) {
        //calculate constant used in weighted distribution
        double ya = -m;
        double yb = 1 - m;
        double a = -1;
        double b = 1;
        int j = 0;
        while (b - a > e) {
            double x12 = (a + b) / 2;
            double xf = (b * ya - a * yb) / (ya - yb);
            double d = Math.min(0.1 * (b - a) * (b - a), Math.abs(x12 - xf));
            double s = Math.signum(x12 - xf);
            double xt = xf + d * s;
            double pk = Math.min(e * Math.pow(2, Math.ceil(-Math.log(e) / Math.log(2)) + 1 - j) - (b - a) / 2, Math.abs(xt - x12));
            double xitp = x12 - s * pk;
            double yitp = (xitp / (1 - xitp * xitp) - 1 + Math.exp(-xitp / (1 - xitp * xitp))) / (1 - Math.exp(-xitp / (1 - xitp * xitp))) / (xitp / (1 - xitp * xitp)) - m;
            if (!Double.isFinite(yitp)) {
                if (Math.abs(xitp - 1) <= Math.abs(xitp) && Math.abs(xitp - 1) <= Math.abs(xitp + 1)) {
                    yitp = 1 - m;
                } else if (Math.abs(xitp - 1) >= Math.abs(xitp) && Math.abs(xitp) <= Math.abs(xitp + 1)) {
                    yitp = 0.5 - m;
                } else {
                    yitp = -m;
                }
            }
            if (yitp > 0) {
                b = xitp;
                yb = yitp;
            } else if (yitp < 0) {
                a = xitp;
                ya = yitp;
            } else {
                a = xitp;
                b = xitp;
                ya = yitp;
                yb = yitp;
            }
            j++;
        }
        return (a + b) / 2 / (1 - (a + b) / 2 * (a + b) / 2);
    }
    public static double probability(int a, int b, int c, int d, double x, double y) {
        double[] wavg = new double[2];
        double[] avg = new double[2];
        double con = getC(x, 0.00001);
        //calculate sum and weighted sum of 20000 samples
        for (int i = 0; i < 20000; i++) {
            //get random probability from distribution
            double seed = Math.log((Math.exp(con) - 1) * Math.random() + 1) / con;
            if (!Double.isFinite(seed)) {
                seed = Math.random();
                if (con < -1) {
                    seed = x;
                }
            }
            //random number used to sample from Markov chain distribution
            double[] first = flexponent(new double[]{seed, 0}, new double[]{a, 0});
            double[] second = flexponent(flAdd(new double[]{1, 0}, new double[]{-seed, 0}), new double[]{b + c, 0});
            double[] third = flexponent(flAdd(flAdd(new double[]{y, 0}, new double[]{-x, 0}), flMultiply(new double[]{x, 0}, new double[]{seed, 0})), new double[]{d, 0});
            wavg = flAdd(wavg, flMultiply(flMultiply(flMultiply(new double[]{seed, 0}, first), second), third));
            avg = flAdd(avg, flMultiply(flMultiply(first, second), third));
        }
        //calculate weighted average
        double[] sol = flDivide(wavg, avg);
        if (Double.isNaN(sol[0] * Math.pow(2, sol[1]))) {
            return x;
        }
        return sol[0] * Math.pow(2, sol[1]);
    }
    public static void calculate(Object[] temp) {
        //this function builds the Markov chain tree
        if ((double)(((double[])temp[2])[0]) > (double)(((double[])temp[2])[1])) {
            int d = temp[0] == null ? 0 : ((int[])((Object[])temp[0])[2])[0];
            int c = temp[0] == null ? 0 : ((int[])((Object[])temp[0])[2])[1];
            int b = temp[1] == null ? 0 : ((int[])((Object[])temp[1])[2])[0];
            int a = temp[1] == null ? 0 : ((int[])((Object[])temp[1])[2])[1];
            double y = ((double[])temp[2])[0];
            double x = ((double[])temp[2])[1];
            double p = probability(a, b, c, d, x, y);
            if (temp[0] == null) {
                temp[0] = new Object[3];
                ((Object[])temp[0])[2] = new double[2];
            }
            if (temp[1] == null) {
                temp[1] = new Object[3];
                ((Object[])temp[1])[2] = new double[2];
            }
            ((Object[])temp[0])[2] = new double[2];
            double[] t = {1 - x / y * (1 - p), x / y * (1 - p)};
            ((Object[])temp[0])[2] = t;
            ((Object[])temp[1])[2] = new double[2];
            double[] s = {1 - p, p};
            ((Object[])temp[1])[2] = s;
        } else {
            int a = temp[0] == null ? 0 : ((int[])((Object[])temp[0])[2])[0];
            int b = temp[0] == null ? 0 : ((int[])((Object[])temp[0])[2])[1];
            int c = temp[1] == null ? 0 : ((int[])((Object[])temp[1])[2])[0];
            int d = temp[1] == null ? 0 : ((int[])((Object[])temp[1])[2])[1];
            double x = ((double[])temp[2])[0];
            double y = ((double[])temp[2])[1];
            double p = probability(a, b, c, d, x, y);
            if (temp[0] == null) {
                temp[0] = new Object[3];
                ((Object[])temp[0])[2] = new double[2];
            }
            if (temp[1] == null) {
                temp[1] = new Object[3];
                ((Object[])temp[1])[2] = new double[2];
            }
            ((Object[])temp[0])[2] = new double[2];
            double[] t = {p, 1 - p};
            ((Object[])temp[0])[2] = t;
            ((Object[])temp[1])[2] = new double[2];
            double[] s = {x / y * (1 - p), 1 - x / y * (1 - p)};
            ((Object[])temp[1])[2] = s;
        }
        if ((Object[])((Object[])temp[0])[0] != null || (Object[])((Object[])temp[0])[1] != null) {
            calculate((Object[])temp[0]);
        }
        if ((Object[])((Object[])temp[1])[0] != null || (Object[])((Object[])temp[1])[1] != null) {
            calculate((Object[])temp[1]);
        }
    }
}
