package com.guilei.coordswitch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * 专门针对于2纬坐标转换，3纬的不适宜
 *
 * @author GuiLei
 *
 */
public class CommonMatrixSwitch
{

    public static double[] calculateMatrixMulThread(final double[] f1, final double[] f2,
                                           final double[] f3, final double[] s1, final double[] s2, final double[] s3)
    {
        final double lenthS2 = distance(s1, s2);
        ThreadPoolExecutor executor = new ThreadPoolExecutor(4, 4, 60,
                TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
        double minDeviation = Double.MAX_VALUE;
        double rotateS2 = 0;
        double rotateS3 = 0;

        final double[] tmp2 = new double[3];
        final List<double[]> results = new ArrayList<double[]>();

        final int loop = 14400;
        final double interval = 40;
        for (int i = 0; i < loop; i++)
        {
            final int index = i;

            executor.execute(new Runnable()
            {
                @Override
                public void run()
                {
                    double minDeviation = Double.MAX_VALUE;
                    double rotateS2 = 0;
                    double rotateS3 = 0;

                    double tmpAngleS2;
                    double tmpAngleS3;
                    double tmpDevation;
                    tmpAngleS2 = index / interval;
                    tmp2[0] = s1[0] + lenthS2 * Math.cos(Math.toRadians(tmpAngleS2));
                    tmp2[1] = s1[1] + lenthS2 * Math.sin(Math.toRadians(tmpAngleS2));
                    tmp2[2] = 1;
                    for (int j = 0; j < loop; j++)
                    {
                        tmpAngleS3 = j / interval;
                        tmpDevation = getDeviation(f1, f2, f3, s1, tmp2, s2, s3,
                                tmpAngleS3);
                        if (minDeviation > tmpDevation)
                        {
                            minDeviation = tmpDevation;
                            rotateS2 = tmpAngleS2;
                            rotateS3 = tmpAngleS3;
                        }
                    }
                    results.add(new double[]{rotateS2, rotateS3, minDeviation});
                }
            });
            do
            {
                try
                {
                    executor.awaitTermination(2, TimeUnit.MILLISECONDS);
                }
                catch (InterruptedException e)
                {
                    e.printStackTrace();
                }
            }
            while (executor.getActiveCount() >= 4);
        }
        executor.shutdown();
        try
        {
            boolean end = true;
            do
            {
                end = !executor.awaitTermination(2, TimeUnit.MILLISECONDS);
            }
            while (end);
        }
        catch (InterruptedException e)
        {
            e.printStackTrace();
        }
        for (double[] ds : results)
        {
            if(minDeviation > ds[2])
            {
                minDeviation = ds[2];
                rotateS2 = ds[0];
                rotateS3 = ds[1];
            }
        }

        tmp2[0] = s1[0] + lenthS2 * Math.cos(Math.toRadians(rotateS2));
        tmp2[1] = s1[1] + lenthS2 * Math.sin(Math.toRadians(rotateS2));
        tmp2[2] = 1;

        double[] matrix = getMatrix(f1, f2, s1, tmp2, rotateS3);
        return matrix;
    }
    /**
     *
     * @param f1
     * @param f2
     * @param f3
     * @param s1
     * @param s2
     * @param s3
     * @return
     */
    public static double[] calculateMatrix(double[] f1, double[] f2,
                                               double[] f3, double[] s1, double[] s2, double[] s3)
    {
        double lenthS2 = distance(s1, s2);

        double minDeviation = Double.MAX_VALUE;
        double rotateS2 = 0;
        double rotateS3 = 0;

        double tmpAngleS2;
        double tmpAngleS3;
        double tmpDevation;
        double[] tmp2 = new double[3];
        for (int i = 0; i < 7200; i++)
        {
            tmpAngleS2 = i * 1d / 20d;
            tmp2[0] = s1[0] + lenthS2 * Math.cos(Math.toRadians(tmpAngleS2));
            tmp2[1] = s1[1] + lenthS2 * Math.sin(Math.toRadians(tmpAngleS2));
            tmp2[2] = 1;
            for (int j = 0; j < 7200; j++)
            {
                tmpAngleS3 = j * 1d / 20d;
                tmpDevation = getDeviation(f1, f2, f3, s1, tmp2, s2, s3,
                        tmpAngleS3);
                if (minDeviation > tmpDevation)
                {
                    minDeviation = tmpDevation;
                    rotateS2 = tmpAngleS2;
                    rotateS3 = tmpAngleS3;
                }
            }
        }
        tmp2[0] = s1[0] + lenthS2 * Math.cos(Math.toRadians(rotateS2));
        tmp2[1] = s1[1] + lenthS2 * Math.sin(Math.toRadians(rotateS2));
        tmp2[2] = 1;
        double[] matrix = getMatrix(f1, f2, s1, tmp2, rotateS3);
        return matrix;
    }

    public static double getDeviation(double[] ff1, double[] ff2, double[] ff3,
                                      double[] ss1, double[] tmp2, double[] ss2, double[] ss3,
                                      double rotate)
    {
        double[] f1 = copy(ff1);
        double[] f2 = copy(ff2);
        double[] s1 = copy(ss1);
        double[] s2 = copy(tmp2);
        double[] s3 = copy(ss3);
        double[] scale;
        double[] translate;
        // 第1步，缩放
        double sx = (s2[0] - s1[0]) / ((f2[0] - f1[0]) == 0 ? 1: (f2[0] - f1[0]));
        double sy = (s2[1] - s1[1]) / ((f2[1] - f1[1]) == 0 ? 1: (f2[1] - f1[1]));

        scale = new double[] { sx, 0, 0, 0, sy, 0, 0, 0, 1 };
        f1[0] *= sx;
        f2[0] *= sx;
        f1[1] *= sy;
        f2[1] *= sy;

        double tx = s1[0] - f1[0];
        double ty = s1[1] - f1[1];

        // 第2步，平移
        translate = new double[] { 1, 0, 0, 0, 1, 0, tx, ty, 1 };
        f1[0] = f1[0] + tx;
        f2[0] = f2[0] + tx;
        f1[1] = f1[1] + ty;
        f2[1] = f2[1] + ty;
        double[] result = multMatrix(scale, translate);

        // 以定位点为中心，回归到0点
        double[] toZero = new double[] { 1, 0, 0, 0, 1, 0, -s1[0], -s1[1], 1 };

        f1[0] += -s1[0];
        f2[0] += -s1[0];
        f1[1] += -s1[1];
        f2[1] += -s1[1];

        // 新坐标系整体沿着逆时针方向旋转rotate度
        double angle = Math.toRadians(rotate);
        double[] rotate2 = new double[] { Math.cos(angle), Math.sin(angle), 0,
                -Math.sin(angle), Math.cos(angle), 0, 0, 0, 1 };
        f1 = matrixMultPoint3d(rotate2, f1);
        f2 = matrixMultPoint3d(rotate2, f2);
        f1[0] += s1[0];
        f2[0] += s1[0];
        f1[1] += s1[1];
        f2[1] += s1[1];

        double[] translate2 = new double[] { 1, 0, 0, 0, 1, 0, s1[0], s1[1], 1 };

        result = multMatrix(result, toZero);
        result = multMatrix(result, rotate2);
        result = multMatrix(result, translate2);
        double[] dr2 = matrixMultPoint3d(result, ff2);
        double[] dr3 = matrixMultPoint3d(result, ff3);
        double newDistanceS2 = distance(dr2, ss2);
        double newDistanceS3 = distance(dr3, s3);
        return newDistanceS2 + newDistanceS3;
    }

    public static double[] getMatrix(double[] ff1, double[] ff2, double[] s1,
                                     double[] tmp2, double rotate)
    {
        double[] f1 = copy(ff1);
        double[] f2 = copy(ff2);
        double[] s2 = copy(tmp2);
        double[] scale;
        double[] translate;
        // 第1步，缩放
        double sx = (s2[0] - s1[0]) / ((f2[0] - f1[0]) == 0 ? 1 : (f2[0] - f1[0]));
        double sy = (s2[1] - s1[1]) / ((f2[1] - f1[1]) == 0 ? 1 : (f2[1] - f1[1]));

        scale = new double[] { sx, 0, 0, 0, sy, 0, 0, 0, 1 };
        f1[0] *= sx;
        f2[0] *= sx;
        f1[1] *= sy;
        f2[1] *= sy;

        double tx = s1[0] - f1[0];
        double ty = s1[1] - f1[1];

        // 第2步，平移
        translate = new double[] { 1, 0, 0, 0, 1, 0, tx, ty, 1 };
        f1[0] = f1[0] + tx;
        f2[0] = f2[0] + tx;
        f1[1] = f1[1] + ty;
        f2[1] = f2[1] + ty;

        double[] result = multMatrix(scale, translate);

        // 以定位点为中心，回归到0点
        double[] toZero = new double[] { 1, 0, 0, 0, 1, 0, -s1[0], -s1[1], 1 };

        f1[0] += -s1[0];
        f2[0] += -s1[0];
        f1[1] += -s1[1];
        f2[1] += -s1[1];

        // 新坐标系整体沿着逆时针方向旋转rotate度
        double angle = Math.toRadians(rotate);// 66.66431416570411
        double[] rotate2 = new double[] { Math.cos(angle), Math.sin(angle), 0,
                -Math.sin(angle), Math.cos(angle), 0, 0, 0, 1 };
        f1 = matrixMultPoint3d(rotate2, f1);
        f2 = matrixMultPoint3d(rotate2, f2);
        f1[0] += s1[0];
        f2[0] += s1[0];
        f1[1] += s1[1];
        f2[1] += s1[1];

        // 重置当初的矩阵
        double[] translate2 = new double[] { 1, 0, 0, 0, 1, 0, s1[0], s1[1], 1 };

        // 最终获取转换矩阵
        result = multMatrix(result, toZero);
        result = multMatrix(result, rotate2);
        result = multMatrix(result, translate2);
        return result;
    }

    private static double[] copy(double[] f1)
    {
        double[] f2 = new double[f1.length];
        System.arraycopy(f1, 0, f2, 0, f1.length);
        return f2;
    }

    private static double distance(double[] p1, double[] p2)
    {
        double value = 0;
        for (int i = 0; i < p2.length; i++)
        {
            value += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        }
        return Math.sqrt(value);
    }

    public static double[] add(double[] p1, double[] p2)
    {
        double[] result = new double[p1.length];
        for (int i = 0; i < result.length; i++)
        {
            result[i] = p1[i] + p2[i];

        }
        return result;
    }

    public static double[] sub(double[] p1, double[] p2)
    {
        double[] result = new double[p1.length];
        for (int i = 0; i < result.length; i++)
        {
            result[i] = p1[i] - p2[i];

        }
        return result;
    }

    public static double[] matrixMultPoint3d(double[] m, double[] v)
    {
        double[] res = new double[] { v[0] * m[0] + v[1] * m[3] + v[2] * m[6],
                v[0] * m[1] + v[1] * m[4] + v[2] * m[7],
                v[0] * m[2] + v[1] * m[5] + v[2] * m[8] };
        return res;
    }

    public static double[] multMatrix(double[] m1, double[] m2)
    {
        return new double[] { m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6],
                m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7],
                m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8],

                m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6],
                m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7],
                m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8],

                m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6],
                m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7],
                m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8] };
    }

    public static void main(String[] args) {
        double[] A = {0, 0, 1};
        double[] B = {0, 1, 1};
        double[] C = {1, 1, 1};
        double[] D = {1, 0, 1};

        double[] AN = {5, 5, 1};
        double[] BN = {6.414, 6.414, 1};
        double[] CN = {7.828,   5, 1};
        double[] DN = {6.414, 3.586, 1};

        double[] matrix = calculateMatrix( A, B, C, AN, BN, CN);
        System.out.println(Arrays.toString(matrix));
        double[] DP = matrixMultPoint3d(matrix, D);
        System.out.println(Arrays.toString(DP));
    }
}
