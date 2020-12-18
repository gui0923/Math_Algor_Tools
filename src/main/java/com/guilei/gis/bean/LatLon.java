package com.guilei.gis.bean;


public class LatLon
{
    public static final LatLon ZERO = new LatLon(Angle.ZERO, Angle.ZERO);

    public Angle latitude = new Angle();
    public Angle longitude = new Angle();

    public LatLon()
    {
//		System.out.println("LatLon.LatLon()");
    }

    public LatLon(double latitude, double longitude)
    {
        set(latitude, longitude);
    }

    /**
     * Constructs a new <code>LatLon</code> from two angles. Neither angle may be null.
     * @param latitude latitude
     * @param longitude longitude
     * @throws IllegalArgumentException if <code>latitude</code> or <code>longitude</code> is null
     */
    public LatLon(Angle latitude, Angle longitude)
    {
        set(latitude, longitude);
    }

    public LatLon(LatLon latLon)
    {
        set(latLon);
    }

    /**
     * Factor method for obtaining a new <code>LatLon</code> from two angles expressed in radians.
     * @param latitude in radians
     * @param longitude in radians
     * @return a new <code>LatLon</code> from the given angles, which are expressed as radians
     */
    public static LatLon fromRadians(double latitude, double longitude)
    {
        return new LatLon(Math.toDegrees(latitude), Math.toDegrees(longitude));
    }

    /**
     * Factory method for obtaining a new <code>LatLon</code> from two angles expressed in degrees.
     * @param latitude in degrees
     * @param longitude in degrees
     * @return a new <code>LatLon</code> from the given angles, which are expressed as degrees
     */
    public static LatLon fromDegrees(double latitude, double longitude)
    {
        return new LatLon(latitude, longitude);
    }

    public LatLon set(double latitude, double longitude)
    {
        this.latitude.setDegrees(latitude);
        this.longitude.setDegrees(longitude);
        return this;
    }

    public LatLon set(Angle latitude, Angle longitude)
    {
        this.latitude.set(latitude);
        this.longitude.set(longitude);
        return this;
    }

    public LatLon set(LatLon latLon)
    {
        this.latitude = latLon.latitude;
        this.longitude = latLon.longitude;
        return this;
    }

    /**
     * Obtains the latitude of this <code>LatLon</code>.
     * @return this <code>LatLon</code>'s latitude
     */
    public final Angle getLatitude()
    {
        return this.latitude;
    }

    /**
     * Obtains the longitude of this <code>LatLon</code>.
     * @return this <code>LatLon</code>'s longitude
     */
    public final Angle getLongitude()
    {
        return this.longitude;
    }

    /**
     * Returns an array of this object's latitude and longitude in degrees.
     * @return the array of latitude and longitude, arranged in that order.
     */
    public double[] asDegreesArray()
    {
        return new double[] { this.getLatitude().degrees, this.getLongitude().degrees };
    }

    /**
     * Returns an array of this object's latitude and longitude in radians.
     * @return the array of latitude and longitude, arranged in that order.
     */
    public double[] asRadiansArray()
    {
        return new double[] { this.getLatitude().radians, this.getLongitude().radians };
    }


    @Override
    public String toString()
    {
        String las = String.format("Lat %7.4f\u00B0", this.getLatitude().getDegrees());
        String los = String.format("Lon %7.4f\u00B0", this.getLongitude().getDegrees());
        return "(" + las + ", " + los + ")";
    }

    @Override
    public boolean equals(Object o)
    {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LatLon latLon = (LatLon) o;

        if (!latitude.equals(latLon.latitude)) return false;
        // noinspection RedundantIfStatement
        if (!longitude.equals(latLon.longitude)) return false;

        return true;
    }

    public static boolean equals(LatLon a, LatLon b)
    {
        return a.getLatitude().equals(b.getLatitude()) && a.getLongitude().equals(b.getLongitude());
    }

    @Override
    public int hashCode()
    {
        int result;
        result = latitude.hashCode();
        result = 29 * result + longitude.hashCode();
        return result;
    }

    /**
     * Compute the forward azimuth between two positions
     * @param p1 first position
     * @param p2 second position
     * @param equatorialRadius the equatorial radius of the globe in meters
     * @param polarRadius the polar radius of the globe in meters
     * @return the azimuth
     */
    public static Angle ellipsoidalForwardAzimuth(LatLon p1, LatLon p2, double equatorialRadius,
                                                  double polarRadius)
    {
        // TODO: What if polar radius is larger than equatorial radius?
        // Calculate flattening
        final double f = (equatorialRadius - polarRadius) / equatorialRadius; // flattening

        // Calculate reduced latitudes and related sines/cosines
        final double U1 = Math.atan((1.0 - f) * Math.tan(p1.latitude.radians));
        final double cU1 = Math.cos(U1);
        final double sU1 = Math.sin(U1);

        final double U2 = Math.atan((1.0 - f) * Math.tan(p2.latitude.radians));
        final double cU2 = Math.cos(U2);
        final double sU2 = Math.sin(U2);

        // Calculate difference in longitude
        final double L = p2.longitude.subtract(p1.longitude).radians;

        // Vincenty's Formula for Forward Azimuth
        // iterate until change in lambda is negligible (e.g. 1e-12 ~= 0.06mm)
        // first approximation
        double lambda = L;
        double sLambda = Math.sin(lambda);
        double cLambda = Math.cos(lambda);

        // dummy value to ensure
        double lambda_prev = Double.MAX_VALUE;
        int count = 0;
        while (Math.abs(lambda - lambda_prev) > 1e-12 && count++ < 100)
        {
            // Store old lambda
            lambda_prev = lambda;
            // Calculate new lambda
            double sSigma = Math
                    .sqrt(Math.pow(cU2 * sLambda, 2) + Math.pow(cU1 * sU2 - sU1 * cU2 * cLambda, 2));
            double cSigma = sU1 * sU2 + cU1 * cU2 * cLambda;
            double sigma = Math.atan2(sSigma, cSigma);
            double sAlpha = cU1 * cU2 * sLambda / sSigma;
            double cAlpha2 = 1 - sAlpha * sAlpha; // trig identity
            // As cAlpha2 approaches zeros, set cSigmam2 to zero to converge on a solution
            double cSigmam2;
            if (Math.abs(cAlpha2) < 1e-6)
            {
                cSigmam2 = 0;
            } else
            {
                cSigmam2 = cSigma - 2 * sU1 * sU2 / cAlpha2;
            }
            double c = f / 16 * cAlpha2 * (4 + f * (4 - 3 * cAlpha2));

            lambda = L + (1 - c) * f * sAlpha
                    * (sigma + c * sSigma * (cSigmam2 + c * cSigma * (-1 + 2 * cSigmam2)));
            sLambda = Math.sin(lambda);
            cLambda = Math.cos(lambda);
        }

        return Angle.fromRadians(Math.atan2(cU2 * sLambda, cU1 * sU2 - sU1 * cU2 * cLambda));
    }

    public static double ellipsoidalDistance(LatLon p1, LatLon p2, double equatorialRadius,
                                             double polarRadius)
    {
        // TODO: I think there is a non-iterative way to calculate the distance. Find it and compare with this one.
        // TODO: What if polar radius is larger than equatorial radius?
        final double F = (equatorialRadius - polarRadius) / equatorialRadius; // flattening = 1.0 / 298.257223563;
        final double R = 1.0 - F;
        final double EPS = 0.5E-13;
        double GLAT1 = p1.getLatitude().radians;
        double GLAT2 = p2.getLatitude().radians;
        double TU1 = R * Math.sin(GLAT1) / Math.cos(GLAT1);
        double TU2 = R * Math.sin(GLAT2) / Math.cos(GLAT2);
        double CU1 = 1. / Math.sqrt(TU1 * TU1 + 1.);
        double SU1 = CU1 * TU1;
        double CU2 = 1. / Math.sqrt(TU2 * TU2 + 1.);
        double S = CU1 * CU2;
        double BAZ = S * TU2;
        double FAZ = BAZ * TU1;
        double GLON1 = p1.getLongitude().radians;
        double GLON2 = p2.getLongitude().radians;
        double X = GLON2 - GLON1;
        double D, SX, CX, SY, CY, Y, SA, C2A, CZ, E, C;
        do
        {
            SX = Math.sin(X);
            CX = Math.cos(X);
            TU1 = CU2 * SX;
            TU2 = BAZ - SU1 * CU2 * CX;
            SY = Math.sqrt(TU1 * TU1 + TU2 * TU2);
            CY = S * CX + FAZ;
            Y = Math.atan2(SY, CY);
            SA = S * SX / SY;
            C2A = -SA * SA + 1.;
            CZ = FAZ + FAZ;
            if (C2A > 0.)
            {
                CZ = -CZ / C2A + CY;
            }
            E = CZ * CZ * 2. - 1.;
            C = ((-3. * C2A + 4.) * F + 4.) * C2A * F / 16.;
            D = X;
            X = ((E * CY * C + CZ) * SY * C + Y) * SA;
            X = (1. - C) * X * F + GLON2 - GLON1;
            // IF(DABS(D-X).GT.EPS) GO TO 100
        } while (Math.abs(D - X) > EPS);

        // FAZ = Math.atan2(TU1, TU2);
        // BAZ = Math.atan2(CU1 * SX, BAZ * CX - SU1 * CU2) + Math.PI;
        X = Math.sqrt((1. / R / R - 1.) * C2A + 1.) + 1.;
        X = (X - 2.) / X;
        C = 1. - X;
        C = (X * X / 4. + 1.) / C;
        D = (0.375 * X * X - 1.) * X;
        X = E * CY;
        S = 1. - E - E;
        S = ((((SY * SY * 4. - 3.) * S * CZ * D / 6. - X) * D / 4. + CZ) * SY * D + Y) * C
                * equatorialRadius * R;

        return S;
    }

    public static int computeHash(Angle latitude, Angle longitude)
    {
        int result;
        result = latitude.hashCode();
        result = 29 * result + longitude.hashCode();
        return result;
    }
}
