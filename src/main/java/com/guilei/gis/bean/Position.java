package com.guilei.gis.bean;

public class Position extends LatLon
{
    public static final Position ZERO = new Position(Angle.ZERO, Angle.ZERO, 0d);

    public double elevation;

    public Position()
    {
//		System.out.println("Position.Position()");
//		return;
    }

    public Position(Angle latitude, Angle longitude, double elevation)
    {
        super(latitude, longitude);
        this.elevation = elevation;
    }

    public Position(LatLon latLon, double elevation)
    {
        super(latLon);
        this.elevation = elevation;
    }

    public static Position fromRadians(double latitude, double longitude, double elevation)
    {
        return new Position().setRadians(latitude, longitude, elevation);
    }

    public Position setRadians(double latitude, double longitude, double elevation)
    {
        return setDegrees(Angle.RADIANS_TO_DEGREES * latitude, Angle.RADIANS_TO_DEGREES * longitude,
                elevation);
    }

    public static Position fromDegrees(double latitude, double longitude, double elevation)
    {
        return new Position().setDegrees(latitude, longitude, elevation);
    }

    public Position setDegrees(double latitude, double longitude, double elevation)
    {
        set(latitude, longitude);
        setElevation(elevation);
        return this;
    }

    public static Position fromDegrees(double latitude, double longitude)
    {
        return new Position().setDegrees(latitude, longitude);
    }

    public Position setDegrees(double latitude, double longitude)
    {
        return set(Angle.fromDegrees(latitude), Angle.fromDegrees(longitude), 0);
    }

    public Position set(Angle latitude, Angle longitude, double elevation)
    {
        super.set(latitude, longitude);
        this.elevation = elevation;
        return this;
    }

    public Position set(LatLon latLon, double elevation)
    {
        super.set(latLon);
        this.elevation = elevation;
        return this;
    }

    /**
     * Obtains the elevation of this position
     * @return this position's elevation
     */
    public double getElevation()
    {
        return this.elevation;
    }

    /**
     * Obtains the elevation of this position
     * @return this position's elevation
     */
    public double getAltitude()
    {
        return this.elevation;
    }

    public Position add(Position that)
    {
        Angle lat = Angle.normalizedLatitude(this.latitude.add(that.latitude));
        Angle lon = Angle.normalizedLongitude(this.longitude.add(that.longitude));

        return new Position(lat, lon, this.elevation + that.elevation);
    }

    public Position subtract(Position that)
    {
        Angle lat = Angle.normalizedLatitude(this.latitude.subtract(that.latitude));
        Angle lon = Angle.normalizedLongitude(this.longitude.subtract(that.longitude));

        return new Position(lat, lon, this.elevation - that.elevation);
    }


    @Override
    public boolean equals(Object o)
    {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        Position position = (Position) o;

        // noinspection RedundantIfStatement
        if (Double.compare(position.elevation, elevation) != 0) return false;

        return true;
    }

    @Override
    public int hashCode()
    {
        int result = super.hashCode();
        long temp;
        temp = elevation != +0.0d ? Double.doubleToLongBits(elevation) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    public String toString()
    {
        return "(" + this.latitude.toString() + ", " + this.longitude.toString() + ", " + this.elevation
                + ")";
    }

    public void setElevation(double elevation)
    {
        this.elevation = elevation;
    }

    public void set(double latitude, double longitude, double elevation)
    {
        super.set(latitude, longitude);
        this.elevation = elevation;
    }

}
