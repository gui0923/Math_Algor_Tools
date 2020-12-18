package com.guilei.gis.coord;


import com.guilei.gis.bean.Angle;
import com.guilei.gis.bean.Position;

/**
 * 坐标转换
 * 国测局坐标，又名火星坐标，国内的地图强制使用的就是火星坐标系，
 * 其中常用的腾讯搜搜地图，阿里云地图，高德地图使用的就是火星坐标系
 *
 * WGS84坐标系，世界标准坐标系，World Wind, Google以及美国的其他地图，使用的就是这个
 *
 * 百度坐标系，百度在国测局坐标系的基础之上进行的二次加密
 *
 * @author 桂磊
 *
 */
public class CoordTransFormUtils
{
    private static final double X_PI = Math.PI * 3000.0 / 180.0;
    private static final double PI = Math.PI;
    private static final double A = 6378245.0;
    private static final double EE = 0.00669342162296594323;// # 扁率

    public static Position gcj02tobd09(Position gcjPosition)
    {
        // 国测局坐标转换成百度坐标
        double lon = gcjPosition.getLongitude().degrees;
        double lat = gcjPosition.getLatitude().degrees;
        double z = Math.sqrt(lon * lon + lat * lat) + 0.00002
                * Math.sin(lat * X_PI);
        double theta = Math.atan2(lat, lon) + 0.000003 * Math.cos(lon * X_PI);
        double bd_lon = z * Math.cos(theta) + 0.0065;
        double bd_lat = z * Math.sin(theta) + 0.006;
        // return new Position(Position.fromDegrees(bd_lat, bd_lon));
        return new Position(Angle.fromDegrees(bd_lat),
                Angle.fromDegrees(bd_lon), 0);
    }

    public static Position bd09togcj02(Position bdPosition)
    {
        // 百度坐标往国测局坐标转换
        double bd_lon = bdPosition.getLongitude().degrees;
        double bd_lat = bdPosition.getLatitude().degrees;
        double x = bd_lon - 0.0065;
        double y = bd_lat - 0.006;
        double z = Math.sqrt(x * x + y * y) - 0.00002 * Math.sin(y * X_PI);
        double theta = Math.atan2(y, x) - 0.000003 * Math.cos(x * X_PI);
        double gg_lon = z * Math.cos(theta);
        double gg_lat = z * Math.sin(theta);
        // return new Position(Position.fromDegrees(gg_lat, gg_lon));
        return new Position(Angle.fromDegrees(gg_lat),
                Angle.fromDegrees(gg_lon), 0);
    }

    public static Position wgs84togcj02(Position wgsPosition)
    {
        //WGS84坐标转换成国测局坐标
        double lng = wgsPosition.getLongitude().degrees;
        double lat = wgsPosition.getLatitude().degrees;
        double mgLat;
        double mgLon;
        if (outOfChina(lat, lng))
        {
            mgLat = lat;
            mgLon = lng;
            return wgsPosition;
        }
        double dlat = transformLat(lng - 105.0, lat - 35.0);
        double dlng = transformLon(lng - 105.0, lat - 35.0);
        double radlat = lat / 180.0 * PI;
        double magic = Math.sin(radlat);
        magic = 1 - EE * magic * magic;
        double sqrtmagic = Math.sqrt(magic);
        dlat = (dlat * 180.0) / ((A * (1 - EE)) / (magic * sqrtmagic) * PI);
        dlng = (dlng * 180.0) / (A / sqrtmagic * Math.cos(radlat) * PI);
        mgLat = lat + dlat;
        mgLon = lng + dlng;
        // return new Position(Position.fromDegrees(mgLat, mgLon));
        return new Position(Angle.fromDegrees(mgLat), Angle.fromDegrees(mgLon),
                0);
    }

    public static Position gcj02towgs84(Position gcjPosition)
    {
        //国测局坐标向WGS84坐标进行转换
        double lng = gcjPosition.getLongitude().degrees;
        double lat = gcjPosition.getLatitude().degrees;
        if (outOfChina(lat, lng))
        {
            return gcjPosition;
        }
        double dlat = transformLat(lng - 105.0, lat - 35.0);
        double dlng = transformLon(lng - 105.0, lat - 35.0);
        double radlat = lat / 180.0 * PI;
        double magic = Math.sin(radlat);
        magic = 1 - EE * magic * magic;
        double sqrtmagic = Math.sqrt(magic);
        dlat = (dlat * 180.0) / ((A * (1 - EE)) / (magic * sqrtmagic) * PI);
        dlng = (dlng * 180.0) / (A / sqrtmagic * Math.cos(radlat) * PI);
        double mglat = lat + dlat;
        double mglng = lng + dlng;
        // return new Position(Position.fromDegrees(lng * 2 - mglng, lat * 2 -
        // mglat));
        return new Position(Angle.fromDegrees(lat * 2 - mglat),
                Angle.fromDegrees(lng * 2 - mglng),
                0);
    }

    public static double[] gcj02towgs84(double lat, double lon)
    {
        //国测局坐标向WGS84坐标进行转换
        if (outOfChina(lat, lon))
        {
            return new double[]{lat, lon};
        }
        double dlat = transformLat(lon - 105.0, lat - 35.0);
        double dlng = transformLon(lon - 105.0, lat - 35.0);
        double radlat = lat / 180.0 * PI;
        double magic = Math.sin(radlat);
        magic = 1 - EE * magic * magic;
        double sqrtmagic = Math.sqrt(magic);
        dlat = (dlat * 180.0) / ((A * (1 - EE)) / (magic * sqrtmagic) * PI);
        dlng = (dlng * 180.0) / (A / sqrtmagic * Math.cos(radlat) * PI);
        double mglat = lat + dlat;
        double mglng = lon + dlng;
        // return new Position(Position.fromDegrees(lng * 2 - mglng, lat * 2 -
        // mglat));
        return new double[]{lat * 2 - mglat, lon * 2 - mglng};
    }

    public static Position transform(Position wgPosition)
    {
        // wgs84 转 国测局
        double wgLat = wgPosition.getLatitude().degrees;
        double wgLon = wgPosition.getLongitude().degrees;
        double mgLat;
        double mgLon;
        if (outOfChina(wgLat, wgLon))
        {
            mgLat = wgLat;
            mgLon = wgLon;
            return wgPosition;
        }
        double dLat = transformLat(wgLon - 105.0, wgLat - 35.0);
        double dLon = transformLon(wgLon - 105.0, wgLat - 35.0);
        double radLat = wgLat / 180.0 * PI;
        double magic = Math.sin(radLat);
        magic = 1 - EE * magic * magic;
        double sqrtMagic = Math.sqrt(magic);
        dLat = (dLat * 180.0) / ((A * (1 - EE)) / (magic * sqrtMagic) * PI);
        dLon = (dLon * 180.0) / (A / sqrtMagic * Math.cos(radLat) * PI);
        mgLat = wgLat + dLat;
        mgLon = wgLon + dLon;
        // return new Position(Position.fromDegrees(mgLat, mgLon));
        return new Position(Angle.fromDegrees(mgLat), Angle.fromDegrees(mgLon),
                0);
    }

    public static boolean outOfChina(double lat, double lon)
    {
        //判断是否在国内，如果在国内才需要进行转换，在国外就不需要了，不过话说，高德地图和百度地图也不支持国外的
        if (lon < 72.004 || lon > 137.8347)
            return true;
        if (lat < 0.8293 || lat > 55.8271)
            return true;
        return false;
    }

    public static double transformLat(double x, double y)
    {
        //转换
        double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y
                + 0.2 * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(y * PI) + 40.0 * Math.sin(y / 3.0 * PI)) * 2.0 / 3.0;
        ret += (160.0 * Math.sin(y / 12.0 * PI) + 320 * Math.sin(y * PI / 30.0)) * 2.0 / 3.0;
        return ret;
    }

    public static double transformLon(double x, double y)
    {
        double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1
                * Math.sqrt(Math.abs(x));
        ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0;
        ret += (20.0 * Math.sin(x * PI) + 40.0 * Math.sin(x / 3.0 * PI)) * 2.0 / 3.0;
        ret += (150.0 * Math.sin(x / 12.0 * PI) + 300.0 * Math.sin(x / 30.0
                * PI)) * 2.0 / 3.0;
        return ret;
    }

    public static Position wgs84tobd09(Position position)
    {
        //WGS84坐标转百度坐标
        return gcj02tobd09(wgs84togcj02(position));
    }

    public static Position bd09towgs84(Position position)
    {
        //百度坐标转WGS84坐标
        return gcj02towgs84(bd09togcj02(position));
    }

    public static void main(String[] args)
    {
        Position[] loPositions = new Position[2];
        loPositions[0] = // new Position(Position.fromDegrees(29.575429778924,
                // 114.21892734521));
                new Position(Angle.fromDegrees(29), Angle.fromDegrees(114), 0);
        loPositions[1] = // new Position(Position.fromDegrees(29.575429778924,
                // 114.21892734521));
                new Position(Angle.fromDegrees(30), Angle.fromDegrees(115), 0);
        for (Position position : loPositions)
        {
            System.out.println(gcj02tobd09(wgs84togcj02(position)));
        }

        loPositions[0] = // new Position(Position.fromDegrees(29.575429778924,
                // 114.21892734521));
                new Position(Angle.fromDegrees(29.00335000762543), Angle.fromDegrees(114.01186648118536), 0);
        loPositions[1] = // new Position(Position.fromDegrees(29.575429778924,
                // 114.21892734521));
                new Position(Angle.fromDegrees(30.003217737932722), Angle.fromDegrees(115.01180010728888), 0);
        for (Position position : loPositions)
        {
            System.out.println(gcj02towgs84(bd09togcj02(position)));
        }
    }
}
