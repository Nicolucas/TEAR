
#ifndef __Lib_SDFgeometry_h__
#define __Lib_SDFgeometry_h__

    struct GeometryParams {
    double angle;
    double radius; 
    };

    void horizontal_sdf(double coor[], struct GeometryParams GeoParamList, double *phi);
    void horizontal_grad_sdf(double coor[], struct GeometryParams GeoParamList, double grad[]);
    void tilted_sdf(double coor[], struct GeometryParams GeoParamList, double *phi);
    void tilted_grad_sdf(double coor[], struct GeometryParams GeoParamList, double grad[]);
    void (*sdf_func[])(double coor[], struct GeometryParams GeoParamList, double *phi);
    void (*sdf_grad_func[])(double coor[], struct GeometryParams GeoParamList, double grad[]);
    void PrintCoordSDFvalues(double coor[], double grad[], double phi);
    void UnitTest_SDFgeometry();
#endif