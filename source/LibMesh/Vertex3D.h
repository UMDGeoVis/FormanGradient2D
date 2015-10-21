#ifndef _VERTEX3D_H
#define	_VERTEX3D_H

#include <math.h>
///An inner-class, extending Vertex2D, representing a vertex in a tetrahedral mesh
class Vertex3D
{
public:
    ///A constructor method
    Vertex3D();
    ///A constructor method
    Vertex3D(const Vertex3D& orig);
    ///A constructor method
    Vertex3D(float x, float y, float z);
    ///A destructor method
    virtual ~Vertex3D();


    ///Getter/Setter for coordinates (x,y,z) and for the scalar value
    float getX() const;
    float getY() const;
    float getZ() const;
    float getF() const;
    //
    void setX(float x);
    void setY(float y);
    void setZ(float z);
    void setF(float f);

    ///Getter/Setter for partial VT* relation
    int VTstar() const;
    void VTstar(int vtstar);


    ///operators
    friend bool operator== (const Vertex3D& p, const Vertex3D &q);
    friend bool operator!= (const Vertex3D& p, const Vertex3D &q);
    void operator+= (const Vertex3D& p);
    void operator/= (const float& val);
    bool operator<(const Vertex3D& p) const;


    ///Following functions do work on x,y,z coordinates only (scalar field is ignored)
    //this function returns the norm of vector vec-v
    float norma(Vertex3D& v){return(sqrt(((v.getX()-x)*(v.getX()-x))+((v.getY()-y)*(v.getY()-y))+((v.getZ()-z)*(v.getZ()-z))));}

    //this function returns the scalar products between vectors v1-vec and v2-vec
    float prodscal(Vertex3D& v1,Vertex3D& v2){return(((v1.getX()-x)*(v2.getX()-x))+((v1.getY()-y)*(v2.getY()-y))+((v1.getZ()-z)*(v2.getZ()-z)));}

    float distance(Vertex3D& v1){return sqrt(pow((v1.getX())-(x),2) + pow((v1.getY())-(y),2) + pow((v1.getZ())-(z),2));}

    //norm of this vector
    float norma() {return sqrt((x)*(x)+(y)*(y)+(z)*(z));}

    //scalar products between vectors vec and v1
    float prodscal(Vertex3D& v1){return(((v1.getX())*(x))+((v1.getY())*(y))+((v1.getZ())*(z)));}

private:
    ///A protected variable representing the x coordinate of the point
    float x;
    ///A protected variable representing the y coordinate of the point
    float y;
    ///A protected variable representing the z coordinate of the point
    float z;
    ///A protected variable representing the scalar value
    float field;


    int vtstar;
};

#endif	/* _VERTEX3D_H */

