#include "Vertex3D.h"

Vertex3D::Vertex3D() {
    this->x=0;
    this->y=0;
    this->z=0;
    this->vtstar=-1;
}

Vertex3D::Vertex3D(const Vertex3D& orig) {
    this->x=orig.getX();
    this->y=orig.getY();
    this->z=orig.getZ();
    this->field=orig.getF();
    this->vtstar = orig.VTstar();
}

Vertex3D::~Vertex3D() {
}

Vertex3D::Vertex3D(float x, float y, float z){
   this->x = x;
   this->y = y;
   this->z = z;
    this->field = 0;
   this->vtstar = -1;
}

void Vertex3D::setF(float f){
    this->field = f;
}

float Vertex3D::getF() const{
    return this->field;
}

void Vertex3D::setZ(float z){
    this->z = z;
}

float Vertex3D::getZ() const{
    return this->z;
}

bool operator== (const Vertex3D &p, const Vertex3D &q) {
        return ((p.x == q.x) && (p.y == q.y) && (p.z == q.z));
}

bool operator !=(const Vertex3D& p, const Vertex3D& q) {
        return !(p == q);
}

void Vertex3D::operator +=(const Vertex3D &p){
    x += p.getX();
    y += p.getY();
    z += p.getZ();
    field += p.getF();
}

void Vertex3D::operator /=(const float &f){
    x /= f;
    y /= f;
    z /= f;
    field /= f;
}

bool Vertex3D::operator <(const Vertex3D &p) const{

    if(this->x != p.getX())
        return this->x < p.getX();
    if(this->y != p.getY())
        return this->y < p.getY();
    if(this->z != p.getZ())
        return this->z < p.getZ();
    if(this->field != p.getF())
        return this->field < p.getF();

    return false;
}


float Vertex3D::getX() const{
    return this->x;
}

float Vertex3D::getY() const{
    return this->y;
}


void Vertex3D::setX(float x){
    this->x = x;
}

void Vertex3D::setY(float y){
    this->y = y;
}

int Vertex3D::VTstar() const
{
    return this->vtstar;
}

void Vertex3D::VTstar(int vtstar)
{
    this->vtstar = vtstar;
}
