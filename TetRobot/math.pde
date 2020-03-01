float EPSILON = 0.000001;
float SQRT2 = sqrt(2);
float SQRT3 = sqrt(3);
float DegToRad(float deg)
{
  return deg / 180.0 * PI;
}

float RadToDef(float rad)
{
  return rad / PI * 180.0; 
}

float Lerp(float x, float s, float y) {
 return x + s * (y-x); 
}

class Vector
   { 
       float x=0,y=0,z=0; 
   // creation    
   Vector () {}; 
   Vector (float px, float py, float pz) {x = px; y = py; z = pz;};
   Vector Set (float px, float py, float pz) {x = px; y = py; z = pz; return this;};
   Vector Set (Vector V) {x = V.x; y = V.y; z = V.z; return this;}; 

   // measure
   float Magnitude() {return(sqrt(sq(x)+sq(y)+sq(z)));}; 

   //alteration
   Vector Add(Vector V) {x+=V.x; y+=V.y; z+=V.z; return this;};
   Vector Sub(Vector V) {x-=V.x; y-=V.y; z-=V.z; return this;};
   Vector Mul(float f) {x*=f; y*=f; z*=f; return this;};
   Vector Div(float f) {x/=f; y/=f; z/=f; return this;};
   Vector Normalize() {float n=Magnitude(); if (n>EPSILON) {Div(n);}; return this;};
   
   // conversion to old version
   vec ToVec() {return V(x, y, z);};
   } // end class Vector
   
   
class Point
   { 
     float x=0,y=0,z=0; 
   Point () {}; 
   Point (float px, float py, float pz) {x = px; y = py; z = pz; };
   Point Set (float px, float py, float pz) {x = px; y = py; z = pz; return this;}; 
   Point Set (Point P) {x = P.x; y = P.y; z = P.z; return this;};
   Point TranslateTowards(float s, Point P) {x+=s*(P.x-x);  y+=s*(P.y-y); z+=s*(P.z-z);  return this;};  // translate by ratio s towards point P
   Point Translate(Vector V) {x+=V.x; y+=V.y; z+=V.z; return this;};                                     // translate by vector v
   
   // conversion to old version
   pt ToPt() {return P(x,y,z);};
   }
   
// =====  vector functions
Vector Vector() {return new Vector(); };                                                                       // make vector (x,y,z)
Vector Vector(float x, float y, float z) {return new Vector(x,y,z); };                                         // make vector (x,y,z)
Vector Vector(Vector V) {return new Vector(V.x,V.y,V.z); };                                                      // make copy of vector V
Vector Vector(Point P, Point Q) {return new Vector(Q.x - P.x, Q.y - P.y, Q.z - P.z); };                          // make vector PQ (Q-P)
Vector Vector(vec V) {return new Vector(V.x,V.y,V.z); };                                                       // convert from vec
Vector Add(Vector A, Vector B) {return Vector(A.x+B.x,A.y+B.y,A.z+B.z); };                                       // A+B
Vector Add(Vector A, Vector B, Vector C) {return Add(Add(A, B), C); };                                           // A+B+C
Vector Sub(Vector U, Vector V) {return Vector(U.x-V.x,U.y-V.y,U.z-V.z);};                                        // U-V
Vector Mul(Vector A, float s) {return new Vector(s*A.x,s*A.y,s*A.z); };                                        // sA
Vector Div(Vector A, float s) {return new Vector(A.x / s,A.y / s,A.z / s); };                                           // (1/s)*A
Vector Lerp(Vector A, float s, Vector B) {return new Vector(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };          // (1-s)A+sB
Vector Normalize(Vector V) {float n = V.Magnitude(); if (n<EPSILON) return Vector(0,0,0); else return Mul(V, 1./n);};            // V / |V| unit of
Vector Cross(Vector U, Vector V) {return Vector( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                            // UxV cross product (normal to both)
float Dot(Vector U, Vector V) {return U.x*V.x + U.y*V.y + U.z*V.z;};                                                               // dot product
float Mixed(Vector U, Vector V, Vector W) {return Dot(U,Cross(V,W));}                                          // U*(VxW)  mixed product (signed scalar)
Vector GetOneNormal(Vector V) {                                                                              // some vector orthogonal to V
  if(abs(V.z)<=min(abs(V.x),abs(V.y))) return Vector(-V.y,V.x,0); 
  if(abs(V.x)<=min(abs(V.z),abs(V.y))) return Vector(0,-V.z,V.y);
  return Vector(V.z,0,-V.x);
  }
  
  
// ===== point functions

Point Point() {return new Point(); };                                                                       // make point (x,y,z)
Point Point(float x, float y, float z) {return new Point(x,y,z); };                                         // make point (x,y,z)
Point Point(Point P) {return new Point(P.x,P.y,P.z); };                                                      // make copy of point
Point Point(pt P) {return new Point(P.x,P.y,P.z); };                                                       // convert from pt
Point Average(Point A, Point B) {return new Point((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); };  // (A+B)/2
Point Average(Point A, Point B, Point C) {return new Point((A.x+B.x+C.x)/3.0,(A.y+B.y+C.y)/3.0,(A.z+B.z+C.z)/3.0); };  // (A+B+C)/3
Point Average(Point A, Point B, Point C, Point D) {return Average(Average(A,B),Average(C,D)); };                       // (A+B+C+D)/4
Point WeightedSum(float a, Point A, float b, Point B) {return Point(a*A.x+b*B.x, a*A.y+b*B.y, a*A.z+b*B.z);};                        // aA+bB 
Point WeightedSum(float a, Point A, float b, Point B, float c, Point C) {return WeightedSum(1.0,WeightedSum(a,A,b,B),c,C);};         // aA+bB+cC
Point Lerp(Point A, float s, Point B) {return new Point(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };          // (1-s)A+sB
Point Translate(Point A, Vector V) {return Point(A.x+V.x, A.y+V.y, A.z+V.z);};                                                // translate by ratio s towards point P
Point Translate(Point A, float s, Vector V) {return Point(A.x+s*V.x, A.y+s*V.y, A.z+s*V.z);};
Point TranslateTowards(Point A, float s, Point P) {return Point(A.x+s*(P.x-A.x), A.y+s*(P.y-A.y), A.z+s*(P.z-A.z));};         // translate by vector V

// ===== mouse TODO
//pt Mouse() {return P(mouseX,mouseY,0);};                                          // current mouse location
//pt Pmouse() {return P(pmouseX,pmouseY,0);};
//vec MouseDrag() {return V(mouseX-pmouseX,mouseY-pmouseY,0);};                     // vector representing recent mouse displacement
//pt ScreenCenter() {return P(width/2,height/2);}                                                        //  point in center of  canvas

// ===== measures
float MagnitudeSqr(Vector V) {return sq(V.x)+sq(V.y)+sq(V.z);};                                                   // V*V    norm squared
float Magnitude(Vector V) {return sqrt(MagnitudeSqr(V));};                                                        // V*V    norm squared
float Distance(Vector P, Vector Q){return sqrt(sq(Q.x-P.x)+sq(Q.y-P.y)+sq(Q.z-P.z)); };                             // ||AB|| distance
boolean IsParallel (Vector U, Vector V) {return Magnitude(Cross(U,V))<Magnitude(U)*Magnitude(V)*EPSILON; }          // true if U and V are almost parallel
float Angle(Vector U, Vector V) {return acos(Dot(U,V)/Magnitude(V)/Magnitude(U)); };                                // angle(U,V) positive  (in 0,PI)
boolean IsCounterClockwise(Vector U, Vector V, Vector W) {return Mixed(U,V,W)>0; };                                   // U*(VxW)>0  U,V,W are clockwise in 3D

float Distance(Point P, Point Q){return sqrt(sq(Q.x-P.x)+sq(Q.y-P.y)+sq(Q.z-P.z)); };                             // ||AB|| distance
float TriangleArea(Point A, Point B, Point C) {return Magnitude(Cross(Vector(A,B), Vector(A,C))) * 0.5;}               // (positive) area of triangle in 3D
float TetVolume(Point A, Point B, Point C, Point D) {return Mixed(Vector(A,B), Vector(A,C), Vector(A,D)) / 6.0;};         // (signed) volume of tetrahedron
boolean IsCounterClockwise(Point A, Point B, Point C, Point D) {return TetVolume(A,B,C,D)>0; };                     // tet is oriented so that A sees B, C, D clockwise 
//TODO
//boolean projectsBetween(pt P, pt A, pt B) {return dot(V(A,P),V(A,B))>0 && dot(V(B,P),V(B,A))>0 ; };
//float disToLine(pt P, pt A, pt B) {return det3(U(A,B),V(A,P)); };
//pt ProjectionOnLine(pt P, pt A, pt B) {return P(A,dot(V(A,B),V(A,P))/dot(V(A,B),V(A,B)),V(A,B));}

// ===== rotations 
Point ProjectionOnLine(Point P, Point A,Point B) {return Translate(A, Dot(Vector(A,B),Vector(A,P))/Dot(Vector(A,B),Vector(A,B)),Vector(A,B));}


Point RotatePointAroundPointInPlane(Point P, float a, Vector I, Vector J, Point G)
{
  float x=Dot(Vector(G,P), I), y=Dot(Vector(G,P),J);
  float c=cos(a), s=sin(a);
  return Translate(P, Add(Mul(I, x*c-x-y*s), Mul(J, x*s+y*c-y))); 
}; // Rotated P by a around G in plane whose orthogonal basis are (I,J)

Point RotatePointAroundAxis(Point P, float r, Point Q, Vector V) // rotate P around QV
{
  Point PP = ProjectionOnLine(P, Q, Translate(Q,V));
  
  return RotatePointAroundPointInPlane(P, r, Normalize(Vector(P,PP)), Normalize(Cross(Vector(P,PP),V)), PP);
}

Vector RotateVecInPlane(Vector V, float a, Vector I, Vector J)
{
  float x=Dot(V,I), y=Dot(V,J); 
  float c=cos(a), s=sin(a); 
  return Add(V, Mul(I, x*c-x-y*s), Mul(J, x*s+y*c-y));
}; // Rotated V by a parallel to plane (I,J) (around the norm of plane)

Vector RotateAroundOrigin(Vector V, float a)
{
  float dx=V.x, dy=V.y, c=cos(a), s=sin(a); 
  return Vector(c*dx+s*dy,-s*dx+c*dy,V.z); 
};  // Q rotated by angle a around the origin (around z axis)


Point RotateAroundOrigin(Point Q, float a)
{
  float dx=Q.x, dy=Q.y, c=cos(a), s=sin(a); 
  return Point(c*dx+s*dy,-s*dx+c*dy,Q.z); 
};  // Q rotated by angle a around the origin (around z axis)

Point RotatePointAroundPoint(Point Q, float a, Point C)
{
  float dx=Q.x-C.x, dy=Q.y-C.y, c=cos(a), s=sin(a); 
  return Point(C.x+c*dx+s*dy, C.y-s*dx+c*dy, Q.z); 
};  // Q rotated by angle a around point P (around z axis)

//TODO
/*
pt R(pt Q, pt C, pt P, pt R)  // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
   {
   vec I0=U(C,P), I1=U(C,R), V=V(C,Q); 
   float c=d(I0,I1), s=sqrt(1.-sq(c)); 
                                       if(abs(s)<0.00001) return Q; // singular cAW
   vec J0=V(1./s,I1,-c/s,I0);  
   vec J1=V(-s,I0,c,J0);  
   float x=d(V,I0), y=d(V,J0);
   return P(Q,x,M(I1,I0),y,M(J1,J0)); 
   } 
*/

// ===== rending functions
void normal(Vector V) {normal(V.x,V.y,V.z);};                                     // changes current normal vector for subsequent smooth shading
void vertex(Point P) {vertex(P.x,P.y,P.z);};                                      // vertex for shading or drawing
void vTextured(Point P, float u, float v) {vertex(P.x,P.y,P.z,u,v);};                          // vertex with texture coordinates
void drawEdge(Point P, Point Q) {line(Q.x,Q.y,Q.z,P.x,P.y,P.z); };                       // render edge (P,Q)
void drawEdge(Point P, Vector V) {line(P.x,P.y,P.z,P.x+V.x,P.y+V.y,P.z+V.z);}
void drawEdge(Point P, float d , Vector V) {line(P.x,P.y,P.z,P.x+d*V.x,P.y+d*V.y,P.z+d*V.z); };
void drawTri(Point A, Point B, Point C) {beginShape(); vertex(A);vertex(B); vertex(C); endShape(CLOSE);};                      // render Triangle(A,B,C)
void drawQuad(Point A, Point B, Point C, Point D) {beginShape(); vertex(A); vertex(B); vertex(C); vertex(D); endShape(CLOSE);};    // render Quad(A,B,C,D)
void drawSphere(Point P, float r){pushMatrix(); translate(P.x,P.y,P.z); sphere(r); popMatrix();};                          // render sphere of radius r and center P

// draw more complex shapes, we use functions in primitives.pde, maybe we can make our own version later 
void drawConeSection(Point P, Point Q, float p, float q) { // surface
  coneSection(P.ToPt(), Q.ToPt(), p, q);
  }
void drawCaplet(Point A, float a, Point B, float b) { // cone section surface that is tangent to Sphere(A,a) and to Sphere(B,b)
  caplet(A.ToPt(), a, B.ToPt(), b);
  } 
//void drawString(Vector P,
//void showCoordinate(Vector P, float s, Vector I, Vector J, Vector K) {noStroke(); fill(yellow); show(P,5); stroke(red); show(P,s,I); stroke(green); show(P,s,J); stroke(blue); show(P,s,K); }; // render coordinate system
//void show(pt P, String s) {text(s, P.x, P.y, P.z); }; // prints string s in 3D at P
//void show(pt P, String s, vec D) {text(s, P.x+D.x, P.y+D.y, P.z+D.z);  }; // prints string s in 3D at P+D (offset vector)
//void showShadow(pt P, float r) {pushMatrix(); translate(P.x,P.y,0); scale(1,1,0.01); sphere(r); popMatrix();}      // render shadow on the floot of sphere of radius r and center P

// TODO
//===== SUBDIVISION
/*
pt B(pt A, pt B, pt C, float s) {return( P(P(B,s/4.,A),0.5,P(B,s/4.,C))); };                          // returns a tucked B towards its neighbors
pt F(pt A, pt B, pt C, pt D, float s) {return( P(P(A,1.+(1.-s)/8.,B) ,0.5, P(D,1.+(1.-s)/8.,C))); };    // returns a bulged mid-edge point 

//===== PARABOLA
void drawParabolaInHat(pt A, pt B, pt C, int rec) {
   if (rec==0) { show(A,B); show(B,C); } //if (rec==0) { beam(A,B,rt); beam(B,C,rt); } 
   else { 
     float w = (d(A,B)+d(B,C))/2;
     float l = d(A,C)/2;
     float t = l/(w+l);
     pt L = L(A,t,B);
     pt R = L(C,t,B);
     pt M = P(L,R);
     drawParabolaInHat(A,L, M,rec-1); drawParabolaInHat(M,R, C,rec-1); 
     };
   };
   */
   
boolean Equal(Point A, Point B)
{
  if(Vector(A,B).Magnitude() < 1) return true;
  else return false;
}
   
//subdivision
Point[] RefineTuck(Point[] A, int N)
{
  Point[] Temp = new Point[Maxlen];
  Point[] A2 = new Point[Maxlen];
  
  for(int i = 0; i < N; i++) A2[i] = A[i];
  
  // refine
  for(int i = 0; i < N-1; i++)
  {
    Temp[2*i] = A2[i];
    Temp[2*i+1] = Average(A2[i], A2[i+1]);
  }
  Temp[2*(N-1)] = A2[N-1];
  for(int i = 1; i < 2*N-1; i++) A2[i] = Temp[i];
  
  // tuck
  for(int j = 0; j < degree / 2; j++)
  {
    Temp[0] = A2[0];
    for(int i = 1; i < 2*(N-1); i++)
    {
      Temp[i] = Average(Average(A2[i-1],A2[i]), Average(A2[i], A2[i+1]));
    }
    Temp[2*(N-1)] = A2[2*(N-1)];
    // copy
    for(int i = 1; i < 2*N-1; i++) A2[i] = Temp[i];
  }
  
  return A2;
}

Point[][] RefineTuck(Point[][] outerframe, int N)
{
  Point[][] Temp = new Point[4][Maxlen];
  Point[][] outerframe2 = new Point[4][Maxlen];
  //copy
  for(int i = 0; i < N; i++) for(int j = 0; j < 4; j++) 
    outerframe2[j][i] = outerframe[j][i];
    
  
  for(int refine = 0; refine < refineCounter; refine++)
  {
    // refine : linear
    for(int j = 0; j < 4; j++)
    {
      for(int i = 0; i < N-1; i++) 
      {
        Temp[j][2*i] = outerframe2[j][i];
        //Temp[j][2*i+1] = Average(outerframe2[j][i], outerframe2[j][i+1]); // to do: Temp[2*i+1] = S(out2[i], 0.5, out2[i+1])
        Temp[j][2*i+1] = Screw(outerframe2, j, i, 0.5, false);
      }
      Temp[j][2*(N-1)] = outerframe2[j][N-1];
    }
    // copyback
    for(int j = 0; j < 4; j++) for(int i = 1; i < 2*N-1; i++)  
      outerframe2[j][i] = Temp[j][i]; // outerframe2[j] = Temp[j]
    
    // tuck
    for(int d = 0; d < degree / 2; d++)
    {
      for(int j = 0; j < 4; j++)
      {
        Temp[j][0] = outerframe2[j][0];
        for(int i = 1; i < 2*(N-1); i++)
        {
          Point[][] TempFrame = new Point[4][2];
          for(int k = 0; k < 4; k++)
          {
            TempFrame[k][0] = Screw(outerframe2, k, i-1, 0.5, false);
            TempFrame[k][1] = Screw(outerframe2, k, i, 0.5, false);
          }
          Temp[j][i] = Screw(TempFrame, j, 0, 0.5, false);
          //Temp[j][i] = Average(Average(outerframe2[j][i-1],outerframe2[j][i]), Average(outerframe2[j][i], outerframe2[j][i+1]));// todo: apply S
          //Temp[j][i] = outerframe2[j][i];
        }
        Temp[j][2*(N-1)] = outerframe2[j][2*(N-1)];
      }
      // copy back
      for(int j = 0; j < 4; j++) for(int i = 1; i < 2*N-1; i++) 
        outerframe2[j][i] = Temp[j][i];
    }
    
    N = 2*N-1;
  }
  
  
  return outerframe2;
}

Point Screw(Point[][] A, int j, int i, float t, boolean draw) // The screw S(A[j][i], t, A[j][i+1])
{
  Vector S;
  float r, dist;
  int condition = 0;
  
  Point O0 = Average(A[0][i], A[1][i], A[2][i], A[3][i]);
  Point O1 = Average(A[0][i+1], A[1][i+1], A[2][i+1], A[3][i+1]);
  Point Q0 = Average(A[1][i], A[2][i]); // Q = (B + C) / 2
  Point Q1 = Average(A[1][i+1], A[2][i+1]);
  // I = OD, J = QA, K = QB(QC)
  Vector I = Sub(Normalize(Vector(O1,A[0][i+1])), Normalize(Vector(O0,A[0][i])));
  Vector J = Sub(Normalize(Vector(Q1,A[3][i+1])), Normalize(Vector(Q0,A[3][i])));
  Vector K = Sub(Normalize(Vector(Q1,A[1][i+1])), Normalize(Vector(Q0,A[1][i])));
  // rotations of basis sysbols
  float ii = Cross(I,J).Magnitude(), jj = Cross(J,K).Magnitude(), kk = Cross(K,I).Magnitude();
  if(ii > jj && ii > kk)
  {
    S = Normalize(Cross(I, J));
    r = Angle(Cross(S,Normalize(Vector(O0,A[0][i]))), Cross(S,Normalize(Vector(O1,A[0][i+1]))));
    condition = 1;
  }
  else if(jj > ii && jj > kk)
  {
    S = Normalize(Cross(J, K));
    r = Angle(Cross(S,Normalize(Vector(Q0,A[3][i]))), Cross(S,Normalize(Vector(Q1,A[3][i+1]))));
    condition = 2;
  }
  else
  {
    S = Normalize(Cross(K, I));
    r = Angle(Cross(S,Normalize(Vector(Q0,A[1][i]))), Cross(S,Normalize(Vector(Q1,A[1][i+1]))));
    condition = 3;
  }
  
  dist = Dot(Vector(O0,O1), S);
  
  // P = (O1+O0)/2 + (SXO)/2tan(r/2)
  Point P = Translate(Average(O0,O1), Div(Cross(Vector(O0,O1),S),2*tan(r/2)));
  Point P2 = Translate(Average(O0,O1), Div(Cross(S,Vector(O0,O1)),2*tan(r/2)));
  // End 1 and 2?
  Point End = RotatePointAroundAxis(Translate(A[j][i], Mul(S, dist)), r, P, S);
  Point End2 = RotatePointAroundAxis(Translate(A[j][i], Mul(Mul(S,-1), dist)), r, P2, Mul(S,-1));
  
  //if(!Equal(End, A[j][i+1])) {S = Mul(S, -1); P = P2; condition = 4;}
  if(Vector(End, A[j][i+1]).Magnitude() > Vector(End2, A[j][i+1]).Magnitude()) {S = Mul(S, -1); P = P2; condition = 4;}
  
  End = RotatePointAroundAxis(Translate(A[j][i], Mul(S, dist)), r, P, S);
  //print(Vector(End, A[j][i+1]).Magnitude(), "\n");
  if(!Equal(End, A[j][i+1])) print("Screw Motion Failed with magnitude =", Vector(End, A[j][i+1]).Magnitude(), "\n");
  
  if(draw) // test
  {
    print(condition, dist, "\n");
    print(i, "\n");
    //t = 1;
    
    //if(End.x == A[j][i+1].x && End.y == A[j][i+1].y && End.z == A[j][i+1].z) print("YES\n");
    //else print("NO\n");
    //fill(white); drawSphere(O0, 10); 
    //fill(yellow); drawSphere(O1, 10); 
    //fill(red); for(int aa = 0; aa < 4; aa++) drawSphere(A[aa][i], 10); 
    //fill(blue); for(int aa = 0; aa < 4; aa++) drawSphere(A[aa][i+1], 10);
    //fill(black); drawSphere(End, 10); 
    //fill(green); for(int aa = 0; aa < 4; aa++) drawSphere(RotatePointAroundAxis(Translate(A[aa][i], Mul(S, t*dist)), t*r, P, S), 10);
    
    fill(red); drawSphere(A[j][i], 10);
    fill(blue); drawSphere(End, 10);
    fill(green); drawSphere(RotatePointAroundAxis(Translate(A[j][i], Mul(S, t*dist)), t*r, P, S), 10);
    
    stroke(black); drawEdge(P, 1000, S); drawEdge(P, 1000, Mul(S,-1));
  }
  
  return RotatePointAroundAxis(Translate(A[j][i], Mul(S, t*dist)), t*r, P, S);
}
