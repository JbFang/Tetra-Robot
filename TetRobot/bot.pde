Vector Up = Vector(0,0,1);
TETROBOT BOT = new TETROBOT();
TETROBOT GHOST = new TETROBOT();

float getRotateAngle(float t)
{
  return (PI - acos(1.0/3.0)) * (1-t);
}

class TETROBOT // class for manipulaitng TETROBOTs
  {
  Point[] hips = new Point[4];     // hip positions
  Point[] outerFrame = new Point[4];
  Point[] feet = new Point[4];
  Point[] knees = new Point[4];
  Vector[] ikDir = new Vector[4];
  Point A,B,C,D, Anext, Bnext, Cnext, Dnext, OldA, OldB, OldC, OldD;
  float l1, l2, h, H, h0;
  float ellipse_a, ellipse_b;
  float footMaxHeight;
  float edgeLength;
  float bodyLength;
  float bodyRatio;
  float currentTime=0;
  color[] col = new color[4];
  int a=0, b=1, c=2, d=3; // I use these to remember which hip is which
  int old_a=0, old_b=1, old_c=2, old_d=3;
  Boolean 
    showTet=false,
    showTiles=false,
    showBody=true,
    showLegs=false,
    showArrow=true,
    isRolling=false;
    
  Point initialPosCenter = Point();
  Point initialPosFront = Point(); // point A
  
  IntList historyCommands = new IntList();
  int currentCommand = -1;
   
  TETROBOT() {}

  void resetBot(float r) // puts the bot in the center of the terrain
    {
    Point X = Point(); // center of circle
    A = Translate(X, Vector(-r,0,0));
    B = RotatePointAroundPoint(A, DegToRad(120), X);
    C = RotatePointAroundPoint(B, DegToRad(120), X);
    edgeLength = Distance(A,B);
    bodyLength = edgeLength * bodyRatio;
    
    outerFrame[0].Set(A);
    outerFrame[1].Set(B);
    outerFrame[2].Set(C);
    feet[0]=Point(A);
    feet[1]=Point(B);
    feet[2]=Point(C);
    Point O = Average(A,B,C); // centroid of triangle
    
    //ikDir[0] = Vector(O,A).Normalize();
    //ikDir[1] = Vector(O,B).Normalize();
    //ikDir[2] = Vector(O,C).Normalize();
    //ikDir[3] = Vector(O,A).Normalize();
    
    initialPosCenter = Point(O);
    initialPosFront = Point(A);
    
    
    D = Point(O).Translate(Mul(Up, edgeLength*sqrt(2./3)));
    
    h0 = edgeLength * sqrt(6)/12 + bodyLength * sqrt(6)/4;
    l1 = h0*1.5;
    l2 = h0;
    
    feet[3] = Translate(A, Mul(Vector(O,Average(B,C)), 6));
    
    h = edgeLength * sqrt(6) / 3;
    H = h0 + l1 + l2 -h;
    
    ellipse_a = edgeLength * 2 * SQRT3 / 3;
    ellipse_b = h0+l1+l2;
    
    
    
    placeBot(A,B,C);
    outerFrame[3].Set(D);
    feet[3].Set(D);
    
    
    a=0; b=1; c=2; d=3;
    historyCommands.clear();
    currentTime = 0;
    currentCommand = -1;
    
    OldA = A;
    OldB = B;
    OldC = C;
    OldD = D;
    }

  void declareBot(float ratio) 
    {
    for (int i=0; i<4; i++) outerFrame[i]=Point();
    
    col[0]=orchid; // STUDENTS: Use your personalized colors
    col[1]=blue;
    col[2]=red;
    col[3]=gold;
    
    bodyRatio = ratio;
    }
    
   int genCmdFromPoint(Point M)
   {
     Point O = Average(A,B,C);
     M.z = O.z;
     float left, right, op;
     Vector OM = Vector(O,M);
     left = Dot(OM, Normalize(Vector(O, Average(A,B))));
     right = Dot(OM, Normalize(Vector(O, Average(A,C))));
     op = Dot(OM, Normalize(Vector(O, Average(B,C))));
     float edge = Distance(A,B);
     if (left > right && left > op && left > 0.5*edge)
     {
       return L_CMD;
     }
     else if (right > left && right > op && right > 0.5*edge)
     {
       return R_CMD;
     }
     else if (op > left && op > right && op > 0.5*edge)
     {
       return O_CMD;
     }
     return -1;
   }


   void placeBot(Point AA, Point BB, Point CC) // puts the bot tet so that one of its facs is the triangle ABC
    {
      Point O = Average(AA,BB,CC); // centroid of triangle
      float e=(Distance(AA,BB)+Distance(BB,CC)+Distance(CC,AA))/3; // average edge length
      Point DD = Point(O).Translate(Mul(Up, e*sqrt(2./3)));
      
      Anext=Point(AA);
      Bnext=Point(BB);
      Cnext=Point(CC);
      Dnext=Point(DD);
    }
    
    void finishRoll()
    {
      A = Anext;
      B = Bnext;
      C = Cnext;
      D = Dnext;
    }
    
    void rollBot(float t)
    {
      currentTime = t;
      Point Center = Point();
      
      if (currentCommand <= O_CMD) {
        Point Q = Average(Bnext, Cnext);
        Point AA,BB,CC,DD;
        
        feet[b].Set(outerFrame[b]);
        feet[c].Set(outerFrame[c]);
        Center = Average(outerFrame[a],outerFrame[b],outerFrame[c],outerFrame[d]);
        // Screw Motion
        if(ScrewMotion){
          
          if(InterpolationMethod1)
          {
            AA = RotatePointAroundPointInPlane(Anext, getRotateAngle(0), Normalize(Vector(Q,Anext)), Up, Q);
            DD = RotatePointAroundPointInPlane(Dnext, getRotateAngle(0), Normalize(Vector(Q,Anext)), Up, Q);
            Point O0 = Average(AA,Bnext,Cnext,DD), O1 = Average(Anext,Bnext,Cnext,Dnext);
            // I: OD(OAA), J: QA(QDD), K: 0
            Vector I = Sub(Normalize(Vector(O1,Anext)), Normalize(Vector(O0,AA)));
            Vector J = Sub(Normalize(Vector(Q,Dnext)), Normalize(Vector(Q,DD)));
            Vector S = Normalize(Cross(I, J));
            float r = Angle(Cross(S,Normalize(Vector(O0,AA))), Cross(S,Normalize(Vector(O1,Anext))));
            float dist = Dot(Vector(O0,O1), S);
            // P = (O1+O0)/2 + (SXO)/2tan(r/2)
            Point P = Translate(Average(O0,O1), Div(Cross(Vector(O0,O1),S),2*tan(r/2)));
            if(showscrew)
            {
              fill(green); drawSphere(P, 10);
              stroke(black); drawEdge(P, 1000, S); drawEdge(P, 1000, Mul(S,-1));
            }
            
            outerFrame[a] = RotatePointAroundAxis(Translate(AA, Mul(S, t*dist)), t*r, P, S);
            outerFrame[b] = RotatePointAroundAxis(Translate(Bnext, Mul(S, t*dist)), t*r, P, S);
            outerFrame[c] = RotatePointAroundAxis(Translate(Cnext, Mul(S, t*dist)), t*r, P, S);
            outerFrame[d] = RotatePointAroundAxis(Translate(DD, Mul(S, t*dist)), t*r, P, S);
          }
          else 
          {
            // Old Frame: AA OldB OldC DD  new Frame: nAA Bnext Cnext nDD
            // relationship: DD <-> nAA
            
            AA = RotatePointAroundPointInPlane(OldA, getRotateAngle(0.5), Normalize(Vector(Average(OldB,OldC),OldA)), Up, Average(OldB,OldC));
            DD = RotatePointAroundPointInPlane(OldD, getRotateAngle(0.5), Normalize(Vector(Average(OldB,OldC),OldA)), Up, Average(OldB,OldC));
            Point nAA = RotatePointAroundPointInPlane(Anext, getRotateAngle(0.5), Normalize(Vector(Q,Anext)), Up, Q);
            Point nDD = RotatePointAroundPointInPlane(Dnext, getRotateAngle(0.5), Normalize(Vector(Q,Anext)), Up, Q);
            Point O0 = Average(AA,OldB,OldC,DD), O1 = Average(nAA,Bnext,Cnext,nDD);
            // I: ODD(OnAA) J: QAA -> ?
            Vector I = Sub(Normalize(Vector(O1,nAA)), Normalize(Vector(O0,DD)));
            Vector J = Vector();
            if(currentCommand == L_CMD) // AA -> Cnext  OldB -> Bnext OldC -> nDD
              J = Sub(Normalize(Vector(Average(Bnext,nDD),Cnext)), Normalize(Vector(Average(OldB,OldC),AA)));
            else//if(currentCommand == R_CMD) // AA -> Bnext OldB -> nDD OldC -> Cnext
              J = Sub(Normalize(Vector(Average(Cnext,nDD),Bnext)), Normalize(Vector(Average(OldB,OldC),AA)));
              
            Vector S = Normalize(Cross(J, I));
            float r = Angle(Cross(S,Normalize(Vector(O0,DD))), Cross(S,Normalize(Vector(O1,nAA))));
            float dist = Dot(Vector(O0,O1), S);
            Point P = Translate(Average(O0,O1), Div(Cross(Vector(O0,O1),S),2*tan(r/2)));
            if(showscrew)
            {
              fill(green); drawSphere(P, 10);
              stroke(black); drawEdge(P, 1000, S); drawEdge(P, 1000, Mul(S,-1));
            }
            
            outerFrame[a] = RotatePointAroundAxis(Translate(DD, Mul(S, t*dist)), t*r, P, S);
            if(currentCommand == L_CMD) // AA -> Cnext  OldB -> Bnext OldC -> nDD
            {
              outerFrame[c] = RotatePointAroundAxis(Translate(AA, Mul(S, t*dist)), t*r, P, S);
              outerFrame[b] = RotatePointAroundAxis(Translate(OldB, Mul(S, t*dist)), t*r, P, S);
              outerFrame[d] = RotatePointAroundAxis(Translate(OldC, Mul(S, t*dist)), t*r, P, S);
            }
            else //if(currentCommand == R_CMD) // AA -> Bnext OldB -> nDD OldC -> Cnext
            {
              outerFrame[b] = RotatePointAroundAxis(Translate(AA, Mul(S, t*dist)), t*r, P, S);
              outerFrame[d] = RotatePointAroundAxis(Translate(OldB, Mul(S, t*dist)), t*r, P, S);
              outerFrame[c] = RotatePointAroundAxis(Translate(OldC, Mul(S, t*dist)), t*r, P, S);
            }
          }
          
        }
        else
        {
          outerFrame[a] = RotatePointAroundPointInPlane(Anext, getRotateAngle(t), Normalize(Vector(Q,Anext)), Up, Q); ///////
          outerFrame[d] = RotatePointAroundPointInPlane(Dnext, getRotateAngle(t), Normalize(Vector(Q,Anext)), Up, Q); ////// x*Qd
          outerFrame[b] = Bnext;
          outerFrame[c] = Cnext;
        }
        
        if(t > 0 && t <= 0.5){ // landing
          feet[a] = RotatePointAroundPointInPlane(Anext, getRotateAngle(2*t), Normalize(Vector(Q,Anext)), Up, Q);
          feet[a] = Translate(feet[a], Mul(Normalize(Vector(Center,feet[a])), H*(1-2*t)));
          //feet[a].z = feet[a].z*(h0+l1+l2)/h;
        }
        if(t > 0.5 && t < 1){ // lifting
          feet[d] = RotatePointAroundPointInPlane(Dnext, getRotateAngle(2*t-1), Normalize(Vector(Q,Anext)), Up, Q);
          //feet[d].z = feet[d].z*(h0+l1+l2)/h;
          feet[d] = Translate(feet[d], Mul(Normalize(Vector(Center,feet[d])), H*(2*t-1)));
        }
      }
      
     
      else {
       // dancing mode
       if (currentCommand == L_DANCE) {
         //print('d');
         outerFrame[c].Set(Cnext);
         outerFrame[a] = RotatePointAroundPoint(Anext, (1.0-t) * (PI / 3), Cnext);
         outerFrame[b] = RotatePointAroundPoint(Bnext, (1.0-t) * (PI / 3), Cnext);
         outerFrame[d] = RotatePointAroundPoint(Dnext, (1.0-t) * (PI / 3), Cnext);
         Center = Average(outerFrame[a],outerFrame[b],outerFrame[c],outerFrame[d]);
         
         feet[c].Set(outerFrame[c]);
         float t1 = 2.0 / 3.0;
         float t2 = 1.0 / 3.0;
         if (t > 0 && t <= t1) {
          feet[a].Set(Lerp(Bnext, t/t1, Anext));
          feet[a].z = -500 * t * (t-t1);
         }
         else {
           feet[a].Set(Anext);
           Point P = Average(Bnext, Cnext);
           Point Last = P.Translate(Vector(Anext, P));
           feet[b].Set(Lerp(Last, (t-t1)/t2, Bnext));
           feet[b].z = -500 * (t-t1) * (t-1);
         }
         //if (t > t2 && t <= 1) {
           
       }
       else if (currentCommand == R_DANCE) {
         outerFrame[b].Set(Bnext);
         outerFrame[a] = RotatePointAroundPoint(Anext, -(1.0-t) * (PI / 3), Bnext);
         outerFrame[c] = RotatePointAroundPoint(Cnext, -(1.0-t) * (PI / 3), Bnext);
         outerFrame[d] = RotatePointAroundPoint(Dnext, -(1.0-t) * (PI / 3), Bnext);
         Center = Average(outerFrame[a],outerFrame[b],outerFrame[c],outerFrame[d]);
         
         feet[b].Set(outerFrame[b]);
         float t1 = 2.0 / 3.0;
         float t2 = 1.0 / 3.0;
         if (t > 0 && t <= t1) {
           feet[a].Set(Lerp(Cnext, t/t1, Anext));
            feet[a].z = -500 * t * (t-t1);
         }
         else {
           feet[a].Set(Anext);
           Point P = Average(Bnext, Cnext);
           Point Last = P.Translate(Vector(Anext, P));
           feet[c].Set(Lerp(Last, (t-t1)/t2, Cnext));
           feet[c].z = -500 * (t-t1) * (t-1);
         }
         
         
       }
       else if (currentCommand == T_DANCE) {
         Point O = Average(Anext,Bnext,Cnext);
         if (t<0.5) t = sin(PI*t)*0.5;
          else t = sin(PI*(t-0.5))*0.5+0.5;
         if (t < 1.0/3) {
           outerFrame[a] = RotatePointAroundPoint(Anext, (1.0-1.5*t) * (2*PI/3), O);
           outerFrame[b] = RotatePointAroundPoint(Bnext, (1.0-1.5*t) * (2*PI/3), O);
           outerFrame[c] = RotatePointAroundPoint(Cnext, (1.0-1.5*t) * (2*PI/3), O);
           
         }
         else if (t > 2.0/3) {
           outerFrame[a] = RotatePointAroundPoint(Anext, (1.5-1.5*t) * (2*PI/3), O);
           outerFrame[b] = RotatePointAroundPoint(Bnext, (1.5-1.5*t) * (2*PI/3), O);
           outerFrame[c] = RotatePointAroundPoint(Cnext, (1.5-1.5*t) * (2*PI/3), O);
         }
         outerFrame[d].Set(D);
         Center = Average(outerFrame[a],outerFrame[b],outerFrame[c],outerFrame[d]);
         
         Point A1, B1, C1;
         A1 = RotatePointAroundPoint(Anext, PI/3, O);
         B1 = RotatePointAroundPoint(Bnext, PI/3, O);
         C1 = RotatePointAroundPoint(Cnext, PI/3, O);
         
         float t1 = 1.0/6;
         float h = 2000;
         if (t > 0 && t <= t1) {
           feet[a].Set(Lerp(Bnext, t/t1, A1));
           feet[a].z = -h * t * (t-t1);
         }
         else if (t <= 2*t1) {
           feet[b].Set(Lerp(Cnext, (t-t1)/t1, B1));
           feet[b].z = -h * (t-t1) * (t-t1*2);
         }
         else if (t <= 3*t1) {
           feet[c].Set(Lerp(Anext, (t-2*t1) / t1, C1));
           feet[c].z = -h * (t-t1*2) * (t-t1*3);
         }
         else if (t <= 4*t1) {
           feet[a].Set(Lerp(A1, (t-3*t1) / t1, Anext));
           feet[a].z = -h * (t-t1*3) * (t-t1*4);
         }
         else if (t <= 5*t1) {
           feet[b].Set(Lerp(B1, (t-4*t1) / t1, Bnext));
           feet[b].z = -h * (t-t1*4) * (t-t1*5);
         }
         else {
           feet[c].Set(Lerp(C1, (t-5*t1) / t1, Cnext));
           feet[c].z = -h * (t-t1*5) * (t-1);
         }
         
       }
      }
      
      for (int i = 0; i<4; i++)
        hips[i] = WeightedSum(1.0-bodyRatio, Center, bodyRatio, outerFrame[i]);
        
      if (currentCommand > 2) {
        for (int i = 0; i<4; i++) {
          hips[i].Translate(Mul(Up, 20*cos(2*PI*t)));
          hips[i].Translate(Mul(Vector(sin(4*PI*t), cos(4*PI*t), 0), 10*cos(2*PI*t)));
        }
        
         Point Top = Translate(hips[d], Mul(Up, l1+l2));
         feet[d].Set(Top);
         feet[d].z *= 0.8;
         float r = 50;
         feet[d].x = Center.x + r * sin(2*PI*t);
         feet[d].y = Center.y + r * cos(2*PI*t);
      }
      
    }
    
    void setCommand(int cmd)
    {
      if (cmd < 0 || cmd >= MAX_CMD) return;
      currentCommand = cmd;
      historyCommands.append(cmd);
      int new_a=a, new_b=b, new_c=c, new_d=d;
      
      old_a = a; old_b = b; old_c = c; old_d = d;
      OldA = A; OldB = B; OldC = C; OldD = D;
      
      if (cmd == L_CMD)
      {
        new_a = d;
        new_b = b;
        new_c = a;
        new_d = c;
        
        Point P = Average(A,B);
        Vector dir = Vector(C,P);
        C = A;
        A = P.Translate(dir);
        placeBot(A,B,C);
      }
      else if (cmd == R_CMD)
      {
        new_a = d;
        new_b = a;
        new_c = c;
        new_d = b;
        
        Point P = Average(A,C);
        Vector dir = Vector(B,P);
        B = A;
        A = P.Translate(dir);
        placeBot(A,B,C);
      }
      else if (cmd == O_CMD)
      {
        new_a = d;
        new_b = c;
        new_c = b;
        new_d = a;
        
        Point P = Average(B,C);
        Vector dir = Vector(A,P);
        A = P.Translate(dir);
        placeBot(A,C,B);
        
        B.Set(outerFrame[1]);
        C.Set(outerFrame[2]);
      }
      else if (cmd == L_DANCE)
      {
        print (cmd);
        new_a = b;
        new_b = c;
        new_c = a;
        new_d = d;
        
        Point P = Average(A,B);
        Vector dir = Vector(C,P);
        C.Set(A);
        A = P.Translate(dir);
        placeBot(A,B,C);
      }
      else if (cmd == R_DANCE) 
      {
        print (cmd);
        new_a = c;
        new_b = a;
        new_c = b;
        new_d = d;
        
        Point P = Average(A,C);
        Vector dir = Vector(B,P);
        B.Set(A);
        A = P.Translate(dir);
        placeBot(A,B,C);
      }
      else if (cmd == T_DANCE) 
      {
        print (cmd);
        //new_a = c;
        //new_b = a;
        //new_c = b;
        new_a = a;
        new_b = b;
        new_c = c;
        new_d = d;
        placeBot(C,A,B);
        A.Set(Anext);
        B.Set(Bnext);
        C.Set(Cnext);
      }
      a = new_a;
      b = new_b;
      c = new_c;
      d = new_d;
      return;
    }
    
    void drawHistoryTiles(float w, float r, color edgeColor)
    {
      if (!showTiles) return;
        int i = 0;
        Point A, B, C;
        A = Point(initialPosFront);
        B = RotatePointAroundPoint(initialPosFront, DegToRad(120), initialPosCenter);
        C = RotatePointAroundPoint(initialPosFront, DegToRad(240), initialPosCenter);
        int idx_a = 0, idx_b = 1, idx_c = 2, idx_d = 3;
        int new_a, new_b, new_c, new_d;
        fill(col[0]); drawSphere(A, r);
        fill(col[1]); drawSphere(B, r);
        fill(col[2]); drawSphere(C, r);
        fill(edgeColor); drawCaplet(A, w, B, w); drawCaplet(C, w, B, w); drawCaplet(A, w, C, w);
        Point new_A;
        while (historyCommands.size() > i)
        {
          int cmd = historyCommands.get(i);
          if (cmd == L_CMD)
          {
            new_a = idx_d;
            new_b = idx_b;
            new_c = idx_a;
            new_d = idx_c;
            
            Point P = Average(A,B);
            Vector dir = Vector(C,P);
            C.Set(A);
            
            new_A = P.Translate(dir);
            fill(edgeColor);
            drawCaplet(new_A, w, B, w);
            drawCaplet(new_A, w, A, w);
            A.Set(new_A);
            
          }
          else if (cmd == R_CMD)
          {
            new_a = idx_d;
            new_b = idx_a;
            new_c = idx_c;
            new_d = idx_b;
            
            Point P = Average(A,C);
            Vector dir = Vector(B,P);
            B.Set(A);
            
            new_A = P.Translate(dir);
            fill(edgeColor);
            drawCaplet(new_A, w, A, w);
            drawCaplet(new_A, w, C, w);
            A.Set(new_A);
          }
          else // O_CMD
          {
            new_a = idx_d;
            new_b = idx_c;
            new_c = idx_b;
            new_d = idx_a;
            
            Point P = Average(B,C);
            Vector dir = Vector(A,P);
            
            new_A = P.Translate(dir);
            fill(edgeColor);
            drawCaplet(new_A, w, B, w);
            drawCaplet(new_A, w, C, w);
            A.Set(new_A);
            
            P.Set(C);
            C.Set(B);
            B.Set(P);
          }
          
          fill(col[new_a]); drawSphere(A,r);
          
          idx_a = new_a;
          idx_b = new_b;
          idx_c = new_c;
          idx_d = new_d;
          
          i++;
        }
    }
    


 Point getCentroid() {return Average(outerFrame[0],outerFrame[1],outerFrame[2],outerFrame[3]);} // average of the 4 vertices
 
 Point FindKnee(Point Hip, Point Foot, float l1, float l2)
  {
    float d = Distance(Hip, Foot);
    if (1.0001*d > (l1+l2))
    {
      return WeightedSum(l2/d, Hip, l1/d, Foot);
    }
    else
    {
      Vector AB = Vector(Hip,Foot);
      AB.Normalize();
      float x = (d*d+l1*l1-l2*l2)/d/2;
      Point M = Translate(Hip, Mul(AB, x));
      Point Center = Average(outerFrame[0],outerFrame[1],outerFrame[2],outerFrame[3]);
      Vector n = Cross(Up, Vector(Center, Foot));
      Vector disp = Cross(n,AB);
      disp.Normalize();
      if (disp.z < 0) disp.Mul(-1);
       Point Knee = M.Translate(disp.Mul(sqrt(l1*l1 - x*x)));
     return Knee;
    }
  }
 
 void drawBot (float w, float r)
    {
    noStroke();
    
    if (showTet)
    {
    for(int i=0; i<4; i++) {fill(col[i]);  drawSphere(outerFrame[i], r);} // shows tet vertices as balls
      fill(green); for (int v=0; v<3; v++)  for (int u=v+1; u<4; u++) drawCaplet(outerFrame[v],w,outerFrame[u],w); // shows tet edges as cylinders
    }
      
      if (showBody) {
          for(int i=0; i<4; i++) {fill(col[i]);  drawSphere(hips[i],r*0.5); }//drawSphere(feet[i], r);}
          
        fill(col[d]); drawTri(hips[a], hips[b], hips[c]);
        fill(col[a]); drawTri(hips[d], hips[b], hips[c]);
        fill(col[b]); drawTri(hips[d], hips[c], hips[a]);
        fill(col[c]); drawTri(hips[d], hips[a], hips[b]);
      }
      
      if (showLegs && showBody)
      {
        fill(grey);
        for (int i=0;i<4; i++) {
          knees[i] = FindKnee(hips[i], feet[i], l1, l2);
          drawSphere(knees[i],r);
          fill(grey);
          drawCaplet(knees[i],w,hips[i],w);
          drawCaplet(knees[i],w,feet[i],w);
        }
        
      }
      
      //Point X = Point(hips[1]).Translate(Vector(hips[1], Average(hips[0], hips[2])).Mul(0.2));
      //Point Y = Point(hips[2]).Translate(Vector(hips[2], Average(hips[0], hips[1])).Mul(0.2));
      //pt X = P(G[b],0.2,V(G[b],P(G[a],G[c])));    
      //pt Y = P(G[c],0.2,V(G[c],P(G[b],G[a])));     
      //fill(red); if(showArrow) arrow(X.ToPt(),Y.ToPt(),5); // shows arrow in the wedge of the rotation axis
    }

     
  } // END TETROBOT CLASS

 
