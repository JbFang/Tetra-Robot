void keyPressed() 
  {
  if(key=='!') snapPicture();
  if(key=='~') filming=!filming;
  if(key=='?') scribeText=!scribeText;
  if(key=='v') tracking=!tracking;
  //if(key=='-') smooth=!smooth; 
  if(key=='+') {degree = 5; degreeChange = true;}
  if(key=='-') {degree = 3; degreeChange = true;}
  if(key=='>') maxf/=2;
  if(key=='<') maxf*=2;
  if(key=='.') cmdIdx++;
  if(key==',') cmdIdx=max(0,cmdIdx-1);
  if(key=='r') {BOT.resetBot(rBase); S=""; cmdIdx=0; f=0; t=0; }
  if(key=='i') {InterpolationMethod1 = !InterpolationMethod1;}
  if(key=='l') {BOT.showLegs=!BOT.showLegs;}
  if(key=='#') {BOT.showTiles=!BOT.showTiles;}
  if(key=='b') {BOT.showBody=!BOT.showBody;}
  if(key=='s') {ScrewMotion = !ScrewMotion;}
  if(key=='t') {showTrajectory = !showTrajectory;}
  if(key=='d') {subdivision = !subdivision; showscrew = !showscrew;}
  if(key=='c') {showscrew = !showscrew;}
  if(key=='^') {BOT.showTet=!BOT.showTet;}
  //if(key=='A') showArrow=!showArrow;
  //if(key=='g') {rolling=true; animating=true; f=0; t=0; BOT.resetBot(rBase);} //BOT.RollBot();
  //if(key=='a') {animating=!animating; }// toggle animation
  
  if(key=='L') {S=S+"L";if(cmdIdx==S.length()-1)f=maxf;}
  if(key=='R') {S=S+"R";if(cmdIdx==S.length()-1)f=maxf;}
  if(key=='O') {S=S+"O";if(cmdIdx==S.length()-1)f=maxf;}
  //if(key=='Z') {S=S.substring(0,S.length()-1); cmdIdx=S.length();}
  if(key=='C') {S=getClipboard(); BOT.resetBot(rBase); cmdIdx=0; f=0; t=0; }
  if(key=='S') setClipboard(S);
  if(key=='*') {currentMode = 1-currentMode;BOT.resetBot(rBase); S=""; cmdIdx=0; f=0; t=0; };
  }

void mouseWheel(MouseEvent event) 
  {
  dz -= event.getAmount(); 
  change=true;
  }

void mousePressed() 
  {
  change=true;
  }
  
void mouseMoved() 
  {
  //if (!keyPressed) 
  if (keyPressed && key==' ') {rx-=PI*(mouseY-pmouseY)/height; ry+=PI*(mouseX-pmouseX)/width;};
  if (keyPressed && key=='`') dz+=(float)(mouseY-pmouseY); // approach view (same as wheel)
  change=true;
  }
  
void mouseDragged() 
  {
  if (keyPressed && key=='t')  // move focus point on plane
    {
    if(center) F.sub(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  if (keyPressed && key=='T')  // move focus point vertically
    {
    if(center) F.sub(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  change=true;
  }  

// **** Header, footer, help text on canvas
void displayHeader()  // Displays title and authors face on screen
    {
    scribeHeader(title,0); scribeHeaderRight(name); 
    fill(white); 
    
    image(face1, width-140,25,150,180); 
    
    image(face2, width-140,200,150,180); 
    }
void displayFooter()  // Displays help text at the bottom
    {
    scribeFooter(guide,1); 
    scribeFooter(menu,0); 
    }

String title ="CS6491-2019-P1: Tetrabot", name ="Jiewen Wang, Jiabin Fang", // STUDENT: PUT YOUR NAMES HERE !!!
       menu="?:help, !:picture, ~:(start/stop)filming, space:rotate, `/wheel:closer, t/T:target, v:tracking, </>:slower/faster",
       guide="s:screwmotion, d:subdivision, +/-:subdivision degree, t:trajectory, b:body, ^:tets, c:screwaxis"; // user's guide
       //guide="C/S:read/save S from/to cliboard, L/O/R:append, Z:undo, g:go, r:reset, i:interactive, ,/.:forward/backward a:animation, b:body, l:legs, ^:tets, #:tiles"; // user's guide
