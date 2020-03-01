// 6491-2019-P1 
// Base-code: Jarek ROSSIGNAC
// Student 1: ??? ????
// Student 2: ??? ????
import processing.pdf.*;    // to save screen shots as PDFs, does not always work: accuracy problems, stops drawing or messes up some curves !!!
import java.awt.Toolkit;
import java.awt.datatransfer.*;
//  ******************* Basecode for P2 ***********************

int SEQUENCE_MODE = 0;
int INTERACTIVE_MODE = 1;
int currentMode = SEQUENCE_MODE;

int L_CMD = 0;
int R_CMD = 1;
int O_CMD = 2;
int L_DANCE = 3;
int R_DANCE = 4;
int T_DANCE = 5;
int MAX_CMD = 6;

Boolean 
  tracking=true, 
  smooth=false, 
  PickedFocus=false, 
  center=true, 
  track=false, 
  ScrewMotion = true, 
  InterpolationMethod1 = true,
  subdivision = true,
  showscrew = false,
  degreeChange = false,
  showTrajectory = false;
float 
  rBase=200, // radius of base triangle
  t=0, 
  s=0, 
  e=1, 
  bodyRatio=0.25;
int
  f=30, maxf=f, level=4, method=5, cmdIdx=0, Maxlen = 100*2*2*2*2*2*2, FrameCount = 0, frame = 0,
  degree = 3, refineCounter = 3, Frames = 0;

String  //S="ttttttrrrrrrttttlllllltttt"; 
  S="RLRLRRRRLRLRLLLLRLRLRRRRLRLRLLLLRLRLR";
// S="LRLRRRLRLRRRLRLRRRLRLRLRLRRRLRLRRRLRLRRRLRLR"; 
//S="LRLRLRRLRLRRLRLRLRLLRLRLL"; 
// S="LLLRRRLLLRRRORORORLLLRRRLLL";
float defectAngle=0;

PImage face1, face2; // picture of author's face, should be: data/pic.jpg in sketch folder

// subdivision
Point[] outerframeA = new Point[Maxlen];
Point[] outerframeB = new Point[Maxlen];
Point[] outerframeC = new Point[Maxlen];
Point[] outerframeD = new Point[Maxlen];
// for records
Point[][] outerframe = new Point[4][Maxlen];
Point[][] Outerframe = new Point[4][Maxlen];

void setup() {
  size(1000, 1000, P3D); // P3D means that we will do 3D graphics
  //size(600, 600, P3D); // P3D means that we will do 3D graphics
  face1 = loadImage("data/jiewen.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  face2 = loadImage("data/jiabin.jpg");
  textureMode(NORMAL);          
  noSmooth();
  frameRate(60);
  //GHOST.declareBot();
  BOT.declareBot(bodyRatio);
  BOT.resetBot(rBase);
  GHOST.declareBot(bodyRatio);
  GHOST.resetBot(rBase);
  _LookAtPt.reset(BOT.getCentroid().ToPt(), 10);


  if (currentMode == INTERACTIVE_MODE) BOT.isRolling = false;
  else if (currentMode == SEQUENCE_MODE && S.length() > 0) BOT.isRolling = true;  
  // pre-computing
  while (true)
  {
    f++; // advance frame counter
    t = 1.*(f-1)/maxf;
    if (smooth) {
      if (t<0.5) t = sin(PI*t)*0.5;
      else t = sin(PI*(t-0.5))*0.5+0.5;
    }
    if (f==maxf) GHOST.finishRoll();
    if (GHOST.isRolling) GHOST.rollBot(t);
    if (f>maxf) // if end of step
    {
      f=0;
      if (cmdIdx<S.length())
      {
        GHOST.isRolling = true;
        char command = S.charAt(cmdIdx++);
        if (command == 'L')
        {
          GHOST.setCommand(L_CMD);
        } else if (command == 'R')
        {
          GHOST.setCommand(R_CMD);
        } else if (command == 'O' || command == 'o')
        {
          GHOST.setCommand(O_CMD);
        } else if (command == 'l')
        {
          frameRate(60);
          GHOST.setCommand(L_DANCE);
        } else if (command == 'r')
        {
          frameRate(60);
          GHOST.setCommand(R_DANCE);
        } else if (command == 't')
        {
          GHOST.setCommand(T_DANCE);
          frameRate(30);
        }
      } else break;
    }
    // recording
    if(t == 0 || t == 0.5) 
    {
      //outerframea[FrameCount] = GHOST.outerFrame[0];
      //outerframeb[FrameCount] = GHOST.outerFrame[1];
      //outerframec[FrameCount] = GHOST.outerFrame[2];
      //outerframed[FrameCount] = GHOST.outerFrame[3];
      outerframe[0][FrameCount] = GHOST.outerFrame[0];
      outerframe[1][FrameCount] = GHOST.outerFrame[1];
      outerframe[2][FrameCount] = GHOST.outerFrame[2];
      outerframe[3][FrameCount] = GHOST.outerFrame[3];
      
      FrameCount ++;
    }
    
  }
  // subdivision
  Frames = FrameCount;
  Outerframe = RefineTuck(outerframe, Frames);
  for(int refine = 0; refine < refineCounter; refine++) FrameCount = FrameCount * 2 - 1;
  
  //outerframeA = outerframe[0];
  //outerframeB = outerframe[1];
  //outerframeC = outerframe[2];
  //outerframeD = outerframe[3];
  outerframeA = Outerframe[0];
  outerframeB = Outerframe[1];
  outerframeC = Outerframe[2];
  outerframeD = Outerframe[3];
  
  cmdIdx = 0;
  f = maxf;
  
}

void draw() {
  background(255);
  hint(ENABLE_DEPTH_TEST); 
  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
  setView();  // see pick tab
  showFloor(); // draws dance floor as yellow mat
  //doPick(); // sets Of and axes for 3D GUI (see pick Tab)

  // prepare for recording video
  //if (millis() < 1000) {
  //  popMatrix(); 
  //  hint(DISABLE_DEPTH_TEST);
  //  return;
  //}
  
  if(degreeChange){
    FrameCount = Frames;
    Outerframe = RefineTuck(outerframe, Frames);
    for(int refine = 0; refine < refineCounter; refine++) FrameCount = FrameCount * 2 - 1;
   
    outerframeA = Outerframe[0];
    outerframeB = Outerframe[1];
    outerframeC = Outerframe[2];
    outerframeD = Outerframe[3];
    
    print("subdivision degree changed.\n");
    degreeChange = false;
  }
  
  if(showTrajectory)
  {
    fill(green);
    for(int i = 0; i < FrameCount; i++) drawSphere(outerframeA[i], 5); 
    fill(red);
    for(int i = 0; i < Frames; i++) drawSphere(outerframe[0][i], 5);
    BOT.drawBot(5, 15);
    
    popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
    return;
    
  }
  
  if (currentMode == SEQUENCE_MODE)  
  {
    f++; // advance frame counter
    t = 1.*f / maxf;//t = 1.*(f - 1)/maxf;
    if (smooth) {
      if (t<0.5) t = sin(PI*t)*0.5;
      else t = sin(PI*(t-0.5))*0.5+0.5;
    }
    if (f==maxf) BOT.finishRoll();
    if (BOT.isRolling) BOT.rollBot(t);
    if (f>maxf) // if end of step
    {
      f=0;
      if (cmdIdx<S.length())
      {
        BOT.isRolling = true;
        char command = S.charAt(cmdIdx++);
        if (command == 'L')
        {
          BOT.setCommand(L_CMD);
        } else if (command == 'R')
        {
          BOT.setCommand(R_CMD);
        } else if (command == 'O' || command == 'o')
        {
          BOT.setCommand(O_CMD);
        } else if (command == 'l')
        {
          frameRate(60);
          BOT.setCommand(L_DANCE);
        } else if (command == 'r')
        {
          frameRate(60);
          BOT.setCommand(R_DANCE);
        } else if (command == 't')
        {
          BOT.setCommand(T_DANCE);
          frameRate(30);
        }
      } else BOT.isRolling=false;
    }

  } else if (currentMode == INTERACTIVE_MODE)  
  {
    f++; // advance frame counter
    if (f>maxf) // if end of step
    {
      f=0;
      // generate command from mouse
      int cmd = BOT.genCmdFromPoint(Point(Of));
      if (cmd >= 0) BOT.setCommand(cmd);
    }
    
  }

  // subdivision
  if(subdivision)
  {
    BOT.outerFrame[0] = outerframeA[frame];
    BOT.outerFrame[1] = outerframeB[frame];
    BOT.outerFrame[2] = outerframeC[frame];
    BOT.outerFrame[3] = outerframeD[frame];
    Point Center = Average(outerframeA[frame], outerframeB[frame], outerframeC[frame], outerframeD[frame]);
    for (int i = 0; i<4; i++)
        BOT.hips[i] = WeightedSum(1.0-bodyRatio, Center, bodyRatio, BOT.outerFrame[i]);
  }
  
  // advance frame
  frame++;
  if(frame >= FrameCount) frame = 0; 
  
  if (tracking&&!mousePressed)   F =_LookAtPt.move(BOT.getCentroid().ToPt()); // F = BOT.BotCentroid();
  
  
  // draw
  BOT.drawBot(5, 15);
  BOT.drawHistoryTiles(3, 10, cyan);



  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas
  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible
  scribeHeader("S="+S+", path="+S.substring(0, min(cmdIdx, S.length()))+", index="+cmdIdx+", t="+nf(t, 1, 3), 1);

  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if (mousePressed) {
    stroke(cyan); 
    strokeWeight(3); 
    noFill(); 
    ellipse(mouseX, mouseY, 20, 20); 
    strokeWeight(1);
  }
  if (keyPressed) { 
    stroke(red); 
    fill(white); 
    ellipse(mouseX+14, mouseY+20, 26, 26); 
    fill(red); 
    text(key, mouseX-5+14, mouseY+4+20); 
    strokeWeight(1);
  }
  if (scribeText) {
    fill(black); 
    displayHeader();
  } // dispalys header on canvas, including my face
  if (scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if (filming) saveFrame("FRAMES/F"+nf(frameCounter++, 4)+".tif");  // save next frame to make a movie
  //change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  //change=true;
}
